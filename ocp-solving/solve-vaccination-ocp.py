import casadi as ca
import casadi.tools as cat
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import networkx as nx
import sys

sys.path.append('.')
from COVIDVaccinationOCP import COVIDVaccinationOCP, rk4_mob, rhs_py_total
from ItalySetup import ItalySetup
from scipy.integrate import solve_ivp

s = ItalySetup(problem_size=13)

ocp_params = {'N': 52,  # Number of control time interval
              'T': 52 * 7,  # N * 7
              'n_int_steps': 10,
              'scaling': 1  # Scaling on population, to work on similar quantities
              }

scaling = ocp_params['scaling']
N = ocp_params['N']
T = ocp_params['T']
n_int_steps = ocp_params['n_int_steps']

# Deprecated params, used with MARIO
model_params = {'sigma': 0.7, 'beta': 0.2, 'mu_b': 0.03615616438, 'gamma': 0.5, 'theta': 3.4476623459780395e-7,
                'lam': 800, 'mu': 4.53e-5, 'rho_v': 1 / (5 * 365), 'rho': 1 / (8 * 365), 'alpha': 0.004002739726,
                'r': 1, 'm': 0.3}

obj_params = {
    'scale_ell': 1e0,
    'scale_If': 1e2,
    'scale_v': 1e-6
}

# Build optimal control problem, can be long...
ocp = COVIDVaccinationOCP(
    N=ocp_params['N'],
    T=ocp_params['T'],
    n_int_steps=ocp_params['n_int_steps'],
    scaling=ocp_params['scaling'],
    setup=s,
    model_params=model_params,
    obj_params=obj_params,
    plot_iterates=True,
    optimize=False
)

# Set initial conditions as constraints
for name in ocp.ic.keys():
    for i in range(s.nnodes):
        ocp.arg['lbx']['x', i, 0, name] = ocp.arg['ubx']['x', i, 0, name] = ocp.ic[name][i]

# Upper bound for the vaccination rate
ocp.arg['ubx']['u', :, :, 'v'] = 15000 * scaling  # TODO: zero if not optimize, duh
# Upper bound of the vaccination rate for a specific place.
# ocp.arg['ubx']['u',s.ind2name.index('Verona'),:,'v'] = 10000  * scaling

# Scaling coef of the different objectives.
ocp.arg['p']['p', 'scale_v'] = 1e-4 / scaling
ocp.arg['p']['p', 'scale_ell'] = 1e2
ocp.arg['p']['p', 'scale_If'] = 1e6

# Changing parameters on the fly
# ocp.arg['p']['p','lam'] = 6.8

# No action for (approximately) the first week
# ocp.arg['lbx']['u',:,:int(7/(T/N)) ,'v'] = 0.
# ocp.arg['ubx']['u',:,:int(7/(T/N)) ,'v'] = 0.


init = ocp.Vars(0)
for t in range(len(s.ind2name)):
    for name in ocp.states.keys():
        if name not in ['B', 'I']:
            for k in range(N + 1):
                init['x', t, k, name] = ocp.ic[name][t]
# init['u'] = 0.


x0 = []
for i in range(len(s.ind2name)):
    for name in ocp.states.keys():
        x0 += [ocp.ic[name][i]]
x0 = np.squeeze(ca.veccat(*x0))
pnum = np.array(ocp.arg['p']['p'])

init = ocp.Vars(0)
# This init with the highest vaccination
for i in range(len(s.ind2name)):
    k_max = int(s.pop_node[i] * scaling / ocp.arg['ubx']['u', i, -1, 'v'] / T * N - 1e-8)
    init['u', i, :k_max, 'v'] = ocp.arg['ubx']['u', i, :k_max, 'v']

sol0 = {}
dt = T / N / n_int_steps
for i in range(len(s.ind2name)):
    for name in ocp.states.keys():
        init['x', i, 0, name] = ocp.ic[name][i]
sim = np.array([x0]).T

# Rough integration with fixed mobility? as a starting point
for k in range(N):
    uk = np.squeeze(ca.veccat(*init['u', :, k, 'v']))
    covk = np.squeeze(ca.veccat(*ocp.arg['p']['cov', :, k]))

    mob = []
    for i in range(s.nnodes):
        # init['u',i,k,'mob'] = sum( s.mobility[i,j]*sim[j*4+2,-1]/(1+sim[j*4+2,-1]) for j in range(s.nnodes) )
        mob.append(
            float(sum(s.mobility[i, j] * sim[j * 4 + 2, -1] / (1 + sim[j * 4 + 2, -1]) for j in range(s.nnodes))))

    # sol0 = solve_ivp(lambda t,y: rhs_py_total(t, y, uk, covk, pnum, s.nnodes, s.mobility),[0,T/N], S[:,-1])
    # sol0 = solve_ivp(lambda t,y: rhs_py_total_mob(t, y, uk, covk, pnum, s.nnodes, s.mobility, mob),[0,T/N], S[:,-1])
    sol_ = rk4_mob(dt, sim[:, -1], uk, covk, pnum, s.nnodes, s.mobility, n_int_steps, mob, s.pop_node * scaling)
    sol0['y'] = np.array([sol_]).T
    sim = np.append(sim, np.array([sol0['y'][:, -1]]).T, axis=1)
    for i in range(len(s.ind2name)):
        init['x', i, k + 1] = np.array([sol0['y'][i * 4:(i + 1) * 4, -1]])

[fnum, gnum] = ocp.nlpFun(init, ocp.arg['p'])
gnum = ocp.g(gnum)

ocp.arg['x0'] = init
# ocp.arg['x0'] = opt

ocp.arg['ubg']['vaccines'] = 3e6 * scaling

# ocp.arg['lbx']['u',:,:,'v'] = init['u',:,:,'v']
# ocp.arg['ubx']['u',:,:,'v'] = init['u',:,:,'v']

# ------------------------------------- SOLVING OCP --------------------------------------
opt = ocp.solveOCP()
# ------------------------------------- SOLVING OCP --------------------------------------

print(ocp.costTerms(opt, ocp.arg['p']))

[fnum, gnum] = ocp.nlpFun(ocp.opt, ocp.arg['p'])
print(f"""
Vaccines stockpile: 
    {float(ocp.arg['ubg']['vaccines'] / scaling):010f} total.
    {float(ocp.g(gnum)['vaccines'] / scaling):010f} spent.
    {float((ocp.arg['ubg']['vaccines'] - ocp.g(gnum)['vaccines']) / scaling):010f} left.""")

# ocp.arg['p']['p','scale_v'] = 1e-10
# ocp.arg['x0'] = opt
# ocp.arg['lam_g0'] = ocp.lam_g
# ocp.arg['lam_x0'] = ocp.lam_x
# opt = ocp.solveOCP()
# print(ocp.costTerms(opt,ocp.arg['p']))

pnum = np.array(ocp.arg['p']['p'])
x0 = np.squeeze(ca.veccat(*opt['x', :, 0]))

dt = T / N / n_int_steps
sim = np.array([x0]).T
sim1 = np.array([x0]).T
sim2 = np.array([x0]).T
for k in range(N):
    uk = np.squeeze(ca.veccat(*opt['u', :, k, 'v']))
    covk = np.squeeze(ca.veccat(*ocp.arg['p']['cov', :, k]))

    mob = []
    for i in range(s.nnodes):
        mob.append(
            float(sum(s.mobility[i, j] * sim2[j * 4 + 2, -1] / (1 + sim2[j * 4 + 2, -1]) for j in range(s.nnodes))))
    # x_ = rk4(dt,sim1[:,-1],uk, covk, pnum, s.nnodes, s.mobility,n_int_steps,s.pop_node*scaling)
    sol0 = solve_ivp(lambda t, y: rhs_py_total(t, y, uk, covk, pnum, s.nnodes, s.mobility, s.pop_node * scaling),
                     [0, T / N], sim[:, -1], rtol=1e-8, atol=1e-8)
    # sol1 = solve_ivp(lambda t,y: rhs_py_total_mob(t, y, uk, covk, pnum, s.nnodes,
    # s.mobility, mob, s.pop_node*scaling),[0,T/N], sim2[:,-1],rtol=1e-8,atol=1e-8)

    sim = np.append(sim, np.array([sol0['y'][:, -1]]).T, axis=1)
    # sim1 = np.append( sim1, np.array([x_]).T, axis=1)
    # sim2 = np.append( sim2, np.array([sol1['y'][:,-1]]).T, axis=1)
# (sol0['y'][:,-1]-np.squeeze(ca.veccat(*opt['x',:,k+1])))/abs(np.squeeze(ca.veccat(*opt['x',:,k+1])))


T = ocp.T
N = ocp.N


def plot_traj(name, T, N, opt, sim, scaling, pop_node):
    dT = T / N
    for i in [list(s.ind2name).index(name)]:
        fig, axes = plt.subplots(2, 3, num=2 + i + 1000)
        for ax in axes.flat: ax.cla()
        S = pop_node[i] * scaling - (
                np.array(ca.veccat(*opt['x', i, :, 'I'])) + np.array(ca.veccat(*opt['x', i, :, 'R'])) + np.array(
            ca.veccat(*opt['x', i, :, 'V'])))
        axes[0, 0].plot(np.arange(0, T + dT, dT), S / scaling, label="S", color='orange', linewidth=1)
        axes[0, 1].plot(np.arange(0, T + dT, dT), np.array(ca.veccat(*opt['x', i, :, 'I'])) / scaling, label="I",
                        color='darkred', linewidth=2)
        axes[0, 2].plot(np.arange(0, T + dT, dT), np.array(ca.veccat(*opt['x', i, :, 'R'])) / scaling, label="R",
                        color='red', linewidth=1)
        axes[1, 0].plot(np.arange(0, T + dT, dT), np.array(ca.veccat(*opt['x', i, :, 'B'])), label="B", color='blue',
                        linewidth=1)
        axes[1, 1].step(np.arange(0, T + dT, dT),
                        np.array(ca.veccat(ca.veccat(*opt['u', i, :, 'v']), opt['u', i, -1, 'v'])) / scaling, 'k',
                        label=r"$\nu(t)$", where="post")
        # axes[1,2].step(np.arange(0, T+dT, dT),np.array(arg['p']['p',:,i]+[arg['p']['p',-1,i]]),'b',label=r"$J(t)$",where="post")
        axes[1, 2].plot(np.arange(0, T + dT, dT), np.array(ca.veccat(*opt['x', i, :, 'V'])) / scaling, label="V",
                        color='orange', linewidth=2)

        Ssim = pop_node[i] * scaling - (sim[i * 4 + 0, :] + sim[i * 4 + 1, :] + sim[i * 4 + 3, :])
        axes[0, 0].plot(np.arange(0, T + dT, dT), Ssim / scaling, label="Sm", color='orange', linewidth=1,
                        linestyle=':')
        axes[0, 1].plot(np.arange(0, T + dT, dT), sim[i * 4 + 0, :] / scaling, label="Im", color='darkred', linewidth=1,
                        linestyle=':')
        axes[0, 2].plot(np.arange(0, T + dT, dT), sim[i * 4 + 1, :] / scaling, label="Rm", color='red', linewidth=1,
                        linestyle=':')
        axes[1, 0].plot(np.arange(0, T + dT, dT), sim[i * 4 + 2, :], label="Bm", color='blue', linewidth=1,
                        linestyle=':')
        axes[1, 2].plot(np.arange(0, T + dT, dT), sim[i * 4 + 3, :] / scaling, label="Vm", color='orange', linewidth=1,
                        linestyle=':')
        for ax in axes.flat:
            ax.legend(loc="upper right")
            ax.grid()
        fig.suptitle(s.ind2name[i])
        # fig.tight_layout()


plot_traj('Torino', T, N, opt, sim, scaling, s.pop_node)
plot_traj('Genova', T, N, opt, sim, scaling, s.pop_node)
# plot_traj('MOMBASA',T,N,opt,sim,scaling,s.pop_node)

plt.draw()

# print(ocp.arg['p']['p','scale_v']*sum( y**2 for x in opt['u',:,:,'v'] for y in x ))
# print(ocp.arg['p']['p','scale_ell']*sum( y for x in opt['x',:,:-1,'I'] for y in x ))
# print(ocp.arg['p']['p','scale_If']*sum( x for x in opt['x',:,-1,'I'] ))


# def mobility_graph(mobility, ind2name, pos_node, pop_node, opt):

#     G = nx.Graph()
#     G.position = {}
#     G.population = {}
#     G.comp = {}

#     for i, node in enumerate(ind2name):

#         G.add_node(node)

#         G.position[node] = (pos_node[i,0], pos_node[i,1])
#         G.population[node] = pop_node[i]
#         G.comp[node] = (sum(np.array(ca.veccat(ca.veccat(*opt['u',i,:,'v']),opt['u',i,-1,'v']))))
#         #BUG: sum over time means depends on control interval scale 1/7

#         for j, connection in enumerate(mobility[i]):
#             G.add_edge(node, ind2name[j], weight=connection)

#     return G

# G = mobility_graph(s.mobility, ind2name, pos_node, s.pop_node, ocp.opt)

# print("digraph has %d nodes with %d edges" % (nx.number_of_nodes(G), nx.number_of_edges(G)))

# # draw with matplotlib/pylab
# plt.figure(figsize=(10, 10))
# # with nodes colored by degree sized by population

# nx.draw(G, 
#         G.position, 
#         node_size=1000./max(s.pop_node) * np.array([G.population[v] for v in G]),
#         #node_color=[float(G.degree(v)) for v in G],
#         #node_color=[G.population[v] for v in G],
#         node_color=[G.comp[v][0]/(G.population[v] * scaling) for v in G],
#         width = 10* np.array([a['weight'] for u,v,a in G.edges(data=True)]),
#         edge_color= np.array([a['weight'] for u,v,a in G.edges(data=True)]),
#         edge_cmap = mpl.cm.viridis,
#         node_cmap = mpl.cm.viridis,
#         with_labels=True)

#     # scale the axes equally
# plt.xlim(-5000, 500)
# plt.ylim(-2000, 3500)

# cmap= mpl.cm.viridis
# vmin = min([G.comp[v][0]/(G.population[v] * scaling) for v in G])
# vmax = max([G.comp[v][0]/(G.population[v] * scaling) for v in G])
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin = vmin, vmax=vmax))
# sm._A = []
# plt.colorbar(sm)

# plt.draw()


# def mobility_graph(mobility, ind2name, pos_node, pop_node, comp):

#     G = nx.Graph()
#     G.position = {}
#     G.population = {}
#     G.comp = {}

#     for i, node in enumerate(ind2name):

#         G.add_node(node)

#         G.position[node] = (pos_node[i,0], pos_node[i,1])
#         G.population[node] = pop_node[i]
#         G.comp[node] = (ocp.ic['S'][i],ocp.ic['I'][i],ocp.ic['R'][i])

#         for j, connection in enumerate(mobility[i]):
#             G.add_edge(node, ind2name[j], weight=connection)

#     return G

# G = mobility_graph(mob_cdr, ind2name, pos_node, s.pop_node, ocp.ic)

# print("digraph has %d nodes with %d edges" % (nx.number_of_nodes(G), nx.number_of_edges(G)))

# # draw with matplotlib/pylab
# plt.figure(figsize=(10, 10))
# # with nodes colored by degree sized by population

# nx.draw(G, 
#         G.position, 
#         node_size=1000/max(s.pop_node) * np.array([G.population[v] for v in G]),
#         #node_color=[float(G.degree(v)) for v in G],
#         node_color=[G.population[v] for v in G],
#         #node_color=[G.comp[v][1] for v in G],
#         width = 40* np.array([a['weight'] for u,v,a in G.edges(data=True)]),
#         edge_color=10* np.array([a['weight'] for u,v,a in G.edges(data=True)]),
#         edge_cmap = mpl.cm.viridis,
#         #with_labels=False,
#         with_labels=True)
#     # scale the axes equally
# #plt.xlim(-5000, 500)
# #plt.ylim(-2000, 3500)

# plt.draw()


# G = mobility_graph(mob_gravity, ind2name, pos_node, s.pop_node, ocp.ic)

# print("digraph has %d nodes with %d edges" % (nx.number_of_nodes(G), nx.number_of_edges(G)))


# # draw with matplotlib/pylab
# plt.figure(figsize=(10, 10))
# # with nodes colored by degree sized by population

# nx.draw(G, 
#         G.position, 
#         node_size=1000/max(s.pop_node) * np.array([G.population[v] for v in G]),
#         #node_color=[float(G.degree(v)) for v in G],
#         node_color=[G.population[v] for v in G],
#         #node_color=[G.comp[v][1] for v in G],
#         width = 40* np.array([a['weight'] for u,v,a in G.edges(data=True)]),
#         edge_color=10* np.array([a['weight'] for u,v,a in G.edges(data=True)]),
#         edge_cmap = mpl.cm.viridis,
#         with_labels=False)
#         #with_labels=True)
#     # scale the axes equally
# #plt.xlim(-5000, 500)
# #plt.ylim(-2000, 3500)

# #plt.draw()"
