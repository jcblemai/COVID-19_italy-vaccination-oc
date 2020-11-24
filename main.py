import copy

import casadi as ca
import matlab.engine
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ItalySetup import ItalySetup
from covidOCP import utils
import networkx

eng = matlab.engine.start_matlab()
eng.cd('geography-paper-master/', nargout=0)

solve_every = 10  # days
model_size = 107  # node s

# Horizon for each problem
N = 30
T = 31

setup = ItalySetup(model_size)

n_int_steps = 1
freq = '1D'  # Unsupported to change this
model_days = pd.date_range(setup.start_date, setup.end_date, freq=freq)
model_step = pd.date_range(setup.start_date, setup.end_date, freq=freq)
mobintime = setup.mobility_ts.resample(freq).mean()
# N = len(model_step) - 1
# T = len(model_days)


eng.run('single_build.m', nargout=0)

nx = 9
states_list = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']
states = states_list
S, E, P, I, A, Q, H, R, V = np.arange(nx)

integ_matlab = np.array(eng.eval('x'))

p_dict, mobfrac, mobmat, betaratiointime, x0 = utils.get_parameters_from_matlab(eng, setup, model_size, model_days,freq)

dt = T / N / n_int_steps

obj_params = {
    'scale_ell': 1e0,
    'scale_If': 0,
    'scale_v': 1e-6
}
betaP0 = p_dict['betaP0']
epsilonA = p_dict['epsilonA']
epsilonI = p_dict['epsilonI']
r = p_dict['r']

p_dict.pop('betaP0')
p_dict.pop('r')
p_dict.pop('epsilonA')
p_dict.pop('epsilonI')
C = mobmat

model_params = copy.copy(p_dict)
ind2name = setup.ind2name

pop_node = setup.pop_node

M = setup.nnodes
mobility = setup.mobility

c = setup.mobility

scale_ell = obj_params['scale_ell']
scale_If = obj_params['scale_If']
scale_v = obj_params['scale_v']
pnum = utils.params_to_vector(model_params)

pnum.append(scale_ell)
pnum.append(scale_If)
pnum.append(scale_v)


# arg['ubx']['u', :, :, 'v']  = 0
# If restart:
# for i, name in enumerate(states.keys()):
#    for k in range(N + 1):
#        for nd in range(M):
#            init['x', nd, k, name] = opt['x',nd, k, name]
# for k in range(N):
#    for nd in range(M):
#        init['u', nd, k, 'v'] = opt['u',nd,k,'v']


def build_graph(setup, opt, mobmat):
    G = networkx.Graph()
    G.position = {}
    G.population = {}
    G.comp = {}
    G.epi = {}
    setup.shp['vacc'] = np.nan
    setup.shp['Rend'] = np.nan
    for i, node in enumerate(setup.ind2name):
        G.add_node(node)
        G.position[node] = (setup.pos_node[i, 0], setup.pos_node[i, 1])
        G.population[node] = setup.pop_node[i]
        # G.comp[node] = (ocp.ic['S'][i], ocp.ic['I'][i],ocp.ic['R'][i])
        try:
            G.epi[node] = {'vacc': sum(np.array(ca.veccat(ca.veccat(*opt['u', i, :, 'v']))))[0],
                           'Rend': float(opt['x', i, -1, 'R'])}
            setup.shp.loc[i, 'vacc'] = sum(np.array(ca.veccat(ca.veccat(*opt['u', i, :, 'v']))))[0]
            setup.shp.loc[i, 'Rend'] = float(opt['x', i, -1, 'R'])
        except NameError as e:
            # print(f'epi data failed, {e}')
            G.epi[node] = {'vacc': np.nan,
                           'Rend': np.nan}
            setup.shp.loc[i, 'vacc'] = np.nan
            setup.shp.loc[i, 'Rend'] = np.nan

        setup.shp.loc[i, 'population'] = setup.pop_node[i]  # overwrite
        for j, connection in enumerate(mobmat[i]):
            if connection != 0:
                G.add_edge(node, setup.ind2name[j], weight=connection)
    return G


# G.number_of_edges()


def plot_graph(G):
    fig, ax = plt.subplots(1, 1, figsize=(20, 20))

    # noinspection PyUnresolvedReferences
    networkx.draw(G,
                  G.position,
                  node_size=1000 / max(setup.pop_node) * np.array([G.population[v] for v in G]),
                  # node_color=[float(G.degree(v)) for v in G],
                  # node_color=[G.population[v] for v in G],
                  node_color=[G.epi[v]['vacc'] / G.population[v] for v in G],
                  width=200 * np.array([max(a['weight'], 0.001) for u, v, a in G.edges(data=True)]),
                  edge_color=10 * np.array([a['weight'] for u, v, a in G.edges(data=True)]),
                  edge_cmap=mpl.cm.viridis,
                  ax=ax,
                  with_labels=False
                  )

    #     # scale the axes equally
    # plt.xlim(min(setup.pos_node[:,0]) - 100000, max(setup.pos_node[:,0])+ 100000)
    # plt.ylim(min(setup.pos_node[:,1]) - 100000, max(setup.pos_node[:,1])+ 100000)

    # setup.shp.plot(ax = ax, column='' cmap='OrRd', facecolor="none", edgecolor="black")

    setup.shp.boundary.plot(ax=ax, edgecolor="black", linewidth=.11)

    plt.draw()


def plot_cloropleth():
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    setup.shp.plot(ax=ax, column='Rend', cmap='OrRd')  # ,  edgecolor="black") #facecolor="none",


def plotscatter():
    import seaborn as sns
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    plt.scatter(setup.shp['vacc'] / setup.shp['population'], setup.shp['Rend'] / setup.shp['population'],
                c=setup.shp['population'])
    ax.set_xlabel("prop. vaccinated")
    ax.set_ylabel("prop. recovered")
    ax.set_xlim(0)
    ax.set_ylim(0, 0.0002)

    sns.scatterplot(setup.shp['vacc'], setup.shp['population'] * 100, hue=setup.shp['population'])
