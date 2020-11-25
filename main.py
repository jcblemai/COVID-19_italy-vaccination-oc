import copy
import casadi as ca
import matlab.engine
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx
from ItalySetup import ItalySetup
from covidOCP import ocp_utils
from covidOCP import COVIDVaccinationOCP

nx = 9
states_names = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']

eng = matlab.engine.start_matlab()
eng.cd('geography-paper-master/', nargout=0)

model_size = 107  # nodes

# Horizon for each problem

T = 31
N = T - 1
n_int_steps = 1

setup = ItalySetup(model_size)
M = setup.nnodes

freq = '1D'  # Unsupported to change this
model_days = pd.date_range(setup.start_date, setup.end_date, freq=freq)

eng.run('single_build.m', nargout=0)
integ_matlab = np.array(eng.eval('x'))


class OCParameters:
    def __init__(self):
        self.mobintime = setup.mobility_ts.resample(freq).mean()
        p_dict, self.mobfrac, self.mobmat, self.betaratiointime, self.x0 = ocp_utils.get_parameters_from_matlab(eng,
                                                                                                                setup,
                                                                                                                model_size,
                                                                                                                model_days,
                                                                                                                freq)

        self.params_structural = {'betaP0': p_dict['betaP0'],
                                  'epsilonA': p_dict['epsilonA'],
                                  'epsilonI': p_dict['epsilonI'],
                                  'r': p_dict['r']}

        p_dict.pop('betaP0')
        p_dict.pop('r')
        p_dict.pop('epsilonA')
        p_dict.pop('epsilonI')

        self.model_params = copy.copy(p_dict)

        self.hyper_params = {
            'scale_ell': 1e0,
            'scale_If': 0,
            'scale_v': 1e-6
        }

        mob_prun = 0.0006
        self.mobmat_pr = self.prune_mobility(mob_prun)

    def prune_mobility(self, mob_prun):
        mobK = self.mobintime.to_numpy().T[:, 0]
        betaR = self.betaratiointime.to_numpy().T[:, 0]
        C = self.params_structural['r'] * self.mobfrac.flatten() * mobK * self.mobmat
        np.fill_diagonal(C, 1 - C.sum(axis=1) + C.diagonal())
        print(f'pruning {C[C < mob_prun].size} non-diagonal mobility elements of {C.size - M}.')
        C[C < mob_prun] = 0  # Prune elements
        mobmat_pr = np.copy(self.mobmat)

        mobmat_pr[C == 0] = 0
        print(f'nnz before: {np.count_nonzero(self.mobmat)}, after: {np.count_nonzero(mobmat_pr)}')

        return mobmat_pr

    def get_pvector(self):
        pvector = ocp_utils.params_to_vector(self.model_params)
        pvector.append(ocp_utils.params_to_vector(self.hyper_params))
        pvector_names = list(self.model_params.keys()) + list(self.hyper_params.keys())
        return pvector, pvector_names


p = OCParameters()

ocp = COVIDVaccinationOCP.COVIDVaccinationOCP(N=N, T=T, n_int_steps=n_int_steps,
                                              s=setup, p=p)


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
