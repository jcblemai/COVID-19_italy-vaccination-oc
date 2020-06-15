import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime

plt.ion()


class ItalySetup:
    def __init__(self, problem_size='full'):  # small (M=10),  medium or large (M=68)

        geodata = pd.read_csv('italy-data/geodata.csv')

        if problem_size == 'full':
            problem_size = len(geodata)

        keep = geodata['name'][:problem_size].to_list()

        geodata = geodata[geodata['name'].isin(keep)].reset_index(drop=True)

        ind2name = geodata.name.to_list()
        nnodes = len(geodata)
        mobility = np.loadtxt('italy-data/mobility.txt')
        mobility = mobility[:problem_size, :problem_size]
        pop_node = geodata.population.to_numpy()
        pos_node = np.zeros((nnodes, 2))

        ic = {'I': np.zeros(nnodes),
              'R': np.zeros(nnodes),
              'B': np.zeros(nnodes),
              'V': np.zeros(nnodes)}
        seeding_place = 'Genova'
        ic['I'][ind2name.index(seeding_place)] = 0.05 * pop_node[ind2name.index(seeding_place)]

        to_plot = ['Torino', 'Genova']

        self.ind2name = ind2name
        self.mobility = mobility
        self.ic = ic

        self.ind_to_plot = [ind2name.index(place) for place in to_plot]
        self.pop_node = pop_node
        self.pos_node = pos_node
        self.nnodes = len(ind2name)

        # FIX: This is to be set manualy in the R code.
        self.start_date = datetime.date(2017, 12, 31)  # fix lentgh
        self.end_date = datetime.date(2018, 12, 31)
        self.end_sim = datetime.date(2019, 12, 31)

        print(f'Created setup with {self.nnodes} nodes.')


if __name__ == '__main__':
    # Variables passed to R
    setup = ItalySetup(10)
    ind2name = setup.ind2name
    mobility = setup.mobility
    ic = setup.ic
    ind_to_plot = setup.ind_to_plot
    pop_node = setup.pop_node
    pos_node = setup.pos_node
    nnodes = setup.nnodes

# print(ocp.arg['p']['p','scale_v']*sum( y**2 for x in opt['u',:,:,'v'] for y in x ))
# print(ocp.arg['p']['p','scale_ell']*( y for x in opt['x',:,:-1,'I'] for y in x ))
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

# G = mobility_graph(mobility, ind2name, pos_node, pop_node, ocp.opt)

# print("digraph has %d nodes with %d edges" % (nx.number_of_nodes(G), nx.number_of_edges(G)))

# # draw with matplotlib/pylab
# plt.figure(figsize=(10, 10))
# # with nodes colored by degree sized by population

# nx.draw(G, 
#         G.position, 
#         node_size=1000./max(pop_node) * np.array([G.population[v] for v in G]),
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
#         G.comp[node] = (ic['S'][i],ic['I'][i],ic['R'][i])

#         for j, connection in enumerate(mobility[i]):
#             G.add_edge(node, ind2name[j], weight=connection)

#     return G

# G = mobility_graph(mob_cdr, ind2name, pos_node, pop_node, ic)

# print("digraph has %d nodes with %d edges" % (nx.number_of_nodes(G), nx.number_of_edges(G)))

# # draw with matplotlib/pylab
# plt.figure(figsize=(10, 10))
# # with nodes colored by degree sized by population

# nx.draw(G, 
#         G.position, 
#         node_size=1000/max(pop_node) * np.array([G.population[v] for v in G]),
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


# G = mobility_graph(mob_gravity, ind2name, pos_node, pop_node, ic)

# print("digraph has %d nodes with %d edges" % (nx.number_of_nodes(G), nx.number_of_edges(G)))


# # draw with matplotlib/pylab
# plt.figure(figsize=(10, 10))
# # with nodes colored by degree sized by population

# nx.draw(G, 
#         G.position, 
#         node_size=1000/max(pop_node) * np.array([G.population[v] for v in G]),
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
