import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import fiona
import geopandas as gpd

plt.ion()


class ItalySetup:
    def __init__(self, nnodes='full', ndays=45, when='future'):  # small (M=10),  medium or large (M=68)

        geodata = pd.read_csv('italy-data/geodata.csv')

        if nnodes == 'full':
            nnodes = len(geodata)

        keep = geodata['name'][:nnodes].to_list()

        geodata = geodata[geodata['name'].isin(keep)].reset_index(drop=True)

        ind2name = geodata.name.to_list()
        mobility = np.loadtxt('italy-data/mobility.txt')
        mobility = mobility[:nnodes, :nnodes]
        pop_node = geodata.population.to_numpy()

        self.ind2name = ind2name
        self.mobility = mobility

        self.pop_node = pop_node
        self.nnodes = len(ind2name)

        self.ndays = ndays

        if when == 'past':
            self.start_date = datetime.date(2020, 1, 20)  # fix lentgh
        elif when == 'future':
            self.start_date = datetime.date(2021, 1, 4)  # fix lentgh
        self.end_date = self.start_date + datetime.timedelta(days=ndays-1)

        self.model_days = pd.date_range(self.start_date, self.end_date, freq='1D')

        mobility_ts = pd.read_csv('italy-data/mobility_ts.csv', header=None, index_col=0, parse_dates=True)
        mobility_ts.columns = np.arange(len(mobility_ts.columns))
        mobility_ts.index = mobility_ts.index.rename('date')
        self.mobility_ts = mobility_ts.iloc[:, :nnodes]
        if when == 'future':
            self.mobility_ts = pd.DataFrame(1, columns=np.arange(len(self.mobility_ts.columns)),
                                            index=self.model_days)

        self.mobintime = self.mobility_ts.resample('1D').mean()
        self.mobintime_arr = self.mobintime.to_numpy().T

        shp = gpd.read_file('italy-data/shp/ProvCM01012019_g_WGS84.shp').sort_values('COD_PROV')

        self.shp = shp.merge(geodata)
        self.pos_node = np.zeros((nnodes, 2))
        self.pos_node[:, 0] = self.shp.centroid.x
        self.pos_node[:, 1] = self.shp.centroid.y

        # Load age stratified data:
        agestrat = pd.read_csv('italy-data/province_age_data.csv', header=1)
        agestrat['Totale'] = agestrat['Totale Maschi'] + agestrat['Totale Femmine']
        agestrat = agestrat[agestrat['Età'] != 'Totale']
        agestrat['Category'] = ''
        agestrat.loc[agestrat['Età'].astype(int) < 16, 'Category'] = '0-16'
        agestrat.loc[agestrat['Età'].astype(int) > 65, 'Category'] = '65+'
        agestrat.loc[agestrat['Category'] == '', 'Category'] = '16-65'
        agestrat = agestrat[['Category', 'Codice provincia', 'Totale', 'Provincia']]
        agestrat = agestrat.groupby(['Category', 'Codice provincia', 'Provincia']).sum().reset_index().sort_values(
            'Codice provincia').reset_index()
        agestrat = agestrat.pivot(index='Codice provincia', columns='Category', values='Totale')
        self.pop_node_ag = agestrat.to_numpy()

        print(f'Loaded Italy Setup with {self.nnodes} nodes.')



if __name__ == '__main__':
    # Variables passed to R
    s = ItalySetup(13)
    ind2name = s.ind2name
    mobility = s.mobility
    ic = s.ic
    ind_to_plot = s.ind_to_plot
    pop_node = s.pop_node
    pos_node = s.pos_node
    nnodes = s.nnodes

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
