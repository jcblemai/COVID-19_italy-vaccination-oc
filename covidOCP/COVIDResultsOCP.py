


class OCPplots:
    """ plot both from csv and from program"""
    def __init__(self, setup, result = None):
        self.setup = setup




    def
        fig, axes = plt.subplots(5, 2, figsize=(10, 10))
        for i, st in enumerate(states_names):
            for k in range(M):
                axes.flat[i].plot(np.array(ca.veccat(*opt['x', k, :til, st])), lw=2, ls='--')
                # if st != 'V':
                #    axes.flat[i].plot(np.array(integ_matlab.T[k+107*i,:til].T),
                #           lw = .5)
                axes.flat[i].set_title(st)
                axes.flat[-1].step(np.arange(len(np.array(ca.veccat(ca.veccat(*opt['u', k, :til, 'v']))))),
                                   np.array(ca.veccat(ca.veccat(*opt['u', k, :til, 'v'])))
                               )


        def plot_node(self, node):
            fig, axes = plt.subplots(2, 5, figsize=(20, 10))
            fig.patch.set_facecolor('white')

            til = self.T

            for i, st in enumerate(states_names):
                axes.flat[i].plot(np.array(ca.veccat(*self.opt['x', node, :til, st])),
                                  linestyle=':', lw=4, color='r')
                # if st != 'V':
                # axes.flat[i].plot(np.array(integ_matlab.T[node + 107 * i, :til].T),
                #                  linestyle='-', lw=2, color='k')

                axes.flat[i].set_title(st);

            axes.flat[-1].step(np.array(ca.veccat(ca.veccat(*self.opt['u', node, :til, 'v']))),  # ,opt['u',node,-1,'v'])),
                               'k', label=r"$\nu(t)$")


        def plot_all(self, opt):
            til = self.T
            fig, axes = plt.subplots(5, 2, figsize=(10, 10))
            for i, st in enumerate(states_names):
                for k in range(self.M):
                    axes.flat[i].plot(np.array(ca.veccat(*opt['x', k, :til, st])), lw=2, ls='--')
                    # if st != 'V':
                    #    axes.flat[i].plot(np.array(integ_matlab.T[k+107*i,:til].T),
                    #           lw = .5)
                    axes.flat[i].set_title(st)
                    axes.flat[-1].step(np.arange(len(np.array(ca.veccat(ca.veccat(*opt['u', k, :til, 'v']))))),
                                       np.array(ca.veccat(ca.veccat(*opt['u', k, :til, 'v'])))
                                       )


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
