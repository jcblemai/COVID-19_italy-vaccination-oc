import casadi as ca
import casadi.tools as cat
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm
from timeit import default_timer as timer

if "Agg" not in mpl.get_backend():
    mpl.interactive(True)
plt.ion()

nx = 9
states_names = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']


class PlotIterates(ca.Callback):
    def __init__(self, name, nx, ng, np, ind_to_plot, T, N, V, ind2name, mobility, pos_node, pop_node, opts={}):
        ca.Callback.__init__(self)

        self.ng = ng
        self.np = np
        self.T = T
        self.N = N
        self.V = V
        self.dT = T / N
        self.ind_to_plot = ind_to_plot
        self.ind2name = ind2name
        self.sol = 0
        self.mobility = mobility
        self.pos_node = pos_node
        self.pop_node = pop_node
        self.nx = nx

        self.fig = []
        self.axes = []
        for k, i in enumerate(self.ind_to_plot):
            figk, axesk = plt.subplots(2, 5, num=2 + i)
            for ax in axesk.flat:
                ax.autoscale(enable=True, axis='both')
            self.fig.append(figk)
            self.axes.append(axesk)

        figk, axesk = plt.subplots(1, 1, num=100)
        self.fig.append(figk)
        self.axes.append(axesk)
        plt.draw()
        plt.show()
        plt.pause(0.01)
        self.construct(name, opts)

    def get_n_in(self):
        return ca.nlpsol_n_out()

    def get_n_out(self):
        return 1

    def get_name_in(self, i):
        return ca.nlpsol_out(i)

    def get_name_out(self, i):
        return "ret"

    def get_sparsity_in(self, i):
        n = ca.nlpsol_out(i)
        if n == 'f':
            return ca.Sparsity.scalar()
        elif n in ('x', 'lam_x'):
            return ca.Sparsity.dense(self.nx)
        elif n in ('g', 'lam_g'):
            return ca.Sparsity.dense(self.ng)
        elif n in ('p', 'lam_p'):
            return ca.Sparsity.dense(self.np)
        else:
            return ca.Sparsity(0, 0)

    def eval(self, arg):
        darg = {}
        for (i, s) in enumerate(ca.nlpsol_out()):
            darg[s] = arg[i]

        sol = self.V(darg['x'])
        self.sol = sol
        if hasattr(self, 'lines'):
            for k, i in enumerate(self.ind_to_plot):
                for isn, sn in enumerate(states_names):
                    self.lines[k][isn][0].set_ydata(np.array(ca.veccat(*sol['x', i, :, sn])))

                self.lines[k][-1][0].set_ydata(
                    np.array(ca.veccat(ca.veccat(*sol['u', i, :, 'v']), sol['u', i, -1, 'v'])))
                for j, ax in enumerate(self.axes[k].flat):
                    ax.relim()
                    ax.autoscale_view()
                    # ax.update_datalim(self.lines[k*5+j][0].get_xydata(), updatex=True, updatey=True)
                # plt.draw()
                self.fig[k].canvas.draw_idle()
                plt.pause(0.1)
        else:
            self.lines = [None] * len(list(self.ind_to_plot))
            for k, i in enumerate(self.ind_to_plot):
                self.lines[k] = []
                for ax in self.axes[k].flat: ax.cla()
                T = self.T
                dT = self.dT
                for isn, sn in enumerate(states_names):
                    self.lines[k].append(self.axes[k].flat[isn].plot(np.arange(0, T + dT, dT),
                                                                     np.array(ca.veccat(*sol['x', i, :, sn])),
                                                                     label=sn, linewidth=2))

                self.lines[k].append(self.axes[k].flat[-1].step(np.arange(0, T + dT, dT),
                                                                np.array(ca.veccat(ca.veccat(*sol['u', i, :, 'v']),
                                                                                   sol['u', i, -1, 'v'])),
                                                                'k', label=r"$\nu(t)$",
                                                                where="post"))

                for ax in self.axes[k].flat:
                    # ax.legend(loc="upper right")
                    ax.grid('on')
                    ax.autoscale(enable=True, axis='both')

                try:
                    self.fig[k].suptitle(self.ind2name[i])
                except:
                    self.fig[k].suptitle(i)

                self.fig[k].canvas.draw_idle()
                plt.pause(0.1)

        if False:
            plt.figure(100)
            self.fig[-1].clear()
            G = mobility_graph(self.mobility, self.ind2name, self.pos_node, self.pop_node, sol, self.N)
            networkx.draw(G,
                          G.position,
                          node_size=0.1 / scaling * np.array([G.inf[v] for v in G]),
                          node_color=np.array([G.vac[v] / G.pop[v] for v in G]),
                          width=5 * np.array([a['weight'] for u, v, a in G.edges(data=True)]),
                          edge_color=np.array([a['exchange'] for u, v, a in G.edges(data=True)]),
                          node_cmap=mpl.cm.viridis,
                          edge_cmap=mpl.cm.viridis,
                          with_labels=False)

            plt.draw()
            plt.pause(0.1)

        return [0]


def rhs_py(t, x, u, cov, p, mob, pop_node, p_foi):
    S, E, P, I, A, Q, H, R, V = x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]

    deltaE, deltaP, sigma, eta, gammaI, gammaA, gammaQ, gammaH, alphaI, alphaH, zeta, gammaV = p[0], p[1], p[2], p[3], \
                                                                                               p[4], p[5], p[6], p[7], \
                                                                                               p[8], p[9], p[10], p[11]

    Cii, betaP0, betaR, epsilonA, epsilonI = p_foi[0], p_foi[1], p_foi[2], p_foi[3], p_foi[4]

    foi_ii = Cii * ((Cii * (betaP0 * betaR * (P + epsilonA * A)) + epsilonI * betaP0 * betaR * I) / (
            Cii * (S + E + P + R + A + V) + I + 1e-10))

    # v = u[0]
    foi = mob / 1e5 + foi_ii
    rhs = [None] * nx
    vaccrate = 0
    rhs[0] = -(foi + vaccrate) * S + gammaV * V  # S
    rhs[1] = foi * S - deltaE * E  # E
    rhs[2] = deltaE * E - deltaP * P  # P
    rhs[3] = sigma * deltaP * P - (eta + gammaI + alphaI) * I  # I
    rhs[4] = (1 - sigma) * deltaP * P - gammaA * A  # A
    rhs[5] = zeta * eta * I - gammaQ * Q  # Q
    rhs[6] = (1 - zeta) * eta * I - (gammaH + alphaH) * H  # H
    rhs[7] = gammaI * I + gammaA * A + gammaH * H + gammaQ * Q  # R
    rhs[8] = vaccrate * S - gammaV * V  # V

    rhs_ell = [None] * 3
    rhs_ell[0] = gammaH * H  # recovered from the hospital
    rhs_ell[1] = alphaH * H  # total death
    rhs_ell[2] = (1 - zeta) * eta * I  # cumulative hospitalized cases

    return rhs, rhs_ell


def frhs_integrate(y, p, foi, pop_node, p_foi=[0] * 5):
    y, ell = rhs_py(t=0, x=y, u=0, cov=0, p=p, mob=foi, pop_node=pop_node, p_foi=p_foi)
    return np.array(y), ell[1]


def rk4_integrate(y, pvector, mob, pop_node, p_foi, dt):
    k1, k1ell = frhs_integrate(y, foi=mob, p=pvector, pop_node=pop_node, p_foi=p_foi)
    k2, k2ell = frhs_integrate(y + dt / 2 * k1, foi=mob, p=pvector, pop_node=pop_node, p_foi=p_foi)
    k3, k3ell = frhs_integrate(y + dt / 2 * k2, foi=mob, p=pvector, pop_node=pop_node, p_foi=p_foi)
    k4, k4ell = frhs_integrate(y + dt * k3, foi=mob, p=pvector, pop_node=pop_node, p_foi=p_foi)
    x_next = y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    # No need for because we sum it week by week few lines below.
    ell_next = dt / 6 * (k1ell + 2 * k2ell + 2 * k3ell + k4ell)

    return x_next, ell_next


def euler_integrate(y, pvector, mob, pop_node, p_foi, dt):
    a, b = frhs_integrate(y, foi=mob, p=pvector, pop_node=pop_node, p_foi=p_foi)
    x_next = y + dt * a
    ell_next = dt * b
    return x_next, ell_next


def integrate(N, setup, parameters, controls, n_rk4_steps=10, method='rk4', save_to=None):
    M = setup.nnodes

    pvector, pvector_names = parameters.get_pvector()

    dt = (N + 1) / N / n_rk4_steps

    S, E, P, I, A, Q, H, R, V = np.arange(nx)

    y = np.zeros((M, N + 1, nx))
    yell = np.zeros((M, N + 1))
    mob = np.zeros((M, N + 1))

    for cp, name in enumerate(states_names):
        for i in range(M):
            y[i, 0, cp] = parameters.x0[i * nx + cp]

    print(f"===> Integrating for {save_to}""")
    for k in tqdm(range(N)):
        mobK = parameters.mobintime_arr[:, k]
        betaR = parameters.betaratiointime_arr[:, k]
        C = parameters.params_structural['r'] * parameters.mobfrac.flatten() * mobK * parameters.mobmat_pr
        np.fill_diagonal(C, 1 - C.sum(axis=1) + C.diagonal())
        C_foi = np.copy(C.diagonal())
        np.fill_diagonal(C, np.zeros_like(C.diagonal()))
        # plt.spy(C)
        Sk, Ek, Pk, Rk, Ak, Ik, Vk = y[:, k, S], y[:, k, E], y[:, k, P], y[:, k, R], y[:, k, A], y[:, k, I], y[:, k, V]

        foi_sup = []
        foi_inf = []
        for n in range(M):
            foi_sup.append(parameters.params_structural['betaP0'] * betaR[n] * (
                    Pk[n] + parameters.params_structural['epsilonA'] * Ak[n]))
            foi_inf.append(Sk[n] + Ek[n] + Pk[n] + Rk[n] + Ak[n] + Vk[n])

        foi = []
        for m in range(M):
            foi.append((sum(C[n, m] * foi_sup[n] for n in range(M)) + parameters.params_structural['epsilonI'] *
                        parameters.params_structural['betaP0'] * betaR[m] * Ik[m]) /
                       (sum(C[l, m] * foi_inf[l] for l in range(M)) + Ik[m] + 1e-15))

        for i in range(M):
            mob_ik = sum(C[i, m] * foi[m] for m in range(M)) * 1e5

            mob[i, k] = mob_ik

            x_ = y[i, k, :]

            VacPpl = Sk[i] + Ek[i] + Pk[i] + Ak[i] + Rk[i]
            vaccrate = controls[i, k] / (VacPpl + 1e-10)
            x_[S] -= vaccrate * Sk[i]
            x_[V] += vaccrate * Sk[i]

            p_foi = [C_foi[i], parameters.params_structural['betaP0'], betaR[i],
                     parameters.params_structural['epsilonA'], parameters.params_structural['epsilonI']]

            ell = 0.
            for nt in range(n_rk4_steps):
                if method == 'rk4':
                    x_, ell_ = rk4_integrate(x_, pvector, mob_ik, setup.pop_node[i], p_foi, dt)
                elif method == 'euler':
                    x_, ell_ = euler_integrate(x_, pvector, mob_ik, setup.pop_node[i], p_foi, dt)
            ell += ell_

            y[i, k + 1, :] = x_
            yell[i, k + 1] = ell

    results = pd.DataFrame(columns=['date', 'comp', 'place', 'value', 'placeID'])

    for nd in range(M):
        results = pd.concat(
            [results, pd.DataFrame.from_dict(
                {'value': np.append(controls[nd, :], controls[nd, -1]).ravel(),
                 'date': setup.model_days,
                 'place': setup.ind2name[nd],
                 'placeID': int(nd),
                 'comp': 'vacc'})])
        for i, st in enumerate(states_names):
            results = pd.concat(
                [results, pd.DataFrame.from_dict({'value': y[nd, :, i].ravel(),
                                                  'date': setup.model_days,
                                                  'place': setup.ind2name[nd],
                                                  'placeID': int(nd),
                                                  'comp': st})])
    results['placeID'] = results['placeID'].astype(int)

    if save_to is not None:
        results.to_csv(f'{save_to}.csv', index=False)

    return results, y, yell, mob


class COVIDVaccinationOCP:
    def __init__(self, N, n_int_steps, setup, parameters, integ='rk4', show_steps=True):
        timer_start = timer()
        self.N = N
        self.n_int_steps = n_int_steps
        self.setup = setup
        self.parameters = parameters
        self.M = M = setup.nnodes

        _, pvector_names = parameters.get_pvector()

        dt = (N + 1) / N / n_int_steps

        states = cat.struct_symSX(states_names)
        [S, E, P, I, A, Q, H, R, V] = states[...]

        controls = cat.struct_symSX(['v', 'mob'])
        [v, mob] = controls[...]

        covar = cat.struct_symSX(['mobility_t', 'betaratio_t'])
        [mobility_t, betaratio_t] = covar[...]

        params = cat.struct_symSX(pvector_names)
        [deltaE, deltaP, sigma, eta, gammaI, gammaA, gammaQ, gammaH, alphaI, alphaH, zeta, gammaV,
         scale_ell, scale_If, scale_v] = params[...]

        pop_nodeSX = ca.SX.sym('pop_node')

        p_foiSX = cat.struct_symSX(['Cii', 'betaP', 'betaR', 'epsilonI', 'epsilonA'])

        # The rhs is at time zero, the time is also no used in the equation so that explain
        rhs, rhs_ell = rhs_py(0, states.cat, controls.cat, covar.cat, params.cat, mob, pop_nodeSX, p_foiSX.cat)
        rhs = ca.veccat(*rhs)
        rhs_ell = ca.veccat(*rhs_ell)  # mod

        frhs = ca.Function('frhs', [states, controls, covar, params, pop_nodeSX, p_foiSX],
                           [rhs, rhs_ell[1]])  # scale_ell * rhs_ell[1] + scale_v * v * v])# mod ICI juste ell[1]

        # ---- dynamic constraints --------
        if integ == 'rk4':
            k1, k1ell = frhs(states, controls, covar, params, pop_nodeSX, p_foiSX)
            k2, k2ell = frhs(states + dt / 2 * k1, controls, covar, params, pop_nodeSX, p_foiSX)
            k3, k3ell = frhs(states + dt / 2 * k2, controls, covar, params, pop_nodeSX, p_foiSX)
            k4, k4ell = frhs(states + dt * k3, controls, covar, params, pop_nodeSX, p_foiSX)
            x_next = states + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            # No need for because we sum it week by week few lines below.
            ell_next = dt / 6 * (k1ell + 2 * k2ell + 2 * k3ell + k4ell)
            rk4_step = ca.Function('rk4_step', [states, controls, covar, params, pop_nodeSX, p_foiSX],
                                   [x_next, ell_next])

        elif integ == 'euler':
            a, b = frhs(states, controls, covar, params, pop_nodeSX, p_foiSX)
            x_next = states + dt * a
            ell_next = dt * b
            rk4_step = ca.Function('rk4_step', [states, controls, covar, params, pop_nodeSX, p_foiSX],
                                   [x_next, ell_next])

        x_ = ca.veccat(*states[...])
        u_ = ca.veccat(*controls[...])
        VacPpl = states['S'] + states['E'] + states['P'] + states['A'] + states['R']
        vaccrate = controls['v'] / (VacPpl + 1e-10)
        x_[0] -= vaccrate * states['S']
        x_[8] += vaccrate * states['S']
        ell = 0.
        for k in range(n_int_steps):
            x_, ell_ = rk4_step(x_, u_, covar, params, pop_nodeSX, p_foiSX)
            ell += ell_

        rk4_int = ca.Function('rk4_int', [states, ca.veccat(controls, covar, params, pop_nodeSX, p_foiSX)], [x_, ell],
                              ['x0', 'p'], ['xf', 'qf'])

        # BUG TODO Isn't this a double multiplication by the scale parameter since ell is already multiplied ?
        ell = ca.Function('ell', [states, controls, covar, params, pop_nodeSX, p_foiSX],
                          [scale_ell * ell + scale_v * v * v, scale_ell * ell,
                           scale_v * v * v])  # Very dependent on regularization factor

        self.Vars = cat.struct_symMX([
            (
                cat.entry("x", struct=states, repeat=[M, N + 1]),
                cat.entry("u", struct=controls, repeat=[M, N]),
            ),
        ])
        self.Params = cat.struct_symMX([
            (
                cat.entry("cov", struct=covar, repeat=[M, N]),
                cat.entry("p", struct=params),
            ),
        ])

        f = vaccines = cases = reg = cdot_T = 0

        dyn = [None] * N
        spatial = [None] * N
        Sgeq0 = [None] * N
        print(f"===> Building OCP {M} nodes:""")
        for k in tqdm(range(N)):
            mobK = self.Params['cov', :, k, 'mobility_t']  # mobintime.to_numpy().T[:,k]
            betaR = self.Params['cov', :, k, 'betaratio_t']  # betaratiointime.to_numpy().T[:,k]
            C = parameters.params_structural['r'] * parameters.mobfrac.flatten() * mobK * parameters.mobmat_pr
            np.fill_diagonal(C, 1 - C.sum(axis=1) + C.diagonal())
            C_foi = np.copy(C.diagonal())
            np.fill_diagonal(C, np.zeros_like(C.diagonal()))

            # Should this be k+1 ? to have the foi mobility.
            Sk, Ek, Pk, Rk, Ak, Ik, Vk = ca.veccat(*self.Vars['x', :, k, 'S']), ca.veccat(*self.Vars['x', :, k, 'E']), \
                                         ca.veccat(*self.Vars['x', :, k, 'P']), ca.veccat(*self.Vars['x', :, k, 'R']), \
                                         ca.veccat(*self.Vars['x', :, k, 'A']), ca.veccat(*self.Vars['x', :, k, 'I']), \
                                         ca.veccat(*self.Vars['x', :, k, 'V'])

            foi_sup = []
            foi_inf = []
            for n in range(M):
                foi_sup.append(parameters.params_structural['betaP0'] * betaR[n] * (
                        Pk[n] + parameters.params_structural['epsilonA'] * Ak[n]))
                foi_inf.append(Sk[n] + Ek[n] + Pk[n] + Rk[n] + Ak[n] + Vk[n])
            foi = []
            for m in range(M):
                foi.append((sum(C[n, m] * foi_sup[n] for n in range(M)) + parameters.params_structural['epsilonI'] *
                            parameters.params_structural['betaP0'] * betaR[m] * Ik[m]) /
                           (sum(C[l, m] * foi_inf[l] for l in range(M)) + Ik[m] + 1e-10))

            dyn[k] = []
            spatial[k] = []
            Sgeq0[k] = []
            for i in range(M):
                p_foi = ca.veccat(C_foi[i], parameters.params_structural['betaP0'], betaR[i],
                                  parameters.params_structural['epsilonA'],
                                  parameters.params_structural['epsilonI'])

                [X_, ell_ik] = rk4_int(self.Vars['x', i, k],
                                       ca.veccat(self.Vars['u', i, k],
                                                 self.Params['cov', i, k],
                                                 self.Params['p'],
                                                 self.setup.pop_node[i],
                                                 p_foi))

                dyn[k].append(self.Vars['x', i, k + 1] - X_)
                ell_ik_, cases_ik, reg_ik = ell(self.Vars['x', i, k], self.Vars['u', i, k], self.Params['cov', i, k],
                                                self.Params['p'],
                                                self.setup.pop_node[i],
                                                p_foi)
                f += ell_ik_  # MOD: before  ell_ik_
                cases += cases_ik
                reg += reg_ik
                mob_ik = sum(C[i, m] * foi[m] for m in range(M)) * 1e5

                # spatial, vaccines and dyn are put in g(x),
                # with constraints that spatial and dyn are equal to zero
                # thus imposing the dynamics and coupling.
                spatial[k].append(self.Vars['u', i, k, 'mob'] - mob_ik)
                VacPpl = sum(self.Vars['x', i, k, comp] for comp in ['S', 'E', 'P', 'A', 'R'])
                Sgeq0[k].append(self.Vars['x', i, k, 'S'] - self.Vars['u', i, k, 'v'] / (VacPpl + 1e-10))
                # Number of vaccine spent = num of vaccine rate * 7 (number of days)
                vaccines += self.Vars['u', i, k, 'v'] * (N + 1) / N

        f /= (N + 1)  # Average over interval for cost ^ but not terminal cost

        print('-----> Writing constraints, ...', end='')
        self.g = cat.struct_MX([
            cat.entry("dyn", expr=dyn),
            cat.entry("spatial", expr=spatial),
            cat.entry("vaccines", expr=vaccines),
            cat.entry("Sgeq0", expr=Sgeq0)
        ])

        # costTerms = ca.Function('costTerms', [self.Vars, self.Params], [cases, reg])
        print('DONE')

        print('-----> Building NLP function...', end='')
        tsnlf = timer()
        # NLOP arguments:
        # 'x' : variable names to optimize: here states and control
        # 'p' : params that can be change after
        # 'f' : is the objective function of 'x', that we ought to minize
        # 'g' : is a function that goes lbg < g(x,p) < ubg. If you want equality
        #     constraint then ubg = lbg = the number yoiu want it to be.
        nlp = {'x': self.Vars, 'p': self.Params, 'f': f, 'g': self.g}
        self.nlpFun = ca.Function('nlpFun', [self.Vars, self.Params], [f, self.g])
        print(f'DONE in {timer() - tsnlf:.1f} s')

        # print('-----> Building Jacobian function...', end='')
        # self.nlpJac = self.nlpFun.factory('nlpJac', ['i0', 'i1'], ['jac:o1:i0'])
        # print('DONE')

        print('-----> Building Solver...', end='')
        tsbs = timer()
        options = {'ipopt': {}}
        options['ipopt']["linear_solver"] = "ma86"  # "ma57"  "ma86"
        options['ipopt']["print_level"] = 12
        options['ipopt']["max_iter"] = 50  # prevent of for beeing clogged in a good scenario
        if show_steps:
            self.callback = PlotIterates('plot_iterates', self.Vars.size, self.g.size, self.Params.size, [0, 1], N + 1,
                                         N, self.Vars,
                                         setup.ind2name, parameters.mobmat_pr, setup.pos_node, setup.pop_node)
            options["iteration_callback"] = self.callback
        self.solver = ca.nlpsol('solver', "ipopt", nlp, options)
        print(f'DONE in {timer() - tsbs:.1f} s')

        # To store the solution:
        self.sol = None
        self.opt = None
        self.lam_g = None
        self.lam_x = None
        self.Jgnum = None
        self.gnum = None
        self.arg = {}
        self.scenario_name = 'no_update'
        print(f'Total build time {timer() - timer_start:.1f}')

    def update(self, parameters, max_total_vacc, max_vacc_rate, states_initial, control_initial, scenario_name='test'):
        # This initialize
        lbg = self.g(0)
        ubg = self.g(0)

        lbx = self.Vars(-np.inf)
        ubx = self.Vars(np.inf)

        ubg['vaccines'] = max_total_vacc  # 2000 * (T * .6) * M  # 8e6 #*M
        lbg['vaccines'] = -np.inf

        ubg['Sgeq0'] = np.inf

        lbx['u', :, :, 'v'] = 0.
        for k in range(self.N):
            for nd in range(self.M):
                ubx['u', nd, k, 'v'] = max_vacc_rate[nd, k]

        # Set initial conditions as constraints
        for cp, name in enumerate(states_names):
            for i in range(self.M):
                lbx['x', i, 0, name] = ubx['x', i, 0, name] = parameters.x0[i * nx + cp]

        # ----> Starting value of the optimizer
        init = self.Vars(0)
        for i, name in enumerate(states_names):
            for k in range(self.N + 1):
                for nd in range(self.M):
                    init['x', nd, k, name] = states_initial[nd, k, i]

        for k in range(self.N):
            for nd in range(self.M):
                init['u', nd, k] = control_initial[nd, k]

        self.arg = {'lbg': lbg,
                    'ubg': ubg,
                    'lbx': lbx,
                    'ubx': ubx,
                    'x0': init,
                    'p': self.Params()}

        self.arg['p']['p'], _ = self.parameters.get_pvector()

        for nd in range(self.M):
            for k in range(self.N):
                self.arg['p']['cov', nd, k, 'mobility_t'] = self.parameters.mobintime_arr[nd, k]
                self.arg['p']['cov', nd, k, 'betaratio_t'] = self.parameters.betaratiointime_arr[nd, k]

        self.scenario_name = scenario_name

    def solveOCP(self, save=True):
        print(f">>> Solving OCP {self.scenario_name} nodes:""")
        timer_start = timer()
        self.sol = self.solver(**self.arg)
        self.opt = self.Vars(self.sol['x'])
        self.lam_g = self.g(self.sol['lam_g'])
        self.lam_x = self.Vars(self.sol['lam_x'])
        [fnum, gnum] = self.nlpFun(self.opt, self.arg['p'])
        # self.Jgnum = self.nlpJac(self.opt, self.arg['p'])
        self.gnum = self.g(gnum)  # 2times ?
        print(f"""
        Vaccines stockpile: 
            {float(self.arg['ubg']['vaccines']):010f} total.
            {float(self.g(self.gnum)['vaccines']):010f} spent.
            {float((self.arg['ubg']['vaccines'] - self.g(self.gnum)['vaccines'])):010f} left.""")

        print(f'Total Solving done in {timer() - timer_start:.1f}s')
        if save:
            self.saveOCP()

    def saveOCP(self):
        results = pd.DataFrame(columns=['date', 'comp', 'place', 'value'])

        for nd in range(self.M):
            results = pd.concat(
                [results, pd.DataFrame.from_dict(
                    {'value': np.array(
                        ca.veccat(ca.veccat(*self.opt['u', nd, :, 'v']), self.opt['u', nd, -1, 'v'])).ravel(),
                     'date': self.setup.model_days,
                     'place': self.setup.ind2name[nd],
                     'placeID': int(nd),
                     'comp': 'vacc'})])
            for i, st in enumerate(states_names):
                results = pd.concat(
                    [results, pd.DataFrame.from_dict({'value': np.array(ca.veccat(*self.opt['x', nd, :, st])).ravel(),
                                                      'date': self.setup.model_days,
                                                      'place': self.setup.ind2name[nd],
                                                      'placeID': int(nd),
                                                      'comp': st})])
        results['placeID'] = results['placeID'].astype(int)
        results.to_csv(f'{self.scenario_name}.csv', index=False)
