import casadi as ca
import casadi.tools as cat
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm
import pickle
from timeit import default_timer as timer

# if "Agg" not in mpl.get_backend():
#    mpl.interactive(True)
# plt.ion()

nx = 13
states_names = ['S1', 'S2', 'S3', 'S4', 'S5', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']
mob_scaling = 1e7


def rhs_py(t, x, u, cov, p, mob, pop_node, p_foi):
    S1, S2, S3, S4, S5, E, P, I, A, Q, H, R, V = x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12]

    deltaE, deltaP, sigma, eta, gammaI, gammaA, gammaQ, gammaH, alphaI, alphaH, zeta, gammaV = p[0], p[1], p[2], p[3], \
                                                                                               p[4], p[5], p[6], p[7], \
                                                                                               p[8], p[9], p[10], p[11]

    pop_node1, pop_node2, pop_node3, pop_node4, pop_node5 = pop_node[0], pop_node[1], pop_node[2], pop_node[3], pop_node[4]
    pop_node = pop_node1 + pop_node2 + pop_node3 + pop_node4 + pop_node5

    Cii, betaP0, betaR, epsilonA, epsilonI = p_foi[0], p_foi[1], p_foi[2], p_foi[3], p_foi[4]

    S = S1+S2+S3+S4+S5

    foi_ii = Cii * ((Cii * (betaP0 * betaR * (P + epsilonA * A)) + epsilonI * betaP0 * betaR * I) / (
            Cii * (S + E + P + R + A + V) + I + 1e-10))

    # v = u[0]
    foi = mob / mob_scaling + foi_ii
    rhs = [None] * nx
    vaccrate = 0
    #rhs[0] = -(foi + vaccrate) * S #+ gammaV * V  # S
    rhs[5] = foi * S - deltaE * E  # E
    rhs[6] = deltaE * E - deltaP * P  # P
    rhs[7] = sigma * deltaP * P - (eta + gammaI) * I  # I alphaI
    rhs[8] = (1 - sigma) * deltaP * P - gammaA * A  # A
    rhs[9] = zeta * eta * I - gammaQ * Q  # Q
    rhs[10] = (1 - zeta) * eta * I - (gammaH + alphaH) * H  # H
    rhs[11] = gammaI * I + gammaA * A + gammaH * H + gammaQ * Q  # R
    rhs[12] = vaccrate * S #- gammaV * V  # V

    p1 = S1/S * 12.13/100 * pop_node/pop_node1
    p2 = S2/S * 24.23/100 * pop_node/pop_node2
    p3 = S3/S * 33.92/100 * pop_node/pop_node3
    p4 = S4/S * 19.58/100 * pop_node/pop_node4
    p5 = S5/S * 10.14/100 * pop_node/pop_node5

    p = p1+p2+p3+p4+p5
    p1 /= p
    p2 /= p
    p3 /= p
    p4 /= p
    p5 /= p

    E1 = foi * S * p1
    E2 = foi * S * p2
    E3 = foi * S * p3
    E4 = foi * S * p4
    E5 = foi * S * p5
    rhs[0] = -E1
    rhs[1] = -E2
    rhs[2] = -E3
    rhs[3] = -E4
    rhs[4] = -E5


    rhs_ell = [None] * 2
    rhs_ell[0] = foi * S  # recovered from the hospital
    #rhs_ell[1] = foi * S #alphaH * H  # total death   # ONLY THIS ONE IS USED
    rhs_ell[1] = 0.01/100*E1+0.04/100*E2+0.43/100*E3+6.06/100*E4+20.83/100*E5

    return rhs, rhs_ell


def frhs_integrate(y, p, foi, pop_node, p_foi=[0] * 5):
    y, ell = rhs_py(t=0, x=y, u=0, cov=0, p=p, mob=foi, pop_node=pop_node, p_foi=p_foi)
    return np.array(y), np.array(ell)


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
    print(f"===> Integrating for {save_to}")
    M = setup.nnodes

    pvector, pvector_names = parameters.get_pvector()

    dt = (N + 1) / N / n_rk4_steps

    S1, S2, S3, S4, S5, E, P, I, A, Q, H, R, V = np.arange(nx)

    y = np.zeros((M, N + 1, nx))
    yell_death = np.zeros((M, N + 1))
    yell_infection = np.zeros((M, N + 1))
    mob = np.zeros((M, N + 1))

    for cp, name in enumerate(states_names):
        for i in range(M):
            if 'S' in name:
                y[i, 0, cp] = parameters.x0_Sagpost[i, cp]
            else:
                y[i, 0, cp] = parameters.x0[i, cp - 4]

    for k in tqdm(range(N)):
        mobK = setup.mobintime_arr[:, k]
        betaR = parameters.betaratiointime_arr[:, k]
        C = parameters.params_structural['r'] * parameters.mobfrac.flatten() * mobK * parameters.mobmat_pr
        np.fill_diagonal(C, 1 - C.sum(axis=1) + C.diagonal())
        C_foi = np.copy(C.diagonal())
        np.fill_diagonal(C, np.zeros_like(C.diagonal()))
        # plt.spy(C)
        Sk, Ek, Pk, Rk, Ak, Ik, Vk = y[:, k, S1]+y[:, k, S2]+y[:, k, S3]+y[:, k, S4]+y[:, k, S5], y[:, k, E], y[:, k, P], y[:, k, R], y[:, k, A], y[:, k, I], y[:, k, V]
        S1k,S2k,S3k,S4k,S5k = y[:, k, S1], y[:, k, S2], y[:, k, S3], y[:, k, S4], y[:, k, S5]

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
            mob_ik = sum(C[i, m] * foi[m] for m in range(M)) * mob_scaling

            mob[i, k] = mob_ik

            x_ = np.copy(y[i, k, :])

            VacPpl = Sk[i] + Ek[i] + Pk[i] + Ak[i] + Rk[i]
            vaccrate = controls[i, k] / (VacPpl + 1e-10)
            x_[S1] -= vaccrate * S1k[i]
            x_[S2] -= vaccrate * S2k[i]
            x_[S3] -= vaccrate * S3k[i]
            x_[S4] -= vaccrate * S4k[i]
            x_[S5] -= vaccrate * S5k[i]
            x_[E] -= vaccrate * Ek[i]
            x_[P] -= vaccrate * Pk[i]
            x_[A] -= vaccrate * Ak[i]
            x_[R] -= vaccrate * Rk[i]
            x_[V] += controls[i, k]

            p_foi = [C_foi[i], parameters.params_structural['betaP0'], betaR[i],
                     parameters.params_structural['epsilonA'], parameters.params_structural['epsilonI']]

            ell = 0.
            for nt in range(n_rk4_steps):
                if method == 'rk4':
                    x_, ell_ = rk4_integrate(x_, pvector, mob_ik, setup.pop_node_agpost[i], p_foi, dt)
                elif method == 'euler':
                    x_, ell_ = euler_integrate(x_, pvector, mob_ik, setup.pop_node_agpost[i], p_foi, dt)
                ell += ell_

            y[i, k + 1, :] = x_
            yell_infection[i, k + 1] = ell[0]
            yell_death[i, k + 1] = ell[1]

    results = pd.DataFrame(columns=['date', 'comp', 'place', 'value', 'placeID'])

    for nd in range(M):
        results = pd.concat(
            [results,
             pd.DataFrame.from_dict(
                {'value': np.append(controls[nd, :], controls[nd, -1]).ravel(),
                 'date': setup.model_days,
                 'place': setup.ind2name[nd],
                 'placeID': int(nd),
                 'comp': 'vacc'}),
             pd.DataFrame.from_dict(
                 {'value': yell_infection[nd],
                  'date': setup.model_days,
                  'place': setup.ind2name[nd],
                  'placeID': int(nd),
                  'comp': 'yell_infection'}),
             pd.DataFrame.from_dict(
                 {'value': yell_death[nd],
                  'date': setup.model_days,
                  'place': setup.ind2name[nd],
                  'placeID': int(nd),
                  'comp': 'yell_death'})
             ])
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

    return results, y, yell_infection, yell_death, mob


class COVIDVaccinationOCPagpost:
    def __init__(self, N, n_int_steps, setup, parameters, integ='rk4', show_steps=True, objective='death'):
        timer_start = timer()
        self.N = N
        self.n_int_steps = n_int_steps
        self.setup = setup
        self.parameters = parameters
        self.M = M = setup.nnodes

        _, pvector_names = parameters.get_pvector()

        dt = (N + 1) / N / n_int_steps

        if objective == 'death':
            ellID = 1
        elif objective == 'infection':
            ellID = 0

        states = cat.struct_symSX(states_names)
        [S1, S2, S3, S4, S5, E, P, I, A, Q, H, R, V] = states[...]

        controls = cat.struct_symSX(['v', 'mob'])
        [v, mob] = controls[...]

        covar = cat.struct_symSX(['mobility_t', 'betaratio_t'])
        [mobility_t, betaratio_t] = covar[...]

        params = cat.struct_symSX(pvector_names)
        [deltaE, deltaP, sigma, eta, gammaI, gammaA, gammaQ, gammaH, alphaI, alphaH, zeta, gammaV,
         scale_ell, scale_If, scale_v] = params[...]

        pop_nodeSX = cat.struct_symSX([f'popnode_{i}' for i in np.arange(1, 5+1)])

        p_foiSX = cat.struct_symSX(['Cii', 'betaP', 'betaR', 'epsilonI', 'epsilonA'])

        # The rhs is at time zero, the time is also no used in the equation so that explain
        rhs, rhs_ell = rhs_py(0, states.cat, controls.cat, covar.cat, params.cat, mob, pop_nodeSX.cat, p_foiSX.cat)
        rhs = ca.veccat(*rhs)
        rhs_ell = ca.veccat(*rhs_ell)  # mod

        frhs = ca.Function('frhs', [states, controls, covar, params, pop_nodeSX, p_foiSX],
                           [rhs, rhs_ell[ellID]])  # scale_ell * rhs_ell[1] + scale_v * v * v])# mod ICI juste ell[1]

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
        VacPpl = states['S1'] + states['S2'] + states['S3'] + states['S4'] + states['S5'] + states['E'] + states['P'] + states['A'] + states['R']
        vaccrate = controls['v'] / VacPpl #  + 1e-10)
        x_[0] -= vaccrate * states['S1']
        x_[1] -= vaccrate * states['S2']
        x_[2] -= vaccrate * states['S3']
        x_[3] -= vaccrate * states['S4']
        x_[4] -= vaccrate * states['S5']
        x_[5] -= vaccrate * states['E']
        x_[6] -= vaccrate * states['P']
        x_[7] -= vaccrate * states['A']
        x_[8] -= vaccrate * states['R']
        x_[9] += controls['v']

        ell = 0.
        for k in range(n_int_steps):
            x_, ell_ = rk4_step(x_, u_, covar, params, pop_nodeSX, p_foiSX)
            ell += ell_

        rk4_int = ca.Function('rk4_int', [states, ca.veccat(controls, covar, params, pop_nodeSX, p_foiSX)], [x_, ell],
                              ['x0', 'p'], ['xf', 'qf'])

        # cat.dotdraw(x_, figsize=(10, 10))

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

        f = cases = reg = cdot_T = 0

        dyn = [None] * N
        spatial = [None] * N
        Sgeq0 = [None] * N
        vaccines = [None] * N
        print(f"===> Building OCP {M} nodes:""")
        for k in tqdm(range(N)):
            mobK = self.Params['cov', :, k, 'mobility_t']  # mobintime.to_numpy().T[:,k]
            betaR = self.Params['cov', :, k, 'betaratio_t']  # betaratiointime.to_numpy().T[:,k]
            C = parameters.params_structural['r'] * parameters.mobfrac.flatten() * mobK * parameters.mobmat_pr
            np.fill_diagonal(C, 1 - C.sum(axis=1) + C.diagonal())
            C_foi = np.copy(C.diagonal())
            np.fill_diagonal(C, np.zeros_like(C.diagonal()))

            # Should this be k+1 ? to have the foi mobility.
            Sk, Ek, Pk, Rk, Ak, Ik, Vk = ca.veccat(*self.Vars['x', :, k, 'S1'])+ca.veccat(*self.Vars['x', :, k, 'S2'])+ca.veccat(*self.Vars['x', :, k, 'S3'])+ca.veccat(*self.Vars['x', :, k, 'S4']) + ca.veccat(*self.Vars['x', :, k, 'S5']),  \
                                         ca.veccat(*self.Vars['x', :, k, 'E']), \
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
                                                 self.setup.pop_node_agpost[i],
                                                 p_foi))

                dyn[k].append(self.Vars['x', i, k + 1] - X_)
                ell_ik_, cases_ik, reg_ik = ell(self.Vars['x', i, k], self.Vars['u', i, k], self.Params['cov', i, k],
                                                self.Params['p'],
                                                self.setup.pop_node_agpost[i],
                                                p_foi)
                f += ell_ik_  # MOD: before  ell_ik_
                cases += cases_ik
                reg += reg_ik
                mob_ik = sum(C[i, m] * foi[m] for m in range(M)) * mob_scaling

                # spatial, vaccines and dyn are put in g(x),
                # with constraints that spatial and dyn are equal to zero
                # thus imposing the dynamics and coupling.
                spatial[k].append(self.Vars['u', i, k, 'mob'] - mob_ik)
                VacPpl = sum(self.Vars['x', i, k, comp] for comp in ['S1', 'S2' 'S3', 'S4', 'S5', 'E', 'P', 'A', 'R'])
                # Sgeq0[k].append(self.Vars['x', i, k, 'S'] - self.Vars['u', i, k, 'v'] / (VacPpl + 1e-10))
                Sgeq0[k].append(VacPpl - self.Vars['u', i, k, 'v'])
                # Number of vaccine spent = num of vaccine rate * 7 (number of days)
                #vaccines[k] = vaccines[k] + self.Vars['u', i, k, 'v'] #* (N + 1) / N

            vaccines[k] = sum([sum(self.Vars['u', :, j, 'v']) for j in range(k+1)])

        f /= (N + 1)  # Average over interval for cost ^ but not terminal cost

        # cat.dotdraw(mob_ik, figsize=(10, 10))
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
        # options['ipopt']["print_level"] = 12
        options['ipopt']["max_iter"] = 10000
        options['ipopt']["print_info_string"] = "yes"
        #options['ipopt']["print_level"] = 12
        #options['ipopt']["ma86_print_level"] = 12  # Never uncomment, makes it write Gb of outputs
        #options['ipopt']["derivative_test"] = "second-order"
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

    def update(self, parameters, stockpile_national_constraint, maxvaccrate_regional, states_initial, control_initial, mob_initial,
               scenario_name='test'):
        # This initialize
        lbg = self.g(0)
        ubg = self.g(0)

        lbx = self.Vars(-np.inf)
        ubx = self.Vars(np.inf)

        if isinstance(stockpile_national_constraint, (int, float, complex)):
            ubg['vaccines'] = stockpile_national_constraint  # 2000 * (T * .6) * M  # 8e6 #*M
            lbg['vaccines'] = -np.inf
        else:
            for k in range(self.N):
                ubg['vaccines', k] = stockpile_national_constraint[k]
                lbg['vaccines', k] = -np.inf

        ubg['Sgeq0'] = np.inf

        lbx['u', :, :, 'v'] = 0.
        for k in range(self.N):
            for nd in range(self.M):
                ubx['u', nd, k, 'v'] = maxvaccrate_regional[nd, k]

        # Set initial conditions as constraints
        for cp, name in enumerate(states_names):
            for i in range(self.M):
                if 'S' in name:
                    lbx['x', i, 0, name] = ubx['x', i, 0, name] = parameters.x0_Sagpost[i, cp]
                else:
                    lbx['x', i, 0, name] = ubx['x', i, 0, name] = parameters.x0[i, cp-4]

        # ----> Starting value of the optimizer
        init = self.Vars(0)
        for i, name in enumerate(states_names):
            for k in range(self.N + 1):
                for nd in range(self.M):
                    init['x', nd, k, name] = states_initial[nd, k, i]

        for k in range(self.N):
            for nd in range(self.M):
                init['u', nd, k, 'v'] = control_initial[nd, k]
                init['u', nd, k, 'mob'] = mob_initial[nd, k]

        self.arg = {'lbg': lbg,
                    'ubg': ubg,
                    'lbx': lbx,
                    'ubx': ubx,
                    'x0': init,
                    'p': self.Params()}

        self.arg['p']['p'], _ = self.parameters.get_pvector()

        for nd in range(self.M):
            for k in range(self.N):
                self.arg['p']['cov', nd, k, 'mobility_t'] = self.setup.mobintime_arr[nd, k]
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
        try:
            print(f"""
            Vaccines stockpile: 
                {float(self.arg['ubg']['vaccines']):.1f} total.
                {float(self.g(self.gnum)['vaccines']):.1f} spent.
                {float((self.arg['ubg']['vaccines'] - self.g(self.gnum)['vaccines'])):.1f} left.""")
        except TypeError:
            pass

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

        with open(f'{self.scenario_name}_lamg.pkl', 'wb') as out:
            pickle.dump(self.lam_g, out, pickle.HIGHEST_PROTOCOL)
        with open(f'{self.scenario_name}_lamx.pkl', 'wb') as out:
            pickle.dump(self.lam_x, out, pickle.HIGHEST_PROTOCOL)
