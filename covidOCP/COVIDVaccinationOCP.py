import casadi as ca
import casadi.tools as cat
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if "Agg" not in mpl.get_backend():
    mpl.interactive(True)

plt.ion()

nx = 9
states_names = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']


def rhs_py(t, x, u, cov, p, mob, pop_node):
    S, E, P, I, A, Q, H, R, V = x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]

    deltaE, deltaP, sigma, eta, gammaI, gammaA, gammaQ, gammaH, alphaI, alphaH, zeta, gammaV = p[0], p[1], p[2], p[3], \
                                                                                               p[4], p[5], p[6], p[7], \
                                                                                               p[8], p[9], p[10], p[11]

    # v = u[0]
    foi = mob
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


def mobility_graph(mobility, ind2name, pos_node, pop_node, opt, N):
    pass


class COVIDVaccinationOCP:
    def __init__(self, N, n_int_steps, setup, parameters, integ='euler'):
        self.N = N
        self.n_int_steps = n_int_steps
        self.setup = setup
        self.parameters = parameters
        self.M = M = setup.nnodes

        _, pvector_names = parameters.get_pvector()

        dt = (N+1) / N / n_int_steps

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

        # The rhs is at time zero, the time is also no used in the equation so that explain
        rhs, rhs_ell = rhs_py(0, states.cat, controls.cat, covar.cat, params.cat, mob, pop_nodeSX)
        rhs = ca.veccat(*rhs)
        rhs_ell = ca.veccat(*rhs_ell)  # mod

        frhs = ca.Function('frhs', [states, controls, covar, params, pop_nodeSX],
                           [rhs, rhs_ell[1]])  # scale_ell * rhs_ell[1] + scale_v * v * v])# mod ICI juste ell[1]

        # ---- dynamic constraints --------
        if integ == 'rk4':
            k1, k1ell = frhs(states, controls, covar, params, pop_nodeSX)
            k2, k2ell = frhs(states + dt / 2 * k1, controls, covar, params, pop_nodeSX)
            k3, k3ell = frhs(states + dt / 2 * k2, controls, covar, params, pop_nodeSX)
            k4, k4ell = frhs(states + dt * k3, controls, covar, params, pop_nodeSX)
            x_next = states + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            # No need for because we sum it week by week few lines below.
            ell_next = dt / 6 * (k1ell + 2 * k2ell + 2 * k3ell + k4ell)
            rk4_step = ca.Function('rk4_step', [states, controls, covar, params, pop_nodeSX], [x_next, ell_next])

        elif integ == 'euler':
            a, b = frhs(states, controls, covar, params, pop_nodeSX)
            x_next = states + dt * a
            ell_next = dt * b
            rk4_step = ca.Function('rk4_step', [states, controls, covar, params, pop_nodeSX], [x_next, ell_next])

        x_ = ca.veccat(*states[...])
        u_ = ca.veccat(*controls[...])
        VacPpl = states['S'] + states['E'] + states['P'] + states['A'] + states['R']
        vaccrate = controls['v'] / (VacPpl + 1e-10)
        x_[0] -= vaccrate * states['S']
        x_[8] += vaccrate * states['S']
        ell = 0.
        for k in range(n_int_steps):
            x_, ell_ = rk4_step(x_, u_, covar, params, pop_nodeSX)
            ell += ell_

        rk4_int = ca.Function('rk4_int', [states, ca.veccat(controls, covar, params, pop_nodeSX)], [x_, ell],
                              ['x0', 'p'], ['xf', 'qf'])

        # BUG TODO Isn't this a double multiplication by the scale parameter since ell is already multiplied ?
        ell = ca.Function('ell', [states, controls, covar, params, pop_nodeSX],
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

        for k in range(N):
            mobK = self.Params['cov', :, k, 'mobility_t']  # mobintime.to_numpy().T[:,k]
            betaR = self.Params['cov', :, k, 'betaratio_t']  # betaratiointime.to_numpy().T[:,k]
            C = parameters.params_structural['r'] * parameters.mobfrac.flatten() * mobK * parameters.mobmat_pr
            np.fill_diagonal(C, 1 - C.sum(axis=1) + C.diagonal())

            Sk, Ek, Pk, Rk, Ak, Ik = ca.veccat(*self.Vars['x', :, k, 'S']), ca.veccat(*self.Vars['x', :, k, 'E']), \
                                     ca.veccat(*self.Vars['x', :, k, 'P']), ca.veccat(*self.Vars['x', :, k, 'R']), \
                                     ca.veccat(*self.Vars['x', :, k, 'A']), ca.veccat(*self.Vars['x', :, k, 'I'])
            foi_sup = []
            foi_inf = []
            for n in range(M):
                foi_sup.append(parameters.params_structural['betaP0'] * betaR[n] * (
                        Pk[n] + parameters.params_structural['epsilonA'] * Ak[n]))
                foi_inf.append(Sk[n] + Ek[n] + Pk[n] + Rk[n] + Ak[n])
            foi = []
            for m in range(M):
                foi.append((sum(C[n, m] * foi_sup[n] for n in range(M)) + parameters.params_structural['epsilonI'] *
                            parameters.params_structural['betaP0'] * betaR[m] * Ik[m]) /
                           (sum(C[l, m] * foi_inf[l] for l in range(M)) + Ik[m]))

            print(f'{k}:', end='')
            dyn[k] = []
            spatial[k] = []
            Sgeq0[k] = []
            for i in range(M):
                [X_, ell_ik] = rk4_int(self.Vars['x', i, k],
                                       ca.veccat(self.Vars['u', i, k],
                                                 self.Params['cov', i, k],
                                                 self.Params['p'],
                                                 self.setup.pop_node[i]))

                dyn[k].append(self.Vars['x', i, k + 1] - X_)
                ell_ik_, cases_ik, reg_ik = ell(self.Vars['x', i, k], self.Vars['u', i, k], self.Params['cov', i, k],
                                                self.Params['p'],
                                                self.setup.pop_node[i])
                f += ell_ik_  # MOD: before  ell_ik_
                cases += cases_ik
                reg += reg_ik
                mob_ik = sum(C[i, m] * foi[m] for m in range(M))

                # spatial, vaccines and dyn are put in g(x),
                # with constraints that spatial and dyn are equal to zero
                # thus imposing the dynamics and coupling.
                spatial[k].append(self.Vars['u', i, k, 'mob'] - mob_ik)
                VacPpl = sum(self.Vars['x', i, k, comp] for comp in ['S', 'E', 'P', 'A', 'R'])
                Sgeq0[k].append(self.Vars['x', i, k, 'S'] - self.Vars['u', i, k, 'v'] / (VacPpl + 1e-10))
                # Number of vaccine spent = num of vaccine rate * 7 (number of days)
                vaccines += self.Vars['u', i, k, 'v'] * (N+1) / N

        f /= (N+1)  # Average over interval for cost ^ but not terminal cost

        print('-> Writing constraints, ...', end='')
        self.g = cat.struct_MX([
            cat.entry("dyn", expr=dyn),
            cat.entry("spatial", expr=spatial),
            cat.entry("vaccines", expr=vaccines),
            cat.entry("Sgeq0", expr=Sgeq0)
        ])

        # costTerms = ca.Function('costTerms', [self.Vars, self.Params], [cases, reg])
        print('DONE')

        print('-> Building NLP function...', end='')
        # NLOP arguments:
        # 'x' : variable names to optimize: here states and control
        # 'p' : params that can be change after
        # 'f' : is the objective function of 'x', that we ought to minize
        # 'g' : is a function that goes lbg < g(x,p) < ubg. If you want equality
        #     constraint then ubg = lbg = the number yoiu want it to be.
        nlp = {'x': self.Vars, 'p': self.Params, 'f': f, 'g': self.g}
        self.nlpFun = ca.Function('nlpFun', [self.Vars, self.Params], [f, self.g])
        print('DONE')

        print('-> Building Jacobian function...', end='')
        self.nlpJac = self.nlpFun.factory('nlpJac', ['i0', 'i1'], ['jac:o1:i0'])
        print('DONE')

        print('-> Building Solver...', end='')
        options = {'ipopt': {}}
        options['ipopt']["linear_solver"] = "ma86"  # "ma57"  "ma86"
        self.solver = ca.nlpsol('solver', "ipopt", nlp, options)
        print('DONE')

        # To store the solution:
        self.sol = None
        self.opt = None
        self.lam_g = None
        self.lam_x = None
        self.Jgnum = None
        self.gnum = None
        self.arg = {}
        self.scenario_name = 'no_update'

    def update(self, parameters, max_total_vacc, max_vacc_rate, states_initial, control_initial, scenario_name = 'test'):
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

    def solveOCP(self, save = True):
        self.sol = self.solver(**self.arg)
        self.opt = self.Vars(self.sol['x'])
        self.lam_g = self.g(self.sol['lam_g'])
        self.lam_x = self.Vars(self.sol['lam_x'])
        [fnum, gnum] = self.nlpFun(self.opt, self.arg['p'])
        self.Jgnum = self.nlpJac(self.opt, self.arg['p'])
        self.gnum = self.g(gnum)  # 2times ?
        print(f"""
        Vaccines stockpile: 
            {float(self.arg['ubg']['vaccines']):010f} total.
            {float(self.g(gnum)['vaccines']):010f} spent.
            {float((self.arg['ubg']['vaccines'] - self.g(self.gnum)['vaccines'])):010f} left.""")

        if save:
            self.saveOCP()

    def saveOCP(self):
        results = pd.DataFrame(columns=['date', 'comp', 'place', 'value'])

        for nd in range(self.M):
            results = pd.concat(
                [results, pd.DataFrame.from_dict(
                    {'value': np.array(ca.veccat(ca.veccat(*self.opt['u', nd, :, 'v']), self.opt['u', nd, -1, 'v'])).ravel(),
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
        results.to_csv(f'df_{self.scenario_name}.csv', index=False)

