import casadi as ca
import casadi.tools as cat
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm
import pickle
import copy
from timeit import default_timer as timer

# if "Agg" not in mpl.get_backend():
#    mpl.interactive(True)
# plt.ion()

nx = 9
states_names = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']
mob_scaling = 1e7
nc = 3  # Number of age classes
ages_names = ['Y', 'M', 'O']

ag_trans_mult = {'Y': 1,  'M': 1, 'O': 1}
ag_death_mult = {'Y': 1, 'M': 1, 'O': 1}


def rhs_py(t, x, u, cov, p, mob, pop_node, p_foi):
    S, E, P, I, A, Q, H, R, V = np.arange(nx)

    S, E, P, I, A, Q, H, R, V = [x[i * nx + S] for i in range(nc)], [x[i * nx + E] for i in range(nc)], \
                                [x[i * nx + P] for i in range(nc)], [x[i * nx + I] for i in range(nc)], \
                                [x[i * nx + A] for i in range(nc)], [x[i * nx + Q] for i in range(nc)], \
                                [x[i * nx + H] for i in range(nc)], [x[i * nx + R] for i in range(nc)], \
                                [x[i * nx + V] for i in range(nc)]

    deltaE, deltaP, sigma, eta, gammaI, gammaA, gammaQ, gammaH, alphaI, alphaH, zeta, gammaV = p[0], p[1], p[2], p[3], \
                                                                                               p[4], p[5], p[6], p[7], \
                                                                                               p[8], p[9], p[10], p[11]

    Cii, betaP0, betaR, epsilonA, epsilonI = p_foi[0], p_foi[1], p_foi[2], p_foi[3], p_foi[4]

    Asum, Psum, Isum = 0, 0, 0
    for ag_id, ag in enumerate(ages_names):
        Asum += ag_trans_mult[ag] * A[ag_id]
        Psum += ag_trans_mult[ag] * P[ag_id]
        Isum += ag_trans_mult[ag] * I[ag_id]

    foi_ii = Cii * ((Cii * (betaP0 * betaR * (Psum + epsilonA * Asum)) + epsilonI * betaP0 * betaR * Isum) / (
            Cii * (sum(S + E + P + R + A + V)) + sum(I) + 1e-10))

    foi = mob / mob_scaling + foi_ii
    rhs = [None] * nx * nc
    rhs_ell = [None] * 2

    # v = u[0]
    for ag_id, ag in enumerate(ages_names):
        vaccrate = 0
        rhs[0 + nx * ag_id] = -(foi + vaccrate) * S[ag_id] #+ gammaV * V[ag_id]  # S
        rhs[1 + nx * ag_id] = foi * S[ag_id] - deltaE * E[ag_id]  # E
        rhs[2 + nx * ag_id] = deltaE * E[ag_id] - deltaP * P[ag_id]  # P
        rhs[3 + nx * ag_id] = sigma * deltaP * P[ag_id] - (eta + gammaI) * I[ag_id]  # I
        rhs[4 + nx * ag_id] = (1 - sigma) * deltaP * P[ag_id] - gammaA * A[ag_id]  # A
        rhs[5 + nx * ag_id] = zeta * eta * I[ag_id] - gammaQ * Q[ag_id]  # Q
        rhs[6 + nx * ag_id] = (1 - zeta) * eta * I[ag_id] - (gammaH + alphaH * ag_death_mult[ag]) * H[ag_id]  # H
        rhs[7 + nx * ag_id] = gammaI * I[ag_id] + gammaA * A[ag_id] + gammaH * H[ag_id] + gammaQ * Q[ag_id]  # R
        rhs[8 + nx * ag_id] = vaccrate * S[ag_id] #- gammaV * V[ag_id]  # V

    rhs_ell[0] = sum(alphaH * ag_death_mult[ag] * H[agid] for agid, ag in enumerate(ages_names))  # total death
    rhs_ell[1] = sum(foi * S[agid] for agid, ag in enumerate(ages_names))  # infection

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
    print(f"===> Integrating for {save_to}""")
    M = setup.nnodes

    pvector, pvector_names = parameters.get_pvector()

    dt = (N + 1) / N / n_rk4_steps

    S, E, P, I, A, Q, H, R, V = np.arange(nx)

    y = np.zeros((M, N + 1, nc, nx))
    yell_death = np.zeros((M, N + 1))
    yell_infection = np.zeros((M, N + 1))
    mob = np.zeros((M, N + 1))

    for cp, name in enumerate(states_names):
        for i in range(M):
            y[i, 0, :, cp] = parameters.x0_ag[i, :, cp]

    for k in tqdm(range(N)):
        mobK = setup.mobintime_arr[:, k]
        betaR = parameters.betaratiointime_arr[:, k]
        C = parameters.params_structural['r'] * parameters.mobfrac.flatten() * mobK * parameters.mobmat_pr
        np.fill_diagonal(C, 1 - C.sum(axis=1) + C.diagonal())
        C_foi = np.copy(C.diagonal())
        np.fill_diagonal(C, np.zeros_like(C.diagonal()))
        # plt.spy(C)
        Sk, Ek, Pk, Rk, Ak, Ik, Vk = [y[:, k, i, S] for i in range(nc)], [y[:, k, i, E] for i in range(nc)], \
                                     [y[:, k, i, P] for i in range(nc)], [y[:, k, i, R] for i in range(nc)], \
                                     [y[:, k, i, A] for i in range(nc)], [y[:, k, i, I] for i in range(nc)], \
                                     [y[:, k, i, V] for i in range(nc)]

        foi_sup = []
        foi_inf = []
        for n in range(M):
            Aksum, Pksum = 0, 0
            for ag_id, ag in enumerate(['Y', 'M']):
                Aksum += ag_trans_mult[ag] * Ak[ag_id][n]
                Pksum += ag_trans_mult[ag] * Pk[ag_id][n]
            foi_sup.append(parameters.params_structural['betaP0'] * betaR[n] * (
                    Pksum + parameters.params_structural['epsilonA'] * Aksum))
            nfoi_inf = 0
            for ag_id, ag in enumerate(['Y', 'M']):
                nfoi_inf += Sk[ag_id][n] + Ek[ag_id][n] + Pk[ag_id][n] + Rk[ag_id][n] + Ak[ag_id][n] + Vk[ag_id][n]
            foi_inf.append(nfoi_inf)

        foi = []
        for m in range(M):
            Iksum, IksumVanilla = 0, 0
            for ag_id, ag in enumerate(['Y', 'M']):
                Iksum += ag_trans_mult[ag] * Ik[ag_id][n]
                IksumVanilla += Ik[ag_id][n]
            foi.append((sum(C[n, m] * foi_sup[n] for n in range(M)) + parameters.params_structural['epsilonI'] *
                        parameters.params_structural['betaP0'] * betaR[m] * Iksum) /
                       (sum(C[l, m] * foi_inf[l] for l in range(M)) + IksumVanilla + 1e-10))

        for i in range(M):
            mob_ik = sum(C[i, m] * foi[m] for m in range(M)) * mob_scaling

            mob[i, k] = mob_ik

            x_ = np.copy(y[i, k])
            for ag_id, ag in enumerate(ages_names):
                VacPpl = Sk[ag_id][i] + Ek[ag_id][i] + Pk[ag_id][i] + Ak[ag_id][i] + Rk[ag_id][i]
                vaccrate = controls[i, k, ag_id] / (VacPpl + 1e-10)
                x_[ag_id][S] -= vaccrate * Sk[ag_id][i]
                x_[ag_id][E] -= vaccrate * Ek[ag_id][i]
                x_[ag_id][P] -= vaccrate * Pk[ag_id][i]
                x_[ag_id][A] -= vaccrate * Ak[ag_id][i]
                x_[ag_id][R] -= vaccrate * Rk[ag_id][i]
                x_[ag_id][V] += controls[i, k, ag_id] #* Sk[ag_id][i]

            p_foi = [C_foi[i], parameters.params_structural['betaP0'], betaR[i],
                     parameters.params_structural['epsilonA'], parameters.params_structural['epsilonI']]

            x_ = x_.reshape(nx*nc)
            ell = 0.
            for nt in range(n_rk4_steps):
                if method == 'rk4':
                    x_, ell_ = rk4_integrate(x_, pvector, mob_ik, setup.pop_node[i], p_foi, dt)
                elif method == 'euler':
                    x_, ell_ = euler_integrate(x_, pvector, mob_ik, setup.pop_node[i], p_foi, dt)
            ell += ell_

            y[i, k + 1] = x_.reshape((nc, nx))
            yell_death[i, k + 1] = ell[0]
            yell_infection[i, k + 1] = ell[1]

    results = pd.DataFrame(columns=['date', 'comp', 'place', 'cat', 'value', 'placeID'])

    for nd in range(M):
        results = pd.concat([results, pd.DataFrame.from_dict(
                    {'value': yell_infection[nd],
                     'date': setup.model_days,
                     'place': setup.ind2name[nd],
                     'cat': 'all',
                     'placeID': int(nd),
                     'comp': 'yell_infection'}),
                pd.DataFrame.from_dict(
                    {'value': yell_death[nd],
                     'date': setup.model_days,
                     'place': setup.ind2name[nd],
                     'cat': 'all',
                     'placeID': int(nd),
                     'comp': 'yell_death'})
                ])
        for ag_id, ag in enumerate(ages_names):
            results = pd.concat(
                [results, pd.DataFrame.from_dict(
                    {'value': np.append(controls[nd, :, ag_id], controls[nd, -1, ag_id]).ravel(),
                     'date': setup.model_days,
                     'place': setup.ind2name[nd],
                     'cat':ages_names[ag_id],
                     'placeID': int(nd),
                     'comp': 'vacc'})
                 ])
            for i, st in enumerate(states_names):
                results = pd.concat(
                    [results, pd.DataFrame.from_dict({'value': y[nd, :, ag_id, i].ravel(),
                                                      'date': setup.model_days,
                                                      'place': setup.ind2name[nd],
                                                      'cat': ages_names[ag_id],
                                                      'placeID': int(nd),
                                                      'comp': st})])
    results['placeID'] = results['placeID'].astype(int)

    if save_to is not None:
        results.to_csv(f'{save_to}.csv', index=False)

    return results, y, yell_death, yell_infection, mob


class COVIDVaccinationOCP:
    def __init__(self, N, n_int_steps, setup, parameters, integ='rk4', objective='death'):
        timer_start = timer()
        self.N = N
        self.n_int_steps = n_int_steps
        self.setup = setup
        self.parameters = parameters
        self.M = M = setup.nnodes

        if objective == 'death':
            ellID = 0
        elif objective == 'infection':
            ellID = 1

        _, pvector_names = parameters.get_pvector()

        dt = (N + 1) / N / n_int_steps

        statesin = cat.struct_symSX(states_names)
        # [S, E, P, I, A, Q, H, R, V] = states[...]

        states = cat.struct_symSX([
            cat.entry("Y", struct=statesin),
            cat.entry("M", struct=statesin),
            cat.entry("O", struct=statesin)]
        )

        # controls = cat.struct_symSX(['v', 'mob'])
        # [_, mob] = controls[...]

        controls = cat.struct_symSX(['vY', 'vM', 'vO', 'mob'])
        [vY, vM, vO, mob] = controls[...]

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
                           [rhs, rhs_ell[ellID]])  # CHANGE HERE FOR DEATH -> INFECTION

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

        # x_ = states.cat
        x_ = ca.veccat(*states[...])
        u_ = ca.veccat(*controls[...])
        self.x_1 = copy.copy(x_)
        for ag_id, ag in enumerate(ages_names):
            VacPpl = states[f'{ag}', 'S'] + states[f'{ag}', 'E'] + states[f'{ag}', 'P'] + states[f'{ag}', 'A'] + states[f'{ag}', 'R']

            vaccrate = controls[f'v{ag}'] / VacPpl # + 1e-10)

            x_[0 + nx * ag_id] -= vaccrate * states[f'{ag}', 'S']
            x_[1 + nx * ag_id] -= vaccrate * states[f'{ag}', 'E']
            x_[2 + nx * ag_id] -= vaccrate * states[f'{ag}', 'P']
            x_[4 + nx * ag_id] -= vaccrate * states[f'{ag}', 'A']
            x_[7 + nx * ag_id] -= vaccrate * states[f'{ag}', 'R']
            x_[8 + nx * ag_id] += controls[f'v{ag}'] #* states[f'{ag}', 'S']

        self.x_2 = copy.copy(x_)
        ell = 0.
        for k in range(n_int_steps):
            x_, ell_ = rk4_step(x_, u_, covar, params, pop_nodeSX, p_foiSX)
            ell += ell_

        rk4_int = ca.Function('rk4_int', [states, ca.veccat(controls, covar, params, pop_nodeSX, p_foiSX)], [x_, ell],
                              ['x0', 'p'], ['xf', 'qf'])

        ell = ca.Function('ell', [states, controls, covar, params, pop_nodeSX, p_foiSX],
                          [scale_ell * ell + scale_v * (vO * vO + vM * vM + vY * vY), scale_ell * ell,
                           scale_v * (vO * vO + vM * vM + vY * vY)])  # Very dependent on regularization factor

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
        rate = [None] * N
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
            Sk, Ek, Pk, Rk, Ak, Ik, Vk = [ca.veccat(*self.Vars['x', :, k, ag, 'S']) for ag in ages_names], \
                                         [ca.veccat(*self.Vars['x', :, k, ag, 'E']) for ag in ages_names], \
                                         [ca.veccat(*self.Vars['x', :, k, ag, 'P']) for ag in ages_names], \
                                         [ca.veccat(*self.Vars['x', :, k, ag, 'R']) for ag in ages_names], \
                                         [ca.veccat(*self.Vars['x', :, k, ag, 'A']) for ag in ages_names], \
                                         [ca.veccat(*self.Vars['x', :, k, ag, 'I']) for ag in ages_names], \
                                         [ca.veccat(*self.Vars['x', :, k, ag, 'V']) for ag in ages_names]

            foi_sup = []
            foi_inf = []
            for n in range(M):
                Aksum, Pksum = 0, 0
                for ag_id, ag in enumerate(['Y', 'M']):
                    Aksum += ag_trans_mult[ag] * Ak[ag_id][n]
                    Pksum += ag_trans_mult[ag] * Pk[ag_id][n]
                foi_sup.append(parameters.params_structural['betaP0'] * betaR[n] * (
                        Pksum + parameters.params_structural['epsilonA'] * Aksum))
                nfoi_inf = 0
                for ag_id, ag in enumerate(['Y', 'M']):
                    nfoi_inf += Sk[ag_id][n] + Ek[ag_id][n] + Pk[ag_id][n] + Rk[ag_id][n] + Ak[ag_id][n] + Vk[ag_id][n]
                foi_inf.append(nfoi_inf)

            foi = []
            for m in range(M):
                Iksum, IksumVanilla = 0, 0
                for ag_id, ag in enumerate(['Y', 'M']):
                    Iksum += ag_trans_mult[ag] * Ik[ag_id][n]
                    IksumVanilla += Ik[ag_id][n]
                foi.append((sum(C[n, m] * foi_sup[n] for n in range(M)) + parameters.params_structural['epsilonI'] *
                            parameters.params_structural['betaP0'] * betaR[m] * Iksum) /
                           (sum(C[l, m] * foi_inf[l] for l in range(M)) + IksumVanilla + 1e-10))
            dyn[k] = []
            spatial[k] = []
            rate[k] = []
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
                mob_ik = sum(C[i, m] * foi[m] for m in range(M)) * mob_scaling

                # spatial, vaccines and dyn are put in g(x),
                # with constraints that spatial and dyn are equal to zero
                # thus imposing the dynamics and coupling.
                spatial[k].append(self.Vars['u', i, k, 'mob'] - mob_ik)
                #rate[k].append(self.Vars['u', i, k, 'vAll'] - sum(self.Vars['u', i, k, f'v{ag}'] for ag in ages_names))
                rate[k].append(sum(self.Vars['u', i, k, f'v{ag}'] for ag in ages_names))
                Sgeq0[k].append([sum(self.Vars['x', i, k, ag, comp] for comp in ['S', 'E', 'P', 'A', 'R']) - self.Vars['u', i, k, f'v{ag}'] for ag in ages_names])
                # Number of vaccine spent = num of vaccine rate * 7 (number of days)
                #vaccines += sum(self.Vars['u', i, k, f'v{ag}'] for ag in ages_names) * (N + 1) / N

            subsum = 0
            for j in range(k+1):
                for ag in ages_names:
                    subsum += sum(self.Vars['u', :, j, f'v{ag}'])
            vaccines[k] = subsum

        f /= (N + 1)  # Average over interval for cost ^ but not terminal cost

        # cat.dotdraw(mob_ik, figsize=(10, 10))
        print('-----> Writing constraints, ...', end='')
        self.g = cat.struct_MX([
            cat.entry("dyn", expr=dyn),
            cat.entry("spatial", expr=spatial),
            cat.entry("rate", expr=rate),
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
        print(f'DONE in {timer() - tsnlf:.1f} s')  # Check here

        # print('-----> Building Jacobian function...', end='')
        # self.nlpJac = self.nlpFun.factory('nlpJac', ['i0', 'i1'], ['jac:o1:i0'])
        # print('DONE')

        print('-----> Building Solver...', end='')
        tsbs = timer()
        options = {'ipopt': {}}
        options['ipopt']["linear_solver"] = "ma86"  # "ma57"  "ma86"
        # options['ipopt']["print_level"] = 12
        options['ipopt']["max_iter"] = 10000  # prevent of for beeing clogged in a good scenario
        options['ipopt']["print_info_string"] = "yes"
        #options['ipopt']["derivative_test"] = "second-order"

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

        for ag in ages_names:
            lbx['u', :, :, f'v{ag}'] = 0

        for k in range(self.N):
            for nd in range(self.M):
                ubg['rate', k, nd] = maxvaccrate_regional[nd, k]
                lbg['rate', k, nd] = -np.inf


        # Set initial conditions as constraints
        for cp, name in enumerate(states_names):
            for i in range(self.M):
                for ag_id, ag in enumerate(ages_names):
                    lbx['x', i, 0, ag, name] = ubx['x', i, 0, ag, name] = parameters.x0_ag[i, ag_id, cp]

        # ----> Starting value of the optimizer
        init = self.Vars(0)
        for i, name in enumerate(states_names):
            for k in range(self.N + 1):
                for nd in range(self.M):
                    for ag_id, ag in enumerate(ages_names):
                        init['x', nd, k, ag, name] = states_initial[nd, k, ag_id, i]

        for k in range(self.N):
            for nd in range(self.M):
                for ag_id, ag in enumerate(ages_names):
                    init['u', nd, k, f'v{ag}'] = control_initial[nd, k, ag_id]
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
        results = pd.DataFrame(columns=['date', 'comp', 'place', 'cat', 'value', 'placeID'])

        for nd in range(self.M):
            for ag_id, ag in enumerate(ages_names):
                results = pd.concat(
                    [results, pd.DataFrame.from_dict(
                        {'value': np.array(
                            ca.veccat(ca.veccat(*self.opt['u', nd, :, f'v{ag}']), self.opt['u', nd, -1, f'v{ag}'])).ravel(),
                         'date': self.setup.model_days,
                         'place': self.setup.ind2name[nd],
                         'cat': ages_names[ag_id],
                         'placeID': int(nd),
                         'comp': 'vacc'})])
                for i, st in enumerate(states_names):
                    results = pd.concat(
                        [results, pd.DataFrame.from_dict({'value': np.array(ca.veccat(*self.opt['x', nd, :, ag, st])).ravel(),
                                                          'date': self.setup.model_days,
                                                          'place': self.setup.ind2name[nd],
                                                          'cat': ages_names[ag_id],
                                                          'placeID': int(nd),
                                                          'comp': st})])
        results['placeID'] = results['placeID'].astype(int)
        results.to_csv(f'{self.scenario_name}.csv', index=False)