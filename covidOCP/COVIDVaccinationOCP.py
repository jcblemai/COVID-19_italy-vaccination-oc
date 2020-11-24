
import casadi as ca
import casadi.tools as cat
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import networkx as nx
import sys, copy
from scipy.integrate import solve_ivp
from utils import *
import networkx
import geopandas as gpd

if "Agg" not in mpl.get_backend():
    mpl.interactive(True)

plt.ion()


def rhs_py(t, x, u, cov, p, mob, pop_node):
    S, E, P, I, A, Q, H, R, V = x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]

    deltaE, deltaP,sigma,eta,gammaI,gammaA,gammaQ,gammaH,alphaI,alphaH, zeta,  gammaV =     p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11]
   
    v = u[0]
    foi = mob
    rhs = [None] * nx
    vaccrate = v / (S + E + P + A + Q + H + R+1)  # NOT I NOT H ?
    vaccrate = v / ca.sqrt((S+E+P+A+Q+H+R)**2+10)   # 10 (1e7 ?1e3 ?)is to tune depending on scale, approx abs. value (brings quantities non negative by vaccinating)
    vaccrate = 0
    rhs[0] = -(foi + vaccrate) * S + gammaV * V                 # S
    rhs[1] = foi * S - deltaE * E;                              # E 
    rhs[2] = deltaE * E - deltaP * P;                           # P
    rhs[3] = sigma * deltaP * P - (eta + gammaI + alphaI) * I;  # I
    rhs[4] = (1 - sigma) * deltaP * P - gammaA * A;             # A
    rhs[5] = zeta * eta * I - gammaQ * Q;                       # Q
    rhs[6] = (1 - zeta) * eta * I - (gammaH + alphaH) * H;      # H
    rhs[7] = gammaI * I + gammaA * A + gammaH * H + gammaQ * Q;  # R
    rhs[8] = vaccrate * S - gammaV * V                           # V
    rhs_ell = [None] * 3
    rhs_ell[0] = gammaH * H;  # recovered from the hospital
    rhs_ell[1] = alphaH * H;  # total death
    rhs_ell[2] = (1 - zeta) * eta * I;  # cumulative hospitalized cases

    return rhs, rhs_ell





def rhs_py_total_mob(t, x, u, covar, p, M, c, mob, pop_node):
    """ 
        Same as rhs_py_total, however the mobility as to be provided as an array. It allows to better 
        reproduce what the ocp is doing. Works well.
    """
    nu = 1
    nc = 1
    X = []
    U = []
    C = []
    for i in range(M):  # Uses flat vectors
        X.append(x[i * nx:(i + 1) * nx])
        U.append(u[i * nu:(i + 1) * nu])
        C.append(covar[i * nc:(i + 1) * nc])

    rhs = np.array([])
    for i in range(M):
        rhs = np.append(rhs, rhs_py(t, X[i], U[i], C[i], p, mob[i], pop_node[i])[0])
    return rhs


def rk4_step(dt, states, controls, covar, params, M, c, pop_node):
    """ A step of RK4 with the right mobility, should work well"""
    k1 = rhs_py_total(0, states, controls, covar, params, M, c, pop_node)
    k2 = rhs_py_total(0, states + dt / 2 * k1, controls, covar, params, M, c, pop_node)
    k3 = rhs_py_total(0, states + dt / 2 * k2, controls, covar, params, M, c, pop_node)
    k4 = rhs_py_total(0, states + dt * k3, controls, covar, params, M, c, pop_node)

    return states + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)






def rk4_mob(dt, states, controls, covar, params, M, c, n_int_steps, mob, pop_node):
    x_ = states
    for k in range(n_int_steps):  # Shouldn't mob be calculated here ? no because it's for 1 integ interval
        x_ = rk4_step_mob(dt, x_, controls, covar, params, M, c, mob, pop_node)

    return x_


class PlotIterates(ca.Callback):
    def __init__(self, name, nx, ng, np, ind_to_plot, T, N, V, ind2name, mobility, pos_node, pop_node, scaling,
                 opts={}):
        ca.Callback.__init__(self)

        self.nx = nx
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
        self.scaling = scaling

        self.fig = []
        self.axes = []
        for k, i in enumerate(self.ind_to_plot):
            figk, axesk = plt.subplots(2, 3, num=2 + i)
            for ax in axesk.flat: ax.autoscale(enable=True, axis='both')
            self.fig.append(figk)
            self.axes.append(axesk)

        figk, axesk = plt.subplots(1, 1, num=100)
        self.fig.append(figk)
        self.axes.append(axesk)

        # plt.figure(1)

        # plt.title('Title')
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

        scaling = self.scaling
        pop_node = self.pop_node
        darg = {}
        for (i, s) in enumerate(ca.nlpsol_out()): darg[s] = arg[i]

        sol = self.V(darg['x'])
        self.sol = sol
        if hasattr(self, 'lines'):
            for k, i in enumerate(self.ind_to_plot):
                S = pop_node[i] - (np.array(ca.veccat(*sol['x', i, :, 'I'])) + np.array(
                    ca.veccat(*sol['x', i, :, 'R'])) + np.array(ca.veccat(*sol['x', i, :, 'V'])))
                self.lines[k][0][0].set_ydata(S / scaling)
                self.lines[k][1][0].set_ydata(np.array(ca.veccat(*sol['x', i, :, 'I'])) / scaling)
                self.lines[k][2][0].set_ydata(np.array(ca.veccat(*sol['x', i, :, 'R'])) / scaling)
                self.lines[k][3][0].set_ydata(np.array(ca.veccat(*sol['x', i, :, 'B'])))
                self.lines[k][4][0].set_ydata(
                    np.array(ca.veccat(ca.veccat(*sol['u', i, :, 'v']), sol['u', i, -1, 'v'])) / scaling)
                self.lines[k][5][0].set_ydata(np.array(ca.veccat(*sol['x', i, :, 'V'])) / scaling)
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
                S = pop_node[i] - (np.array(ca.veccat(*sol['x', i, :, 'I'])) + np.array(
                    ca.veccat(*sol['x', i, :, 'R'])) + np.array(ca.veccat(*sol['x', i, :, 'V'])))
                self.lines[k].append(
                    self.axes[k][0, 0].plot(np.arange(0, T + dT, dT), S / scaling, label="S", color='orange',
                                            linewidth=1))
                self.lines[k].append(self.axes[k][0, 1].plot(np.arange(0, T + dT, dT),
                                                             np.array(ca.veccat(*sol['x', i, :, 'I'])) / scaling,
                                                             label="I", color='darkred', linewidth=2))
                self.lines[k].append(self.axes[k][0, 2].plot(np.arange(0, T + dT, dT),
                                                             np.array(ca.veccat(*sol['x', i, :, 'R'])) / scaling,
                                                             label="R", color='red', linewidth=1))
                self.lines[k].append(
                    self.axes[k][1, 0].plot(np.arange(0, T + dT, dT), np.array(ca.veccat(*sol['x', i, :, 'B'])),
                                            label="B", color='blue', linewidth=1))
                self.lines[k].append(self.axes[k][1, 1].step(np.arange(0, T + dT, dT), np.array(
                    ca.veccat(ca.veccat(*sol['u', i, :, 'v']), sol['u', i, -1, 'v'])) / scaling, 'k', label=r"$\nu(t)$",
                                                             where="post"))
                self.lines[k].append(self.axes[k][1, 2].plot(np.arange(0, T + dT, dT),
                                                             np.array(ca.veccat(*sol['x', i, :, 'V'])) / scaling,
                                                             label="V", color='orange', linewidth=2))

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

        if True:
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


def mobility_graph(mobility, ind2name, pos_node, pop_node, opt, N):



class COVIDVaccinationOCP:
    def __init__(self, N, T, n_int_steps, scaling, setup, plot_iterates=False, optimize=True):
        self.N = N
        self.T = T
        self.n_int_steps = n_int_steps

        states = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']
        states = cat.struct_symSX(states)
        [S, E, P, I, A, Q, H, R, V] = states[...]

        controls = cat.struct_symSX(['v', 'mob'])
        [v, mob] = controls[...]

        covar = cat.struct_symSX(['mobility_t','betaratio_t'])
        [mobility_t, betaratio_t] = covar[...]

        params = cat.struct_symSX(list(model_params.keys()) + ['scale_ell', 'scale_If', 'scale_v'])
        [deltaE, deltaP, sigma, eta, gammaI, gammaA, gammaQ, gammaH, alphaI, alphaH, zeta, gammaV,
        scale_ell, scale_If, scale_v] = params[...]

        pop_nodeSX = ca.SX.sym('pop_node')


        # The rhs is at time zero, the time is also no used in the equation so that exoplain
        rhs, rhs_ell = rhs_py(0, states.cat, controls.cat, covar.cat, params.cat, mob, pop_nodeSX)
        rhs = ca.veccat(*rhs)
        rhs_ell = ca.veccat(*rhs_ell)  # mod




        frhs = ca.Function('frhs', [states, controls, covar, params, pop_nodeSX],
                                [rhs, rhs_ell[1]])#scale_ell * rhs_ell[1] + scale_v * v * v])# mod ICI juste ell[1]



        dt = T / N / n_int_steps  # length of an integration interval
        # ---- dynamic constraints --------
        k1, k1ell = frhs(states, controls, covar, params, pop_nodeSX)
        k2, k2ell = frhs(states + dt / 2 * k1, controls, covar, params, pop_nodeSX)
        k3, k3ell = frhs(states + dt / 2 * k2, controls, covar, params, pop_nodeSX)
        k4, k4ell = frhs(states + dt * k3, controls, covar, params, pop_nodeSX)
        x_next = states + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ell_next = dt / 6 * (
                k1ell + 2 * k2ell + 2 * k3ell + k4ell)  # No need for because we sum it week by week few lines below.

        rk4_step = ca.Function('rk4_step', [states, controls, covar, params, pop_nodeSX], [x_next, ell_next])



        # Overwrite last cell with eurler.
        dt = T / N / n_int_steps  # length of an integration interval
        a, b = frhs(states, controls, covar, params, pop_nodeSX)
        x_next = states + dt*a
        ell_next = dt*b
        rk4_step = ca.Function('rk4_step', [states, controls, covar, params, pop_nodeSX], [x_next, ell_next])



        x_ = ca.veccat(*states[...])
        u_ = ca.veccat(*controls[...])
        ell = 0.
        VacPpl = states['S'] + states['E'] + states['P'] + states['A'] + states['R']
        vaccrate = controls['v']/(VacPpl+1e-10)
        x_[0] -= vaccrate * states['S']
        x_[8] += vaccrate * states['S']
        for k in range(n_int_steps):
            x_, ell_ = rk4_step(x_, u_,covar, params, pop_nodeSX)
            ell += ell_



        rk4_int = ca.Function('rk4_int', [states, ca.veccat(controls, covar, params, pop_nodeSX)], [x_, ell],
                            ['x0', 'p'], ['xf', 'qf'])

        # BUG TODO Isn't this a double multiplication by the scale parameter since ell is already multiplied ?
        ell = ca.Function('ell', [states, controls, covar, params, pop_nodeSX],
                        [scale_ell * ell + scale_v * v * v, scale_ell * ell,
                        scale_v * v * v])  # Very dependent on regularization factor
        #ell = ca.Function('ell', [states, controls, covar, params, pop_nodeSX],
        #                  [ell, scale_ell * ell,
        #                   scale_v * v * v])  # Very dependent on regularization factor




        Vars = cat.struct_symMX([
            (
                cat.entry("x", struct=states, repeat=[M, N + 1]),
                cat.entry("u", struct=controls, repeat=[M, N]),
            ),
        ])

        Params = cat.struct_symMX([
            (
                cat.entry("cov", struct=covar, repeat=[M, N]), 
                cat.entry("p", struct=params),
            ),
        ])


        lbx = Vars(-np.inf)
        ubx = Vars(np.inf)

        f = 0
        vaccines = 0
        cases = 0
        reg = 0
        cdot_T = 0
        dyn = [None] * N
        spatial = [None] * N

        Sgeq0 = [None] * N

        mob_prun = 0.0006

        # Pruning mobility once:
        k=0
        mobility_history = []
        mobK = mobintime.to_numpy().T[:,k]
        betaR =  betaratiointime.to_numpy().T[:,k]
        C = r*mobfrac.flatten()*mobK*mobmat
        np.fill_diagonal(C,1-C.sum(axis=1)+ C.diagonal())
        print(f'pruning {C[C<mob_prun].size} non-diagonal mobility elements of {C.size-M}.')  
        C[C<mob_prun] = 0 # Prune elements
        mobility_history.append(C)

        mobmat_pr = np.copy(mobmat)


        mobmat_pr[mobility_history[0] == 0] = 0
        print(f'nnz before: {np.count_nonzero(mobmat)}, after: {np.count_nonzero(mobmat_pr)}')

        for k in range(N):
            mobK = Params['cov',:,k,'mobility_t'] #mobintime.to_numpy().T[:,k]
            betaR =  Params['cov',:,k,'mobility_t']# betaratiointime.to_numpy().T[:,k]
            C = r*mobfrac.flatten()*mobK*mobmat_pr
            np.fill_diagonal(C,1-C.sum(axis=1)+ C.diagonal())
            
            dyn[k] = []
            spatial[k] = []
            Sgeq0[k] = []
            Sk, Ek, Pk, Rk, Ak, Ik = ca.veccat(*Vars['x', :, k, 'S']), ca.veccat(*Vars['x', :, k, 'E']), 
                                     ca.veccat(*Vars['x', :, k, 'P']), ca.veccat(*Vars['x', :, k, 'R']), 
                                    ca.veccat(*Vars['x', :, k, 'A']), ca.veccat(*Vars['x', :, k, 'I'])

            #if k == 0 or k == N-1:
            #    print(f'pruning {C[C<mob_prun].size} non-diagonal mobility elements of {C.size-M}.')  
            #C[C<mob_prun] = 0 # Prune elements
            #mobility_history.append(C)
            
            foi_sup = []
            foi_inf = []
            for n in range(M):
                foi_sup.append(betaP0*betaR[n]*(Pk[n]+epsilonA*Ak[n]))
                foi_inf.append(Sk[n]+Ek[n]+Pk[n]+Rk[n]+Ak[n])
            
            foi = []
            for m in range(M):
                foi.append((sum(C[n, m] * foi_sup[n] for n in range(M)) + epsilonI*betaP0*betaR[m]*Ik[m]) /
                        (sum(C[l, m] * foi_inf[l] for l in range(M)) + Ik[m]))
                

            print(f'{k}:', end='')
            for i in range(M):
                [X_, ell_ik] = rk4_int(Vars['x', i, k],
                                    ca.veccat(Vars['u', i, k], Params['cov', i, k], Params['p'], pop_node[i]))

                dyn[k].append(Vars['x', i, k + 1] - X_)
                ell_ik_, cases_ik, reg_ik = ell(Vars['x', i, k], Vars['u', i, k], Params['cov', i, k], Params['p'], pop_node[i])
                f += ell_ik_ #MOD: before  ell_ik_
                cases += cases_ik
                reg += reg_ik
                #mob_ik = sum(C[i, m] * (
                #            (sum(C[n, m] * (betaP0*betaR[n]*(Pk[n]+epsilonA*Ak[n])) for n in range(M)) + 
                #             epsilonI*betaP0*betaR[m]*Ik[m]) /
                #            (sum(C[l, m]*(Sk[l]+Ek[l]+Pk[l]+Rk[l]+Ak[l]) for l in range(M)) + Ik[m]))
                #             for m in range(M))
                mob_ik = sum(C[i, m] * foi[m] for m in range(M))
                
                spatial[k].append(Vars['u', i, k, 'mob'] - mob_ik) # spatial, vaccines and dyn are put in g(x), with constraints that spatial and dyn are equal to zero
                VacPpl = sum(Vars['x', i, k,name] for name in ['S','E','P','A','R'])
                Sgeq0[k].append( Vars['x', i, k,'S'] - Vars['u', i, k, 'v']/(VacPpl+1e-10) )
                # thus imposing the dynamics and coupling.
                vaccines += Vars[ 'u', i, k, 'v'] * T / N  # Number of vaccine spent = num of vaccine rate * 7 (number of days)

        f /= T  # Average over interval for cost ^ but not terminal cost

    def build(self):

        print('Writing constraints, ...', end='')
        g = cat.struct_MX([
            cat.entry("dyn", expr=dyn),
            cat.entry("spatial", expr=spatial),
            cat.entry("vaccines", expr=vaccines),
            cat.entry("Sgeq0", expr=Sgeq0)
        ])

        costTerms = ca.Function('costTerms', [Vars, Params], [cases, reg])


        # This initialize
        lbg = g(0)
        ubg = g(0)


        ubg['vaccines'] = 2000*(T*.6)*M #8e6 #*M
        lbg['vaccines'] = -np.inf

        ubg['Sgeq0'] = np.inf

        optimize = 0
        lbx['u', :, :, 'v'] = 0.
        ubx['u', :, :, 'v'] = 2000 * optimize  # = 0 if we don't want to optimize
        # ubx['u',:,:,'v'] = 0
        #ubx['u', :, :1, 'v'] = 0.


        # Set initial conditions as constraints
        for cp, name in enumerate(states.keys()):
            for i in range(M):
                lbx['x', i, 0, name] = ubx['x', i, 0, name] = x0[i*nx+cp]

        print('DONE')
        print('Building NLP function...', end='')

        # NLOP arguments:
        # 'x' : variable names to optimize: here states and control
        # 'p' : params that can be change after
        # 'f' is the objective function of 'x', that we ought to minize
        # 'g' is a function that goes lbg < g(x,p) < ubg. If you want equality constraint then ubg = lbg = the number yoiu want it to be.
        nlp = {'x': Vars, 'p': Params, 'f': f, 'g': g}
        nlpFun = ca.Function('nlpFun', [Vars, Params], [f, g])
        print('DONE')
        print('Building Jacobian function...', end='')
        nlpJac = nlpFun.factory('nlpJac', ['i0', 'i1'], ['jac:o1:i0'])

        print('DONE')
        print('Building Solver...', end='')

        options = {}
        options['ipopt'] = {}
        #options['ipopt']["linear_solver"] = "ma57"
        options['ipopt']["linear_solver"] = "ma86"
        # options['ipopt']["linear_solver"] = "ma86"
        # options['ipopt']["linear_solver"] = "ma97"
        # options['ipopt']['bound_frac'] = 1e-4
        # options['ipopt']['bound_push'] = 1e-4
        # options['ipopt']['slack_bound_frac'] = 1e-4
        # options['ipopt']['slack_bound_push'] = 1e-6
        # options['ipopt']['hessian_constant'] = 'yes'
        # options['ipopt']["tol"] = 1e-8

        solver = ca.nlpsol('solver', "ipopt", nlp, options)

        # options = {}
        # options['ls_step_factor'] = 0.5
        # options['pr_tol'] = 1e-8
        # options['du_tol'] = 1e-8
        # options['qp_solver'] = 'nlpsol'
        # options['qp_solver_options']  = {'nlpsol':'ipopt','nlpsol_options':{'print_time':False,'ipopt.print_level':0,'ipopt.hessian_constant':'yes','ipopt.linear_solver':'ma97','ipopt.tol':1e-16}}
        # # options['regularization'] = 'regularize_reduced'
        # options['regularization'] = 'regularize_reduced_diag'
        # options['regularization_tol'] = 1e-10
        # options['linesearch'] = 'filter'
        # options["iteration_callback"] = plotIterates
        # self.solver = SQPmethod(nlp,options)

        print('DONE')


    def update(self):
        init = Vars(0)
        for i, name in enumerate(states.keys()):
            for k in range(N + 1):
                # init['x',:,k,name] = sol0.y[i,k]
                for nd in range(M):
                    init['x', nd, k, name] = integ_matlab.T[nd+107*i, k].T # p['x0']  # self.ic[name][i]

        init['u'] = 0.

        arg = {}
        arg['lbg'] = lbg
        arg['ubg'] = ubg

        arg['lbx'] = lbx
        arg['ubx'] = ubx
        arg['x0'] = init
        arg['p'] = Params()
        # for i in range(M):

        #    for k in range(N):
        #        arg['p']['cov',i,k] = rainfall_norm[i,k]**params(pnum)['r']
        # arg['p']['cov',i,k] = np.random.random()**params(pnum)['r']

        arg['p']['p'] = pnum

        for i in range(M):
            for k in range(N):
                arg['p']['cov',i,k,'mobility_t'] =  mobintime.to_numpy().T[i,k]
                arg['p']['cov',i,k,'betaratio_t'] = betaratiointime.to_numpy().T[i,k]

        

        pnum[-3] = 1e5     # scale_ell
        pnum[-2] = 0       # scale_If
        pnum[-1] = 1e-10   # scale_V
        pnum[-4] = 1/(9*30)  #gammaV 
        arg['p']['p'] = pnum
        arg['x0'] = init




    def solveOCP(self):
        sol = solver(**arg)
        opt = Vars(sol['x'])
        lam_g = g(sol['lam_g'])
        lam_x = Vars(sol['lam_x'])
        [fnum,gnum]=nlpFun(opt,arg['p'])
        Jgnum=nlpJac(opt, arg['p'])
        gnum = g(gnum)   # 2times ?
        print(f"""
        Vaccines stockpile: 
            {float(arg['ubg']['vaccines']):010f} total.
            {float(g(gnum)['vaccines']):010f} spent.
            {float((arg['ubg']['vaccines']  - g(gnum)['vaccines'])):010f} left.""")

    def plot_node(self, node):
        fig, axes = plt.subplots(2,5, figsize = (20,10))
        fig.patch.set_facecolor('white')

        node = 1
        til =  T

        for i, st in enumerate(states_list):
            axes.flat[i].plot(np.array(ca.veccat(*opt['x',node,:til,st])), 
                        linestyle=':', lw = 4, color='r')
            if st != 'V':
                axes.flat[i].plot(np.array(integ_matlab.T[node+107*i,:til].T), 
                        linestyle='-', lw = 2, color='k')

            axes.flat[i].set_title(st);

        axes.flat[-1].step(np.array(ca.veccat(ca.veccat(*opt['u',node,:til,'v']))),#,opt['u',node,-1,'v'])),
                        'k',label=r"$\nu(t)$")


    def plot_all(self):
        fig, axes = plt.subplots(5,2, figsize = (10,10))


        for i, st in enumerate(states_list):
            for k in range(M):
                axes.flat[i].plot(np.array(ca.veccat(*opt['x',k,:til,st])), 
                        lw = 2, ls = '--')
                #if st != 'V':
                #    axes.flat[i].plot(np.array(integ_matlab.T[k+107*i,:til].T), 
                #           lw = .5)
                axes.flat[i].set_title(st);
                axes.flat[-1].step(np.arange(len(np.array(ca.veccat(ca.veccat(*opt['u',k,:til,'v']))))),
                                                np.array(ca.veccat(ca.veccat(*opt['u',k,:til,'v'])))
                                )

