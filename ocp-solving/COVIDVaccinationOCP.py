#!/usr/bin/env python
# coding: utf-8

import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import casadi as ca
import casadi.tools as cat
import time
import copy
import networkx
from scipy.integrate import solve_ivp

# from SQPmethod import *

if "Agg" not in mpl.get_backend():
    mpl.interactive(True)

plt.ion()


def rhs_py(t, x, u, cov, p, mob, pop_node):
    """ 
       The one that work ^^
    """
    I, R, V = x[0], x[1], x[2], x[3]
    v = u[0]
    J = cov[0]
    sigma = p[0]
    beta = p[1]
    mu_b = p[2]
    gamma = p[3]
    theta = p[4]
    lam = p[5]
    mu = p[6]
    rho_v = p[7]
    rho = p[8]
    alpha = p[9]
    r = p[10]
    m = p[11]

    foi = beta * ((1 - m) * B / (1 + B) + m * mob)
    foi *= ca.atan(1e5 * foi) / np.pi + 0.5
    # foi += 5e-5
    # foi = beta * mob      # Mobility term now has everything
    rhs = [None] * 4

    S = pop_node - R - V - I
    # S *= ca.atan(1e4*S)/np.pi+0.5
    # S += 3.2e-5

    dI = sigma * foi * S - (gamma + mu) * I
    dI *= ca.atan(1e4 * (I)) / np.pi + 0.5
    # dI += 5e-4

    dR = (1 - sigma) * foi * S + gamma * I - (
                rho + mu + v / (pop_node - V - I + 1 + 1e2 * ca.exp(-1e2 * (pop_node - V - I - 1)))) * R
    dR *= ca.atan(1e4 * (R)) / np.pi + 0.5
    # dR += 5e-4

    dB = -mu_b * B + theta * I * (1. + lam * J)
    dB *= ca.atan(1e4 * B) / np.pi + 0.5
    dB += 5e-4

    dV = v * (pop_node - V - I) / (pop_node - V - I + 1 + 1e2 * ca.exp(-1e2 * (pop_node - V - I - 1))) - (
                rho_v + mu) * V
    dV *= ca.atan(1e4 * (V)) / np.pi + 0.5
    # dV += 5e-4

    rhs[0] = dI
    rhs[1] = dR
    rhs[2] = -mu_b * B + theta * I * (1. + lam * J)
    rhs[3] = dV

    rhs_ell = sigma * foi * S

    return rhs, rhs_ell


def rhs_py_good(t, x, u, cov, p, mob, pop_node):
    """ 
       The one that work ^^
    """
    I, R, B, V = x[0], x[1], x[2], x[3]
    v = u[0]
    J = cov[0]
    sigma = p[0]
    beta = p[1]
    mu_b = p[2]
    gamma = p[3]
    theta = p[4]
    lam = p[5]
    mu = p[6]
    rho_v = p[7]
    rho = p[8]
    alpha = p[9]
    r = p[10]
    m = p[11]

    foi = beta * ((1 - m) * B / (1 + B) + m * mob)
    foi *= ca.atan(1e7 * foi) / np.pi + 0.5
    foi += 4e-8
    # foi = beta * mob      # Mobility term now has everything
    rhs = [None] * 4
    S = pop_node - R - V - I
    rhs[0] = sigma * foi * S - (gamma + mu) * I
    rhs[1] = (1 - sigma) * foi * S + gamma * I - (
                rho + mu + v / (pop_node - V - I + 1 + 1e4 * ca.exp(-1e4 * (pop_node - V - I - 1)))) * R
    rhs[2] = -mu_b * B + theta * I * (1. + lam * J)
    rhs[3] = v * (pop_node - V - I) / (pop_node - V - I + 1 + 1e4 * ca.exp(-1e4 * (pop_node - V - I - 1))) - (
                rho_v + mu) * V

    rhs_ell = sigma * foi * S

    return rhs, rhs_ell


def rhs_py_notricks(t, x, u, cov, p, mob, pop_node):
    """ 
        same  as rhs_py without any tricks for stability. Is not used.
    """
    I, R, B, V = x[0], x[1], x[2], x[3]
    v = u[0]
    J = cov[0]
    sigma = p[0]
    beta = p[1]
    mu_b = p[2]
    gamma = p[3]
    theta = p[4]
    lam = p[5]
    mu = p[6]
    rho_v = p[7]
    rho = p[8]
    alpha = p[9]
    r = p[10]
    m = p[11]

    foi = beta * ((1 - m) * B / (1 + B) + m * mob)

    rhs = [None] * 4
    S = pop_node - R - V - I
    rhs[0] = sigma * foi * S - (gamma + mu) * I
    rhs[1] = (1 - sigma) * foi * S + gamma * I - (rho + mu + v / (pop_node - V - I + 1)) * R
    rhs[2] = -mu_b * B + theta * I * (1. + lam * J)
    rhs[3] = v * (pop_node - V - I) / (pop_node - V - I + 1) - (rho_v + mu) * V

    rhs_ell = sigma * foi * S

    return rhs, rhs_ell


def rhs_py_total(t, x, u, covar, p, M, c, pop_node):
    """ 
        Give the complete rhs of the equation, plugin-in the mobility as it is usually implemented. 
        I think it works well
    """
    nx = 4
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
        mob_i = sum(c[i, j] * X[j][2] / (1 + X[j][2]) for j in range(M))
        rhs = np.append(rhs, rhs_py(t, X[i], U[i], C[i], p, mob_i, pop_node[i])[0])
    return rhs


def rhs_py_total_mob(t, x, u, covar, p, M, c, mob, pop_node):
    """ 
        Same as rhs_py_total, however the mobility as to be provided as an array. It allows to better 
        reproduce what the ocp is doing. Works well.
    """
    nx = 4
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


def rk4_step_mob(dt, states, controls, covar, params, M, c, mob, pop_node):
    """ A step of RK4 with the mobility provided, should work well"""
    k1 = rhs_py_total_mob(0, states, controls, covar, params, M, c, mob, pop_node)
    k2 = rhs_py_total_mob(0, states + dt / 2 * k1, controls, covar, params, M, c, mob, pop_node)
    k3 = rhs_py_total_mob(0, states + dt / 2 * k2, controls, covar, params, M, c, mob, pop_node)
    k4 = rhs_py_total_mob(0, states + dt * k3, controls, covar, params, M, c, mob, pop_node)

    return states + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def rk4(dt, states, controls, covar, params, M, c, n_int_steps, pop_node):
    x_ = states
    for k in range(n_int_steps):
        x_ = rk4_step(dt, x_, controls, covar, params, M, c, pop_node)

    return x_


def rk4_mob(dt, states, controls, covar, params, M, c, n_int_steps, mob, pop_node):
    x_ = states
    for k in range(n_int_steps):  # Shouldn't mob be calculated here ? no because it's for 1 integ interval
        x_ = rk4_step_mob(dt, x_, controls, covar, params, M, c, mob, pop_node)

    return x_


class plot_iterates(ca.Callback):
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
    G = networkx.Graph()
    G.position = {}
    G.pop = {}
    G.vac = {}
    G.inf = {}

    for i, node in enumerate(ind2name):

        G.add_node(node)

        G.position[node] = (pos_node[i, 0], pos_node[i, 1])
        G.pop[node] = pop_node[i]
        G.vac[node] = max(np.array(ca.veccat(*opt['x', i, :, 'V'])))[0]
        G.inf[node] = sum(np.array(ca.veccat(*opt['x', i, :, 'I'])))[0]

        for j, connection in enumerate(mobility[i]):
            G.add_edge(node, ind2name[j], weight=connection,
                       exchange=sum(mobility[i, j] * opt['x', j, k, 'B'] / (1 + opt['x', j, k, 'B']) for k in range(N)))

    return G


class COVIDVaccinationOCP:
    def __init__(self, N, T, n_int_steps, scaling, setup, model_params, obj_params, plot_iterates_flag=False,
                 optimize=True):
        self.N = N
        self.T = T
        self.n_int_steps = n_int_steps
        self.scaling = scaling
        self.pos_node = setup.pos_node

        ind2name = setup.ind2name
        ind_to_plot = setup.ind_to_plot
        pop_node = setup.pop_node

        # First thing to do is to apply the scaling:

        self.pop_node = pop_node * scaling
        M = setup.nnodes
        mobility = setup.mobility
        ic = copy.deepcopy(setup.ic)
        ic['I'] = ic['I'] * scaling
        ic['R'] = ic['R'] * scaling
        ic['V'] = ic['V'] * scaling

        self.ic = ic  # scaled by scaling

        print(f'Building OCP with {M} nodes')

        c = setup.mobility

        sigma = model_params['sigma']
        beta = model_params['beta']
        mu_b = model_params['mu_b']
        gamma = model_params['gamma']
        theta = model_params['theta'] / scaling
        lam = model_params['lam']
        mu = model_params['mu']
        rho_v = model_params['rho_v']
        rho = model_params['rho']
        alpha = model_params['alpha']
        r = model_params['r']
        m = model_params['m']

        scale_ell = obj_params['scale_ell']
        scale_If = obj_params['scale_If']
        scale_v = obj_params['scale_v']

        pnum = [sigma, beta, mu_b, gamma, theta, lam, mu, rho_v, rho, alpha, r, m, scale_ell, scale_If, scale_v]

        # ---- decision variables ---------
        states = cat.struct_symSX(['I', 'R', 'V'])
        [I, R, V] = states[...]

        controls = cat.struct_symSX(['v', 'mob'])
        [v, mob] = controls[...]

        covar = cat.struct_symSX(['J'])
        [J] = covar[...]

        params = cat.struct_symSX(
            ['sigma', 'beta', 'mu_b', 'gamma', 'theta', 'lam', 'mu', 'rho_v', 'rho', 'alpha', 'r', 'm', 'scale_ell',
             'scale_If', 'scale_v'])
        [sigma, beta, mu_b, gamma, theta, lam, mu, rho_v, rho, alpha, r, m, scale_ell, scale_If, scale_v] = params[...]

        pop_nodeSX = ca.SX.sym('pop_node')
        # mob = ca.SX.sym('mob')

        self.states = states
        self.controls = controls
        self.covar = covar
        self.params = params

        # The rhs is at time zero, the time is also no used in the equation so that exoplain
        rhs, rhs_ell = rhs_py(0, states.cat, controls.cat, covar.cat, params.cat, mob, pop_nodeSX)
        rhs = ca.veccat(*rhs)

        # t = np.linspace(0, T, N+1)
        # sol0 = solve_ivp(lambda t,y: rhs_py(t, y, [0,0],[0], pnum),
        #  [0,T], [ic[name][0] for name in states.keys()],t_eval=t)

        frhs = ca.Function('frhs', [states, controls, covar, params, pop_nodeSX],
                           [rhs, scale_ell * rhs_ell + scale_v * v * v])

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

        x_ = ca.veccat(*states[...])
        ell = 0.
        # x_eul = x_ + dt*k1
        for k in range(n_int_steps):
            x_, ell_ = rk4_step(x_, controls, covar, params, pop_nodeSX)
            ell += ell_

        rk4_int = ca.Function('rk4_int', [states, ca.veccat(controls, covar, params, pop_nodeSX)], [x_, ell_],
                              ['x0', 'p'], ['xf', 'qf'])
        # rk4 = ca.Function('rk4', [states, controls], [x_eul])

        # dae = {'x':states, 'p':ca.veccat(controls,covar,params), 'ode':rhs, 'quad':scale_ell*rhs_ell + scale_v*v*v}
        # # opts = {'tf':T/N, 'collocation_scheme':'radau', 'interpolation_order':2, 'number_of_finite_elements':5}
        # # F = ca.integrator('F', 'collocation', dae, opts)
        # opts = {'tf':T/N, 'number_of_finite_elements':n_int_steps}
        # F = ca.integrator('F', 'rk', dae, opts)

        # BUG TODO Isn't this a double multiliation by the scale parameter since ell is already multiplied ?
        ell = ca.Function('ell', [states, controls, covar, params, pop_nodeSX],
                          [scale_ell * ell + scale_v * v * v, scale_ell * ell,
                           scale_v * v * v])  # Very dependent on regularization factor

        Vars = cat.struct_symMX([
            (
                cat.entry("x", struct=states, repeat=[M, N + 1]),
                cat.entry("u", struct=controls, repeat=[M, N]),
            ),
        ])

        Params = cat.struct_symMX([
            (
                cat.entry("cov", struct=covar, repeat=[M, N]),  # TODO BUG this make the rainfall weekly !
                cat.entry("p", struct=params),
            ),
        ])

        self.Vars = Vars
        self.Params = Params
        # This initialize
        lbx = Vars(-np.inf)
        ubx = Vars(np.inf)

        f = 0
        vaccines = 0
        cases = 0
        reg = 0
        cdot_T = 0
        dyn = [None] * M
        spatial = [None] * M
        for i in range(M):
            dyn[i] = []
            spatial[i] = []
            for k in range(N):
                # Fk = F(x0=Vars['x',i,k],p=ca.veccat(Vars['u',i,k],Params['cov',i,k],Params['p']))
                # X_ = Fk['xf']
                # ell_ik = Fk['qf']
                [X_, ell_ik] = rk4_int(Vars['x', i, k],
                                       ca.veccat(Vars['u', i, k], Params['cov', i, k], Params['p'], self.pop_node[i]))
                dyn[i].append(Vars['x', i, k + 1] - X_)
                ell_ik_, cases_ik, reg_ik = ell(Vars['x', i, k], Vars['u', i, k], Params['cov', i, k], Params['p'],
                                                self.pop_node[i])
                f += ell_ik
                cases += cases_ik
                reg += reg_ik

                mob_ik = sum(c[i, j] * Vars['x', j, k, 'B'] / (1 + Vars['x', j, k, 'B']) for j in range(M))
                spatial[i].append(Vars[
                                      'u', i, k, 'mob'] - mob_ik)  # spatial, vaccines and dyn are put in g(x), with constraints that spatial and dyn are equal to zero
                # thus imposing the dynamics and coupling.
                vaccines += Vars[
                                'u', i, k, 'v'] * T / N  # Number of vaccine spent = num of vaccine rate * 7 (number of days)

        f /= T  # Average over interval for cost ^ but not terminal cost .

        # Final time objective:
        for i in range(M):
            Bi1Bi = Vars['x', i, -1, 'B'] / (1 + Vars['x', i, -1, 'B'])
            mob_i = sum(c[i, j] * Vars['x', j, -1, 'B'] / (1 + Vars['x', j, -1, 'B']) for j in range(M))
            foi = Params['p', 'beta'] * ((1 - Params['p', 'm']) * Bi1Bi + Params['p', 'm'] * mob_i)
            # R0 = Params['p','theta']*Params['p','beta']*Vars['x',i,-1,'S']
            # R0 /= Params['p','mu_b']*( Params['p','gamma'] + Params['p','mu'] + Params['p','alpha'] )
            # f += Params['p','scale_If']*Vars['x',i,-1,'I']
            f += Params['p', 'scale_If'] * Params['p', 'sigma'] * (
                        self.pop_node[i] - sum(Vars['x', i, -1, name] for name in ['R', 'I', 'V'])) * foi
            cdot_T += Params['p', 'scale_If'] * Params['p', 'sigma'] * (
                        self.pop_node[i] - sum(Vars['x', i, -1, name] for name in ['R', 'I', 'V'])) * foi

        g = cat.struct_MX([
            cat.entry("dyn", expr=dyn),
            cat.entry("spatial", expr=spatial),
            cat.entry("vaccines", expr=vaccines),
        ])

        costTerms = ca.Function('costTerms', [Vars, Params], [cases, reg, cdot_T])
        self.costTerms = costTerms

        # This initialize
        lbg = g(0)
        ubg = g(0)

        self.g = g

        ubg['vaccines'] = 100000. * M
        lbg['vaccines'] = -np.inf

        lbx['u', :, :, 'v'] = 0.
        ubx['u', :, :, 'v'] = 15000 * optimize  # = 0 if we don't want to optimize
        # ubx['u',:,:,'v'] = 0
        ubx['u', :, :1, 'v'] = 0.

        # Set initial conditions as constraints
        for name in states.keys():
            for i in range(M):
                lbx['x', i, 0, name] = ubx['x', i, 0, name] = ic[name][i]

        # NLOP arguments:
        # 'x' : variable names to optimize: here states and control
        # 'p' : params that can be change after
        # 'f' is the objective function of 'x', that we ought to minize
        # 'g' is a function that goes lbg < g(x,p) < ubg. If you want equality constraint then ubg = lbg = the number yoiu want it to be.
        nlp = {'x': Vars, 'p': Params, 'f': f, 'g': g}
        self.nlpFun = ca.Function('nlpFun', [Vars, Params], [f, g])
        self.nlpJac = self.nlpFun.factory('nlpJac', ['i0', 'i1'], ['jac:o1:i0'])

        plotIterates = plot_iterates('plot_iterates', Vars.size, g.size, Params.size, ind_to_plot, T, N, Vars, ind2name,
                                     mobility, self.pos_node, self.pop_node, self.scaling)
        self.plotIterates = plotIterates

        options = {}
        options['ipopt'] = {}
        options['ipopt']["linear_solver"] = "ma57"
        # options['ipopt']["linear_solver"] = "ma86"
        # options['ipopt']["linear_solver"] = "ma97"
        # options['ipopt']['bound_frac'] = 1e-4
        # options['ipopt']['bound_push'] = 1e-4
        # options['ipopt']['slack_bound_frac'] = 1e-4
        # options['ipopt']['slack_bound_push'] = 1e-6
        # options['ipopt']['hessian_constant'] = 'yes'
        # options['ipopt']["tol"] = 1e-8
        if plot_iterates_flag == True:
            options["iteration_callback"] = plotIterates
        self.solver = ca.nlpsol('solver', "ipopt", nlp, options)

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

        init = Vars(0)
        for i, name in enumerate(states.keys()):
            for k in range(N + 1):
                # init['x',:,k,name] = sol0.y[i,k]
                init['x', :, k, name] = ic[name][i]
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

        self.arg = arg

        print('Done builing OCP, ready for use.')

    def solveOCP(self):

        sol = self.solver(**self.arg)
        self.sol = sol

        opt = self.Vars(sol['x'])
        lam_g = self.g(sol['lam_g'])
        lam_x = self.Vars(sol['lam_x'])
        [fnum, gnum] = self.nlpFun(opt, self.arg['p'])
        Jgnum = self.nlpJac(opt, self.arg['p'])
        gnum = self.g(gnum)

        self.opt = opt
        self.lam_g = lam_g
        self.lam_x = lam_x

        return opt
        # [fnum1,gnum1] = nlpFun(init)
        # [fnum2,gnum2] = nlpFun(opt)

        # plt.ion()

        # dT = dt*n_int_steps

        # # Plot solution from the solver:
        # #fig, axes = plt.subplots(2, 3,num=10)
        # #for ax in axes.flat: ax.cla()
        # #axes[0,0].plot(sol0.t,sol0.y[0],label="S", color = 'orange', linestyle = '--', linewidth=1)
        # #axes[0,1].plot(sol0.t,sol0.y[1],label="I", color = 'darkred', linewidth=2)
        # #axes[0,2].plot(sol0.t,sol0.y[2],label="R", color = 'red', linestyle = '-.', linewidth=1)
        # #axes[1,0].plot(sol0.t,sol0.y[3],label="B", color = 'blue', linestyle = '-.', linewidth=1)
        # #axes[1,2].plot(sol0.t,sol0.y[4],label="V", color = 'orange', linewidth=2)
        # #for ax in axes.flat:
        # #   ax.legend(loc="upper right")
        # #   ax.grid()
        # #fig.tight_layout()

        # for i in range(M):
        #     fig, axes = plt.subplots(2, 3,num=2+i)
        #     for ax in axes.flat: ax.cla()
        #     axes[0,0].plot(np.arange(0, T+dT, dT),np.array(ca.veccat(*opt['x',i,:,'S'])),label="S", color = 'orange', linestyle = '--', linewidth=1)
        #     axes[0,1].plot(np.arange(0, T+dT, dT),np.array(ca.veccat(*opt['x',i,:,'I'])),label="I", color = 'darkred', linewidth=2)
        #     axes[0,2].plot(np.arange(0, T+dT, dT),np.array(ca.veccat(*opt['x',i,:,'R'])),label="R", color = 'red', linestyle = '-.', linewidth=1)
        #     axes[1,0].plot(np.arange(0, T+dT, dT),np.array(ca.veccat(*opt['x',i,:,'B'])),label="B", color = 'blue', linestyle = '-.', linewidth=1)
        #     axes[1,1].step(np.arange(0, T+dT, dT),np.array(ca.veccat(ca.veccat(*opt['u',i,:,'v']),opt['u',i,-1,'v'])),'k',label=r"$\nu(t)$",where="post")
        #     #axes[1,2].step(np.arange(0, T+dT, dT),np.array(arg['p']['p',:,i]+[arg['p']['p',-1,i]]),'b',label=r"$J(t)$",where="post")
        #     axes[1,2].plot(np.arange(0, T+dT, dT),np.array(ca.veccat(*opt['x',i,:,'V'])),label="V", color = 'orange', linewidth=2)
        #     for ax in axes.flat:
        #         ax.legend(loc="upper right")
        #         ax.grid()
        #     fig.suptitle(ind2name[i])
        #     fig.tight_layout()

        # plt.figure(3+i)
        # plt.clf()
        # plt.spy(Jgnum)

        # fig, ax = plt.subplots(1, 1,num=1)
        # for ax in axes.flat: ax.cla()

        # axt =  ax.twinx()

        # ax.plot(np.arange(0, T+dT, dT),np.array(ca.veccat(*opt['x',:,i,'S'])),label="S", color = 'orange', linestyle = '--', linewidth=.5)
        # ax.plot(np.arange(0, T+dT, dT),np.array(ca.veccat(*opt['x',:,i,'I'])),label="I", color = 'darkred', linewidth=2)
        # ax.plot(np.arange(0, T+dT, dT),np.array(ca.veccat(*opt['x',:,i,'R'])),label="R", color = 'red', linestyle = '-.', linewidth=.5)
        # ax.plot(np.arange(0, T+dT, dT),1e5*np.array(ca.veccat(*opt['x',:,i,'B'])),label="B", color = 'blue', linestyle = '-.', linewidth=.5)
        # #plt.plot(limit(np.arange(s.N+1)),'r--',label="speed limit")
        # axt.step(np.arange(0, T+dT, dT),np.array(ca.veccat(ca.veccat(*opt['u',:,i,'v']),opt['u',-1,i,'v'])),'k',label=r"$\nu(t)$",where="post")
        # axt.step(np.arange(0, T+dT, dT),np.array(arg['p']['p',:,i]+[arg['p']['p',-1,i]]),'b',label=r"$J(t)$",where="post")
        # ax.legend(loc="upper right")
        # axt.legend(loc="upper left")

        # ax.set_xlabel(r'Time (days)')
        # ax.grid()
        # #ax.set_title('Basic reproduction number & rainfall', fontsize=14)
        # ax.set_ylabel(r'# of individuals')
        # axt.set_ylabel(r'Vaccination rate [vacc./day]')
        # fig.tight_layout()

        # fig, axes = plt.subplots(2, 3,num=3)
        # for ax in axes.flat: ax.cla()
        # axes[0,0].plot(np.arange(0, T+dT, dT),np.array(ca.veccat(*opt['x',:,i,'S'])).squeeze()-sol0.y[0],label="S", color = 'orange', linestyle = '--', linewidth=1)
        # axes[0,1].plot(np.arange(0, T+dT, dT),np.array(ca.veccat(*opt['x',:,i,'I'])).squeeze()-sol0.y[1],label="I", color = 'darkred', linewidth=2)
        # axes[0,2].plot(np.arange(0, T+dT, dT),np.array(ca.veccat(*opt['x',:,i,'R'])).squeeze()-sol0.y[2],label="R", color = 'red', linestyle = '-.', linewidth=1)
        # axes[1,0].plot(np.arange(0, T+dT, dT),np.array(ca.veccat(*opt['x',:,i,'B'])).squeeze()-sol0.y[3],label="B", color = 'blue', linestyle = '-.', linewidth=1)
        # axes[1,1].step(np.arange(0, T+dT, dT),np.array(ca.veccat(ca.veccat(*opt['u',:,i,'v']),opt['u',-1,i,'v'])),'k',label=r"$\nu(t)$",where="post")
        # #axes[1,2].step(np.arange(0, T+dT, dT),np.array(arg['p']['p',:,i]+[arg['p']['p',-1,i]]),'b',label=r"$J(t)$",where="post")
        # axes[1,2].plot(np.arange(0, T+dT, dT),np.array(ca.veccat(*opt['x',:,i,'V'])),label="V", color = 'orange', linewidth=2)
        # for ax in axes.flat:
        #     ax.legend(loc="upper right")
        #     ax.grid()
        # fig.tight_layout()

        # ---- post-processing        ------
        # plt.figure()
        # plt.spy(sol.value(cadi.jacobian(opti.g,opti.x)))
        # plt.figure()
        # plt.spy(sol.value(cadi.hessian(opti.f+cadi.dot(opti.lam_g,opti.g),opti.x)[0]))
