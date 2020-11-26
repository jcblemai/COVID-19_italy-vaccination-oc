import numpy as np
import pandas as pd
import copy
import datetime

nx = 9
states_names = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']
S, E, P, I, A, Q, H, R, V = np.arange(nx)


class OCParameters:
    def __init__(self, eng, setup, M, freq='1D'):
        self.M = M
        self.mobintime = setup.mobility_ts.resample(freq).mean()

        matlab_start_date = datetime.date(2020, 1, 20)  # fix lentgh
        matlab_end_date = datetime.date(2020, 7, 1)
        matlab_model_days = pd.date_range(matlab_start_date, matlab_end_date, freq='1D')
        p_dict, self.mobfrac, self.mobmat, self.betaratiointime, self.x0 = get_parameters_from_matlab(eng,
                                                                                                      setup,
                                                                                                      M,
                                                                                                      matlab_model_days,
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
            'scale_ell': 1e5,
            'scale_If': 0,
            'scale_v': 1e-10
        }

        mob_prun = 0.0006
        self.mobmat_pr = self.prune_mobility(mob_prun)

        # Numpy array from dataframes:
        self.mobintime_arr = self.mobintime.to_numpy().T
        self.betaratiointime_arr = self.betaratiointime.to_numpy().T

    def prune_mobility(self, mob_prun):
        mobK = self.mobintime.to_numpy().T[:, 0]
        betaR = self.betaratiointime.to_numpy().T[:, 0]
        C = self.params_structural['r'] * self.mobfrac.flatten() * mobK * self.mobmat
        np.fill_diagonal(C, 1 - C.sum(axis=1) + C.diagonal())
        # print(f'pruning {C[C < mob_prun].size} non-diagonal mobility elements of {C.size - self.M}.')
        C[C < mob_prun] = 0  # Prune elements
        mobmat_pr = np.copy(self.mobmat)

        mobmat_pr[C == 0] = 0
        print(f'nnz before: {np.count_nonzero(self.mobmat)}, after: {np.count_nonzero(mobmat_pr)}')

        return mobmat_pr

    def get_pvector(self):
        pvector = params_to_vector(self.model_params)
        for v in params_to_vector(self.hyper_params):
            pvector.append(v)
        pvector_names = list(self.model_params.keys()) + list(self.hyper_params.keys())
        return pvector, pvector_names


def get_parameters_from_matlab(eng, s, model_size, model_days, freq):
    p = {}  # 1D parameters
    p['deltaE'] = eng.eval('deltaE')
    p['deltaP'] = eng.eval('deltaP')
    p['sigma'] = eng.eval('sigma')
    p['eta'] = eng.eval('eta')
    p['gammaI'] = eng.eval('gammaI')
    p['gammaA'] = eng.eval('gammaA')
    p['gammaQ'] = eng.eval('gammaQ')
    p['gammaH'] = eng.eval('gammaH')
    p['alphaI'] = eng.eval('alphaI')
    p['alphaH'] = eng.eval('alphaH')
    p['zeta'] = eng.eval('V.zeta')
    p['eta'] = eng.eval('eta')
    p['r'] = eng.eval('r')
    p['betaP0'] = eng.eval('betaP0')
    p['epsilonA'] = eng.eval('epsilonI')
    p['epsilonI'] = eng.eval('epsilonA')
    p['gammaV'] = 1 / (9 * 30)
    x0_matlab = np.array(eng.eval('V.x0')).flatten()
    x0 = np.zeros(9 * s.nnodes)
    for i in range(s.nnodes):
        x0[i * nx:(i + 1) * nx] = [x0_matlab[107 * S + i],
                                   x0_matlab[107 * E + i],
                                   x0_matlab[107 * P + i],
                                   x0_matlab[107 * I + i],
                                   x0_matlab[107 * A + i],
                                   x0_matlab[107 * Q + i],
                                   x0_matlab[107 * H + i],
                                   x0_matlab[107 * R + i],
                                   np.zeros_like(x0_matlab[107 * R + i])]

    beta_ratio = np.array(eng.eval('beta_ratio'))[:model_size]
    beta_ratio_ts = pd.DataFrame(beta_ratio.T, index=model_days, columns=np.arange(s.nnodes))
    betaratiointime = beta_ratio_ts.resample(freq).mean()

    mobile_frac = np.array(eng.eval('V.p'))[:model_size].flatten()
    mobility_matrix = np.array(eng.eval('full(V.q)'))[:model_size, :model_size]

    return p, mobile_frac, mobility_matrix, betaratiointime, x0


def params_to_vector(p):
    plist = []
    for key in p.keys():
        plist.append(p[key])
    return plist
