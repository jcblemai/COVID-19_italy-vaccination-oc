import numpy as np
import pandas as pd

nx = 9
states_names = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']
S, E, P, I, A, Q, H, R, V = np.arange(nx)


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
    p['gammaV'] = 1 / (9*30)
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
