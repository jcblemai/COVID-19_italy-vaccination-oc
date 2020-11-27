import numpy as np
from ItalySetup import ItalySetup
from covidOCP import COVIDVaccinationOCP, COVIDParametersOCP
import pickle

matlab = True
save_param = True

# Check if convariates in N or T

nx = 9
states_names = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']

model_size = 'full'  # nodes

# Horizon for each problem
ndays = 'full'
n_int_steps = 1

setup = ItalySetup(model_size, ndays)
M = setup.nnodes
N = len(setup.model_days) - 1

if matlab:
    import matlab.engine

    eng = matlab.engine.start_matlab()
    eng.cd('geography-paper-master/', nargout=0)
    eng.run('single_build.m', nargout=0)
    integ_matlab = np.array(eng.eval('x'))

    matlab_initial = np.zeros((M, N + 1, nx))
    for i, name in enumerate(states_names):
        for k in range(N + 1):
            for nd in range(M):
                matlab_initial[nd, k, i] = integ_matlab.T[nd + 107 * i, k].T

    p = COVIDParametersOCP.OCParameters(eng=eng, setup=setup, M=M)
    if save_param:
        with open('parameters.pkl', 'wb') as out:
            pickle.dump(p, out, pickle.HIGHEST_PROTOCOL)

else:
    with open('parameters.pkl', 'rb') as inp:
        p = pickle.load(inp)

scenario_name = 'opt-test10000'
control_initial = np.zeros((M, N))
max_vacc_rate = np.zeros((M, N))
for k in range(N):
    for nd in range(M):
        max_vacc_rate[nd, k] = 0
        control_initial[nd, k] = 0

initial = np.zeros((M, N + 1, nx))
for i, name in enumerate(states_names):
    for k in range(N + 1):
        for nd in range(M):
            initial[nd, k, i] = 0

ocp = COVIDVaccinationOCP.COVIDVaccinationOCP(N=N, n_int_steps=n_int_steps,
                                              setup=setup, parameters=p)

# p.prune_mobility(-1)
results = COVIDVaccinationOCP.integrate(N,
                                        setup=setup,
                                        parameters=p,
                                        controls=control_initial,
                                        save_to='integ_0',
                                        n_rk4_steps=10)

ocp.update(parameters=p,
           max_total_vacc=np.inf,
           max_vacc_rate=max_vacc_rate,
           states_initial=matlab_initial,
           control_initial=control_initial,
           scenario_name=scenario_name)

ocp.solveOCP()
