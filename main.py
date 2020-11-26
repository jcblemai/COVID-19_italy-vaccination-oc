import copy
import casadi as ca
import matlab.engine
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx
from ItalySetup import ItalySetup
from covidOCP import COVIDVaccinationOCP, COVIDParametersOCP

# Check if convariates in N or T

nx = 9
states_names = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']

eng = matlab.engine.start_matlab()
eng.cd('geography-paper-master/', nargout=0)

model_size = 10  # nodes

# Horizon for each problem
ndays = 30
n_int_steps = 1

setup = ItalySetup(model_size, ndays)
M = setup.nnodes
N = len(setup.model_days) - 1


eng.run('single_build.m', nargout=0)
integ_matlab = np.array(eng.eval('x'))

matlab_initial = np.zeros((M, N + 1, nx))
for i, name in enumerate(states_names):
    for k in range(N + 1):
        for nd in range(M):
            matlab_initial[nd, k, i] = integ_matlab.T[nd + 107 * i, k].T

p = COVIDParametersOCP.OCParameters(eng=eng, setup=setup, M=M)

ocp = COVIDVaccinationOCP.COVIDVaccinationOCP(N=N, n_int_steps=n_int_steps,
                                              setup=setup, parameters=p)

scenario_name = 'test10000'
control_initial = np.zeros((M, N))
max_vacc_rate = np.zeros((M, N))
for k in range(N):
    for nd in range(M):
        max_vacc_rate[nd, k] = 2000
        control_initial[nd, k] = 0

ocp.update(parameters=p,
           max_total_vacc=10000,
           max_vacc_rate=max_vacc_rate,
           states_initial=matlab_initial,
           control_initial=control_initial,
           scenario_name=scenario_name)

ocp.solveOCP()
