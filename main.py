import numpy as np
from ItalySetup import ItalySetup
from covidOCP import COVIDVaccinationOCP, COVIDParametersOCP
import pickle
import matplotlib.pyplot as plt
import click
import sys, os

nx = 9
states_names = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']
outdir = 'model_output/'

os.makedirs(outdir, exist_ok=True)

# All arrays here are (nnodes, ndays, (nx))

nnodes = 107
ndays = 'full'
use_matlab = False
file_prefix = 'test'

if use_matlab:
    save_param = True

n_int_steps = 1

setup = ItalySetup(nnodes, ndays)
M = setup.nnodes
N = len(setup.model_days) - 1

if use_matlab:
    import matlab.engine

    eng = matlab.engine.start_matlab()
    eng.cd('geography-paper-master/', nargout=0)
    eng.run('single_build.m', nargout=0)

    p = COVIDParametersOCP.OCParameters(eng=eng, setup=setup, M=M)
    if save_param:
        with open(f'{outdir}parameters_{nnodes}.pkl', 'wb') as out:
            pickle.dump(p, out, pickle.HIGHEST_PROTOCOL)

else:
    with open(f'{outdir}parameters_{nnodes}.pkl', 'rb') as inp:
        p = pickle.load(inp)

control_initial = np.zeros((M, N))
max_vacc_rate = np.zeros((M, N))
for k in range(N):
    for nd in range(M):
        max_vacc_rate[nd, k] = 200
        control_initial[nd, k] = 0

initial = np.zeros((M, N + 1, nx))
for i, name in enumerate(states_names):
    for k in range(N + 1):
        for nd in range(M):
            initial[nd, k, i] = 0

p.prune_mobility(-1)
results, y, yell, mob = COVIDVaccinationOCP.integrate(N,
                                                      setup=setup,
                                                      parameters=p,
                                                      controls=control_initial,
                                                      save_to=f'{outdir}{file_prefix}-integ{nnodes}',
                                                      n_rk4_steps=10)
plt.figure(figsize = (10,10))
plt.step(np.arange(mob.T.shape[0]), mob.T)

sys.exit(0)

p.prune_mobility()
ocp = COVIDVaccinationOCP.COVIDVaccinationOCP(N=N, n_int_steps=n_int_steps,
                                              setup=setup, parameters=p)

p.prune_mobility(-1)

ocp.update(parameters=p,
           max_total_vacc=1000,
           max_vacc_rate=max_vacc_rate,
           states_initial=p.matlab_initial,
           control_initial=control_initial,
           scenario_name=f'{outdir}{file_prefix}-opt{nnodes}-10000')

ocp.solveOCP()

# @click.command()
# @click.option("-n", "--nnodes", "nnodes", default=10, envvar="OCP_NNODES", help="Spatial model size to run")
# @click.option("-t", "--ndays", "ndays", default='full', envvar="OCP_NDAYS", help="Number of days to run")
# @click.option("--use_matlab", "use_matlab", envvar="OCP_MATLAB", type=bool, default=False, show_default=True,
#              help="whether to use matlab for the current run")
# @click.option("-f", "--file_prefix", "file_prefix", envvar="OCP_PREFIX", type=str, default='test', show_default=True,
#              help="file prefix to add to identify the current set of runs.")
# def run_scenarios(nnodes, ndays, use_matlab, file_prefix):

# if __name__ == '__main__':
#    ocp, p = run_scenarios()
