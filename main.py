import numpy as np
from ItalySetup import ItalySetup
from covidOCP import COVIDVaccinationOCP as COVIDVaccinationOCP
from covidOCP import COVIDParametersOCP
import pickle
import matplotlib.pyplot as plt
import click
import sys, os

nx = 9
states_names = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']
outdir = 'model_output/'


@click.command()
@click.option("-s", "--scenario", "scenario", default='all', help="Index of scenario to run")
@click.option("-n", "--nnodes", "nnodes", default=107, envvar="OCP_NNODES", help="Spatial model size to run")
@click.option("-t", "--ndays", "ndays", default='full', envvar="OCP_NDAYS", help="Number of days to run")
@click.option("--use_matlab", "use_matlab", envvar="OCP_MATLAB", type=bool, default=False, show_default=True,
              help="whether to use matlab for the current run")
@click.option("-f", "--file_prefix", "file_prefix", envvar="OCP_PREFIX", type=str, default='test', show_default=True,
              help="file prefix to add to identify the current set of runs.")
def cli(scenario, nnodes, ndays, use_matlab, file_prefix):
    return scenario, nnodes, ndays, use_matlab, file_prefix


if __name__ == '__main__':
    scenario, nnodes, ndays, use_matlab, file_prefix = cli()
    os.makedirs(outdir, exist_ok=True)

    scenario_specifications = {'ndays': [60, 90, 120, 'full']}

    # All arrays here are (nnodes, ndays, (nx))7
    ocp = True

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
            max_vacc_rate[nd, k] = 0
            control_initial[nd, k] = 0

    initial = np.zeros((M, N + 1, nx))
    for i, name in enumerate(states_names):
        for k in range(N + 1):
            for nd in range(M):
                initial[nd, k, i] = 0

    results, state_initial, yell, mob = COVIDVaccinationOCP.integrate(N,
                                                                      setup=setup,
                                                                      parameters=p,
                                                                      controls=control_initial,
                                                                      save_to=f'{outdir}{file_prefix}-int{nnodes}-r0-m0',
                                                                      method='rk4',
                                                                      n_rk4_steps=n_int_steps)

    if ocp:
        ocp = COVIDVaccinationOCP.COVIDVaccinationOCP(N=N, n_int_steps=n_int_steps,
                                                      setup=setup, parameters=p,
                                                      show_steps=False)

        ocp.update(parameters=p,
                   max_total_vacc=1e6,
                   max_vacc_rate=max_vacc_rate,
                   states_initial=state_initial,  # p.matlab_initial,#
                   control_initial=control_initial,
                   mob_initial=mob,
                   scenario_name=f'{outdir}{file_prefix}-opt{nnodes}-r0-m0')

        # ocp.solveOCP()
    # plt.figure(figsize=(10, 10))
    # plt.step(np.arange(mob.T.shape[0]), mob.T)
    # plt.show()

    scn_maxvacc = [1e6, 4e6, 8e6, 12e6, 16e6, 20e6]

    scn_maxvacc = [16e6]

    # scn_maxvacc = [m*(nnodes/107)*(ndays/160) for m in scn_maxvacc]
    #scn_maxvacc = [int(m * (nnodes / 107)) for m in scn_maxvacc]

    mvr = 8000

    for scn_id, scn_maxvacc in enumerate(scn_maxvacc):

        control_initial = np.zeros((M, N))
        max_vacc_rate = np.zeros((M, N))
        allocated_total = 0
        unvac_nd = np.copy(setup.pop_node)

        for k in range(N):
            for nd in range(M):
                max_vacc_rate[nd, k] = mvr
                if (allocated_total + mvr < scn_maxvacc) and (unvac_nd[nd] - mvr > 0):
                    control_initial[nd, k] = mvr
                    allocated_total += mvr
                    unvac_nd[nd] -= mvr

        results, state_initial, yell, mob = COVIDVaccinationOCP.integrate(N,
                                                                          setup=setup,
                                                                          parameters=p,
                                                                          controls=control_initial,
                                                                          save_to=f'{outdir}{file_prefix}-int{nnodes}-r{mvr}-m{int(scn_maxvacc)}',
                                                                          n_rk4_steps=n_int_steps)
        if ocp:
            ocp.update(parameters=p,
                       max_total_vacc=scn_maxvacc,
                       max_vacc_rate=max_vacc_rate,
                       states_initial=state_initial,
                       control_initial=control_initial,
                       mob_initial=mob,
                       scenario_name=f'{outdir}{file_prefix}-opt{nnodes}-r{mvr}-m{int(scn_maxvacc)}')

            ocp.solveOCP()