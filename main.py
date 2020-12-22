import numpy as np
from ItalySetup import ItalySetup
from covidOCP import COVIDVaccinationOCP as COVIDVaccinationOCP
from covidOCP import COVIDParametersOCP
import pickle
import matplotlib.pyplot as plt
import click
import sys, os
import itertools

nx = 9
states_names = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']
outdir = 'model_output/'
when = 'future'
ocp = True
n_int_steps = 5


@click.command()
@click.option("-s", "--scenario_id", "scn_id", default=0, help="Index of scenario to run")
@click.option("-n", "--nnodes", "nnodes", default=107, envvar="OCP_NNODES", help="Spatial model size to run")
@click.option("-t", "--ndays", "ndays", default='full', envvar="OCP_NDAYS", help="Number of days to run")
@click.option("--use_matlab", "use_matlab", envvar="OCP_MATLAB", type=bool, default=True, show_default=True,
              help="whether to use matlab for the current run")
@click.option("-f", "--file_prefix", "file_prefix", envvar="OCP_PREFIX", type=str, default='',
              show_default=True,
              help="file prefix to add to identify the current set of runs.")
def cli(scn_id, nnodes, ndays, use_matlab, file_prefix):
    return scn_id, nnodes, ndays, use_matlab, file_prefix


def pick_scenario(setup, scn_id):
    scenarios_specs = {
        'vaccpermonthM': [1, 2.5, 7.5],
        'vacctotalM': [10, 15, 20]
    }

    # Compute all permutatios
    keys, values = zip(*scenarios_specs.items())
    permuted_specs = [dict(zip(keys, v)) for v in itertools.product(*values)]

    scn_spec = permuted_specs[scn_id]

    tot_pop = setup.pop_node.sum()
    scenario = {'name': f"FR{scn_spec['vaccpermonthM']}-T{scn_spec['vacctotalM']}",
                'vacctotal': scn_spec['vacctotalM']*1e6,
                'rate_fomula': f"({scn_spec['vaccpermonthM'] / tot_pop / 30}*pop_nd)"
    }

    return scenario


if __name__ == '__main__':
    # standalone_mode: so click doesn't exit, see
    # https://stackoverflow.com/questions/60319832/how-to-continue-execution-of-python-script-after-evaluating-a-click-cli-function
    scn_id, nnodes, ndays, use_matlab, file_prefix = cli(standalone_mode=False)

    os.makedirs(outdir, exist_ok=True)
    ndays = 31

    # All arrays here are (nnodes, ndays, (nx))
    if use_matlab:
        save_param = True

    setup = ItalySetup(nnodes, ndays, when)
    M = setup.nnodes
    N = len(setup.model_days) - 1

    scenario = pick_scenario(setup, scn_id)
    file_prefix = file_prefix + scenario['name']

    print(f"""Running scenario {scn_id}: {scenario['name']}, building setup with
    ndays: {ndays}
    nnodes: {nnodes}
    use_matlab: {use_matlab}
    when?  {when}
    rk_steps: {n_int_steps}
    ---> Saving results to prefix: {file_prefix}""")

    if use_matlab:
        import matlab.engine
        eng = matlab.engine.start_matlab()
        if when == 'past':
            eng.cd('geography-paper-master/', nargout=0)
            eng.run('single_build.m', nargout=0)
        if when == 'future':
            eng.cd('data-assimilation/', nargout=0)
            # The realization with the maximum infected at the end of the 3 months is realization 33.
            # The realization with the median number of infected at the end of the 3 months is realization 24.
            eng.workspace['i'] = 24
            eng.run('minimal_interface.m', nargout=0)

        p = COVIDParametersOCP.OCParameters(eng=eng, setup=setup, M=M, when=when)

        if save_param:
            with open(f'{outdir}parameters_{nnodes}_{when}.pkl', 'wb') as out:
                pickle.dump(p, out, pickle.HIGHEST_PROTOCOL)
    else:
        with open(f'{outdir}parameters_{nnodes}_{when}.pkl', 'rb') as inp:
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

    #scn_maxvacc = [3000 * 107 * N * 2 / 3]

    # scn_maxvacc = [m*(nnodes/107)*(ndays/160) for m in scn_maxvacc]
    # scn_maxvacc = [int(m * (nnodes / 107)) for m in scn_maxvacc]

    mvr = 3000

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
