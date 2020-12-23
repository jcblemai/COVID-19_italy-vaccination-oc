import numpy as np
from ItalySetup import ItalySetup
from covidOCP import COVIDVaccinationOCP as COVIDVaccinationOCP
from covidOCP import COVIDParametersOCP
import pickle
import matplotlib.pyplot as plt
import click
import sys, os
import scipy.interpolate
import itertools

nx = 9
states_names = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']
outdir = 'model_output/'
when = 'future'
optimize = True
n_int_steps = 5
ocp = None


@click.command()
@click.option("-s", "--scenario_id", "scn_ids", default=[0, 1], help="Index of scenario to run")
@click.option("-n", "--nnodes", "nnodes", default=10, envvar="OCP_NNODES", help="Spatial model size to run")
@click.option("-t", "--ndays", "ndays", default=90, envvar="OCP_NDAYS", help="Number of days to run")
@click.option("--use_matlab", "use_matlab", envvar="OCP_MATLAB", type=bool, default=False, show_default=True,
              help="whether to use matlab for the current run")
@click.option("-f", "--file_prefix", "file_prefix", envvar="OCP_PREFIX", type=str, default='',
              show_default=True,
              help="file prefix to add to identify the current set of runs.")
def cli(scn_ids, nnodes, ndays, use_matlab, file_prefix):
    if not isinstance(scn_ids, list):
        scn_ids = [scn_ids]
    return scn_ids, nnodes, ndays, use_matlab, file_prefix


def pick_scenario(setup, scn_id):
    scenarios_specs = {
        'vacctotalM': [10, 15, 20],
        'vaccpermonthM': [1, 2.5, 7.5],
        'epicourse': ['U', 'L']  # 'U'
    }

    # Compute all permutatios
    keys, values = zip(*scenarios_specs.items())
    permuted_specs = [dict(zip(keys, v)) for v in itertools.product(*values)]

    scn_spec = permuted_specs[scn_id]

    tot_pop = setup.pop_node.sum()
    scenario = {'name': f"FR-{scn_spec['epicourse']}-R{scn_spec['vaccpermonthM']}-T{scn_spec['vacctotalM']}",
                'vacctotal': scn_spec['vacctotalM'] * 1e6,
                'rate_fomula': f"({scn_spec['vaccpermonthM']*1e6 / tot_pop / 30}*pop_nd)"
                }
    # Build beta scenarios:
    if scn_spec['epicourse'] == 'C':
        scenario['beta_mult'] = np.ones((setup.nnodes, setup.ndays))
    elif scn_spec['epicourse'] == 'U':
        course = scipy.interpolate.interp1d([0, setup.ndays / 2, setup.ndays], [2, .4, 2], kind='quadratic')
        course = course(np.arange(0, setup.ndays))
        scenario['beta_mult'] = np.ones((setup.nnodes, setup.ndays)) * course
    elif scn_spec['epicourse'] == 'L':
        course = scipy.interpolate.interp1d([0, setup.ndays / 2, setup.ndays], [.45, .45, .45], kind='quadratic')
        course = course(np.arange(0, setup.ndays))
        scenario['beta_mult'] = np.ones((setup.nnodes, setup.ndays)) * course
    return scenario


if __name__ == '__main__':
    # standalone_mode: so click doesn't exit, see
    # https://stackoverflow.com/questions/60319832/how-to-continue-execution-of-python-script-after-evaluating-a-click-cli-function
    scn_ids, nnodes, ndays, use_matlab, file_prefix = cli(standalone_mode=False)

    os.makedirs(outdir, exist_ok=True)

    # All arrays here are (nnodes, ndays, (nx))
    setup = ItalySetup(nnodes, ndays, when)
    M = setup.nnodes
    N = setup.ndays - 1

    if use_matlab:
        p = COVIDParametersOCP.OCParameters(setup=setup, M=M, when=when)
        if True:
            with open(f'{outdir}parameters_{nnodes}_{when}.pkl', 'wb') as out:
                pickle.dump(p, out, pickle.HIGHEST_PROTOCOL)
    else:
        with open(f'{outdir}parameters_{nnodes}_{when}.pkl', 'rb') as inp:
            p = pickle.load(inp)

    for scn_id in scn_ids:
        scenario = pick_scenario(setup, scn_id)
        prefix = file_prefix + scenario['name']

        print(f"""Running scenario {scn_id}: {scenario['name']}, building setup with
        ndays: {ndays}
        nnodes: {nnodes}
        use_matlab: {use_matlab}
        when?  {when}
        rk_steps: {n_int_steps}
        ---> Saving results to prefix: {prefix}""")

        p.apply_epicourse(setup, scenario['beta_mult'])

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
                                                                          save_to=f'{outdir}{prefix}-int{nnodes}-nc',
                                                                          method='rk4',
                                                                          n_rk4_steps=n_int_steps)

        if optimize and ocp is None:
            ocp = COVIDVaccinationOCP.COVIDVaccinationOCP(N=N, n_int_steps=n_int_steps,
                                                          setup=setup, parameters=p,
                                                          show_steps=False)

        control_initial = np.zeros((M, N))
        max_vacc_rate = np.zeros((M, N))
        allocated_total = 0
        unvac_nd = np.copy(setup.pop_node)

        for k in range(N):
            for nd in range(M):
                pop_nd = setup.pop_node[nd]
                max_vacc_rate[nd, k] = eval(scenario['rate_fomula'])
                if (allocated_total + max_vacc_rate[nd, k] < scenario['vacctotal']) and (unvac_nd[nd] - max_vacc_rate[nd, k] > 0):
                    control_initial[nd, k] = max_vacc_rate[nd, k]
                    allocated_total += max_vacc_rate[nd, k]
                    unvac_nd[nd] -= max_vacc_rate[nd, k]

        results, state_initial, yell, mob = COVIDVaccinationOCP.integrate(N,
                                                                          setup=setup,
                                                                          parameters=p,
                                                                          controls=control_initial,
                                                                          save_to=f'{outdir}{prefix}-int{nnodes}',
                                                                          n_rk4_steps=n_int_steps)
        if optimize:
            ocp.update(parameters=p,
                       max_total_vacc=scenario['vacctotal'],
                       max_vacc_rate=max_vacc_rate,
                       states_initial=state_initial,
                       control_initial=control_initial,
                       mob_initial=mob,
                       scenario_name=f'{outdir}{prefix}-opt{nnodes}')

            ocp.solveOCP()
