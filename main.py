import numpy as np
from ItalySetup import ItalySetup
from covidOCP import COVIDVaccinationOCP as COVIDVaccinationOCP
from covidOCP import COVIDAgeStructuredOCP as COVIDAgeStructuredOCP
from covidOCP import COVIDParametersOCP
import pickle
import matplotlib.pyplot as plt
import click
import sys, os
from scenarios_utils import pick_scenario, build_scenario

nx = 9
states_names = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']
outdir = 'model_output/'
os.makedirs(outdir, exist_ok=True)
when = 'future'
optimize = True
n_int_steps = 6
ocp = None
nc = 1


@click.command()
@click.option("-s", "--scenario_id", "scn_ids", default=5, help="Index of scenario to run")
@click.option("-n", "--nnodes", "nnodes", default=107, envvar="OCP_NNODES", help="Spatial model size to run")
@click.option("-t", "--ndays", "ndays", default=15, envvar="OCP_NDAYS", help="Number of days to run")
@click.option("--use_matlab", "use_matlab", envvar="OCP_MATLAB", type=bool, default=False, show_default=True,
              help="whether to use matlab for the current run")
@click.option("-a", "--age_struct", "age_struct", type=bool, default=False, show_default=True,
              help="Whether to use agestructured OCP")
@click.option("-f", "--file_prefix", "file_prefix", envvar="OCP_PREFIX", type=str, default='',
              show_default=True, help="file prefix to add to identify the current set of runs.")
def cli(scn_ids, nnodes, ndays, use_matlab, age_struct, file_prefix):
    if not isinstance(scn_ids, list):
        scn_ids = [int(scn_ids)]
    return scn_ids, nnodes, ndays, use_matlab, age_struct, file_prefix


if __name__ == '__main__':
    # standalone_mode: so click doesn't exit, see
    # https://stackoverflow.com/questions/60319832/how-to-continue-execution-of-python-script-after-evaluating-a-click-cli-function
    scn_ids, nnodes, ndays, use_matlab, age_struct, file_prefix = cli(standalone_mode=False)
    # scn_ids = np.arange(18)

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
        initial = np.zeros((M, N + 1, nx))

        results, state_initial, yell, mob = COVIDVaccinationOCP.integrate(N,
                                                                          setup=setup,
                                                                          parameters=p,
                                                                          controls=control_initial,
                                                                          save_to=f'{outdir}{prefix}-int{nnodes}_{ndays}-nc',
                                                                          method='rk4',
                                                                          n_rk4_steps=n_int_steps)

        if optimize and ocp is None:
            ocp = COVIDVaccinationOCP.COVIDVaccinationOCP(N=N, n_int_steps=n_int_steps,
                                                          setup=setup, parameters=p,
                                                          show_steps=False)

        max_vacc_rate, vacc_total, control_initial = build_scenario(setup, scenario)

        results, state_initial, yell, mob = COVIDVaccinationOCP.integrate(N,
                                                                          setup=setup,
                                                                          parameters=p,
                                                                          controls=control_initial,
                                                                          save_to=f'{outdir}{prefix}-int{nnodes}_{ndays}',
                                                                          n_rk4_steps=n_int_steps)

        if optimize:
            ocp.update(parameters=p,
                       max_total_vacc=vacc_total,
                       max_vacc_rate=max_vacc_rate,
                       states_initial=state_initial,
                       control_initial=control_initial,
                       mob_initial=mob,
                       scenario_name=f'{outdir}{prefix}-opt{nnodes}_{ndays}')

            ocp.solveOCP()