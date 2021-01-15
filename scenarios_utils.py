import scipy.interpolate
import itertools
import datetime
import numpy as np


def pick_scenario(setup, scn_id):
    if setup.nnodes == 107:
        scenarios_specs = {
            # 'vacctotalM': [2, 5, 10, 15, 20],
            'newdoseperweek': [125000, 250000, 479700, 1e6, 2e6],
            'vaccpermonthM': [1, 2, 15],  # ax.set_ylim(0.05, 0.4)
            'epicourse': ['U', 'L']  # 'U'
        }
    elif setup.nnodes == 10:
        scenarios_specs = {
            'newdoseperweek': [125000/10, 250000/10],
            'vaccpermonthM': [125000/5, 1],
            'epicourse': ['U']  # 'U', 'L'
        }

    # Compute all permutatios
    keys, values = zip(*scenarios_specs.items())
    permuted_specs = [dict(zip(keys, v)) for v in itertools.product(*values)]

    scn_spec = permuted_specs[scn_id]
    # check if the scenario is useless:
    # if scn_spec['vaccpermonthM']*setup.ndays/30 < scn_spec['vacctotalM']:
    #    raise ValueError("Scenario is useless")

    tot_pop = setup.pop_node.sum()
    scenario = {'name': f"{scn_spec['epicourse']}-r{scn_spec['vaccpermonthM']}-t{int(scn_spec['newdoseperweek'])}",
                'newdoseperweek': scn_spec['newdoseperweek'],
                'rate_fomula': f"({scn_spec['vaccpermonthM'] * 1e6 / tot_pop / 30}*pop_nd)"
                }
    # Build beta scenarios:
    if scn_spec['epicourse'] == 'C':
        scenario['beta_mult'] = np.ones((setup.nnodes, setup.ndays))
    elif scn_spec['epicourse'] == 'U':
        course = scipy.interpolate.interp1d([0, 50, 100, 1000], [1.3, .7, 1.2, 1], kind='linear')
        course = course(np.arange(0, setup.ndays))
        scenario['beta_mult'] = np.ones((setup.nnodes, setup.ndays)) * course
    elif scn_spec['epicourse'] == 'L':
        course = scipy.interpolate.interp1d([0, 50, 100, 1000], [.75, .85, .45, .45], kind='linear')
        course = course(np.arange(0, setup.ndays))
        scenario['beta_mult'] = np.ones((setup.nnodes, setup.ndays)) * course
    return scenario


def build_scenario(setup, scenario):
    M = setup.nnodes
    N = setup.ndays - 1
    control_initial = np.zeros((M, N))
    maxvaccrate_regional = np.zeros((M, N))
    unvac_nd = np.copy(setup.pop_node)
    stockpile_national = np.zeros(N)

    stockpile = 0
    for k in range(N):
        if (setup.start_date + datetime.timedelta(days=k)).weekday() == 0:  # if monday
            stockpile_national[k] = scenario['newdoseperweek']
        else:
            stockpile_national[k] = 0
        stockpile += stockpile_national[k]
        for nd in range(M):
            pop_nd = setup.pop_node[nd]
            maxvaccrate_regional[nd, k] = eval(scenario['rate_fomula'])
            to_allocate = stockpile * pop_nd / setup.pop_node.sum()
            to_allocate = min(to_allocate, maxvaccrate_regional[nd, k], unvac_nd[nd] - 100)
            control_initial[nd, k] = to_allocate
            stockpile -= to_allocate
            unvac_nd[nd] -= to_allocate

    stockpile_national = np.cumsum(stockpile_national)
    stockpile_national_constraint = np.copy(stockpile_national)
    for k in range(N):
        if (setup.start_date + datetime.timedelta(days=k)).weekday() != 0:  # if NOt monday:
            stockpile_national_constraint[k] = np.inf

    return maxvaccrate_regional, stockpile_national, stockpile_national_constraint, control_initial
