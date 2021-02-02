import scipy.interpolate
import itertools
import datetime
import numpy as np


def pick_scenario(setup, scn_id):
    if setup.nnodes == 107:
        scenarios_specs = {
            # 'vacctotalM': [2, 5, 10, 15, 20],
            'newdoseperweek': [125000, 250000, 479700, 1e6],
            'vaccpermonthM': [15, 150],  # ax.set_ylim(0.05, 0.4)
            'epicourse': ['U', 'L']  # 'U'
        }
    elif setup.nnodes == 10:
        scenarios_specs = {
            'newdoseperweek': [125000/10, 250000/10],
            'vaccpermonthM': [1/10, 15/10],
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
    scenario = {'name': f"{scn_spec['epicourse']}-r{int(scn_spec['vaccpermonthM'])}-t{int(scn_spec['newdoseperweek'])}-id{scn_id}",
                'newdoseperweek': scn_spec['newdoseperweek'],
                'rate_fomula': f"({scn_spec['vaccpermonthM'] * 1e6 / tot_pop / 30}*pop_nd)"
                }
    # Build beta scenarios:
    if scn_spec['epicourse'] == 'C':
        scenario['beta_mult'] = np.ones((setup.nnodes, setup.ndays))
    elif scn_spec['epicourse'] == 'U':
        course = scipy.interpolate.interp1d([0, 50, 100, 1000], [1.2,.65, 1.2, 1], kind='linear')
        course = course(np.arange(0, setup.ndays))
        scenario['beta_mult'] = np.ones((setup.nnodes, setup.ndays)) * course
    elif scn_spec['epicourse'] == 'L':
        course = scipy.interpolate.interp1d([0, 50, 100, 1000], [.65, .8, .85, .75], kind='linear')
        course = course(np.arange(0, setup.ndays))
        scenario['beta_mult'] = np.ones((setup.nnodes, setup.ndays)) * course
    return scenario


def build_scenario(setup, scenario, strategy = None):
    M = setup.nnodes
    N = setup.ndays - 1
    control_initial = np.zeros((M, N))
    maxvaccrate_regional = np.zeros((M, N))
    unvac_nd = np.copy(setup.pop_node)
    delivery_national = np.zeros(N)

    if strategy is None:
        strategy = setup.pop_node

    stockpile = 0
    for k in range(N):
        if (setup.start_date + datetime.timedelta(days=k)).weekday() == 0:  # if monday
            delivery_national[k] = scenario['newdoseperweek']
        else:
            delivery_national[k] = 0
        stockpile += delivery_national[k]
        for nd in range(M):
            pop_nd = setup.pop_node[nd]
            maxvaccrate_regional[nd, k] = eval(scenario['rate_fomula'])
            to_allocate = stockpile * strategy[nd] / strategy.sum()
            to_allocate = min(to_allocate, maxvaccrate_regional[nd, k]*.9, unvac_nd[nd]*.9)
            control_initial[nd, k] = to_allocate
            stockpile -= to_allocate
            unvac_nd[nd] -= to_allocate

    stockpile_national_constraint = np.cumsum(delivery_national)
    for k in range(N):
        if (setup.start_date + datetime.timedelta(days=k)).weekday() != 6:  # if NOt sunday
            stockpile_national_constraint[k] = np.inf

    if stockpile_national_constraint[-1] == np.inf:
        # np.nanmax(stockpile_national_constraint[stockpile_national_constraint != np.inf]) + scenario['newdoseperweek']
        stockpile_national_constraint[-1] = np.cumsum(delivery_national).max()

    return maxvaccrate_regional, delivery_national, stockpile_national_constraint, control_initial
