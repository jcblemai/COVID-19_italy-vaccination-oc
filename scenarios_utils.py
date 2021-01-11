import scipy.interpolate
import itertools
import numpy as np


def pick_scenario(setup, scn_id):
    if setup.nnodes == 107:
        scenarios_specs = {
            #'vacctotalM': [2, 5, 10, 15, 20],
            'newdoseperweek': [250000, 479700, 1e6, 2e6],
            'vaccpermonthM': [1, 2, 15], # ax.set_ylim(0.05, 0.4)
            'epicourse': ['U', 'L']  # 'U'
        }
    elif setup.nnodes == 10:
        scenarios_specs = {
            'newdoseperweek': [10000, 50000],
            'vaccpermonthM': [1],
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
    scenario = {'name': f"mr-{scn_spec['epicourse']}-L{scn_spec['vaccpermonthM']}-S{scn_spec['newdoseperweek']}",
                'newdoseperday': scn_spec['newdoseperweek']/7,
                'rate_fomula': f"({scn_spec['vaccpermonthM'] * 1e6 / tot_pop / 30}*pop_nd)"
                }
    # Build beta scenarios:
    if scn_spec['epicourse'] == 'C':
        scenario['beta_mult'] = np.ones((setup.nnodes, setup.ndays))
    elif scn_spec['epicourse'] == 'U':
        course = scipy.interpolate.interp1d([0, 50, 100], [1.3, .7, 2], kind='quadratic')
        course = course(np.arange(0, setup.ndays))
        scenario['beta_mult'] = np.ones((setup.nnodes, setup.ndays)) * course
    elif scn_spec['epicourse'] == 'L':
        course = scipy.interpolate.interp1d([0, 50, 100], [.75, .85, .45], kind='quadratic')
        course = course(np.arange(0, setup.ndays))
        scenario['beta_mult'] = np.ones((setup.nnodes, setup.ndays)) * course
    return scenario


def build_scenario(setup, scenario):
    M = setup.nnodes
    N = setup.ndays - 1
    control_initial = np.zeros((M, N))
    max_vacc_rate = np.zeros((M, N))
    unvac_nd = np.copy(setup.pop_node)
    vacc_total = np.zeros(N)

    stockpile = 0
    for k in range(N):
        vacc_total[k] = scenario['newdoseperday']
        stockpile += vacc_total[k]
        for nd in range(M):
            pop_nd = setup.pop_node[nd]
            max_vacc_rate[nd, k] = eval(scenario['rate_fomula'])
            to_allocate = stockpile * pop_nd / setup.pop_node.sum()
            to_allocate = min(to_allocate, max_vacc_rate[nd, k], unvac_nd[nd] - 100)
            control_initial[nd, k] = to_allocate
            stockpile -= to_allocate
            unvac_nd[nd] -= to_allocate

    vacc_total = np.cumsum(vacc_total)

    return max_vacc_rate, vacc_total, control_initial
