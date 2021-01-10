import scipy.interpolate
import itertools
import numpy as np


def pick_scenario(setup, scn_id):
    if setup.nnodes == 107:
        scenarios_specs = {
            #'vacctotalM': [2, 5, 10, 15, 20],
            'vacctotalperweek': [250000, 479700, 1e6],
            'vaccpermonthM': [1, 15], # ax.set_ylim(0.05, 0.4)
            'epicourse': ['U', 'L']  # 'U'
        }
    elif setup.nnodes == 10:
        scenarios_specs = {
            'vacctotalM': [.25, .5, 1.5],
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
    scenario = {'name': f"mr-{scn_spec['epicourse']}-R{scn_spec['vaccpermonthM']}-T{scn_spec['vacctotalperweek']}",
                'vacctotalday': scn_spec['vacctotalperweek']/7,
                'rate_fomula': f"({scn_spec['vaccpermonthM'] * 1e6 / tot_pop / 30}*pop_nd)"
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


def build_scenario(setup, scenario):
    M = setup.nnodes
    N = setup.ndays - 1
    control_initial = np.zeros((M, N))
    max_vacc_rate = np.zeros((M, N))
    allocated_total = 0
    unvac_nd = np.copy(setup.pop_node)
    vacc_total = np.zeros(N)

    for k in range(N):
        vacc_total[k] = scenario['vacctotalday']
        for nd in range(M):
            pop_nd = setup.pop_node[nd]
            max_vacc_rate[nd, k] = eval(scenario['rate_fomula'])
            if (allocated_total + max_vacc_rate[nd, k] < scenario['vacctotal']) and (
                    unvac_nd[nd] - max_vacc_rate[nd, k] > 0):
                control_initial[nd, k] = max_vacc_rate[nd, k]
                allocated_total += max_vacc_rate[nd, k]
                unvac_nd[nd] -= max_vacc_rate[nd, k]

    vacc_total = np.cumsum(vacc_total)

    return max_vacc_rate, vacc_total, control_initial
