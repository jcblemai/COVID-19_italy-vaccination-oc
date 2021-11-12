import copy
import numpy as np
from ItalySetup import ItalySetupProvinces
from covidOCP import COVIDVaccinationOCP as COVIDVaccinationOCP
from covidOCP import COVIDAgeStructuredOCP as COVIDAgeStructuredOCP
from covidOCP import COVIDParametersOCP
import pickle
import matplotlib.pyplot as plt
import click
import time
import sys, os
from scenarios_utils import pick_scenario, build_scenario
import pandas as pd
import multiprocessing as mp
import tqdm

# Replace the jupyter notebook that was based on matlab called generate_all_scn.

nx = 9
states_names = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']
when = 'future-mobintime'

# to load the optimal strategy
input_directory = 'helvetios-runs/2021-11-05-107_90'
input_prefix = f'week'
# to output the now files
output_directory = 'model_output/2021-11-09'
output_prefix = f'altstratint'

nnodes = 107  # nodes
ndays_ocp = 90
ndays = 90

#setup_ocp = ItalySetupProvinces(nnodes, ndays_ocp, when)

os.makedirs(f'{output_directory}', exist_ok=True)

class AlternativeStrategy:
    def __init__(self, setup, scenario, decision_variable, alloc_function=None, dv_per_pop=False, require_projection=False, alloc_arr=None):
        self.maxvaccrate_regional, self.delivery_national, self.stockpile_national_constraint, _ = build_scenario(setup, scenario)
        self.M = setup.nnodes
        self.pop_node = setup.pop_node
        self.ind2name = setup.ind2name  # so we get rid of setup when pickled
        self.ndays = setup.ndays
        self.maxvaccrate_regional = self.maxvaccrate_regional[:,0] # stays the same over the course, so take the first one

        # updated states variables
        self.unvaccinated = np.copy(setup.pop_node)
        self.stockpile = 0

        self.decision_variable = decision_variable
        self.require_projection = require_projection

        # The decision variable is per habitant
        self.name = decision_variable.capitalize()
        if self.decision_variable == 'novacc':
            self.name = 'Baseline'
        self.shortname = decision_variable[:3]

        self.dv_per_pop = dv_per_pop
        self.divider = np.ones(self.M)
        if dv_per_pop:
            self.divider = self.pop_node
            self.shortname += '_pp'
            self.name += ' per hab.'

        if alloc_function == 'focused':
            self.alloc_function = self.focused_alloc
            self.name += ' (focused)'
            self.shortname += '_f'
        elif alloc_function == 'proportional':
            self.alloc_function = self.proportional_alloc
            self.name += ' (proportional)'
            self.shortname += '_p'

        self.compute_new_strat = True

        if alloc_arr is not None:
            self.alloc_arr = alloc_arr
            self.compute_new_strat = False
        else:
            self.alloc_arr = np.ones((self.M, self.ndays-1)) * -1 # to be filled

    def focused_alloc(self, decision_df_sorted, nd, nodename):
        return self.maxvaccrate_regional[nd]

    # this can be vectorized...
    def proportional_alloc(self, decision_df_sorted, nd, nodename):
        return self.stockpile * decision_df_sorted.loc[nodename]['value'] / decision_df_sorted['value'].sum()

    def allocate_now(self, decision_variable_array, today_idx):
        # Sort the decision variable dataframe:
        self.stockpile += self.delivery_national[today_idx]
        #optimize when already allocated
        if self.stockpile <= 1:
            return np.zeros(self.M)

        decision_variable_df = pd.DataFrame(decision_variable_array, index=self.ind2name, columns=['value'])
        decision_variable_df.sort_values('value', ascending=False, inplace=True)
        alloc_now = np.zeros(self.M)
        for nodename in decision_variable_df.index:
            nd = self.ind2name.index(nodename)
            to_allocate = self.alloc_function(decision_variable_df, nd, nodename)
            to_allocate = min(to_allocate, self.unvaccinated[nd], self.stockpile, self.maxvaccrate_regional[nd])
            alloc_now[nd] = to_allocate
            self.stockpile        -= to_allocate
            self.unvaccinated[nd] -= to_allocate
            if self.stockpile <= 1:
                return alloc_now
        return alloc_now

    def get_today_allocation(self, today_idx, susceptible=None, incidence=None):
        if self.compute_new_strat:
            self.alloc_arr[:, today_idx] = self.compute_today_allocation(today_idx, susceptible, incidence)

        # return from memory
        return self.alloc_today_from_memory(today_idx)

    def alloc_today_from_memory(self, today_idx):
        return self.alloc_arr[:, today_idx]

    def compute_today_allocation(self, today_idx, susceptible, incidence):
        if 'susceptible' in self.decision_variable:
            return self.allocate_now(susceptible/self.divider, today_idx)
        elif 'incidence' in self.decision_variable:
            return self.allocate_now(incidence/self.divider, today_idx)
        elif 'population' in self.decision_variable:
            return self.allocate_now(self.pop_node/self.divider, today_idx)
        elif 'novacc' in self.decision_variable:
            return np.zeros(self.M)
        elif 'optimal' in self.decision_variable:
            raise ValueError('No you cannot compute today_allocation for optimal, only from memory')
        else:
            raise ValueError(f'impossible to compute allocation from {self.decision_variable}')


def create_all_alt_strategies(setup, scenario_name, scenario):
    # create scenarios
    decisions_variables = ['susceptible', 'population', 'incidence']
    alt_strategies = {}
    for decision_variable in decisions_variables:
        require_projection = False
        if decision_variable == 'incidence':
            require_projection = True
        for alloc_function in ['focused', 'proportional']:
            for dv_per_pop in [True, False]:
                alt_strat = AlternativeStrategy(setup,
                                                scenario,
                                                decision_variable,
                                                alloc_function,
                                                dv_per_pop,
                                                require_projection)
                alt_strategies[alt_strat.shortname] = alt_strat
    # get the optimal strategy
    fname = f"{input_directory}/{input_prefix}-{scenario_name}-opt-{nnodes}_{ndays_ocp}.csv"
    optimal_df = pd.read_csv(fname, index_col='date', parse_dates=True)
    optimal_alloc = optimal_df[optimal_df['comp'] == 'vacc'][['value', 'placeID']].pivot(columns='placeID', values='value').T
    optimal_alloc_array = optimal_alloc.sort_index().to_numpy()
    alt_strat = AlternativeStrategy(setup,
                                    scenario,
                                    'optimal',
                                    alloc_arr=optimal_alloc_array)
    alt_strategies[alt_strat.shortname] = alt_strat

    alt_strat = AlternativeStrategy(setup,
                                    scenario,
                                    'novacc')

    alt_strategies[alt_strat.shortname] = alt_strat

    print(f'generated {len(alt_strategies.keys())} strategies: {list(alt_strategies.keys())} for scenario {scenario_name}')
    return alt_strategies

def worker_one_posterior_realization(post_real, scenario_name, scenario, alt_strategies):
    tic1 = time.time()

    # create object here so not shared:
    with open(f'italy-data/full_posterior/setup_{nnodes}_{when}.pkl', 'rb') as inp:
        setup = pickle.load(inp)
    #setup = ItalySetupProvinces(nnodes, ndays, when)

    print(f"{scenario_name}, {post_real}")

    with open(f'italy-data/full_posterior/parameters_{nnodes}_{when}_{post_real}.pkl', 'rb') as inp:
        p = pickle.load(inp)
    p.apply_epicourse(setup, scenario['beta_mult'])

    all_results = pd.DataFrame(columns=['method_short', 'method', 'infected', 'post_sample', 'doses', 'scenario-beta', 'scenario-rate', 'scenario-tot', 'scenario', 'newdoseperweek'])

    for shortname, strat in alt_strategies.items():
        tic = time.time()
        results, state_initial, yell, = COVIDVaccinationOCP.accurate_integrate(len(setup.model_days) - 1,
                                                                               setup=setup,
                                                                               parameters=p,
                                                                               controls=None,
                                                                               save_to=None,#f'{output_directory}/{output_prefix}-{scenario_name}-{shortname}-{post_real}',
                                                                               only_yell=True,
                                                                               alloc_strat=strat)
        yell_tot = results[results['comp'] == 'yell'].pivot(values='value', columns='place', index='date').sum().sum()
        vacc_tot = results[results['comp'] == 'yell'].pivot(values='value', columns='place', index='date').sum().sum()


        all_results = pd.concat([all_results, pd.DataFrame.from_dict(
            {'method_short': [shortname],
             'method': [strat.name],
             'infected': [yell_tot],
             'post_sample': [post_real],
             'doses': [vacc_tot],
             'scenario-beta': [scenario_name.split('-')[0]],
             'scenario-rate': [scenario_name.split('-')[1]],
             'scenario-tot': [scenario_name.split('-')[2]],
             'scenario': [scenario_name],
             'newdoseperweek': [int(scenario_name.split('-')[2][1:])]
             })])

        print(f"{scenario_name}, {post_real}, {shortname} done in {time.time()-tic} s, with compute a new strat set as {strat.compute_new_strat}")

    print(f"{scenario_name}, {post_real} done in {time.time()-tic1} seconds")
    return all_results


setup_shared = ItalySetupProvinces(nnodes, ndays, when) #shared between thread, don't use everywhere

# Generate posterior
if False:
    for post_real in tqdm.tqdm(np.arange(1, 102+1)):
        p = COVIDParametersOCP.OCParameters(setup=setup_shared, M=M, when=when, posterior_draw=post_real)
        with open(f'italy-data/full_posterior/parameters_{nnodes}_{when}_{post_real}.pkl', 'wb') as out:
            pickle.dump(p, out, pickle.HIGHEST_PROTOCOL)
        with open(f'italy-data/full_posterior/setup_{nnodes}_{when}.pkl', 'wb') as out:
            pickle.dump(setup_shared, out, pickle.HIGHEST_PROTOCOL)
    exit(0)


# Pick the right scenarios
scenarios = {pick_scenario(setup_shared, i)['name']: pick_scenario(setup_shared, i) for i in np.arange(15)}
pick = 'r15-'
scenarios = {k:v for (k,v) in scenarios.items() if pick in k}
print(f'doing {len(scenarios)}: {list(scenarios.keys())}')


pool = mp.Pool(mp.cpu_count())
if __name__ == '__main__':
    all_results = []
    tic = time.time()

    alt_strategies = {}
    for scenario_name, scenario in scenarios.items():
        alt_strategies[scenario_name] = create_all_alt_strategies(setup_shared, scenario_name, scenario)

    print("computing all scenarios on realization 102, the median realization, to construct all the alternative strategies")
    results_scn = pool.starmap(worker_one_posterior_realization, [(102, scenario_name, scenario, alt_strategies[scenario_name]) for scenario_name, scenario in scenarios.items()])
    all_results.append(pd.concat(results_scn))


    for scenario_name, scenario in scenarios.items():
        print(f'>>> Doing scenario {scenario_name}')
        for shortname, strat in alt_strategies[scenario_name].items():
            strat.compute_new_strat = False
        results_scn = pool.starmap(worker_one_posterior_realization,
                        [(post_real, scenario_name, copy.deepcopy(scenario), copy.deepcopy(alt_strategies[scenario_name])) for post_real in np.arange(1, 101+1)])
        all_results.append(pd.concat(results_scn))

    all_results = pd.concat(all_results)
    all_results.to_csv(f'{output_directory}/{output_prefix}-ALL.csv', index=False)

    print(f"Terminating succesfuly in {(time.time() - tic1)/3600} hours")











