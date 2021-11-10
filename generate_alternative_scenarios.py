import numpy as np
from ItalySetup import ItalySetupProvinces
from covidOCP import COVIDVaccinationOCP as COVIDVaccinationOCP
from covidOCP import COVIDAgeStructuredOCP as COVIDAgeStructuredOCP
from covidOCP import COVIDParametersOCP
import pickle
import matplotlib.pyplot as plt
import click
import sys, os
from scenarios_utils import pick_scenario, build_scenario
import pandas as pd
import multiprocessing as mp



# Replace the jupyter notebook that was based on matlab called generate_all_scn.

nx = 9
states_names = ['S', 'E', 'P', 'I', 'A', 'Q', 'H', 'R', 'V']
when = 'future-mobintime'
prefix = f'altstratint'
outdir = 'helvetios-runs/2021-11-05-107_90/'
generated_dir = 'model_output/2021-11-09'

nnodes = 107  # nodes
ndays_ocp = 90
ndays = 90

setup = ItalySetupProvinces(nnodes, ndays, when)
setup_ocp = ItalySetupProvinces(nnodes, ndays_ocp, when)
M = setup.nnodes
N = len(setup.model_days) - 1

with open(f'italy-data/parameters_{nnodes}_{when}.pkl', 'rb') as inp:
    p = pickle.load(inp)

os.makedirs(f'{generated_dir}', exist_ok=True)

scenarios = {pick_scenario(setup, i)['name']: pick_scenario(setup, i) for i in np.arange(15)}
pick = 'r15-'
scenarios = {k:v for (k,v) in scenarios.items() if pick in k}
print(f'doing {len(scenarios)}: {list(scenarios.keys())}')

pool = mp.Pool(mp.cpu_count())

class AlternativeStrategy:
    def __init__(self, setup, scenario, alloc_function, decision_variable, dv_per_pop, require_projection):
        self.maxvaccrate_regional, self.delivery_national, self.stockpile_national_constraint, _ = build_scenario(setup, scenario)
        self.M = setup.nnodes
        self.pop_node = setup.pop_node
        self.maxvaccrate_regional = self.maxvaccrate_regional[:,0] # stays the same over the course, so take the first one

        # updated states variables
        self.unvaccinated = np.copy(setup.pop_node)
        self.stockpile = 0

        self.decision_variable = decision_variable
        self.require_projection = require_projection

        # The decision variable is per habitant
        self.name = decision_variable
        self.shortname = decision_variable[:3]

        self.dv_per_pop = dv_per_pop
        self.divider = np.ones(self.M)
        if dv_per_pop:
            self.divider = self.pop_node
            self.shortname += '_pp'

        if alloc_function == 'focused':
            self.alloc_function = self.focused_alloc
            self.name += ' (focused)'
            self.shortname += '_f'
        elif alloc_function == 'proportional':
            self.alloc_function = self.proportional_alloc
            self.name += ' (proportional)'
            self.shortname += '_p'

    def focused_alloc(self, decision_df_sorted, nd, nodename):
        return self.maxvaccrate_regional[nd]

    def proportional_alloc(self, decision_df_sorted, nd, nodename):
        return self.stockpile * decision_df_sorted.loc[nodename]['value'] /  decision_df_sorted['value'].sum()

    def allocate_now(self, decision_variable_array, today_idx):
        # Sort the decision variable dataframe:
        print(decision_variable_array)
        decision_variable_df = pd.DataFrame(decision_variable_array, index=setup.ind2name, columns=['value'])
        print(decision_variable_df)
        decision_variable_df.sort_values('value', ascending=False, inplace=True)
        print(decision_variable_df)
        alloc_now = np.zeros(self.M)
        self.stockpile += self.delivery_national[today_idx]
        for nodename in decision_variable_df.index:
            nd = setup.ind2name.index(nodename)
            to_allocate = self.alloc_function(decision_variable_df, nd, nodename)
            print(to_allocate.shape, self.maxvaccrate_regional.shape)
            to_allocate = min(to_allocate, self.unvaccinated[nd], self.stockpile, self.maxvaccrate_regional[nd])
            alloc_now[nd] = to_allocate
            self.stockpile        -= to_allocate
            self.unvaccinated[nd] -= to_allocate
        return alloc_now

    def get_allocation(self, today_idx, susceptible, incidence):
        if 'susceptible' in self.decision_variable:
            return self.allocate_now(susceptible/self.divider, today_idx)
        elif 'incidence' in self.decision_variable:
            return self.allocate_now(incidence/self.divider, today_idx)
        elif 'population' in self.decision_variable:
            return self.allocate_now(self.pop_node/self.divider, today_idx)


def create_all_alt_strategies(scenario):
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
                                                alloc_function,
                                                decision_variable,
                                                dv_per_pop,
                                                require_projection)
                alt_strategies[alt_strat.shortname] = alt_strat

    print(f'generated {len(alt_strategies.keys())} strategies: {list(alt_strategies.keys())}')

    return alt_strategies

for scenario_name, scenario in scenarios.items():
    print(f'Doing scenario {scenario_name}')
    p.apply_epicourse(setup, scenario['beta_mult'])
    control_initial = np.zeros((M, N))

    alt_strategies = create_all_alt_strategies(scenario)
    for shortname, strat in alt_strategies.items():
        results, state_initial, yell, = COVIDVaccinationOCP.accurate_integrate(N,
                                                                               setup=setup,
                                                                               parameters=p,
                                                                               controls=None,
                                                                               save_to=f'{outdir}{prefix}-{shortname}-{nnodes}_{ndays}-nc',
                                                                               only_yell = False,
                                                                               alloc_strat = strat)
    exp_accurate_py = results[results['comp'] == 'yell'].pivot(values='value', columns='place', index='date')

    #all_sims = pool.starmap(compute_alt_scenarios,
    #                    [(post_rea, scenario_name, scenario) for post_rea in range(100)])










