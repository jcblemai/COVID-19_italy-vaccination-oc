# Optimal control of the spatial allocation of COVID-19 vaccines: Italy as a case study

This repository contains the code associated with the manuscript **Optimal control of the spatial allocation of COVID-19 vaccines: Italy as a case study**, by _Joseph C. Lemaitre, Damiano Pasetto, Mario Zanon, Enrico Bertuzzo, Lorenzo Mari, Stefano Miccoli, Renato Casagrandi, Marino Gatto and Andrea Rinaldo_ in PLoS Computational Biology.

The full code, along with the whole commit history from the project inception in June 2020 to the publication in June 2022 is made available, along with serveral unfruitful exploratory trials (using collocation methods, age-based optimizations.

## Manuscript abstract

While campaigns of vaccination against SARS-CoV-2 are underway across the world, communities face the challenge of a fair and effective distribution of a limited supply of doses. Current vaccine allocation strategies are based on criteria such as age or risk. In the light of strong spatial heterogeneities in disease history and transmission, we explore spatial allocation strategies as a complement to existing approaches. Given the practical constraints and complex epidemiological dynamics, designing effective vaccination strategies at a country scale is an intricate task.  We propose a novel optimal control framework to derive the best possible vaccine allocation for given disease transmission projections and constraints on vaccine supply and distribution logistics. As a proof-of-concept, we couple our framework with an existing spatially explicit compartmental COVID-19 model tailored to the Italian geographic and epidemiological context. We optimize the vaccine allocation on scenarios of unfolding disease transmission across the 107 provinces of Italy, from January to April 2021. For each scenario, the optimal solution significantly outperforms alternative strategies that prioritize provinces based on incidence, population distribution, or prevalence of susceptibles. Our results suggest that the complex interplay between the mobility network and the spatial heterogeneities implies highly non-trivial prioritization strategies for effective vaccination campaigns. Our work demonstrates the potential of optimal control for complex and heterogeneous epidemiological landscapes at country, and possibly global, scales.

## Figures

![fig1](https://user-images.githubusercontent.com/7485811/172430841-410d4769-268c-44de-8aac-d60dad1872eb.png)
_A) Diagram representing the compartments of the epidemiological model and the possible transitions in a single province. We control the vaccination rate (teal arrows), aiming at minimizing incident infections (pink arrow). Individuals in compartments outside of the yellow block are able to move along the mobility network shown in (B), hence the force of infection in a province is coupled with the dynamics of other connected provinces. To reduce the problem to a tractable size, we only consider the most important connections (red edges) when optimizing, but we use the full network (red and grey edges) to assess our strategies. A discussion on the effect of this simplification is provided in SI. Nodes' size and color display each province's population, and edges' width shows the strength of the coupling between a pair of provinces._

<img width="1416" alt="Screen Shot 2022-06-07 at 18 17 55" src="https://user-images.githubusercontent.com/7485811/172431308-aae3ffd0-9fb1-43a9-af07-6b2e5c5b7df8.png">
_Optimal vaccine allocation for the baseline, pessimistic scenario. (A) Cumulative proportion of vaccine doses administered in the 107 provinces, some of which are highlighted. The local distribution rate is limited by a rate that is proportional to the population. This logistic constraint is visualized here as the maximum slope, equal for every province.
    (B) Stacked cumulative absolute number of vaccines in the 107 provinces of Italy. The national stockpile is shown in black and is replenished every week (on Mondays) with 479'700 doses. We display the name of the provinces with a final allocation of more than 150'000 doses._



## Files Descriptions:

### Optimal control for infectious disease
The file to run to solve the optimal control `main.py`, with the following arguments (support the --help keyword):
```
  -s, --scenario_id INTEGER    Index of scenario to run
  -n, --nnodes INTEGER         Spatial model size to run
  -t, --ndays INTEGER          Number of days to run
  --use_matlab BOOLEAN         whether to use matlab for the current run
                               [default: False]
  -a, --age_struct BOOLEAN     Whether to use agestructured OCP  [default:
                               False]
  -f, --file_prefix TEXT       file prefix to add to identify the current set
                               of runs.  [default: test]
  -d, --output_directory TEXT  Where to write runs  [default: model_output/]
  -o, --optimize BOOLEAN       Whether to optimize  [default: True]
```
Other variations of this file (main-ag with age stratification, main-ag-reg for age-stratification at regional scale, main-eq for optimization with equality constraints, main-agpost for age as a postprocessing output) are not used in the final manuscript.

From this file, `ItalySetup.py` handles the spatial setup (population, mobility, shapefiles), and `scenario_utils.py` all the information about scenarios (id for output filenames, maximum stockpile through time, feasible initial condition).

Folder `batch_scripts/` contains slurm scripts to run both the optimization and the alternative strategies on HPC clusters.

Folder `italy-data/` contains the data of the Italian epidemic.

Folder `matlab/` contains the data assimilation code. This provides the parameters to start the optimization from. Namely, module `OCParameters` in file  `covidOCP/COVIDParametersOCP.py` creates a picklable parameter object that contains every parameter of the setup. This parameter object is either build by running a matlab engine, and running file  `minimal_interface.m` to populate the namespace (with argument `future-mobintime` for the results displayed in the paper with mobility that change in time, `future` for the results in the medRxiv preprint, and `past for the parameters at the start of 2020). It's also possible to choose the posterior draw, with 102 being the index of the median posterior draw used for optimization.

Folder `covidOCP` contains the optimal control problem definition and solving capability. As for `main.py`, several variations exist (with equity, flavors of age-stratification, ...). For the results in the papers, we used file `covidOCP/COVIDVaccinationOCP.py`. It has functions defining the right-hand side of the within province dynamic, code to integrate it with either an euler or a runge-kutta 4 scheme, two functions to solve it without optimization (e.g to get a solution, assess a defined strategy):
- `integrate()` integrate as the optimal control problem does, i.e with fixed daily mobility component, pruned mobility... It is mainly used to provide initial conditions for the optimal control problem.
- `accurate_integrate()` integrate with all bells and whistles, as the data assimilation does. It is used to evaluate strategies.

Class `COVIDVaccinationOCP` builds the optimal control problem from the above component, iteratively for each node and each day. It's quite long! Method `update` allows one to update some part of the problem without re-building it.  Then, `solveOCP` gets the optimal vaccine allocation and `saveOCP` saves it to disk with the right file names.

#### Some history:
Old version contains:
- CVO_implicit: implicit force of infection: 1st order approximation of the force of infection to mitigate
the fact that the relationship between nodes is not taken into account
- Between ocp_trial and ocp_trial-vector (diff: in vector: loop over time then over nodes, and compute foi once for all nodes. Moreover use matrix multiplication over loop (which I later found that it produced non-structural zero, so this change was undone) )

### Results and figure
First, after optimizing the optimal control for one scenarion (running `main.py`) or multiple scenarios with the batch script, you'd need to generate the alternative strategies outcomes to compare them with the optimal strategy.
Scripts:
- `generate_alternative_scenarios.py` create and evaluate all alternative strategies. Each strategy is identified by a name, and it runs in parallel as it's very long. It also evaluates the optimal solution.
- `generate_alternative_scenarios_sensitivity_analysis.py` does the same but instead of evaluating on the posterior, performs a sort of sensitivity analysis where it shuffles all province's dynamics (see the SI).

After that, it's all notebooks :
- `analysis_beta.ipynb` : scenarios of the rate of transmission in time (pessimistic and optimistic) are defined here with some plots of the dynamics with no control (paper Fig 2 and some SI figures)
- `analysis_opt-newFIG.ipynb`: analyse the optimal allocation for one scenario (paper Fig 5,6 and some SI figures)
- `analysis_opt.ipynb`: same as above, but different figures.
- `analysis_scn.ipynb`: compare across the alternative allocation and optimal allocation (paper fig 3 and some SI)
- `analysis_ts.ipynb`: Temporal analysis of the optimal solution (Fig 4 and some SI Figures)
- `compare integration.ipynb`: compare the simplified and full integration scheme, for a SI figure.
- `dev-agestratified.ipynb`: development stage (scratchpad) notebook
- `dev-integration.ipynb`: development stage (scratchpad) notebook
- `dev-ocp_trial-vector.ipynb`: development stage (scratchpad) notebook. my starting version. 
- `dev-ocp_trial-vector2-crtl-out.ipynb`: development stage (scratchpad) notebook. control out of the loop, but no structural zero (so matrix formula)
- `dev-ocp_trial-vector2-crtl-out-struct0.ipynb`: development stage (scratchpad) notebook. same with strutural zero. Most refined version, the final program is based on this.
- `dev-ocp_trial-vector3-JIT.ipynb`: development stage (scratchpad) notebook. just in time compiled trial. Not really happy yet but might be worth considering.
- `generate_all_scn.ipynb`: template for the file `generate_alternative_scenarios.py`
- `post-review-analysis.ipynb`: some more analysis in paper figure.

The remaining notebooks display some exploration of age-stratified optimal allocations.

In order to faciliate reproducing results without spending thousand of $ on compute, we include the raw output of optimal allocation used in the manuscript in folder `helvetios-runs/2021-11-17-107_90/` along with Ipopt outputs. This is the output from launching slurm file `batch_script/medium.run`. As such, all notebooks above should run seamlessly in a local folder. A file is named e.g `week-L-r150-t479700-id15-opt-107_90.csv` where:
- `week` is the run prefix I've choosed for these runs.
- `L` means Low transmission (hence optimistic scenario). vs. `U`
- `r150` is the maximum rate of distribution across italy (logistic constraints). per month in Million
- `t479700` is the weekly stockpile delivery in doses
- `id15` is the scenario id (to use later, summarize all the above
- `opt` means that the allocation is optimal (int when the allocation is the baseline, no vaccine one. It could also be an alternative strategy short name
- `107_90` means the full problem with all 107 provinces and optimizing 90 days ahead.

Ah, and I almost forgot, use the conda environment in file `environment_cross.yml` (or `environment.yml`if you have the same hardware as me).

Feel free to contact me before using this (jo.lemaitresamra[at]gmail.com, I'll be happy to guide you through choices we've made, an to send you the log of the project in time with every exploration.
