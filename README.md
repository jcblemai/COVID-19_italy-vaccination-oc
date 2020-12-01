Old version contains:
- CVO_implicit: implicit force of infection: 1st order approximation of the force of infection to mitigate
the fact that the relationship between nodes is not taken into account
- Between ocp_trial and ocp_trial-vector (diff: in vector: loop over time then over nodes, and compute foi once for all nodes. Moreover use matrix multiplication over loop (which I later found that it produced non-structural zero, so this change was undone) )

Notebooks:

- ocp_trial_vector: my starting version
- ocp_trial_vector-crtl: control out of the loop, but no structural zero (so matrix formula)
- ocp_trial_vector-crtl-struct0: same with strutural zero. Most refined version, the final program 
is based on this.
- ocp_trial-vector3-JIT: just in time trial. Not really happy yet but might be worth considering.
