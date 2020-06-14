# csr-summary-stats
Package for introductory analysis of Spatial Point Processes, including basic utility functions, summary statistics and metrics for 2D point processes. To be primarily used to test the hypothesis of CSR.

## Included Summary Statistics
These functions calculate the associated metric of an observed Point Process `PP`. When `restrict_domain` is a valid argument, the metric is calculated on a subset of the original domain, in a first attempt to account for edge effects.

`F_hat(PP, r, num_sampled_points=1000, plot=False)`

`G_hat(PP, r, plot=False)`

`J_hat(PP, r, num_F_sampled_points=1000, plot=False)`

`K_hat(PP, r, restrict_domain=True, plot=False)`

`L_hat(PP, r, restrict_domain=True, plot=False)`

`O_hat(PP, r, bandwidth=0.1, restrict_domain=True, kernel="BK", plot=False)`

`PC_hat(PP, r, bandwidth=0.1, restrict_domain=True, kernel="BK", plot=False)`

## Simulation
Any of the above summary statistic functions can be simulated by calling the following. Resulting max/min envelopes will then be plotted and can be compared against a given realisation. Currently, only the default arguments of a summary statistic function can be used during simulation.

`simulate_summary_statistic(PP, r, summary_func, n_sims=100, plot=False)`

## Real-Valued Hypothesis Test
`CSR_hypothesis_test_dmin(observed_PP, significance_level, n_sim=500, plot_dist=False)`

## Testing
`pytest test/point_process_utils_test.py`
