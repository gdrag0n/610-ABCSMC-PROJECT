# 610-ABCSMC-PROJECT

# Features of ABC_SMC_and_Rejection_Sampling.Rmd

Models: 
m0: i.i.d. Bernoulli.
m1: Ising-like model.

Algorithms:
ABC SMC with a decreasing tolerance schedule.
Benchmarking with ABC Rejection Sampling.

Visualizations:
Tolerance schedule, ESS, and posterior parameter distribution.

Outputs
Posterior probabilities (P(m=0) and P(m=1)).
Diagnostic plots for analysis.


# Features of ABC_SMC_Model_Selection.Rmd

Models:
m0: i.i.d. Bernoulli model.
m1: Ising-like model.

Features
Simulates data for both models and excludes degenerate datasets.
Uses sufficient statistics (S0 and S1) and a distance metric.
Refines posterior probabilities through particle resampling over multiple populations.
Visualizes true vs. inferred posterior probabilities.
