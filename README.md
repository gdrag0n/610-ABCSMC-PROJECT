# 610-ABCSMC-PROJECT


The file ABC_SMC_and_Rejection_Sampling.Rmd implements ABC SMC and ABC Rejection Sampling for model selection and parameter inference between two models:

m0: i.i.d. Bernoulli.
m1: Ising-like model.

# Features

Algorithms:
ABC SMC with a decreasing tolerance schedule.
Benchmarking with ABC Rejection Sampling.

Visualizations:
Tolerance schedule, ESS, and posterior parameter distribution.

Outputs
Posterior probabilities (P(m=0) and P(m=1)).
Diagnostic plots for analysis.


The file ABC_SMC_Model_Selection.Rmd implements Approximate Bayesian Computation Sequential Monte Carlo (ABC SMC) for model selection and parameter inference. It compares two models:

m0: i.i.d. Bernoulli model.
m1: Ising-like model.

# Features
Simulates data for both models and excludes degenerate datasets.
Uses sufficient statistics (S0 and S1) and a distance metric.
Refines posterior probabilities through particle resampling over multiple populations.
Visualizes true vs. inferred posterior probabilities.
