# 610-ABCSMC-PROJECT

This repository contains R Markdown files implementing Approximate Bayesian Computation Sequential Monte Carlo (ABC SMC) and ABC Rejection Sampling for model selection and parameter inference. It compares two models:

Model m0: i.i.d. Bernoulli model.
Model m1: Ising-like model.

Files
1. ABC_SMC_and_Rejection_Sampling.Rmd

Algorithms:
ABC SMC with a decreasing tolerance schedule.
ABC Rejection Sampling for benchmarking.

Visualizations:

Tolerance schedule.

Effective Sample Size (ESS).

Posterior parameter distribution.

Outputs:
Posterior probabilities (P(m=0) and P(m=1)).
Diagnostic plots for analysis.

2. ABC_SMC_Model_Selection.Rmd

Features:

Simulates data for both models and excludes degenerate datasets.

Uses sufficient statistics (S0 and S1) and a distance metric.

Refines posterior probabilities through particle resampling over multiple populations.

Visualizes true vs. inferred posterior probabilities.

Key Components

Models:
m0: i.i.d. Bernoulli distribution.
m1: Ising-like model with state dependence.

Algorithms:
ABC SMC for iterative posterior refinement.
ABC Rejection Sampling for comparison.
