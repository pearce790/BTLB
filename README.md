# BTL-Binomial MFM Model with Rater Effects

This Github repository contains source code for estimation of the BTL-Binomial model and its variants under rater effects and preference heterogeneity. All functions necessary for random data generation, density calculation, MAP estimation and MFMM estimation of the models proposed and discussed in the paper are provided in the file "source.R".

## Instructions for Replication:

1. To replicate Section 4.1 of the paper, open the "Section 4.1 Simulation Study" folder and run "main.R".
2. To replicate Section 4.2 of the paper and associated Appendix B.2, open the "Section 4.2 Conference" folder and run "main.R".
3. To replicate Section 4.3 of the paper and associated Appendices B.3 and B.4, open the "Section 4.3 Panel Review" folder and run "main.R".

## Note on Run Times

Due to the large number of distinct simulation settings and MCMC iterations included in the analysis, the replication files may be slow to run without parallelization. On a relatively standard but high-powered 2023 MacBook Pro, the replication files listed above will take approximately 2--4 days to run in total without parallelization (1--2 days for replication [1] and [2] and about 45 minutes for replication [3]). However, the code to replicate [1] and [2] may be completely parallelized, permitting each to run in about 1 hour total. Details are provided in the "main.R" file associated with each replication folder, which include slurm files for replication.
