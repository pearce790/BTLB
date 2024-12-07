# BTL-Binomial MFM Model with Rater Effects

This Github repository contains source code for estimation of the BTL-Binomial model and its variants under rater effects and preference heterogeneity. All functions necessary for random data generation, density calculation, MAP estimation and MFMM estimation of the models proposed and discussed in the paper are provided in "20240429_source.R".

Instructions for Replication:

1. To replicate Section 4.1 of the paper, open the "Section 4.1 Simulation Study" folder and run "main.R".
2. To replicate Section 4.2 of the paper and associated Appendix B.2, open the "Section 4.2 Conference" folder and run "main.R".
3. To replicate Section 4.3 of the paper and associated Appendices B.3 and B.4, open the "Section 4.3 Panel Review" folder and run "main.R".

Due to the large number of MCMC iterations and distinct simulation settings included in the analysis, the replication files may be slow to run. On a standard 2023 MacBook Pro with 32GB memory, the replication files listed above took approximately ..., ...., and 45 minutes to run, respectively.