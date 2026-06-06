# Intro
This is the Github repo for the paper "Making Recursive Bayesian Inference Robust."

While Bayesian inference has become increasingly popular with advances in computational resources, its algorithms can be computationally prohibitive and may not scale with large datasets. This has led to growing interest in alternative algorithms, such as approximation methods and variants of Markov chain Monte Carlo. Among these approaches, prior proposal–recursive Bayesian (PP-RB) inference facilitates scalable Bayesian computation by recursively updating the posterior distribution across stages and utilizing parallel computing resources. While the well-known ``degeneracy'' issue in PP-RB has been studied, another limitation that PP-RB can yield incorrect inferences when posterior distributions shift substantially between stages has remained unsolved. To address this, we propose parallel-tempered prior proposal-recursive Bayesian (PPP-RB) inference, which extends PP-RB by leveraging the key idea underlying Metropolis-coupled Markov chain Monte Carlo. We show both theoretically and empirically that PPP-RB targets the true posterior distribution. We illustrate PPP-RB through numerical studies and real data analysis in application to earthquake count data and sea surface salinity in the North Atlantic region. In these applications, we compare PPP-RB with PP-RB and a standard MCMC, demonstrating that PPP-RB is more efficient in terms of effective sample size per elapsed time.

# Data
The earthquake dataset are provided by 
The Green sea turtle dataset are provided by Florida Fish and Wildlife Conservation CommissionFish and Wildlife Research Institute (2024)[1]. The BTS datasets are available as a USGS data release [2].

# Script
## Simulation_study
- lang_model_5_3_check_multiply.stan: Stan code used to fit the model.
- stan_bts_public.R: R code for the analysis.

## Data_analysis
- lang_model_3.stan: Stan code used to fit the model.
- stan_sea_grass_public.R: R code for the analysis.

# Reference
[1] Florida Fish and Wildlife Conservation Commission-Fish and Wildlife Research Institute (2024). Seagrass Florida. https://geodata.myfwc.com/datasets/seagrass-habitat-in-florida/explore?location=27.466248

[2] Fueka, A. B., Nafus, M. G., Bailey, L., Yackel Adams, A. A., and Hooten, M. B. (2022). Exogenous and endogenous factors influence invasive reptile movement at multiple scales, 2018–2019: U.S. Geological Survey data release. https://doi.org/10.5066/P948KRN3
