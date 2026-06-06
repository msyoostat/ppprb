# Intro
This is the Github repo for the paper "Making Recursive Bayesian Inference Robust."

While Bayesian inference has become increasingly popular with advances in computational resources, its algorithms can be computationally prohibitive and may not scale with large datasets. This has led to growing interest in alternative algorithms, such as approximation methods and variants of Markov chain Monte Carlo. Among these approaches, prior proposal–recursive Bayesian (PP-RB) inference facilitates scalable Bayesian computation by recursively updating the posterior distribution across stages and utilizing parallel computing resources. While the well-known ``degeneracy'' issue in PP-RB has been studied, another limitation that PP-RB can yield incorrect inferences when posterior distributions shift substantially between stages has remained unsolved. To address this, we propose parallel-tempered prior proposal-recursive Bayesian (PPP-RB) inference, which extends PP-RB by leveraging the key idea underlying Metropolis-coupled Markov chain Monte Carlo. We show both theoretically and empirically that PPP-RB targets the true posterior distribution. We illustrate PPP-RB through numerical studies and real data analysis in application to earthquake count data and sea surface salinity in the North Atlantic region. In these applications, we compare PPP-RB with PP-RB and a standard MCMC, demonstrating that PPP-RB is more efficient in terms of effective sample size per elapsed time.

# Data
world_series_earthquake_2_2.RData: The earthquake dataset are obtined from United States Geological Survey (2026) [1].
sss_7.RData: SSS dataset are obtained from ESA Sea Surface Salinity Climate Change Initiative [2].


# Script

## Data_analysis
- lang_model_3.stan: Stan code used to fit the model.
- stan_sea_grass_public.R: R code for the analysis.

# Reference
[1] United States Geological Survey (2026). USGS earthquake catalog. https://earthquake.usgs.gov/earthquakes/search/.

[2] Boutin, J., Vergely, J.-L., Reul, N., Catany, R., Jouanno, J., Martin, A., Rouffi, F., Bertino, L., Bonjean, F., Corato, G., Gévaudan, M., Guimbard, S., Khvorostyanov, D., Kolodziejczyk, N., Matthews, M., Olivier, L., Raj, R., Rémy, E., Reverdin, G., Supply,A., Thouvenin-Masson, C., Vialard, J., Sabia, R., and Mecklenburg, S. (2024). ESA Sea Surface Salinity Climate Change Initiative (Sea_Surface_Salinity_cci): Monthly sea surface salinity product for the Northern Hemisphere on a 25km EASE grid, v04.41, for 2010 to 2022.
