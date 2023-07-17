# MatchIt.jl
GitHub Repo: [https://github.com/eohne/MatchIt.jl](https://github.com/eohne/MatchIt.jl)

Performs Propensity Score Matching in pure Julia. The package is somewhat inspired by the R MatchIt package.  
The current implementation allows propsensity scores to be calculated using both Logistic and Probit regressions.  
The package implements one-to-one nearest neighbor matching with and without replacement.  
 - Matching with replacement uses KDTrees from the `NearestNeighbor.jl` package. 
 - Matching without replacement loops through (all) potential matches and is therefore slow.  
The package further allows for matching exactly on some covariates (e.g. match to the closest propensity score conditional on the two observations being from the same Firm and date).

The currently implemented functionality has been crosschecked with the `lalonda.csv` test data from the R MatchIt package and matches the output of the R package.  


## Installation

<!-- The package is NOT registered in the [`General`](https://github.com/JuliaRegistries/General) registry (yet). -->

You can manually install it by: 
