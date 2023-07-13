# Propensity Score Matching
[![Build Status](https://github.com/eohne/PSMatching.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/eohne/PSMatching.jl/actions/workflows/CI.yml)

Probably not the most efficient package but does the job.  

Can use Logit or Probit regressions to calculate propensity scores.  

Allows for some exact matching on some covariates. Amongst the control observation that match exactly the one with the closest propensity score is chosen.  

Uses KDTrees from NearestNeighbors.jl for fast matching.