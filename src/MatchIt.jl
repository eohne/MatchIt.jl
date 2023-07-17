module MatchIt

  # import packages:
    using DataFrames, GLM, PlotlyJS, Statistics,NearestNeighbors
    using HypothesisTests 
    using Random: randperm

    export MatchedIt, matchit,balance_plot
    # Re-exports from Base
    export getfield , show, summary
    # Re-exports from GLM
    export LogitLink, ProbitLink, @formula

    include("Structs.jl")
    include("Propensities.jl")
    include("Matching.jl")

"""
    matchit(df::DataFrame, f::FormulaTerm; link::GLM.Link01=LogitLink(), exact=[], maxDist::Float64=Inf,replacement::Bool=true,order::String="data",tolerance::Float64=-1.)

Performs propensity score matching. Also allows for exact matching on certain variables. 

# Arguments:  
  * df`::DataFrame`
  * f`::FormulaTerm`::The formula used in the GLM regression. (use @formula() to create it)
  * T`::String`: The name of the treatment variable (should be 0 or 1)
  * X`::Vector{String}`: The names of the covariates that are used in the propensity score matching
  * exact`::Vector`: The names of the variables that should be exactly matched. Defaults to no variables (`[]`)
  * link`::GLM.Link01`: The link function used for the regression. Either `LogitLink()` (default) or `ProbitLink()`
  * maxDist: The maximum distance between matched treatment and control observations to be included in the sample (defaults to `Inf`)
  * replacement`::Bool`: Whether the matching should occure with or without replacement. Without replacement is much slower and the result depends on the row order. (defaults to `true`)
  * order`::String`: In which order should the matching without replacement take place (defaults to `"data"`). Allows for:  
    - `"data"`: Follows the row order of the data  
    - `"smallest"`: Matches the smalles distances/propensity scores first  
    - `"largest"`: Matches the largest distances/propensity scores first  
    - `"random"`: Matches in a randomised order  
  * tolerance`::Float64`: In matching with replacement this value defines the minimum distance that is good enough such that the search for better matches is halted. (defaults to `-1` meaning that all control observations are compared to each treatment observation)


# Output:  
Returns a `MatchedIt` struct. It contains:
  * `df`: the original `DataFrame` 
  * `matched`: the matched `DataFrame`
  * `link`: the Link used (either `LogitLink()` or `ProbitLink()`)
  * `f`: the formula used in the matching process
  * `T`: the name of the treatment indicator variable

# Example  
Simple nearest neighbor PSM matching with replacement.  

```julia-repl
julia> m = matchit(input_data, @formula(treat ~ age + educ + race + nodegree + re74 + re75))

A MatchIt object:
Treatmeant variable: treat
Matching formula: treat ~ age + educ + race + nodegree + re74 + re75
Link: LogitLink
Number of observations: 614 (original), 370 (matched)
```
Returning the matched sample:  

```julia-repl
julia> m.matched
370×12 DataFrame
 Row │ Column1  treat  age    educ   race     married  nodegree  re74      re75        ⋯
     │ String7  Int64  Int64  Int64  String7  Int64    Int64     Float64   Float64     ⋯
─────┼──────────────────────────────────────────────────────────────────────────────────
   1 │ PSID368      0     40     11  black          1         1      0.0       0.0     ⋯
   2 │ PSID202      0     20      9  hispan         1         1      0.0    1283.66     
   3 │ PSID293      0     31     12  black          1         0      0.0      42.9677   
   4 │ PSID196      0     18     11  black          0         1      0.0    1367.81     
  ⋮  │    ⋮       ⋮      ⋮      ⋮       ⋮        ⋮        ⋮         ⋮          ⋮       ⋱
 367 │ NSW182       1     25     14  black          1         0  35040.1   11536.6     ⋯
 368 │ NSW183       1     35      9  black          1         1  13602.4   13830.6      
 369 │ NSW184       1     35      8  black          1         1  13732.1   17976.2      
 370 │ NSW185       1     33     11  black          1         1  14660.7   25142.2      
                                                          3 columns and 362 rows omitted
```

"""
function matchit(df::DataFrame, f::FormulaTerm; link::GLM.Link01=LogitLink(), exact=[], maxDist::Float64=Inf,
  replacement::Bool=true,order::String="data",tolerance::Float64=-1.)
  if typeof(exact) <:AbstractString
    exact = [exact]
  end
    data = copy(df)
    get_propensities!(data,f,link)
    T = String(StatsModels.termvars(f.lhs)[1])
    matched =NearestNeighbor(data,T,exact,maxDist,replacement, order,tolerance)
    return MatchedIt(data,matched,link,f,Symbol(T))
end

    include("Summary.jl")
    include("Plotting.jl")    
end