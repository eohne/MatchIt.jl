module MatchIt

  # import packages:
    using DataFrames, GLM, PlotlyJS, Statistics,NearestNeighbors
    using HypothesisTests 
    using Random: randperm

    export MatchedIt, matchit, balance_plot,match_table
    # Re-exports from Base
    export getfield 
    # Re-exports from GLM
    export LogitLink, ProbitLink, @formula

    include("Structs.jl")
    include("Propensities.jl")
    include("Matching.jl")

"""
    matchit(matchit(df::DataFrame, f::FormulaTerm; link::GLM.Link01=LogitLink(), exact, maxDist::Float64=Inf,replacement::Bool=true,order::String="data",tolerance::Float64=-1.))

Performs propensity score matching. Also allows for exact matching on certain variables. 

# Arguments:  
  * df`::DataFrame`
  * f`::FormulaTerm`::The formula used in the GLM regression. (use @formula() to create it)
  * T`::String`: The name of the treatment variable (should be zero or 1)
  * X`::Vector{String}`: The names of the covariates that are used in the propensity score matching
  * exact`::Vector{String}`: The names of the variables that should be exactly matched. Defaults to no variables (`String[]`)
  * link`::GLM.Link01`: The link function used for the regression. Either `LogitLink()` (default) or `ProbitLink()`
  * maxDist: The maximum distance between matched treatment and control observations to be included in the sample (defaults to `Inf`)
  * replacement`::Bool`: Whether the matching should occure with or without replacement. Without replacement is much slower and the result depends on the row order. (defaults to `true`)
  * order`::String`: In which order should the matching without replacement take place (defaults to `"data"`). Allows for:  
    - `"data"`: Follows the row order of the data  
    - `"smallest"`: Matches the smalles distances/propensity scores first  
    - `"largest"`: Matches the largest distances/propensity scores first  
    - `"random"`: Matches in a randomised order  
  * tolerance`::Float64`: In matching with replacement it defines the minimum distance that is good enough and no better matches are searched for. (defaults to `-1` meaning that all control observations are compared to each treatment observation)


# Output:  
Returns a `MatchedIt` struct. It contains:
  * `df`: the original `DataFrame` 
  * `matched`: the matched `DataFrame`
  * `link`: the Link used (either `LogitLink()` or `ProbitLink()`)
  * `T`: the name of the treatment indicator variable

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
    MatchedIt(data,matched,link,Symbol(T))
end

    include("Balance Table.jl")
    include("Plotting.jl")    
end
