module PSMatching

    using DataFrames, GLM, PlotlyJS, Statistics,NearestNeighbors
    using HypothesisTests

    export PSM, psmatch, balance_plot,match_table
    # Re-exports from Base
    export getfield 
    # Re-exports from GLM
    export LogitLink, ProbitLink

    include("Structs.jl")
    include("Propensities.jl")
    include("Matching.jl")

"""
    psmatch(df::DataFrame, T::String, X::Vector{String}, exact::Vector{String}=[], link::GLM.Link01=LogitLink(), maxDist=Inf)

Performs propensity score matching. Also allows for exact matching on certain variables. 

# Arguments:  
  * df`::DataFrame`
  * T`::String`: The name of the treatment variable (should be zero or 1)
  * X`::Vector{String}`: The names of the covariates that are used in the propensity score matching
  * exact`::Vector{String}`: The names of the variables that should be exactly matched. Defaults to no variables (`String[]`)
  * link`::GLM.Link01`: The link function used for the regression. Either `LogitLink()` (default) or `ProbitLink()`
  * maxDist: The maximum distance between matched treatment and control observations to be included in the sample (defaults to `Inf`)

# Output:  
Returns a `PSM` struct. It contains:
  * `df`: the original `DataFrame` 
  * `matched`: the matched `DataFrame`
  * `link`: the Link used (either `LogitLink()` or `ProbitLink()`)
  * `T`: the name of the treatment indicator variable

"""
function psmatch(df::DataFrame, T::String, X::Vector{String}, exact::Vector{String}=[], link::GLM.Link01=LogitLink(), maxDist=Inf)
    data = copy(df)
    get_propensities!(data,T,X,link)
    matched =NearestNeighbour(data,T,exact,maxDist)
    PSM(data,matched,link,Symbol(T))
end

    include("Balance Table.jl")
    include("Plotting.jl")    
end