mutable struct PSM
    df::AbstractDataFrame # Entire DataFrame (unmatched)
    matched::AbstractDataFrame # the part that is matched
    link # Type of matching (only NearestNeighbors supported atm)
    T::Symbol # Treatment variable
end

function Base.getindex(x::PSM,idx::Symbol)
    return getfield(x,idx)
end