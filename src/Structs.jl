mutable struct MatchedIt
    df::AbstractDataFrame # Entire DataFrame (unmatched)
    matched::AbstractDataFrame # the part that is matched
    link # Type of matching (only NearestNeighbors supported atm)
    f::FormulaTerm #formula used
    T::Symbol # Treatment variable
    dist_type::String # mahalanobis or glm
end

function Base.getindex(x::MatchedIt,idx::Symbol)
    return getfield(x,idx)
end

function Base.show(io::IO, obj::MatchedIt) 
    print(io::IO,"\nA MatchIt object:\n")
    print(io::IO,"Treatmeant variable: $(string(obj.T))\n")
    print(io::IO,"Matching formula: $(obj.f)\n")
    if isequal(lowercase(obj.dist_type), "glm")
    print(io::IO,"Matched by: $(uppercase(obj.dist_type))\n")
    print(io::IO,"Link: $(string(typeof(obj.link)))\n")
    else
    print(io::IO,"Matched by: $(uppercasefirst(obj.dist_type))\n")
    end
    print(io::IO,"Number of observations: $(nrow(obj.df)) (original), $(nrow(obj.matched)) (matched)\n")
end
