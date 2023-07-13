"""
    NearestNeighbour(data::DataFrame,T::String,exact::Vector{String}=String[],maxDist=Inf)

Performs Nearest Neighbor matching. Allows for exact matches on certain variables before performing the nearest neighbor search.

# Arguments:  
  * data`::DataFrame`: Needs to contain a treatment variable and the propensity score (`:PS`) (optionally the variables that need to be exactly matched)
  * T`::String`: The name of the treatment indicator (0 or 1)
  * exact::`Vector{String}`: a vector of the variable names used for exact matching.
  * maxDist: The maximum allowed distance between nearest NearestNeighbors

# Output:
A `DataFrame` containg the matched sample. Includes all variables from the provided DataFrame plus the distance metric `:dist`
"""
function NearestNeighbour(data::DataFrame,T::String,exact::Vector{String}=String[],maxDist=Inf)
    if isempty(exact)
        data.exact_match .= "1"
    else
       data.exact_match .= string.([string.(data[:,i],"_") for i in exact]...)
    end 

    treated = @view data[isequal.(data[:,T],1),:]
    control = @view data[isequal.(data[:,T],0),:]
    matched = DataFrame()
    no_exact_match = 0
    for exact_m  in unique(data.exact_match)
        c_df = @view control[isequal.(exact_m,control.exact_match),:]
        t_df = @view treated[isequal.(exact_m,treated.exact_match),:]
        if isempty(c_df) || isempty(t_df)
            no_exact_match +=1
            continue
        end
        tree = KDTree(c_df.PS', leafsize = 10)
        idx,dist = nn(tree,t_df.PS')
        out_c = c_df[idx,:]
        out_c.dist = dist
        t_df.dist = dist
        append!(matched, out_c)
        append!(matched,t_df)
    end
    if isequal(no_exact_match, size(unique(data.exact_match),1))
        throw("Zero cases with exact matches!")
    elseif no_exact_match>0
        @warn "Could not find exact matches for $(no_exact_match) groups!"
    end
    matched = matched[matched.dist.<=maxDist,:]
    select!(matched,Not(:exact_match))
    select!(data,Not(:exact_match))
    return matched
end;

