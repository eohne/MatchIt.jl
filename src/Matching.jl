function _nearest_neighbor_matching(treat::AbstractArray{Float64,1}, control::AbstractArray{Float64,1})
    tree = KDTree(control')
    return nn(tree, treat')
end

function _nearest_neighbor_matching(treat::AbstractArray{Float64,1}, control::AbstractArray{Float64,1}, order::String,tolerance::Float64)
    n_treat = size(treat, 1)
    n_control = size(control, 1)

    matching_indices = Int.(zeros(n_treat))
    distances = zeros(n_treat)
    used_indices = Bool.(zeros(n_control))

    if isequal(order, "smallest")
        sort_idx = sortperm(treat)
    elseif isequal(order, "largest")
        sort_idx = sortperm(-treat)
    elseif isequal(order, "random")
        sort_idx = randperm(n_treat)
    else
        sort_idx = collect(1:n_treat)
    end

    matching_indices .= -1
    distances .= Inf 

    for i in 1:n_treat
        if all(used_indices)
            continue
        end
        min_idx = 0
        min_dist = Inf
        
        for j in 1:n_control
            if min_dist<tolerance
                break
            end
            if used_indices[j]
                continue
            end
            dist = abs(control[j] - treat[sort_idx[i]])
            if dist < min_dist
                min_dist = dist
                min_idx = j
            end
        end
        matching_indices[sort_idx[i]] = min_idx
        distances[sort_idx[i]] = min_dist
        used_indices[min_idx] = true
    end
    return matching_indices, distances
end


"""
    NearestNeighbour(data::DataFrame,T::String,exact=String[],maxDist::Float64=Inf, replacement::Bool=true,order::String="data", tolerance::Float64=-1.)

Performs Nearest Neighbor matching. Allows for exact matches on certain variables before performing the nearest neighbor search.

# Arguments:  
  * data`::DataFrame`: Needs to contain a treatment variable and the propensity score (`:Dist`) (optionally the variables that need to be exactly matched)
  * T`::String`: The name of the treatment indicator (0 or 1)
  * exact::`Vector{String}`: A vector of the variable names used for exact matching.
  * maxDist`::Float64`: The maximum allowed distance between nearest neighbors (defaults to `Inf`)
  * replacement`::Bool`: Whether the matching should occure with or without replacement. Without replacement is much slower and the result depends on the row order. (defaults to `true`)
  * order`::String`: In which order should the matching without replacement take place (defaults to `"data"`). Allows for:  
    - `"data"`: Follows the row order of the data  
    - `"smallest"`: Matches the smalles distances/propensity scores first  
    - `"largest"`: Matches the largest distances/propensity scores first  
    - `"random"`: Matches in a randomised order  
  * tolerance`::Float64`: In matching with replacement it defines the minimum distance that is good enough and no better matches are searched for. (defaults to `-1` meaning that all control observations are compared to each treatment observation)

# Output:
A `DataFrame` containg the matched sample. Includes all variables from the provided DataFrame plus the distance metric `:Dist`
"""
function NearestNeighbor(data::DataFrame,T::String,exact=String[],maxDist::Float64=Inf, replacement::Bool=true,order::String="data", tolerance::Float64=-1.)
    if isempty(exact)
        treated = @view data[isequal.(data[:,T],1),:]
        control = @view data[isequal.(data[:,T],0),:]
        if replacement
            idx, dist = _nearest_neighbor_matching(treated.Dist,control.Dist)
        else
            idx, dist = _nearest_neighbor_matching(treated.Dist,control.Dist, order, tolerance)
        end        
        out_c = control[idx,:]
        out_c.distance = dist
        treated.distance = dist
        matched = vcat(out_c,treated)
    else
       data.exact_match .= string.([string.(data[:,i],"_") for i in exact]...)
        treated = @view data[isequal.(data[:,T],1),:]
        control = @view data[isequal.(data[:,T],0),:]
        matched = DataFrame()
        no_exact_match = 0
        exact_m = "1_black_"
        for exact_m  in unique(data.exact_match)
            c_df = @view control[isequal.(exact_m,control.exact_match),:]
            t_df = @view treated[isequal.(exact_m,treated.exact_match),:]
            if isempty(c_df) || isempty(t_df)
                no_exact_match +=1
                continue
            end
            
            if replacement
                idx, dist = _nearest_neighbor_matching(t_df.Dist,c_df.Dist)
            else
                idx, dist = _nearest_neighbor_matching(t_df.Dist,c_df.Dist, order, tolerance)
            end
            out_t = t_df[.!isequal.(idx,-1),:]            
            out_c = c_df[idx[idx.!=-1],:]
            out_c.distance = dist[idx.!=-1]
            out_t.distance = dist[idx.!=-1]
            append!(matched, out_c)
            append!(matched,out_t)
        end

        if isequal(no_exact_match, size(unique(data.exact_match),1))
            throw("Zero cases with exact matches!")
        elseif no_exact_match>0
            @warn "Could not find exact matches for $(no_exact_match) groups!"
        end
        select!(matched,Not(:exact_match))
        select!(data,Not(:exact_match))
    end
    matched = matched[matched.distance.<=maxDist,:]
    return matched
end;