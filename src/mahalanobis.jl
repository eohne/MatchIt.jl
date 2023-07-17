function _mahalanobis_replace(treat::Matrix,control::Matrix)
    s_c = size(control,1)
    s_t = size(treat,1)
    min_dist, min_idx = Vector{Float64}(undef,s_t), Vector{Int}(undef,s_t) 
    covariance = cov(control)
    @inbounds for i in 1:s_t
        d = Inf
        idx = 0 
        for j in 1:s_c
            temp = Mahalanobis(covariance,skipchecks=true)(view(treat,i,:),view(control,j,:))
            if iszero(temp)
                d = temp
                idx = j
                break
            elseif temp < d
                d = temp
                idx = j 
            else
                continue
            end
        end
        min_idx[i] = idx
        min_dist[i] = d
    end 
    return min_dist, min_idx
end

function _mahalanobis_no_replace(treat::Matrix,control::Matrix,order::String)
    
    s_c = size(control,1)
    s_t = size(treat,1)
    min_dist, min_idx = Vector{Float64}(undef,s_t), Vector{Int}(undef,s_t) 
    covariance = cov(control)
    used_c = Int.(zeros(s_c))
    if isequal(order, "data")
        sort_idx = collect(1:s_t)
    elseif isequal(order,"random")
        sort_idx = randperm(s_t)
    else
        throw("""order: $order not supported when matching on Mahalanobis Distance. Use one of: data or random""")
    end
    @inbounds for i in 1:s_t
        d = Inf
        idx = 0 
        for j in 1:s_c
            if used_c[j]==1
                continue
            end
            temp = Mahalanobis(covariance,skipchecks=true)(view(treat,i,:),view(control,j,:))
            if iszero(temp)
                d = temp
                idx = j
                break
            elseif temp < d
                d = temp
                idx = j 
            else
                continue
            end
        end
        min_idx[sort_idx[i]] = idx
        used_c[idx] = 1
        min_dist[sort_idx[i]] = d
    end 
    return min_dist, min_idx
end

function mhn_matching(df::DataFrame, f::FormulaTerm, exact::AbstractArray,replacement::Bool=true,order::String="data")
    if isempty(exact)
        treat = df[isequal.(df[:,StatsModels.termvars(f.lhs)[1]], 1),StatsModels.termvars(f.rhs)] 
        control = df[isequal.(df[:,StatsModels.termvars(f.lhs)[1]], 0),StatsModels.termvars(f.rhs)] 
        if replacement
            dist, idx = _mahalanobis_replace(Matrix(treat),Matrix(control))
        else
            dist, idx = _mahalanobis_no_replace(Matrix(treat),Matrix(control),order)
        end
    else
        df.exact_match .= string.([string.(df[:,i],"_") for i in exact]...)
        treat = view(df,isequal.(df[:,StatsModels.termvars(f.lhs)[1]], 1),:) 
        control = view(df, isequal.(df[:,StatsModels.termvars(f.lhs)[1]], 0),:)
        dist = Vector{Float64}(undef,size(treat,1))
        idx = Vector{Int}(undef,size(treat,1))
        idx .=-1
        for i in unique(df.exact_match)
            idx_t  = findall( isequal.(i,treat.exact_match))
            idx_c  = findall( isequal.(i,control.exact_match))
            s_treat = view(treat,isequal.(treat.exact_match,i),StatsModels.termvars(f.rhs)) |> Matrix
            s_control = view(control,isequal.(control.exact_match,i),StatsModels.termvars(f.rhs)) |> Matrix
            if replacement
                temp_dist, temp_idx = _mahalanobis_replace(s_treat,s_control)
            else
                temp_dist, temp_idx = _mahalanobis_no_replace(s_treat,s_control,order)
            end
            dist[idx_t] = temp_dist
            idx[idx_t] =  idx_c[temp_idx]
        end
    end
    treat[!,StatsModels.termvars(f.lhs)[1]] .=1
    matched_treated = treat[.!isequal.(idx,-1),:]
    matched_treated.Dist = dist[.!isequal.(idx,-1)]
    control[!,StatsModels.termvars(f.lhs)[1]] .=0
    matched_control = control[idx[.!isequal.(idx,-1)],:]
    matched_control.Dist = dist[.!isequal.(idx,-1)]
    matched = vcat(matched_treated,matched_control)
    matched = matched[:, vcat(StatsModels.termvars(f),:Dist)]
    return matched
end


