
"""
    get_propensities(data::DataFrame,T::String,X::Vector{String},link::GLM.Link01)

Calculates the Propensity Scores and returns a `DataFrame` with a new column `:PS`.

# Arguments:
  * data`::DataFrame`: The DataFrame containg the data. This should have no missing values. Rows with missing values will be deleted.
  * T`::String`: the name of the treatment variable
  * X`::Vector{String}`: the names of the covariates to be matched on
  * link`::GLM.Link01`: The link function used. One of `LogitLink()` or `ProbitLink()`

# Output
  * `DataFrame` containing the original data without missing values that has a new column (`:PS`) with the propensity scores. 

"""
function get_propensities(data::DataFrame,T::String,X::Vector{String},link::GLM.Link01)
    out = copy(data)
    if any([any(ismissing.(i)) for i in eachcol(data)])
        before = nrow(out)
        dropmissing!(out)
        after = nrow(out)
        @warn "Missing values in data. $(before-after) rows with missing values were dropped!"
    end
    
    if in("PS", names(out))
        @warn "PS variable already exists in DataFrame. The variable will be dropped and replaced with the calculated propensity score!"
    end
    
    fm = term(Symbol(T)) ~ sum(term.(Symbol.(X)))
    out.PS = predict(glm(fm, out, Binomial(), link ))
    return out
end

"""
    get_propensities!(data::DataFrame,T::String,X::Vector{String},link::GLM.Link01)

Calculates the Propensity Scores and returns a `DataFrame` with a new column `:PS`.

# Arguments:
  * data`::DataFrame`: The DataFrame containg the data. This should have no missing values. Rows with missing values will be deleted.
  * T`::String`: the name of the treatment variable
  * X`::Vector{String}`: the names of the covariates to be matched on
  * link`::GLM.Link01`: The link function used. One of `LogitLink()` or `ProbitLink()`

# Output
  * `DataFrame` containing the original data without missing values that has a new column (`:PS`) with the propensity scores. 

"""
function get_propensities!(data::DataFrame,T::String,X::Vector{String},link::GLM.Link01)
    if any([any(ismissing.(i)) for i in eachcol(data)])
        before = nrow(data)
        dropmissing!(data)
        after = nrow(data)
        @warn "Missing values in data. $(before-after) rows with missing values were dropped!"
    end
        
    if in("PS", names(data))
        @warn "PS variable already exists in DataFrame. The variable will be dropped and replaced with the calculated propensity score!"
    end

    fm = term(Symbol(T)) ~ sum(term.(Symbol.(X)))
    data.PS = predict(glm(fm, data, Binomial(), link ))
    return data::DataFrame
end
