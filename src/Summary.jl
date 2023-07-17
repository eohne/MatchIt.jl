function _round_floats(x::Float64,digits::Int)
    return round(x,digits = digits)
end
function _round_floats(x::Any,digits::Int)
    return x
end

function _find_length(m::Matrix)
    s = Int.(zeros(size(m,2)))
    for c in 1:size(m,2)
        for r in 1:size(m,1)
            l = length(string(m[r,c]))
            if l > s[c]
                s[c] = l
            end
        end
    end
    return s
end
function _find_length(v::Vector)
    return length.(string.(v))
end


function _print_balance_data(x::DataFrame)
    m = _round_floats.(Matrix(x),4)
    s_l = _find_length(m)
    n_l = _find_length(names(x))
    s_l = max.(s_l,n_l)
    # "│" add a horizontal line 
    for i in 1:size(x,2)
        print("$(rpad(names(x)[i],s_l[i]))")
        if i ==1
            print("│")
        end
        if i < length(names(x))
            print("    ")
        end
    end
    print("\n")
    print(join(["─" for i in 1:(s_l[1])]),"│",join(["─" for i in 1:4]))
    print(join(["─" for i in 1:(sum(s_l[2:end])+(4*(size(x,2)-2)))]),"\n")    
    for r in 1:size(x,1)
        for c in 1:size(x,2)
            print("$(rpad(string(m[r,c]),s_l[c]))")
            if c ==1
                print("│")
            end
            if c < size(x,2)
                print("    ") 
            end
        end
        print("\n")
    end
end

function _print_balance_observations(m::Matrix)
    m = _round_floats.(m,0)
    s_l = _find_length(m)
    n_l = _find_length([" ","All","Matched","Unmatched"])
    s_l = max.(s_l,n_l)
    # "│" add a horizontal line 
    for i in 1:size(m,2)
        print("$(rpad([" ","All    ","Matched    ","Unmatched"][i],s_l[i]))")
        if i ==1
            print("│")
        end
        if i < 2
            print("    ")
        end
    end
    print("\n")
    print(join(["─" for i in 1:(s_l[1])]),"│",join(["─" for i in 1:4]))
    print(join(["─" for i in 1:(sum(s_l[2:end])+(4*(size(m,2)-2)))]),"\n")    
    for r in 1:size(m,1)
        for c in 1:size(m,2)
            print("$(rpad(string(m[r,c]),s_l[c]))")
            if c ==1
                print("│")
            end
            if c < size(m,2)
                print("    ") 
            end
        end
        print("\n")
    end
end


function _print_summary(unmatched,matched,pre,n_treat,n_control,n_treat_matched,n_control_matched)
    print("\nMatch Summary\n")
    print(join(["═" for i in 1:45]),"\n\n")
    print("Sampel Sizes:\n")
    obs = ["Treat" n_treat n_treat_matched (n_treat - n_treat_matched); "Control" n_control n_control_matched (n_control - n_control_matched)]
    _print_balance_observations(obs)
    print("\n\n")

    if pre
        print("\nSummary of Balance for All Data:\n")
        _print_balance_data(unmatched)
        print("\n\n")
    end

    print("Summary of Balance for Matched Data:\n")
    _print_balance_data(matched)
end



 """
    Base.summary(obj::MatchedIt,test::Bool=false,pre::Bool=false)

Gives a summary output for the matched sample. The output includes a summary of the number of observations matched and the mean of all variables in the original and matched dataframe of the treatment and control group. If `test=true` a p-value of a 2 sample Welch test is reported.

# Arguments:
  * obj::`MatchedIt`: The output from a call to `matchit`
  * test`::Bool`: Whether a difference in mean test (2 sample Welch test) should be performed between treatment and control observations. (defaults to `false`)
  * pre`::Bool`: Whether to also show the output of the means and T-tests for the sample before matching.
  
"""
function Base.summary(obj::MatchedIt,test::Bool=false,pre::Bool=false)
    n_treat = sum(obj.df[:,obj.T])
    n_control = sum(obj.df[:,obj.T].==0)
    n_treat_matched = sum(obj.matched[:,obj.T])
    n_control_matched = sum(obj.matched[:,obj.T].==0)
    matched = select(obj.matched, findall(col -> all(v -> v isa Number, col), eachcol(obj.matched )))
    select!(matched, Not(:dist))
    after_not = Statistics.mean(Matrix(matched[matched[:,obj.T].==0,:]),dims=1)'
    after_treated = Statistics.mean(Matrix(matched[matched[:,obj.T].==1,:]),dims=1)'
    n = names(matched)

    if pre
        df = select(obj.df, findall(col -> all(v -> v isa Number, col), eachcol(obj.df)))
        before_not = Statistics.mean(Matrix(df[df[:,obj.T].==0,:]),dims=1)'
        before_treated = Statistics.mean(Matrix(df[df[:,obj.T].==1,:]),dims=1)'
    end
 
    if test   
        # Welch's T-test
        if pre
            test_unmatched = []
            test_matched = []
            for i in 1:ncol(df)
                temp = missing
                temp2 = missing
                try
                    temp = HypothesisTests.UnequalVarianceTTest(df[df[:,obj.T].==1,i],df[df[:,obj.T].==0,i]) |> HypothesisTests.pvalue
                catch e

                end
                try
                    temp2 = HypothesisTests.UnequalVarianceTTest(matched[matched[:,obj.T].==1,i],matched[matched[:,obj.T].==0,i]) |> HypothesisTests.pvalue
                catch e
                end
                push!(test_unmatched,temp)
                push!(test_matched,temp2)
            end
            unmatched = hcat(n,before_not,before_treated,test_unmatched)
            matched = hcat(n,after_not,after_treated,test_matched)
            unmatched = DataFrame(unmatched, ["Var","Mean1","Mean2","P_Value"])
            matched = DataFrame(matched, ["Var","Mean1","Mean2","P_Value"])
        else
            test_matched = []
            for i in 1:ncol(matched)
                 temp2 = missing
                try
                    temp2 = HypothesisTests.UnequalVarianceTTest(matched[matched[:,obj.T].==1,i],matched[matched[:,obj.T].==0,i]) |> HypothesisTests.pvalue
                catch e
                end
                 push!(test_matched,temp2)
            end

            matched = hcat(n,after_not,after_treated,test_matched)
            matched = DataFrame(matched, ["Var","Mean1","Mean2","P_Value"])
            unmatched = []
        end
    else
        if pre
            unmatched = hcat(n,before_not,before_treated)
            matched = hcat(n,after_not,after_treated)
            unmatched = DataFrame(unmatched, ["Var","Mean1","Mean2"])
            matched = DataFrame(matched, ["Var","Mean1","Mean2"])
        else
            matched = hcat(n,after_not,after_treated)
            matched = DataFrame(matched, ["Var","Mean1","Mean2"])
            unmatched = []
        end
    end
    _print_summary(unmatched,matched,pre,n_treat,n_control,n_treat_matched,n_control_matched)
 end