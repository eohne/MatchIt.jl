function match_table(d::PSM)
    df = select(d.df, findall(col -> all(v -> v isa Number, col), eachcol(d.df)))
    matched = select(d.matched, findall(col -> all(v -> v isa Number, col), eachcol(d.matched )))
    select!(matched, Not(:dist))
    before_not = Statistics.mean(Matrix(df[df[:,d.T].==0,:]),dims=1)'
    before_treated = Statistics.mean(Matrix(df[df[:,d.T].==1,:]),dims=1)'
    after_not = Statistics.mean(Matrix(matched[matched[:,d.T].==0,:]),dims=1)'
    after_treated = Statistics.mean(Matrix(matched[matched[:,d.T].==1,:]),dims=1)'
    n = names(matched)
    
    # Welch's T-test
    test_unmatched = []
    test_matched = []
    for i in 1:ncol(df)
        temp = missing
        temp2 = missing
        try
            temp = HypothesisTests.UnequalVarianceTTest(df[df[:,d.T].==1,i],df[df[:,d.T].==0,i]) |> HypothesisTests.pvalue
        catch e

        end
        try
            temp2 = HypothesisTests.UnequalVarianceTTest(matched[matched[:,d.T].==1,i],matched[matched[:,d.T].==0,i]) |> HypothesisTests.pvalue
        catch e
        end
        push!(test_unmatched,temp)
        push!(test_matched,temp2)
    end
    unmatched = hcat(n,before_not,before_treated,test_unmatched)
    matched = hcat(n,after_not,after_treated,test_matched)
    unmatched = DataFrame(unmatched, ["Var","Mean1","Mean2","P_Value"])
    matched = DataFrame(matched, ["Var","Mean1","Mean2","P_Value"])
    print("\nBefore Matching:\n ")
    show(unmatched)
    print("\n\nAfter Matching:\n")
    show(matched)
    return unmatched, matched
 end