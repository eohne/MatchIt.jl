using PSMatching
using DataFrames, GLM, PlotlyJS
using Test

data = DataFrame(y = Int.(round.(rand(100),digits=0)), x1 = rand(100),x2=rand(100),id=repeat(collect(1:10),10))
@testset "matchit" begin
    ta = matchit(data,@formula(y + x1 + x2),exact = ["id"]);
    @test typeof(ta) <: MatchedIt
    @test nrow(ta.df) == nrow(data)
    @test nrow(ta.matched) >0
    ta2 = match_table(ta);
    @test typeof(ta2) <: Tuple{DataFrame,DataFrame}
    @test typeof(balance_plot(ta)) <: PlotlyJS.SyncPlot
end
