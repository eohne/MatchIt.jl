using PSMatching
using DataFrames, GLM, PlotlyJS
using Test

data = DataFrame(y = Int.(round.(rand(100),digits=0)), x1 = rand(100),x2=rand(100),id=repeat(collect(1:10),10))
@testset "PSM" begin
    ta = PSM(data,"y",["x1","x2"],["id"],LogitLink(),Inf);
    @test typeof(ta) <: PSMatched
    @test nrow(ta.df) == nrow(data)
    @test nrow(ta.matched) >0
    ta2 = match_table(ta);
    @test typeof(ta2) <: Tuple{DataFrame,DataFrame}
    @test typeof(balance_plot(ta)) <: PlotlyJS.SyncPlot
end
