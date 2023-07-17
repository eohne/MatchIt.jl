using MatchIt
using DataFrames, PlotlyJS, CSV
using Test

data = CSV.File("../Example Data/lalonde.csv") |> DataFrame;
@testset "matchit" begin
    ta = matchit(data,@formula(treat ~ age + educ + race + married + nodegree + re74 + re75),exact = ["married"]);
    ta2 = matchit(data,@formula(treat ~ age + educ + race + married + nodegree + re74 + re75));
    ta3 = matchit(data,@formula(treat ~ age + educ + race + married + nodegree + re74 + re75),link=ProbitLink());
    ta4 = matchit(data,@formula(treat ~ age + educ + race + married + nodegree + re74 + re75),exact = ["married"],replacement=false);
    ta5 = matchit(data,@formula(treat ~ age + educ + race + married + nodegree + re74 + re75),exact = [],replacement=false);
    @test typeof(ta) <: MatchedIt
    @test nrow(ta.df) == nrow(data)
    @test nrow(ta.matched) >0

    @test typeof(ta2) <: MatchedIt
    @test nrow(ta2.df) == nrow(data)
    @test nrow(ta2.matched) >0
    
    @test typeof(ta3) <: MatchedIt
    @test nrow(ta3.df) == nrow(data)
    @test nrow(ta3.matched) >0
    
    @test typeof(ta4) <: MatchedIt
    @test nrow(ta4.df) == nrow(data)
    @test nrow(ta4.matched) >0
    
        
    @test typeof(ta5) <: MatchedIt
    @test nrow(ta5.df) == nrow(data)
    @test nrow(ta5.matched) >0

    s = summary(ta,true,true);
    @test isnothing(s)
    @test typeof(balance_plot(ta)) <: PlotlyJS.SyncPlot
end
