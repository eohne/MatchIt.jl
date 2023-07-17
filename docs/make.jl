push!(LOAD_PATH,"../src/")
using MatchIt
using Documenter
makedocs(
         sitename = "MatchIt.jl",
         format = Documenter.HTML(
            # analytics = "G-LFRFQ0X1VF",
         canonical = "https://eohne.github.io/Matchit.jl/dev/"),
         modules  = [MatchIt],
         pages=[
                "Home" => "index.md",
                "Function Documentation" =>[
                    "Matching" =>"Matching.md",
                    "All Functions" =>"AllFunctions.md",
                ],
                "Version Change Log" => "VersionChanges.md"
               ])
deploydocs(;
    repo="github.com/eohne/MatchIt.jl",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#", "dev" => "dev"]
)