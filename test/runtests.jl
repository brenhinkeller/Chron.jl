using Chron
using Test, Statistics, Distributions

const make_plots = get(ENV, "MAKE_PLOTS", false) == "true"
@testset "Utilities" begin include("testUtilities.jl") end
@testset "Eruption / deposition age distributions" begin include("testDist.jl") end
@testset "Strat only" begin 
    include("testStratOnly.jl") 
    include("testStratOnlyGeneral.jl") 
end
@testset "Radiocarbon" begin include("testRadiocarbon.jl") end
@testset "Coupled model" begin include("testCoupled.jl") end
