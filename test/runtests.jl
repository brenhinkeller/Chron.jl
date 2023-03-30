using Chron
using Test, Statistics, StatsBase

const make_plots = get(ENV, "MAKE_PLOTS", false)
@testset "Utilities" begin include("testUtilities.jl") end
@testset "Eruption / deposition age distributions" begin include("testDist.jl") end
@testset "Strat only" begin include("testStratOnly.jl") end
@testset "Radiocarbon" begin include("testRadiocarbon.jl") end
@testset "Coupled model" begin include("testCoupled.jl") end
