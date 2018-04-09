__precompile__()

module Chron

    # Basic statistics and UI resources
    using StatsBase: fit, Histogram
    using ProgressMeter: @showprogress
    using LsqFit: curve_fit
    using KernelDensity: kde

    # Weighted mean, etc
    include("Utilities.jl");
    # Functions for estimating extrema of a finite-range distribution
    include("DistMetropolis.jl");
    # Functions for stratigraphic modelling
    include("StratMetropolis.jl");

    # Higher-level functions for fitting and plotting
    using Plots, Colors;
    gr();
    include("Fitplot.jl");

    # Structs
    export StratAgeData, HiatusData, StratAgeModelConfiguration, StratAgeModel

    # Functions
    export tMinDistMetropolis, BootstrapDistributionKDE, plotRankOrderErrorbar,
        StratMetropolis, StratMetropolisHiatus, checkDistLogLikelihood,
        checkCrystLogLikelihood, crystMinMaxMetropolis,
        gwmean, awmean, drawFromDistribution, fillFromDistribution,
        doubleLinearExponential, doubleLinearExponentialLL, cntr, pctile,
        findclosest, findclosestbelow, findclosestabove

    # Distribution
    export UniformDistribution, TriangularDistribution, HalfNormalDistribution,
        EllisDistribution, MeltsZirconDistribution, MeltsVolcanicZirconDistribution

end # module
