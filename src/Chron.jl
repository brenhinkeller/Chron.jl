# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                   Chron.jl                                    #
#                                                                               #
#       A two-part framework for (1) estimating eruption/deposition age         #
#  distributions from complex mineral age spectra and (2) subsequently building #
#  a stratigraphic age model based on those distributions. Each step relies on  #
#  a Markov-Chain Monte Carlo model.                                            #
#                                                                               #
#    The first model uses an informative prior distribution to estimate the     #
#  times of first (i.e., saturation) and last  mineral crystallization (i.e.,   #
#  eruption/deposition).                                                        #
#                                                                               #
#    The second model uses the estimated (posterior) eruption/deposition ages   #
#  distributions along with the constraint of stratigraphic superposition to    #
#  produce an age-depth model                                                   #
#                                                                               #
#   Last modified by C. Brenhin Keller 2018-04-09                               #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

__precompile__()

module Chron

    # Basic statistics and UI resources
    using StatsBase: fit, Histogram, percentile
    using ProgressMeter: @showprogress
    using LsqFit: curve_fit
    using KernelDensity: kde
    using Interpolations: interpolate, Gridded, Linear
    using IndirectArrays: IndirectArray

    # Weighted mean, etc
    include("Utilities.jl");
    include("Colormaps.jl")
    # Functions for estimating extrema of a finite-range distribution
    include("DistMetropolis.jl");
    # Functions for stratigraphic modelling
    include("StratMetropolis.jl");

    # Higher-level functions for fitting and plotting
    using Plots
    include("Fitplot.jl");

    # Structs
    export StratAgeData, HiatusData, StratAgeModelConfiguration, StratAgeModel

    # Functions
    export tMinDistMetropolis, BootstrapDistributionKDEfromStrat, StratMetropolis,
        StratMetropolisHiatus, StratMetropolisDist, StratMetropolisDistHiatus,
        checkDistLogLikelihood, checkCrystLogLikelihood, crystMinMaxMetropolis,
        drawFromDistribution, fillFromDistribution, doubleLinearExponential,
        doubleLinearExponentialLL, plotRankOrderErrorbar, tMinDistMetropolisLA,
        crystMinMaxMetropolisLA, BootstrapDistributionKDE

    # Utility functions
    export nanmin, nanmax, nanrange, pctile, nanmedian, nanmean, nanstd,
        linterp1s, linterp1, cntr, gwmean, awmean,
        findclosest, findclosestbelow, findclosestabove

    # Colormaps
    export viridis, inferno, plasma, magma, fire

    # Distributions
    export UniformDistribution, TriangularDistribution, HalfNormalDistribution,
        EllisDistribution, MeltsZirconDistribution, MeltsVolcanicZirconDistribution

end # module
