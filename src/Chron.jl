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

    # Weighted mean, etc
    include("Utilities.jl");
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
    export tMinDistMetropolis, BootstrapDistributionKDE, StratMetropolis,
        StratMetropolisHiatus, StratMetropolisDist, StratMetropolisDistHiatus,
        checkDistLogLikelihood, checkCrystLogLikelihood, crystMinMaxMetropolis,
        gwmean, awmean, drawFromDistribution, fillFromDistribution,
        doubleLinearExponential, doubleLinearExponentialLL, cntr, pctile,
        findclosest, findclosestbelow, findclosestabove, plotRankOrderErrorbar,
        tMinDistMetropolisLA, crystMinMaxMetropolisLA

    # Distributions
    export UniformDistribution, TriangularDistribution, HalfNormalDistribution,
        EllisDistribution, MeltsZirconDistribution, MeltsVolcanicZirconDistribution

end # module
