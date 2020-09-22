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

    # Backwards compatibility
    using Compat
    # Forwards compatibility
    VERSION >= v"0.7" && using Statistics, DelimitedFiles, SpecialFunctions

    # AVX vectorziation tools
    using LoopVectorization
    using SIMDPirates: vifelse, verf
    using VectorizationBase: SVec
    import SpecialFunctions.erf
    erf(x::SVec) = verf(x)

    # Basic statistics and UI resources
    using StatsBase: fit, Histogram, percentile
    using ProgressMeter: @showprogress, Progress, update!
    using LsqFit: curve_fit
    using KernelDensity: kde
    using Interpolations

    # Weighted mean, systematic errors, etc
    include("Utilities.jl")
    # Functions for propagating systematic uncertainties
    using Distributions
    include("Systematic.jl")
    # Custom objects for holding Chron age data
    include("Objects.jl")
    # Intcal2013 calibration curve for radiocarbion
    include("Intcal.jl")
    # Functions for estimating extrema of a finite-range distribution
    include("DistMetropolis.jl")
    # Functions for stratigraphic modelling
    include("StratMetropolis.jl")

    # Higher-level functions for fitting and plotting
    using Colors: Colorant, ColorTypes, RGBX, RGB, N0f8
    include("Colormaps.jl")

    using Plots; gr();
    include("Fitplot.jl")

    # Structs
    export ChronAgeData, HiatusData, StratAgeModelConfiguration, StratAgeModel

    # High-level functions
    export  StratMetropolis, StratMetropolisDist, StratMetropolis14C,
        tMinDistMetropolis, metropolis_min!, metropolis_min,
        metropolis_minmax!, metropolis_minmax,
        bilinear_exponential, bilinear_exponential_ll,
        plot_rankorder_errorbar, plot_rankorder_errorbar!,
        screen_outliers, BootstrapCrystDistributionKDE,
        BootstrapCrystDistributionKDEfromStrat


    # Utility functions
    export nanminimum, nanmaximum, nanrange, nanextrema, pctile, nanmedian, nanmean, nanstd,
        linsp, linterp1s, linterp1, cntr, gwmean, awmean,
        normpdf, normpdf_ll, normcdf, norm_quantile, norm_width,
        findclosest, findclosestbelow, findclosestabove,
        draw_from_distribution, draw_from_distribution!,
        linterp_at_index, interpolate_ll, normproduct, normproduct_ll

    # Colormaps
    export viridis, inferno, plasma, magma, fire, water

    # Distributions
    export UniformDistribution, TriangularDistribution, HalfNormalDistribution, TruncatedNormalDistribution,
        EllisDistribution, MeltsZirconDistribution, MeltsVolcanicZirconDistribution

end # module
