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
#    The second model uses the estimated (posterior) eruption/deposition age    #
#  distributions along with the constraint of stratigraphic superposition to    #
#  produce an age-depth model.                                                  #
#                                                                               #
#    DOI: 10.17605/OSF.IO/TQX3F                                                 #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

module Chron

    using Reexport
    @reexport using NaNStatistics
    @reexport using StatGeochemBase
    @reexport using Isoplot
    using Isoplot: val, err

    # Basic statistics and UI resources
    using ProgressMeter: @showprogress, Progress, update!
    using LsqFit: curve_fit
    using KernelDensity: kde
    using DelimitedFiles
    using Distributions
    using QuadGK

    # Custom objects for holding Chron age data
    include("Objects.jl")
    # General utility functions
    include("Utilities.jl")
    # Functions for propagating systematic uncertainties
    include("Systematic.jl")
    # Intcal2013 calibration curve for radiocarbion
    include("Intcal.jl")
    # Functions for estimating extrema of a finite-range distribution
    include("DistMetropolis.jl")
    # Functions for stratigraphic modelling
    include("StratMetropolis.jl")

    using Plots
    include("Fitplot.jl")

    # Structs
    export ChronAgeData, GeneralAgeData, HiatusData
    export SystematicUncertainty, StratAgeModelConfiguration, StratAgeModel

    # High-level functions
    export  StratMetropolis, StratMetropolisDist, StratMetropolis14C,
        tMinDistMetropolis, metropolis_min!, metropolis_min,
        metropolis_minmax!, metropolis_minmax,
        screen_outliers, BootstrapCrystDistributionKDE

    # Custom distributions and related functions
    export strat_ll, BilinearExponential, Radiocarbon
    
    # Dealing with systematic uncertainty
    export add_systematic_uncert_UPb, add_systematic_uncert_ArAr, add_systematic_uncert_UTh, SystematicUncertainty

    # Relative crystallization/closure distributions
    # (mostly already exported by reexporting Isoplot.jl)
    export TruncatedNormalDistribution, EllisDistribution, ArClosureDistribution

    include("Deprecations.jl")

end # module
