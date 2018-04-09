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
#       This file contains ..., while ... contains the ...                      #
#                                                                               #
#   Last modified by C. Brenhin Keller 2018-03-19                               #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## --- Load external resources

    # # Basic statistics and UI resources
    # using StatsBase: fit, Histogram
    # using ProgressMeter: @showprogress
    # using LsqFit: curve_fit
    # using KernelDensity: kde
    #
    # # Weighted mean, etc
    # include("Utilities.jl");
    # # Functions for estimating extrema of a finite-range distribution
    # include("DistMetropolis.jl");
    # # Functions for stratigraphic modelling
    # include("StratMetropolis.jl");
    #
    # # Higher-level functions for fitting and plotting
    # using Plots, Colors;
    # gr();
    # include("Fitplot.jl");


    using Chron


## --- Define sample properties

# # # # # # # # # # # # Enter sample information here! # # # # # # # # # # # # #
    Name   =        ("KJ08-157", "KJ04-75", "KJ09-66", "KJ04-72", "KJ04-70");
    Height =        [     -52.0,      44.0,      54.0,      82.0,      93.0];
    Height_Sigma =  [       3.0,       1.0,       3.0,       3.0,       3.0];
    InputAgeSigmaLevel = 2; # i.e., are the data files 1-sigma or 2-sigma
    DataPath = "examples/DenverUPbExampleData/" # Where are the data files?
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Count the names to see how many samples we have
    nSamples = length(Name);
    # Make an instance of a ChronSection object with
    smpl = StratAgeData(
        Name,  # Sample names
        Height,  # Sample heights
        Height_Sigma,  # Height_sigma
        fill(NaN,nSamples),  # Sample ages
        fill(NaN,nSamples),  # Sample age uncertainty
        fill(NaN,nSamples),  # Sample age 2.5% CI
        fill(NaN,nSamples),  # Sample age 97.5% CI
        zeros(nSamples), # Sidedness (zeros by default, geochron constraints are two-sided)
        fill(NaN,nSamples,nSamples), # Sample age distribution parameters
        DataPath, # Relative path where we can find .csv data files
        InputAgeSigmaLevel,# i.e., are the data files 1-sigma or 2-sigma
    );

## --- Calculate bootstrapped distribution

    # Bootstrap a KDE of the pre-eruptive (or pre-deposition) zircon distribution
    # shape from individual sample datafiles using a KDE of stacked sample data
    BootstrappedDistribution = BootstrapDistributionKDE(smpl);
    plot(BootstrappedDistribution,xlabel="Time (unitless)",ylabel="")

## --- Estimate the eruption age distributions for each sample

# # # # # # # # # # # # Configure distribution model here! # # # # # # # # # # #
  # Number of steps to run in distribution MCMC
  distSteps = 10^6;
  distBurnin = floor(Int,distSteps/2);

  # Choose the form of the prior distribution to use
  # A variety of potentially useful distributions are provided in DistMetropolis.jl
  # Options include UniformDisribution, TriangularDistribution,
  # BootstrappedDistribution, and MeltsVolcanicZirconDistribution
  # or you can define your own.
  dist = BootstrappedDistribution;
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Run MCMC to estimate saturation and eruption/deposition age distributions
    smpl = tMinDistMetropolis(smpl,distSteps,dist);

    # using JLD: @save
    # @save "smpl.jld" smpl

## --- Run stratigraphic model

# # # # # # # # # # # Configure stratigraphic model here! # # # # # # # # # # #
# If you don't know what these do, you can probably leave them as-is
    resolution = 2.0; # Same units as sample height. Smaller is slower!
    (bottom, top) = extrema(smpl.Height);
    bounding = 0.1; # how far away do we place runaway bounds, as a fraction of total section height
    npoints_approx = round(Int,length(bottom:resolution:top) * (1+2*bounding))
    burnin = 20000*npoints_approx;
    nsteps = 25000;
    sieve = round(Int,npoints_approx);
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Pass configuration to a struct
    config = StratAgeModelConfiguration(
        resolution, # Resolution
        burnin, # Lenght of burnin
        nsteps, # Number of regular MC steps
        sieve, # Seive strength
        bounding, # Runaway bounds
    );

    # Run the model
    (mdl, agedist, lldist) = StratMetropolis(smpl, config);

    # Plot results (mean and 95% confidence interval for both model and data)
    hdl = plot([mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(minimum(mdl.Height),0.5,:blue), label="model")
    plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="")
    plot!(hdl, smpl.Age, smpl.Height, xerror=(smpl.Age-smpl.Age_025CI,smpl.Age_975CI-smpl.Age),label="data",seriestype=:scatter,color=:black)
    plot!(hdl, xlabel="Age (Ma)", ylabel="Height (cm)")
    savefig(hdl,"AgeDepthModel.pdf");

## --- Same but with hiata

    # # A type of object to hold data about hiatuses
    # hiatus = HiatusData(
    #     [20.0,], # Height
    #     [0.0,], # Height_Sigma
    #     [0.5,], # Duration
    #     [0.1,], # Duration_Sigma
    # );
    #
    # (mdl, agedist, hiatusdist, lldist) = StratMetropolisHiatus(smpl, hiatus, config);
    # smpl
