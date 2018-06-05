# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                 examples.jl                                   #
#                                                                               #
#     Illustrates the use of the Chron.jl package for eruption/deposition       #
#  age estimation and production of a stratigraphic age model.                  #
#                                                                               #
#     This file uses code cells (denoted by "## ---"). To evaluate a cell in    #
#  the Julia REPL and move to the next cell, the default shortcut in Atom is    #
#  alt-shift-enter.                                                             #
#                                                                               #
#      You may have to adjust the path below which specifies the location of    #
#  the CSV data files for each sample (from examples/DenverUPbExampleData/      #
#  in the Chron.jl repository).                                                 #
#                                                                               #
#   Last modified by C. Brenhin Keller 2018-04-09                               #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## --- Load the Chron package

    using Plots; gr();
    using Chron

## --- Define sample properties

# # # # # # # # # # # # Enter sample information here! # # # # # # # # # # # # #
    # This example data is from Clyde et al. (2016) "Direct high-precision
    # U–Pb geochronology of the end-Cretaceous extinction and calibration of
    # Paleocene astronomical timescales" EPSL 452, 272–280.
    # doi: 10.1016/j.epsl.2016.07.041
    Name   =        ("KJ08-157", "KJ04-75", "KJ09-66", "KJ04-72", "KJ04-70",);
    Height =        [     -52.0,      44.0,      54.0,      82.0,      93.0,];
    Height_Sigma =  [       3.0,       1.0,       3.0,       3.0,       3.0,];
    inputSigmaLevel = 2; # i.e., are the data files 1-sigma or 2-sigma. Integer.
    Path = "examples/DenverUPbExampleData/" # Where are the data files?
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
        fill(NaN,5,nSamples), # Sample age distribution parameters
        Path, # Relative path where we can find .csv data files
        inputSigmaLevel,# i.e., are the data files 1-sigma or 2-sigma
    );

## --- (Optional) Calculate bootstrapped distribution

    # # Bootstrap a KDE of the pre-eruptive (or pre-deposition) zircon distribution
    # # shape from individual sample datafiles using a KDE of stacked sample data
    # BootstrappedDistribution = BootstrapDistributionKDEfromStrat(smpl);
    # plot(BootstrappedDistribution,xlabel="Time (arbitrary units)",ylabel="Probability Density")

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
  dist = TriangularDistribution;
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Run MCMC to estimate saturation and eruption/deposition age distributions
    smpl = tMinDistMetropolis(smpl,distSteps,distBurnin,dist);

    # # (Optional) Save the sample struct for later use
    # using JLD: @save, @load
    # @save "smpl.jld" smpl

## --- Run stratigraphic model

    # # (Optional) Load the saved sample struct
    # @load "smpl.jld" smpl

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

    # Run the stratigraphic MCMC model
    (mdl, agedist, lldist) = StratMetropolisDist(smpl, config);

    # Plot results (mean and 95% confidence interval for both model and data)
    hdl = plot([mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(minimum(mdl.Height),0.5,:blue), label="model")
    plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="")
    plot!(hdl, smpl.Age, smpl.Height, xerror=(smpl.Age-smpl.Age_025CI,smpl.Age_975CI-smpl.Age),label="data",seriestype=:scatter,color=:black)
    plot!(hdl, xlabel="Age (Ma)", ylabel="Height (cm)")
    savefig(hdl,"AgeDepthModel.pdf");
    display(hdl)

    # Interpolate results at KTB (height = 0)
    KTB = linterp1s(mdl.Height,mdl.Age,0)
    KTB_min = linterp1s(mdl.Height,mdl.Age_025CI,0)
    KTB_max = linterp1s(mdl.Height,mdl.Age_975CI,0)
    @printf("Interpolated age: %0.3f +%0.3f/-%0.3f Ma", KTB, KTB_max-KTB, KTB-KTB_min)


## --- Calculate deposition rate binned by age

    # Set bin width and spacing
    binwidth = 0.01; # Myr
    binoverlap = 10;
    ages = minimum(mdl.Age):binwidth/binoverlap:maximum(mdl.Age);
    spacing = binoverlap;

    # Calculate rates for the stratigraphy of each markov chain step
    dhdt_dist = Array{Float64}(length(ages)-spacing,nsteps);
    for i=1:nsteps
        heights = linterp1(reverse(agedist[:,i]),reverse(mdl.Height),ages);
        dhdt_dist[:,i] = (heights[spacing+1:end]-heights[1:end-spacing])./(ages[1:end-spacing]-ages[spacing+1:end]);
    end

    # Find mean and 1-sigma (68%) CI
    dhdt = nanmean(dhdt_dist,dim=2);
    dhdt_50p = nanmedian(dhdt_dist,dim=2);
    dhdt_16p = pctile(dhdt_dist,15.865,dim=2); # Minus 1-sigma (15.865th percentile)
    dhdt_84p = pctile(dhdt_dist,84.135,dim=2); # Plus 1-sigma (84.135th percentile)
    bincenters = ages[1+Int(binoverlap/2):end-Int(binoverlap/2)]

    # Plot results
    hdl = plot(bincenters,dhdt, label="Mean", color=:black, linewidth=2)
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_16p; reverse(dhdt_84p)], fill=(minimum(mdl.Height),0.4,:blue), linecolor=:white, label="68% CI")
    plot!(hdl,bincenters,dhdt_50p, label="Median", color=:grey, linewidth=1)
    plot!(hdl, ylabel = "Depositional Rate (cm / Myr over $(Int(binwidth*1000)) kyr)", xlabel="Age (Ma)", fg_color_legend=:white)
    savefig(hdl,"DepositionRateModel.pdf");
    display(hdl)

## --- (Optional) If your section has hiata / exposure surfaces of known duration, try this:

    # # A type of object to hold data about hiatuses
    # hiatus = HiatusData(
    #     [20.0,35.0], # Height
    #     [0.0,0.0], # Height_Sigma
    #     [0.2,0.3], # Duration -- this is a minimum duration only
    #     [0.05,0.05], # Duration_Sigma (one-sided, minimum only)
    # );
    #
    # # Run the model
    # (mdl, agedist, hiatusdist, lldist) = StratMetropolisDistHiatus(smpl, hiatus, config);
    #
    # # Plot results (mean and 95% confidence interval for both model and data)
    # hdl = plot([mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(minimum(mdl.Height),0.5,:blue), label="model")
    # plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="")
    # plot!(hdl, smpl.Age, smpl.Height, xerror=(smpl.Age-smpl.Age_025CI,smpl.Age_975CI-smpl.Age),label="data",seriestype=:scatter,color=:black)
    # plot!(hdl, xlabel="Age (Ma)", ylabel="Height (cm)")
