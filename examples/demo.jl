# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                   demo.jl                                     #
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

    using Chron
    using Plots; gr(); default(fmt = :svg);

    if VERSION>=v"0.7"
        using Statistics
        using DelimitedFiles
        using SpecialFunctions
    end

## --- Define sample properties
# # # # # # # # # # # # Enter sample information here! # # # # # # # # # # # # #
    # This example data is from Clyde et al. (2016) "Direct high-precision
    # U–Pb geochronology of the end-Cretaceous extinction and calibration of
    # Paleocene astronomical timescales" EPSL 452, 272–280.
    # doi: 10.1016/j.epsl.2016.07.041
    Name   =        ("KJ08-157", "KJ04-75", "KJ09-66", "KJ04-72", "KJ04-70",);
    Height =        [     -52.0,      44.0,      54.0,      82.0,      93.0,];
    Height_Sigma =  [       3.0,       1.0,       3.0,       3.0,       3.0,];
    Path = "examples/DenverUPbExampleData/" # Where are the data files?
    inputSigmaLevel = 2; # i.e., are the data files 1-sigma or 2-sigma. Integer.

    AgeUnit = "Ma" # Unit of measurement for ages and errors in the data files
    HeightUnit = "cm" # Unit of measurement for Height and Height_Sigma

    # Make an instance of a Chron StratAgeData object with the data we've just enterd
    smpl = NewStratAgeData(Name,Path,inputSigmaLevel)
    smpl.Height[:] = Height
    smpl.Height_Sigma[:] = Height_Sigma
    smpl.Age_Sidedness[:] = zeros(length(Name)) # Sidedness (zeros by default, geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## --- (Optional) Calculate bootstrapped distribution

    # Bootstrap a KDE of the pre-eruptive (or pre-deposition) zircon distribution
    # shape from individual sample datafiles using a KDE of stacked sample data
    BootstrappedDistribution = BootstrapCrystDistributionKDEfromStrat(smpl);
    plot(BootstrappedDistribution,xlabel="Time (arbitrary units)",ylabel="Probability Density",label="Bootstrapped distribution",fg_color_legend=:white)

## --- Estimate the eruption age distributions for each sample

# # # # # # # # # # # # Configure distribution model here! # # # # # # # # # # #
  # Number of steps to run in distribution MCMC
  distSteps = 10^6;
  distBurnin = floor(Int,distSteps/100);

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
# If you in doubt, you can probably leave these parameters as-is
    resolution = 1. # Same units as sample height. Smaller is slower!
    (bottom, top) = extrema(smpl.Height);
    bounding = 0.1; # how far away do we place runaway bounds, as a fraction of total section height
    npoints_approx = round(Int,length(bottom:resolution:top) * (1+2*bounding))
    burnin = 10000*npoints_approx;
    nsteps = 15000;
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

    # # Youngest Zircon
    # for i=1:length(smpl.Name)
    #     data = readdlm(string(smpl.Path, smpl.Name[i], ".csv"),',')
    #     Iyz = indmin(data[:,1]);
    #     smpl.Age[i] = minimum(data[Iyz,1]);
    #     smpl.Age_Sigma[i] = minimum(data[Iyz,2]/smpl.inputSigmaLevel);
    # end
    # (mdl, agedist, lldist) = StratMetropolis(smpl, config);

    # #LNWM
    # for i=1:length(smpl.Name)
    #     data = readdlm(string(smpl.Path, smpl.Name[i], ".csv"),',')
    #     sI = sortperm(data[:,1])
    #     # Weighted mean of youngst 3 zircons per sample (assuming there are at least 3 zircons in sample)
    #     Ns = min(size(data,1),3)
    #     (mu, sigma) = awmean(data[sI[1:Ns],1],data[sI[1:Ns],2]./smpl.inputSigmaLevel)
    #     smpl.Age[i] = mu;
    #     smpl.Age_Sigma[i] = sigma;
    # end
    # (mdl, agedist, lldist) = StratMetropolis(smpl, config);


## --- Plot strat model

    # Plot results (mean and 95% confidence interval for both model and data)
    hdl = plot([mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(minimum(mdl.Height),0.5,:blue), label="model")
    plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="")
    plot!(hdl, smpl.Age, smpl.Height, xerror=(smpl.Age-smpl.Age_025CI,smpl.Age_975CI-smpl.Age),label="data",seriestype=:scatter,color=:black)
    plot!(hdl, xlabel="Age ($AgeUnit)", ylabel="Height ($HeightUnit)")
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
    bincenters = ages[1+Int(binoverlap/2):end-Int(binoverlap/2)]
    spacing = binoverlap;

    # Calculate rates for the stratigraphy of each markov chain step
    dhdt_dist = Array{Float64}(undef,length(ages)-binoverlap,nsteps);
    @time for i=1:nsteps
        heights = linterp1(reverse(agedist[:,i]),reverse(mdl.Height),ages);
        dhdt_dist[:,i] = abs.(heights[1:end-spacing] - heights[spacing+1:end]) ./ binwidth;
    end

    # # Exact (added precision is below sampling resolution, so not useful) and very slow
    # @showprogress "Calculating dh/dt..." for i=1:nsteps
    #     for j=1:length(bincenters)
    #         t = (agedist[:,i] .> bincenters[j] - binwidth/2) .& (agedist[:,i] .< bincenters[j] + binwidth/2)
    #         dhdt_dist[j,i] = sum(t) * resolution / binwidth;
    #     end
    # end

    # Find mean and 1-sigma (68%) CI
    dhdt = nanmean(dhdt_dist,dim=2);
    dhdt_50p = nanmedian(dhdt_dist,dim=2);
    dhdt_16p = pctile(dhdt_dist,15.865,dim=2); # Minus 1-sigma (15.865th percentile)
    dhdt_84p = pctile(dhdt_dist,84.135,dim=2); # Plus 1-sigma (84.135th percentile)

    # Plot results
    hdl = plot(bincenters,dhdt, label="Mean", color=:black, linewidth=2)
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_16p; reverse(dhdt_84p)], fill=(minimum(mdl.Height),0.4,:darkred), linealpha=0, label="68% CI")
    plot!(hdl,bincenters,dhdt_50p, label="Median", color=:grey, linewidth=1)
    plot!(hdl, xlabel="Age ($AgeUnit)", ylabel="Depositional Rate ($HeightUnit / $AgeUnit over $binwidth $AgeUnit)", fg_color_legend=:white)
    ylims!(hdl, 0, 2500)
    savefig(hdl,"DepositionRateModel.pdf");
    display(hdl)

## --- Multiple confidence intervals (every 10)

    dhdt_20p = pctile(dhdt_dist,20,dim=2);
    dhdt_80p = pctile(dhdt_dist,80,dim=2);
    dhdt_25p = pctile(dhdt_dist,25,dim=2);
    dhdt_75p = pctile(dhdt_dist,75,dim=2);
    dhdt_30p = pctile(dhdt_dist,30,dim=2);
    dhdt_70p = pctile(dhdt_dist,70,dim=2);
    dhdt_35p = pctile(dhdt_dist,35,dim=2);
    dhdt_65p = pctile(dhdt_dist,65,dim=2);
    dhdt_40p = pctile(dhdt_dist,40,dim=2);
    dhdt_60p = pctile(dhdt_dist,60,dim=2);
    dhdt_45p = pctile(dhdt_dist,45,dim=2);
    dhdt_55p = pctile(dhdt_dist,55,dim=2);

    # Plot results
    hdl = plot(bincenters,dhdt, label="Mean", color=:black, linewidth=2)
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_16p; reverse(dhdt_84p)], fill=(minimum(mdl.Height),0.2,:darkred), linealpha=0, label="68% CI")
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_20p; reverse(dhdt_80p)], fill=(minimum(mdl.Height),0.2,:darkred), linealpha=0, label="")
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_25p; reverse(dhdt_75p)], fill=(minimum(mdl.Height),0.2,:darkred), linealpha=0, label="")
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_30p; reverse(dhdt_70p)], fill=(minimum(mdl.Height),0.2,:darkred), linealpha=0, label="")
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_35p; reverse(dhdt_65p)], fill=(minimum(mdl.Height),0.2,:darkred), linealpha=0, label="")
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_40p; reverse(dhdt_60p)], fill=(minimum(mdl.Height),0.2,:darkred), linealpha=0, label="")
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_45p; reverse(dhdt_55p)], fill=(minimum(mdl.Height),0.2,:darkred), linealpha=0, label="")
    plot!(hdl,bincenters,dhdt_50p, label="Median", color=:grey, linewidth=1)
    plot!(hdl, xlabel="Age ($AgeUnit)", ylabel="Depositional Rate ($HeightUnit / $AgeUnit over $binwidth $AgeUnit)", fg_color_legend=:white)
    savefig(hdl,"DepositionRateModelCI.pdf");
    display(hdl)

## --- Make heatmap

    rateplotmax = 3*maximum(dhdt); # May want to adjust this -- this is just a guess
    using StatsBase: fit, Histogram
    edges = linsp(0, rateplotmax, length(ages)-spacing+1)
    dhdt_im = Array{Float64}(undef,length(ages)-spacing,length(ages)-spacing);
    for i=1:length(ages)-spacing
        dhdt_im[:,i] = fit(Histogram, dhdt_dist[i, .~ isnan.(dhdt_dist[i,:])], edges, closed=:left).weights
    end

    # Rescale image to fit in UInt8 (0-255)
    using IndirectArrays: IndirectArray
    imSc = dhdt_im./pctile(dhdt_im[:],97.5)*256; # Rescale to include 97.5 CI (may need ot adjust)
    imSc[imSc.>255] = 255; # Cut off to fit in Uint8

    # Apply colormap. Available colormaps include viridis, inferno, plasma, fire
    A = IndirectArray(floor.(UInt8,imSc)+1, inferno);

    # Plot image
    img = plot(bincenters,cntr(edges),A,yflip=false,xflip=false, colorbar=:right);
    plot!(img, xlabel="Age ($AgeUnit)", ylabel="Rate ($HeightUnit / $AgeUnit, $binwidth $AgeUnit Bin)")
    savefig(img,"DepositionRateModelHeatmap.pdf");
    display(img)

    # dhdt_im_log = copy(dhdt_im);
    # dhdt_im_log[dhdt_im .>0] = log10.(dhdt_im[dhdt_im .>0])
    # heatmap(bincenters,cntr(edges),dhdt_im_log, xlabel="Age ($AgeUnit)", ylabel="Rate ($HeightUnit / $AgeUnit, $binwidth $AgeUnit Bin)")


## --- Probability that a given interval of stratigraphy was deposited entirely before/after a given time

    # Stratigraphic height and absoltue age/uncert to test
    testHeight = -40.0
    testAge = 66.0
    testAge_sigma = 0.05

    # Find index of nearest model height node
    nearest = indmin((testHeight-mdl.Height).^2);

    # Cycle through each possible age within testAge +/- 5 sigma, with resolution of 1/50 sigma
    test_ages = (testAge-5*testAge_sigma):testAge_sigma/50:(testAge+5*testAge_sigma)
    test_prob_older = Array{Float64}(undef,size(test_ages))
    # Evaluate the probability that model age is older than each test_age at the given strat level
    for i=1:length(test_ages)
        test_prob_older[i] = sum(agedist[nearest,:] .> test_ages[i]) ./ size(agedist,2);
    end

    # Normalized probability for each distance away from testAge between +5sigma and -5sigma
    prob_norm = normpdf(testAge, testAge_sigma, test_ages) ./ sum(normpdf(testAge, testAge_sigma, test_ages));  # SUM = 1

    # Integrate the product
    prob_older = sum(test_prob_older .* prob_norm)
    print("$(prob_older*100) % chance that $(mdl.Height[nearest]) $HeightUnit was deposited before $testAge +/- $testAge_sigma $AgeUnit Gaussian")


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

## --- (Optional) Add systematic uncertainties for U-Pb data

    # # Tracer (ET2535) uncertainty converted from per cent to relative
    # unc_tracer = 0.03/2/100;
    #
    # # U-238 Decay constant and uncertainty, Myr^-1
    # lambda238 = 1.55125e-10 * 1e6;
    # unc_lambda238 = 0.107/2/100; # converted from per cent to relative
    #
    # # Consider only the distribution of ages at model nodes where we have an ash bed
    # age_dist_X = Array{Float64}(undef,length(smpl.Height),size(agedist,2));
    # for i = 1:length(smpl.Height)
    #    closest_model_node = indmin(abs.(mdl.Height-smpl.Height[i]))
    #    age_dist_X[i,:] = agedist[closest_model_node,:];
    # end
    #
    # # Convert ages to 206Pb/238U ratios of the distribution
    # ratio_dist = exp.(age_dist_X.*lambda238)-1;
    #
    # # Add tracer uncertainty
    # ratio_dist_tracerunc = Array{Float64}(undef,size(ratio_dist));
    # for i=1:size(ratio_dist,2)
    #     ratio_dist_tracerunc[:,i] = ratio_dist[:,i].*(1 + unc_tracer*randn());
    # end
    #
    # # Convert 206/238 ratios back to ages, in Ma
    # age_dist_XY = log.(ratio_dist_tracerunc+1)./lambda238;
    #
    # # Add decay constant uncertainty
    # age_dist_XYZ = Array{Float64}(undef,size(ratio_dist));
    # for i=1:size(ratio_dist,2)
    #     age_dist_XYZ[:,i] = log.(ratio_dist_tracerunc[:,i]+1)./(lambda238.*(1 + unc_lambda238.*randn()));
    # end
    #
    # # Calculate the means and 95% confidence intervals for different levels of systematic uncertainties
    #
    # age_dist_X_mean = mean(age_dist_X,2); # Mean age
    # age_dist_X_std =  std(age_dist_X,2); # Standard deviation
    # age_dist_X_median = median(age_dist_X,2); # Median age
    # age_dist_X_025p = pctile(age_dist_X,2.5,dim=2); # 2.5th percentile
    # age_dist_X_975p = pctile(age_dist_X,97.5,dim=2); # 97.5th percentile
    #
    # age_dist_XY_mean = mean(age_dist_XY,2); # Mean age
    # age_dist_XY_std =  std(age_dist_XY,2); # Standard deviation
    # age_dist_XY_median = median(age_dist_XY,2); # Median age
    # age_dist_XY_025p = pctile(age_dist_XY,2.5,dim=2); # 2.5th percentile
    # age_dist_XY_975p = pctile(age_dist_XY,97.5,dim=2); # 97.5th percentile
    #
    # age_dist_XYZ_mean = mean(age_dist_XYZ,2); # Mean age
    # age_dist_XYZ_std =  std(age_dist_XYZ,2); # Standard deviation
    # age_dist_XYZ_median = median(age_dist_XYZ,2); # Median age
    # age_dist_XYZ_025p = pctile(age_dist_XYZ,2.5,dim=2); # 2.5th percentile
    # age_dist_XYZ_975p = pctile(age_dist_XYZ,97.5,dim=2); # 97.5th percentile
    #
    # age_X_95p = [age_dist_X_mean age_dist_X_975p-age_dist_X_mean age_dist_X_mean-age_dist_X_025p];
    # age_XY_95p = [age_dist_XY_mean age_dist_XY_975p-age_dist_XY_mean age_dist_XY_mean-age_dist_XY_025p];
    # age_XYZ_95p = [age_dist_XYZ_mean age_dist_XYZ_975p-age_dist_XYZ_mean age_dist_XYZ_mean-age_dist_XYZ_025p];

## ---
