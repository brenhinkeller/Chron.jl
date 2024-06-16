# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                             Chron1.0Coupled.jl                                #
#                                                                               #
#     Illustrates the use of the Chron.jl package for eruption/deposition       #
#  age estimation and production of a stratigraphic age model.                  #
#                                                                               #
#      You may have to adjust the path below which specifies the location of    #
#  the CSV data files for each sample, depending on what you want to run.       #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## --- Load the Chron package - - - - - - - - - - - - - - - - - - - - - - - - -

    using Chron
    using Plots

## --- Define sample properties - - - - - - - - - - - - - - - - - - - - - - - -

    # This example data is from Clyde et al. (2016) "Direct high-precision
    # U–Pb geochronology of the end-Cretaceous extinction and calibration of
    # Paleocene astronomical timescales" EPSL 452, 272–280.
    # doi: 10.1016/j.epsl.2016.07.041

    nSamples = 5 # The number of samples you have data for
    smpl = ChronAgeData(nSamples)
    smpl.Name      =  ("KJ08-157", "KJ04-75", "KJ09-66", "KJ04-72", "KJ04-70",)
    smpl.Height   .=  [     -52.0,      44.0,      54.0,      82.0,      93.0,]
    smpl.Height_sigma .= [    3.0,       1.0,       3.0,       3.0,       3.0,]
    smpl.Age_Sidedness .= zeros(nSamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
    smpl.Path = joinpath(@__DIR__, "DenverUPbExampleData") # Where are the data files?
    smpl.inputSigmaLevel = 2 # i.e., are the data files 1-sigma or 2-sigma. Integer.
    smpl.Age_Unit = "Ma" # Unit of measurement for ages and errors in the data files
    smpl.Height_Unit = "cm" # Unit of measurement for Height and Height_sigma

    # IMPORTANT: smpl.Height must increase with increasing stratigraphic height
    # -- i.e., stratigraphically younger samples must be more positive. For this
    # reason, it is convenient to represent depths below surface as negative
    # numbers.

    # For each sample in smpl.Name, we must have a csv file at smpl.Path which
    # contains two columns of data, namely:
    # Age, Age sigma
    # where uncertainty (sigma) is absolute uncertainty.
    # For instance, examples/DenverUPbExampleData/KJ08-157.csv contains:
    #
    #   66.12,0.14
    #   66.115,0.048
    #   66.11,0.1
    #   66.11,0.17
    #   66.096,0.056
    #   66.088,0.081
    #   66.085,0.076
    #   66.073,0.084
    #   66.07,0.11
    #   66.055,0.043
    #   66.05,0.16
    #   65.97,0.12

    # If you are using U-Pb data and want Pb-loss-aware eruption estimation,
    # simply provide five columns of data instea of two, corresponding to
    # ²⁰⁷Pb/²³⁵U, ²⁰⁷Pb/²³⁵U sigma, ²⁰⁶Pb/²³⁸U, ²⁰⁶Pb/²³⁸U sigma, correlation


## --- Bootstrap pre-eruptive distribution - - - - - - - - - - - - - - - - - - -

    # Bootstrap a KDE of the pre-eruptive (or pre-depositional) mineral age
    # distribution using a KDE of stacked sample data from each data file
    BootstrappedDistribution = BootstrapCrystDistributionKDE(smpl)
    x = range(0,1,length=length(BootstrappedDistribution))
    h = plot(x, BootstrappedDistribution,
        label="Bootstrapped distribution",
        xlabel="Time (arbitrary units)",
        ylabel="Probability Density",
        framestyle=:box
    )
    savefig(h, joinpath(smpl.Path,"BootstrappedDistribution.pdf"))
    display(h)

## --- Estimate the eruption age distributions for each sample  - - - - - - - -

    # Configure distribution model here
    distSteps = 1*10^6 # Number of steps to run in distribution MCMC
    distBurnin = distSteps÷10 # Number to discard

    # Choose the form of the prior closure/crystallization distribution to use
    dist = BootstrappedDistribution
    ## You might alternatively consider:
    # dist = UniformDistribution              # A reasonable default
    # dist = MeltsVolcanicZirconDistribution  # A single magmatic pulse, truncated by eruption
    # dist = ExponentialDistribution          # Applicable for survivorship processes, potentially including inheritance/dispersion in Ar-Ar dates


    # Run MCMC to estimate saturation and eruption/deposition age distributions
    @time tMinDistMetropolis(smpl,distSteps,distBurnin,dist)

    # This will save rank-order and distribution plots, and print results to a
    # csv file -- you can find them in smpl.Path

## --- Run stratigraphic model - - - - - - - - - - - - - - - - - - - - - - - - -

    # Configure the stratigraphic Monte Carlo model
    config = StratAgeModelConfiguration()
    # If you in doubt, you can probably leave these parameters as-is
    config.resolution = 1.0 # Same units as sample height. Smaller is slower!
    config.bounding = 0.5 # how far away do we place runaway bounds, as a fraction of total section height
    (bottom, top) = extrema(smpl.Height)
    npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
    config.nsteps = 15000 # Number of steps to run in distribution MCMC
    config.burnin = 10000*npoints_approx # Number to discard
    config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Run the stratigraphic MCMC model
    @time (mdl, agedist, lldist) = StratMetropolisDist(smpl, config)


## --- Plot stratigraphic model - - - - - - - - - - - - - - - - - - - - - - - -

    # Plot results (mean and 95% confidence interval for both model and data)
    hdl = plot(framestyle=:box,
        xlabel="Age ($(smpl.Age_Unit))",
        ylabel="Height ($(smpl.Height_Unit))",
    )
    plot!(hdl, [mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(round(Int,minimum(mdl.Height)),0.5,:blue), label="model") # Age-depth model CI
    plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="") # Center line
    t = smpl.Age_Sidedness .== 0 # Two-sided constraints (plot in black)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],smpl.Age_975CI[t]-smpl.Age[t]),label="data",seriestype=:scatter,color=:black)
    t = smpl.Age_Sidedness .== 1 # Minimum ages (plot in cyan)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],zeros(count(t))),label="",seriestype=:scatter,color=:cyan,msc=:cyan)
    any(t) && zip(smpl.Age[t], smpl.Age[t].+nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
    t = smpl.Age_Sidedness .== -1 # Maximum ages (plot in orange)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(zeros(count(t)),smpl.Age_975CI[t]-smpl.Age[t]),label="",seriestype=:scatter,color=:orange,msc=:orange)
    any(t) && zip(smpl.Age[t], smpl.Age[t].-nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
    savefig(hdl,"AgeDepthModel.pdf")
    display(hdl)


## --- Interpolate model age at a specific stratigraphic height - - - - - - - -

    # Stratigraphic height at which to interpolate
    interp_height = 0

    age_at_height = linterp1s(mdl.Height,mdl.Age,interp_height)
    age_at_height_min = linterp1s(mdl.Height,mdl.Age_025CI,interp_height)
    age_at_height_max = linterp1s(mdl.Height,mdl.Age_975CI,interp_height)
    print("Interpolated age at height=$interp_height: $age_at_height +$(age_at_height_max-age_at_height)/-$(age_at_height-age_at_height_min) $(smpl.Age_Unit)")

    # Optional: interpolate full age distribution
    interpolated_distribution = Array{Float64}(undef,size(agedist,2))
    for i=1:size(agedist,2)
        interpolated_distribution[i] = linterp1s(mdl.Height,agedist[:,i],interp_height)
    end
    hdl = histogram(interpolated_distribution, nbins=50, label="", framestyle=:box)
    plot!(hdl, xlabel="Age ($(smpl.Age_Unit)) at height=$interp_height", ylabel="Likelihood (unnormalized)")
    savefig(hdl, "Interpolated age distribution.pdf")
    display(hdl)

## --- Calculate deposition rate binned by age  - - - - - - - - - - - - - - - -

    # Set bin width and spacing
    binwidth = round(nanrange(mdl.Age)/10,sigdigits=1) # Can also set manually, commented out below
    # binwidth = 0.01 # Same units as smpl.Age
    binoverlap = 10

    agebinedges = collect(minimum(mdl.Age):binwidth/binoverlap:maximum(mdl.Age))
    agebincenters = (agebinedges[1:end-binoverlap] + agebinedges[1+binoverlap:end])/2

    # Calculate rates for the stratigraphy of each markov chain step
    dhdt_dist = zeros(length(agebincenters), config.nsteps)
    @time for i=1:config.nsteps
        heights = linterp1(reverse(agedist[:,i]), reverse(mdl.Height), agebinedges, extrapolate=NaN)
        dhdt_dist[:,i] .= (heights[1:end-binoverlap] - heights[binoverlap+1:end]) ./ binwidth
    end

    # Find mean and 1-sigma (68%) CI
    dhdt = nanmean(dhdt_dist,dim=2)
    dhdt_50p = nanmedian(dhdt_dist,dim=2)
    dhdt_16p = nanpctile(dhdt_dist,15.865,dim=2) # Minus 1-sigma (15.865th percentile)
    dhdt_84p = nanpctile(dhdt_dist,84.135,dim=2) # Plus 1-sigma (84.135th percentile)

    # Plot results
    hdl = plot(agebincenters,dhdt, label="Mean", color=:black, linewidth=2)
    plot!(hdl,[agebincenters; reverse(agebincenters)],[dhdt_16p; reverse(dhdt_84p)], fill=(0,0.2,:darkblue), linealpha=0, label="68% CI")
    for lci in 20:5:45
        dhdt_lp = nanpctile(dhdt_dist,lci,dim=2)
        dhdt_up = nanpctile(dhdt_dist,100-lci,dim=2)
        plot!(hdl,[agebincenters; reverse(agebincenters)],[dhdt_lp; reverse(dhdt_up)], fill=(0,0.2,:darkblue), linealpha=0, label="")
    end
    plot!(hdl, agebincenters,dhdt_50p, label="Median", color=:grey, linewidth=1)
    plot!(hdl, xlabel="Age ($(smpl.Age_Unit))", ylabel="Depositional Rate ($(smpl.Height_Unit) / $(smpl.Age_Unit) over $binwidth $(smpl.Age_Unit))", fg_color_legend=:white)
    # savefig(hdl,"DepositionRateModelCI.pdf")
    display(hdl)

## --- Make heatmap - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    rateplotmax = 3*maximum(dhdt) # May want to adjust this -- this is just a guess
    using StatsBase: fit, Histogram
    edges = range(0, rateplotmax, length=length(ages)-spacing+1)
    dhdt_im = Array{Float64}(undef,length(ages)-spacing,length(ages)-spacing)
    for i=1:length(ages)-spacing
        dhdt_im[:,i] .= fit(Histogram, dhdt_dist[i, .~ isnan.(view(dhdt_dist,i,:))], edges, closed=:left).weights
    end

    # Apply colormap. Available colormaps include viridis, inferno, plasma, fire, etc.
    img = imsc(dhdt_im, inferno, nanpctile(dhdt_im, 2.5), nanpctile(dhdt_im, 97.5))

    # Plot image
    h = plot(bincenters, cntr(edges), img, yflip=false, xflip=false, colorbar=:right, framestyle=:box)
    plot!(h, xlabel="Age ($(smpl.Age_Unit))", ylabel="Rate ($(smpl.Height_Unit) / $(smpl.Age_Unit), $binwidth $(smpl.Age_Unit) Bin)")
    xrange = abs(last(bincenters)-first(bincenters))
    yrange = abs(last(edges) - first(edges))
    plot!(h, ylims = extrema(cntr(edges)), size=(600,400), aspectratio=2/3/(yrange/xrange))
    savefig(h,"AccumulationRateModelHeatmap.pdf")
    display(h)

## --- Probability that a given interval of stratigraphy was deposited entirely before/after a given time

    # Stratigraphic height and absoltue age/uncert to test
    testHeight = -40.0
    testAge = 66.0
    testAge_sigma = 0.05

    # Find index of nearest model height node
    nearest = argmin((testHeight .- mdl.Height).^2)

    # Cycle through each possible age within testAge +/- 5 sigma, with resolution of 1/50 sigma
    test_ages = (testAge-5*testAge_sigma):testAge_sigma/50:(testAge+5*testAge_sigma)
    test_prob_older = Array{Float64}(undef,size(test_ages))
    # Evaluate the probability that model age is older than each test_age at the given strat level
    for i=1:length(test_ages)
        test_prob_older[i] = sum(agedist[nearest,:] .> test_ages[i]) ./ size(agedist,2)
    end

    # Normalized probability for each distance away from testAge between +5sigma and -5sigma
    prob_norm = normpdf.(testAge, testAge_sigma, test_ages) ./ sum(normpdf.(testAge, testAge_sigma, test_ages));  # SUM = 1

    # Integrate the product
    prob_older = sum(test_prob_older .* prob_norm)
    print("$(prob_older*100) % chance that $(mdl.Height[nearest]) $(smpl.Height_Unit) was deposited before $testAge +/- $testAge_sigma $(smpl.Age_Unit) Gaussian")


## --- (Optional) If your section has hiata / exposure surfaces of known duration, try this:

    # Data about hiatuses
    nHiatuses = 2 # The number of hiatuses you have data for
    hiatus = HiatusData(nHiatuses) # Struct to hold data
    hiatus.Height         = [-7.0, 35.0 ]
    hiatus.Height_sigma   = [ 0.0,  0.0 ]
    hiatus.Duration       = [ 0.3,  0.3 ]
    hiatus.Duration_sigma = [ 0.05, 0.05]

    # Run the model. Note the additional `hiatus` arguments
    @time (mdl, agedist, hiatusdist, lldist) = StratMetropolisDist(smpl, hiatus, config)

    # Plot results (mean and 95% confidence interval for both model and data)
    hdl = plot([mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(minimum(mdl.Height),0.5,:blue), label="model")
    plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="")
    plot!(hdl, smpl.Age, smpl.Height, xerror=(smpl.Age-smpl.Age_025CI,smpl.Age_975CI-smpl.Age),label="data",seriestype=:scatter,color=:black)
    plot!(hdl, xlabel="Age (Ma)", ylabel="Height (cm)", framestyle=:box)

## --- (Optional) Add systematic uncertainties for U-Pb data

    # # Tracer (ET2535) uncertainty converted from per cent to relative
    # unc_tracer = 0.03/2/100
    #
    # # U-238 Decay constant and uncertainty, Myr^-1
    # lambda238 = 1.55125e-10 * 1e6
    # unc_lambda238 = 0.107/2/100 # converted from per cent to relative
    #
    # # Consider only the distribution of ages at model nodes where we have an ash bed
    # age_dist_X = Array{Float64}(undef,length(smpl.Height),size(agedist,2))
    # for i = 1:length(smpl.Height)
    #    closest_model_node = argmin(abs.(mdl.Height-smpl.Height[i]))
    #    age_dist_X[i,:] = agedist[closest_model_node,:]
    # end
    #
    # # Convert ages to 206Pb/238U ratios of the distribution
    # ratio_dist = exp.(age_dist_X.*lambda238)-1
    #
    # # Add tracer uncertainty
    # ratio_dist_tracerunc = Array{Float64}(undef,size(ratio_dist))
    # for i=1:size(ratio_dist,2)
    #     ratio_dist_tracerunc[:,i] = ratio_dist[:,i].*(1 + unc_tracer*randn())
    # end
    #
    # # Convert 206/238 ratios back to ages, in Ma
    # age_dist_XY = log.(ratio_dist_tracerunc+1)./lambda238
    #
    # # Add decay constant uncertainty
    # age_dist_XYZ = Array{Float64}(undef,size(ratio_dist))
    # for i=1:size(ratio_dist,2)
    #     age_dist_XYZ[:,i] = log.(ratio_dist_tracerunc[:,i]+1)./(lambda238.*(1 + unc_lambda238.*randn()))
    # end
    #
    # # Calculate the means and 95% confidence intervals for different levels of systematic uncertainties
    #
    # age_dist_X_mean = nanmean(age_dist_X,2) # Mean age
    # age_dist_X_std =  nanstd(age_dist_X,2) # Standard deviation
    # age_dist_X_median = nanmedian(age_dist_X,2) # Median age
    # age_dist_X_025p = nanpctile(age_dist_X,2.5,dim=2) # 2.5th percentile
    # age_dist_X_975p = nanpctile(age_dist_X,97.5,dim=2) # 97.5th percentile
    #
    # age_dist_XY_mean = nanmean(age_dist_XY,2) # Mean age
    # age_dist_XY_std =  nanstd(age_dist_XY,2) # Standard deviation
    # age_dist_XY_median = nanmedian(age_dist_XY,2) # Median age
    # age_dist_XY_025p = nanpctile(age_dist_XY,2.5,dim=2) # 2.5th percentile
    # age_dist_XY_975p = nanpctile(age_dist_XY,97.5,dim=2) # 97.5th percentile
    #
    # age_dist_XYZ_mean = nanmean(age_dist_XYZ,2) # Mean age
    # age_dist_XYZ_std =  nanstd(age_dist_XYZ,2) # Standard deviation
    # age_dist_XYZ_median = nanmedian(age_dist_XYZ,2) # Median age
    # age_dist_XYZ_025p = nanpctile(age_dist_XYZ,2.5,dim=2) # 2.5th percentile
    # age_dist_XYZ_975p = nanpctile(age_dist_XYZ,97.5,dim=2) # 97.5th percentile
    #
    # age_X_95p = [age_dist_X_mean age_dist_X_975p-age_dist_X_mean age_dist_X_mean-age_dist_X_025p]
    # age_XY_95p = [age_dist_XY_mean age_dist_XY_975p-age_dist_XY_mean age_dist_XY_mean-age_dist_XY_025p]
    # age_XYZ_95p = [age_dist_XYZ_mean age_dist_XYZ_975p-age_dist_XYZ_mean age_dist_XYZ_mean-age_dist_XYZ_025p]

## ---
