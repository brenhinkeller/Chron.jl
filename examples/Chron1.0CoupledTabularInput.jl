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
    cd(@__DIR__) # Move to script directory

## --- Define sample properties - - - - - - - - - - - - - - - - - - - - - - - -

    # As in Chron1.0Coupled, but here we read in the information from a csv table:
    ds = importdataset("coupled_example_input.csv", importas=:Tuple)

    smpl = ChronAgeData(length(ds.Name))
    smpl.Name      =  Tuple(ds.Name)
    smpl.Height   .=  ds.Height
    smpl.Height_sigma .= ds.Height_sigma
    smpl.Age_Sidedness .= ds.Age_Sidedness              # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
    smpl.Path = first(ds.Path)                          # Where are the data files?
    smpl.inputSigmaLevel = first(ds.inputSigmaLevel)    # i.e., are the data files 1-sigma or 2-sigma. Integer.
    smpl.Age_Unit = first(ds.Age_Unit)                  # Unit of measurement for ages and errors in the data files
    smpl.Height_Unit = first(ds.Height_Unit)            # Unit of measurement for Height and Height_sigma

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
    exportdataset(NamedTuple(mdl), "AgeDepthModel.csv")

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


## --- Interpolate results at a specific heights read from file  - - - - - - - -

    # Stratigraphic heights at which to interpolate
    hds = importdataset("coupled_example_height_query.csv", importas=:Tuple)

    hds = (;
        Height = hds.Height,
        Age = linterp1s(mdl.Height, mdl.Age, hds.Height),
        Age_sigma = linterp1s(mdl.Height, mdl.Age_sigma, hds.Height),
        Age_025CI = linterp1s(mdl.Height, mdl.Age_025CI, hds.Height),
        Age_975CI = linterp1s(mdl.Height, mdl.Age_975CI, hds.Height),
    )
    exportdataset(hds, "coupled_example_height_results.csv")

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
    hdl = plot(
        xlabel="Age ($(smpl.Age_Unit))", 
        ylabel="Depositional Rate ($(smpl.Height_Unit) / $(smpl.Age_Unit) over $binwidth $(smpl.Age_Unit))", 
        fg_color_legend=:white,
        framestyle=:box,
    )
    plot!(hdl, agebincenters,dhdt, label="Mean", color=:black, linewidth=2)
    plot!(hdl,[agebincenters; reverse(agebincenters)],[dhdt_16p; reverse(dhdt_84p)], fill=(0,0.2,:darkblue), linealpha=0, label="68% CI")
    for lci in 20:5:45
        dhdt_lp = nanpctile(dhdt_dist,lci,dim=2)
        dhdt_up = nanpctile(dhdt_dist,100-lci,dim=2)
        plot!(hdl,[agebincenters; reverse(agebincenters)],[dhdt_lp; reverse(dhdt_up)], fill=(0,0.2,:darkblue), linealpha=0, label="")
    end
    plot!(hdl, agebincenters,dhdt_50p, label="Median", color=:grey, linewidth=1)
    savefig(hdl,"DepositionRateModelCI.pdf")
    display(hdl)

## --- Make heatmap - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    rateplotmax = 3*maximum(dhdt) # May want to adjust this -- this is just a guess
    using StatsBase: fit, Histogram
    ratebinedges = range(0, rateplotmax, length=length(agebincenters)+1)
    dhdt_im = zeros(length(agebincenters),length(agebincenters))
    for i=1:length(agebincenters)
        histcounts!(view(dhdt_im, :, i), view(dhdt_dist,i,:), ratebinedges)
    end

    # Apply colormap. Available colormaps include viridis, inferno, plasma, fire, etc.
    img = imsc(dhdt_im, inferno, nanpctile(dhdt_im, 2.5), nanpctile(dhdt_im, 97.5))

    # Plot image
    h = plot(agebincenters, cntr(ratebinedges), img, yflip=false, xflip=false, colorbar=:right, framestyle=:box)
    plot!(h, xlabel="Age ($(smpl.Age_Unit))", ylabel="Rate ($(smpl.Height_Unit) / $(smpl.Age_Unit), $binwidth $(smpl.Age_Unit) Bin)")
    xrange = abs(last(agebincenters)-first(agebincenters))
    yrange = abs(last(ratebinedges) - first(ratebinedges))
    plot!(h, ylims = extrema(cntr(ratebinedges)), size=(600,400), aspectratio=2/3/(yrange/xrange))
    savefig(h,"AccumulationRateModelHeatmap.pdf")
    display(h)

## --- (Optional) If your section has hiata / exposure surfaces of known duration, try this:

    # Data about hiatuses
    nhiatuses = 2 # The number of hiatuses you have data for
    hiatus = HiatusData(nhiatuses) # Struct to hold data
    hiatus.Height         = [-7.0, 35.0 ]
    hiatus.Height_sigma   = [ 0.0,  0.0 ]
    hiatus.Duration       = [ 0.3,  0.3 ]
    hiatus.Duration_sigma = [ 0.05, 0.05]

    # Run the model. Note the additional `hiatus` arguments
    @time (mdl, agedist, hiatusdist, lldist) = StratMetropolisDist(smpl, hiatus, config)
    exportdataset(NamedTuple(mdl), "AgeDepthModel.csv")

    # Plot results (mean and 95% confidence interval for both model and data)
    hdl = plot([mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(minimum(mdl.Height),0.5,:blue), label="model")
    plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="")
    plot!(hdl, smpl.Age, smpl.Height, xerror=(smpl.Age-smpl.Age_025CI,smpl.Age_975CI-smpl.Age),label="data",seriestype=:scatter,color=:black)
    plot!(hdl, xlabel="Age (Ma)", ylabel="Height (cm)", framestyle=:box)
    savefig(hdl,"AgeDepthModel.pdf")
    display(hdl)

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
