# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                          Chron1.0CoupledConcordia.jl                          #
#                                                                               #
#     Illustrates the use of the Chron.jl package for Pb-loss-aware coupled     #
#  eruption/deposition age estimation and production of a stratigraphic age     #
#  model. Equivalent to the plain "coupled" example, except Pb-loss-aware.      #
#                                                                               #
#      You may have to adjust the path below which specifies the location of    #
#  the CSV data files for each sample, depending on what you want to run.       #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## --- Load the Chron package - - - - - - - - - - - - - - - - - - - - - - - - -

    using Chron
    using Plots

## --- Define sample properties - - - - - - - - - - - - - - - - - - - - - - - -

    # This example data is from MacLennan et al. 2020, 10.1126/sciadv.aay6647
    # This is all _exactly the same_ as for normal eruption/deposition ages,
    # except the CSV files in `smpl.Path` contain five columns of raw isotopic data
    nsamples = 3 # The number of samples you have data for
    smpl = ChronAgeData(nsamples)
    smpl.Name      =  ("KR18-04", "KR18-01", "KR18-05")
    smpl.Height   .=  [      0.0,     100.0,     200.0] # Arbitrary example heights
    smpl.Age_Sidedness .= zeros(nsamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
    smpl.Path = joinpath(@__DIR__, "ConcordiaExampleData") # Where are the data files?
    smpl.inputSigmaLevel = 1 # i.e., are the data files 1-sigma or 2-sigma. Integer.
    smpl.Age_Unit = "Ma" # Unit of measurement for ages and errors in the data files
    smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma

    # IMPORTANT: smpl.Height must increase with increasing stratigraphic height
    # -- i.e., stratigraphically younger samples must be more positive. For this
    # reason, it is convenient to represent depths below surface as negative
    # numbers.

    # For each sample in smpl.Name, we must have a csv file at smpl.Path which
    # contains five columns of data, with columns for
    # ²⁰⁷Pb/²³⁵U, ²⁰⁷Pb/²³⁵U sigma, ²⁰⁶Pb/²³⁸U, ²⁰⁶Pb/²³⁸U sigma, correlation
    # Where all uncertainties (sigma) must be absolute (not relative) uncertainties!
    # For instance, examples/ConcordiaExampleData/KR18-01.csv contains:
    #
    # 1.1002,0.00060511,0.123908,0.00001982528,0.333
    # 1.1003,0.0005226425,0.123893,0.000020442345,0.421
    # 1.1,0.0011,0.123874,0.00002353606,0.281
    # 1.1002,0.00060511,0.123845,0.000025388225,0.418
    # 1.1007,0.0005338395,0.123833,0.000025385765,0.534
    # 1.0991,0.001154055,0.123797,0.000031568235,0.298
    # 1.09931,0.0004067447,0.123762,0.00003898503,0.709
    # 1.09947,0.0004617774,0.123752,0.00002598792,0.579
    # 1.0986,0.00093381,0.123738,0.00003650271,0.288
    # 1.09883,0.00047799105,0.123735,0.0000222723,0.506
    # 1.09904,0.000384664,0.123733,0.000021653275,0.404
    # 1.0758,0.0005379,0.121175,0.00002302325,0.427

    # This is the only change you have to make relative to regular "Chron1.0Coupled"
    # Just provide five columns of data instead of two.

## --- (optional) Bootstrap pre-eruptive distribution - - - - - - - - - - - - -

    # Bootstrap a KDE of the pre-eruptive (or pre-depositional) mineral age
    # distribution using a KDE of stacked sample data from each data file
    BootstrappedDistribution = BootstrapCrystDistributionKDE(smpl, tpbloss=0)
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
    dist = HalfNormalDistribution
    ## You might alternatively consider:
    # dist = BootstrappedDistribution
    # dist = UniformDistribution              # A reasonable default
    # dist = MeltsVolcanicZirconDistribution  # A single magmatic pulse, truncated by eruption


    # Run MCMC to estimate saturation and eruption/deposition age distributions
    @time tMinDistMetropolis(smpl,distSteps,distBurnin,dist)

    # This will save Concordia and distribution plots (including distribution of
    # Pb-loss times), and print results to a csv file.
    # You can find all these results in `smpl.Path`


## --- Run stratigraphic model - - - - - - - - - - - - - - - - - - - - - - - - -


    # Configure the stratigraphic Monte Carlo model
    config = StratAgeModelConfiguration()
    # If you in doubt, you can probably leave these parameters as-is
    config.resolution = 10.0 # Same units as sample height. Smaller is slower!
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
        fg_color_legend=:white,
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
    interp_height = 150

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

## --- End of File
