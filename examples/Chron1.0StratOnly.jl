# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                            Chron1.0StratOnly.jl                               #
#                                                                               #
#     Illustrates the use of the Chron.jl package for the production of a       #
#  stratigraphic age-depth model based on any Gaussian age constraints          #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## --- Load required pacages  - - - - - - - - - - - - - - - - - - - - - - - - -

    using Chron
    using Plots

## --- Define sample properties - - - - - - - - - - - - - - - - - - - - - - - -

    # # # # # # # # # # # Enter sample information here! # # # # # # # # # # # #
    # Input the number of samples we wish to model (must match below)
    nsamples = 4
    # Make an instance of a ChronAgeData object for n samples
    smpl = ChronAgeData(nsamples)
    smpl.Name          = ("Sample 1", "Sample 2", "Sample 3", "Sample 4") # Et cetera
    smpl.Age          .= [ 6991,  7088,  7230,  7540,] # Measured ages
    smpl.Age_sigma    .= [   30,    70,    50,    50,] # Measured 1-σ uncertainties
    smpl.Height       .= [ -355,  -380,-397.0,-411.5,] # Depths below surface should be negative
    smpl.Height_sigma .= fill(0.01, nsamples) # Usually assume little or no sample height uncertainty
    smpl.Age_Sidedness .= zeros(nsamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
    smpl.Age_Unit = "Years BP" # Unit of measurement for ages
    smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma

    # IMPORTANT: smpl.Height must increase with increasing stratigraphic height
    # -- i.e., stratigraphically younger samples must be more positive. For this
    # reason, it is convenient to represent depths below surface as negative
    # numbers.
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


## --- Run stratigraphic model - - - - - - - - - - - - - - - - - - - - - - - - -

    # # # # # # # # # # Configure stratigraphic model here! # # # # # # # # # #
    # Configure the stratigraphic Monte Carlo model
    config = StratAgeModelConfiguration()
    # If you in doubt, you can probably leave these parameters as-is
    config.resolution = 0.2 # Same units as sample height. Smaller is slower!
    config.bounding = 0.5 # how far away do we place runaway bounds, as a fraction of total section height. Larger is slower.
    (bottom, top) = extrema(smpl.Height)
    npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
    config.nsteps = 15000 # Number of steps to run in distribution MCMC
    config.burnin = 10000*npoints_approx # Number to discard
    config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Run the stratigraphic MCMC model
    @time (mdl, agedist, lldist) = StratMetropolis(smpl, config)
    exportdataset(NamedTuple(mdl), "AgeDepthModel.csv")

    # Plot results (mean and 95% CI for model / 2-sigma for data)
    hdl = plot(framestyle=:box,
        fg_color_legend=:white,
        xlabel="Age ($(smpl.Age_Unit))",
        ylabel="Height ($(smpl.Height_Unit))",
    )
    plot!(hdl, [mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(round(Int,minimum(mdl.Height)),0.5,:blue), label="model")
    plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="") # Center line
    t = smpl.Age_Sidedness .== 0 # Two-sided constraints (plot in black)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=2*smpl.Age_sigma[t],label="data",seriestype=:scatter,color=:black)
    t = smpl.Age_Sidedness .== 1 # Minimum ages (plot in cyan)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(2*smpl.Age_sigma[t],zeros(count(t))),label="",seriestype=:scatter,color=:cyan,msc=:cyan)
    any(t) && zip(smpl.Age[t], smpl.Age[t].+nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
    t = smpl.Age_Sidedness .== -1 # Maximum ages (plot in orange)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(zeros(count(t)),2*smpl.Age_sigma[t]),label="",seriestype=:scatter,color=:orange,msc=:orange)
    any(t) && zip(smpl.Age[t], smpl.Age[t].-nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
    savefig(hdl,"AgeDepthModel.pdf")
    display(hdl)

## --- Interpolate results at a specific height - - - - - - - - - - - - - - - -

    # Stratigraphic height at which to interpolate
    height = -404

    age_at_height = linterp1s(mdl.Height,mdl.Age,height)
    age_at_height_min = linterp1s(mdl.Height,mdl.Age_025CI,height)
    age_at_height_max = linterp1s(mdl.Height,mdl.Age_975CI,height)
    print("Interpolated age at height=$height: $age_at_height +$(age_at_height_max-age_at_height)/-$(age_at_height-age_at_height_min) $(smpl.Age_Unit)")

    # Optional: interpolate full age distribution
    interpolated_distribution = Array{Float64}(undef,size(agedist,2))
    for i=1:size(agedist,2)
        interpolated_distribution[i] = linterp1s(mdl.Height,agedist[:,i],height)
    end
    hdl = histogram(interpolated_distribution, nbins=50, label="")
    plot!(hdl, xlabel="Age ($(smpl.Age_Unit)) at height=$height", ylabel="Likelihood (unnormalized)")
    savefig(hdl, "Interpolated age distribution.pdf")
    display(hdl)

## --- Calculate deposition rate binned by age - - - - - - - - - - - - - - - - -

    # Set bin width and spacing
    binwidth = round(nanrange(mdl.Age)/10,sigdigits=1) # Can also set manually, commented out below
    # binwidth = 100 # Same units as smpl.Age
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
## --- Optional: Stratigraphic model including hiatuses - - - - - - - - - - - -

    # Data about hiatuses
    nhiatuses = 2 # The number of hiatuses you have data for
    hiatus = HiatusData(nhiatuses) # Struct to hold data
    hiatus.Height         = [-371.5, -405.0 ]
    hiatus.Height_sigma   = [   0.0,    0.0 ]
    hiatus.Duration       = [ 100.0,   123.0]
    hiatus.Duration_sigma = [  30.5,    20.0]

    # Run the model. Note the additional `hiatus` arguments
    @time (mdl, agedist, hiatusdist, lldist) = StratMetropolis(smpl, hiatus, config); sleep(0.5)
    exportdataset(NamedTuple(mdl), "AgeDepthModel.csv")

    # Plot results (mean and 95% confidence interval for both model and data)
    hdl = plot([mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(minimum(mdl.Height),0.5,:blue), label="model")
    plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="", fg_color_legend=:white)
    plot!(hdl, smpl.Age, smpl.Height, xerror=smpl.Age_sigma*2,label="data",seriestype=:scatter,color=:black)
    plot!(hdl, xlabel="Age ($(smpl.Age_Unit))", ylabel="Height ($(smpl.Height_Unit))")
    savefig(hdl, "Interpolated age distribution.pdf")
    display(hdl)

## --- End of File - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
