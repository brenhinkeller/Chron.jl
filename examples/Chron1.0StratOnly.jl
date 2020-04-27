## --- Load required pacages, install Chron if required

    try
        using Chron
    catch
        using Pkg
        Pkg.add(PackageSpec(url="https://github.com/brenhinkeller/Chron.jl"))
        using Chron
    end

    using Statistics, StatsBase, SpecialFunctions
    using Plots; gr();

## --- Define sample properties

    # # # # # # # # # # # Enter sample information here! # # # # # # # # # # # #
    # Input the number of samples we wish to model (must match below)
    nSamples = 4
    # Make an instance of a ChronSection object for nSamples
    smpl = NewChronAgeData(nSamples)
    smpl.Name          = ("Sample 1", "Sample 2", "Sample 3", "Sample 4") # Et cetera
    smpl.Age           = [ 6991,  7088,  7230,  7540,] # Measured ages
    smpl.Age_sigma     = [   30,    70,    50,    50,] # Measured 1-σ uncertainties
    smpl.Height[:]     = [ -355,  -380,-397.0,-411.5,] # Depths below surface should be negative
    smpl.Height_sigma[:]  = fill(0.01, nSamples) # Usually assume little or no sample height uncertainty
    smpl.Age_Sidedness[:] = zeros(nSamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
    smpl.Age_Unit = "Years BP" # Unit of measurement for ages
    smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma

    # IMPORTANT: smpl.Height must increase with increasing stratigraphic height
    # -- i.e., stratigraphically younger samples must be more positive. For this
    # reason, it is convenient to represent depths below surface as negative
    # numbers.
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


## --- Run stratigraphic model

    # # # # # # # # # # Configure stratigraphic model here! # # # # # # # # # #
    # Configure the stratigraphic Monte Carlo model
    config = NewStratAgeModelConfiguration()
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
    (mdl, agedist, lldist) = StratMetropolis(smpl, config)

    # Plot results (mean and 95% confidence interval for both model and data)
    hdl = plot([mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(round(Int,minimum(mdl.Height)),0.5,:blue), label="model")
    plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="")
    plot!(hdl, smpl.Age, smpl.Height, xerror=smpl.Age_sigma*2, label="data",seriestype=:scatter,color=:black)
    plot!(hdl, xlabel="Age ($(smpl.Age_Unit))", ylabel="Height ($(smpl.Height_Unit))")
    savefig(hdl,"AgeDepthModel.pdf");
    display(hdl)

## --- Interpolate results at a specific height

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

## --- Calculate deposition rate binned by age  - - - - - - - - - - - - - - - -

    # Set bin width and spacing
    binwidth = round(nanrange(mdl.Age)/10,sigdigits=1) # Can also set manually, commented out below
    # binwidth = 100 # Same units as smpl.Age
    binoverlap = 10
    ages = collect(minimum(mdl.Age):binwidth/binoverlap:maximum(mdl.Age))
    bincenters = ages[1+Int(binoverlap/2):end-Int(binoverlap/2)]
    spacing = binoverlap

    # Calculate rates for the stratigraphy of each markov chain step
    dhdt_dist = Array{Float64}(undef, length(ages)-binoverlap, config.nsteps)
    @time for i=1:config.nsteps
        heights = linterp1(reverse(agedist[:,i]), reverse(mdl.Height), ages)
        dhdt_dist[:,i] .= abs.(heights[1:end-spacing] - heights[spacing+1:end]) ./ binwidth
    end

    # Find mean and 1-sigma (68%) CI
    dhdt = nanmean(dhdt_dist,dim=2)
    dhdt_50p = nanmedian(dhdt_dist,dim=2)
    dhdt_16p = pctile(dhdt_dist,15.865,dim=2) # Minus 1-sigma (15.865th percentile)
    dhdt_84p = pctile(dhdt_dist,84.135,dim=2) # Plus 1-sigma (84.135th percentile)
    # Other confidence intervals (10:10:50)
    dhdt_20p = pctile(dhdt_dist,20,dim=2)
    dhdt_80p = pctile(dhdt_dist,80,dim=2)
    dhdt_25p = pctile(dhdt_dist,25,dim=2)
    dhdt_75p = pctile(dhdt_dist,75,dim=2)
    dhdt_30p = pctile(dhdt_dist,30,dim=2)
    dhdt_70p = pctile(dhdt_dist,70,dim=2)
    dhdt_35p = pctile(dhdt_dist,35,dim=2)
    dhdt_65p = pctile(dhdt_dist,65,dim=2)
    dhdt_40p = pctile(dhdt_dist,40,dim=2)
    dhdt_60p = pctile(dhdt_dist,60,dim=2)
    dhdt_45p = pctile(dhdt_dist,45,dim=2)
    dhdt_55p = pctile(dhdt_dist,55,dim=2)

    # Plot results
    hdl = plot(bincenters,dhdt, label="Mean", color=:black, linewidth=2)
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_16p; reverse(dhdt_84p)], fill=(0,0.2,:darkblue), linealpha=0, label="68% CI")
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_20p; reverse(dhdt_80p)], fill=(0,0.2,:darkblue), linealpha=0, label="")
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_25p; reverse(dhdt_75p)], fill=(0,0.2,:darkblue), linealpha=0, label="")
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_30p; reverse(dhdt_70p)], fill=(0,0.2,:darkblue), linealpha=0, label="")
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_35p; reverse(dhdt_65p)], fill=(0,0.2,:darkblue), linealpha=0, label="")
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_40p; reverse(dhdt_60p)], fill=(0,0.2,:darkblue), linealpha=0, label="")
    plot!(hdl,[bincenters; reverse(bincenters)],[dhdt_45p; reverse(dhdt_55p)], fill=(0,0.2,:darkblue), linealpha=0, label="")
    plot!(hdl,bincenters,dhdt_50p, label="Median", color=:grey, linewidth=1)
    plot!(hdl, xlabel="Age ($(smpl.Age_Unit))", ylabel="Depositional Rate ($(smpl.Height_Unit) / $(smpl.Age_Unit) over $binwidth $(smpl.Age_Unit))", fg_color_legend=:white)
    # savefig(hdl,"DepositionRateModelCI.pdf")
    display(hdl)

## --- Optional: Stratigraphic model including hiatuses

    # Data about hiatuses
    nHiatuses = 2 # The number of hiatuses you have data for
    hiatus = NewHiatusData(nHiatuses) # Struct to hold data
    hiatus.Height         = [-371.5, -405.0 ]
    hiatus.Height_sigma   = [   0.0,    0.0 ]
    hiatus.Duration       = [ 100.0,   123.0]
    hiatus.Duration_sigma = [  30.5,    20.0]

    # Run the model. Note the additional `hiatus` arguments
    (mdl, agedist, hiatusdist, lldist) = StratMetropolis(smpl, hiatus, config); sleep(0.5)

    # Plot results (mean and 95% confidence interval for both model and data)
    hdl = plot([mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(minimum(mdl.Height),0.5,:blue), label="model")
    plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="", fg_color_legend=:white)
    plot!(hdl, smpl.Age, smpl.Height, xerror=smpl.Age_sigma*2,label="data",seriestype=:scatter,color=:black)
    plot!(hdl, xlabel="Age ($(smpl.Age_Unit))", ylabel="Height ($(smpl.Height_Unit))")

## --- End of File
