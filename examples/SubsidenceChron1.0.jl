## For subsidence modeling in extensional basins

## --- Load required pacages, install Chron if required

    try
        using Chron
    catch
        using Pkg
        Pkg.add(PackageSpec(url="https://github.com/brenhinkeller/Chron.jl"))
        using Chron
    end

    import Pkg;
    Pkg.add(["StatGeochem","Distributions","Plots","StatsBase","ProgressMeter","LsqFit","KernelDensity","Interpolations","SpecialFunctions"])
    using StatGeochem, Distributions, Plots, Statistics, StatsBase, SpecialFunctions
    using StatsBase: fit, Histogram, percentile
    using ProgressMeter: @showprogress, Progress, update!
    using LsqFit: curve_fit
    using KernelDensity: kde
    using Interpolations 
    using Plots; gr();


## --- Define properties of the stratigraphy

    # # # # # # # # # # # Enter stratigraphic information here! # # # # # # # # # # # #
    data_csv = importdataset("Svalbard.csv",',')
    nLayers = length(data_csv["Thickness"])

    strat = NewStratData(nLayers)
    strat.Lithology          = data_csv["Lithology"] 
    strat.Thickness         .= data_csv["Thickness"] 

    # Run the decompaction and backstripping MC model
    @time (Sₜ, Sμ, Sσ, model_strat_heights) = DecompactBackstrip(strat)

    # Plot results - tectonic subsidence in comparison with present day stratigraphic heights
    p1 = plot(1:2232, Sμ, alpha = 1, yflip = true, xflip = true, label = "Tectonic subsidence", color = "blue")
    plot!(p1, 1:2232, reverse((model_strat_heights[2:end])*1000), yflip = true, label = "Present-day thickness", color = "red")
    plot!(p1, 1:2232, Sₜ[:,2:end], alpha = 0.01, label = "", yflip = true, color = "blue", fg_color_legend=:white)


## --- Define properties of age constraints

    # # # # # # # # # # # Enter age constraint (sample) information here! # # # # # # # # # # # #
    # Input the number of samples we wish to model (must match below)
    nSamples = 3
    # Make an instance of a ChronSection object for nSamples
    smpl = NewChronAgeData(nSamples)
    smpl.Name          = ("Sample 1", "Sample 2", "Sample 3") # Et cetera
    smpl.Age          .= [ 791.1,  737.5,  717] # Measured ages
    smpl.Age_sigma    .= [  2.45,    4.8,  0.4] # Measured 1-σ uncertainties
    smpl.Height       .= [ -1283,   -120,    0] # Depths below surface should be negative
    smpl.Height_sigma .= fill(0.01, nSamples) # Usually assume little or no sample height uncertainty
    smpl.Age_Sidedness .= zeros(nSamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
    smpl.Age_Unit = "Ma" # Unit of measurement for ages
    smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma

    # IMPORTANT: smpl.Height must increase with increasing stratigraphic height
    # -- i.e., stratigraphically younger samples must be more positive. For this
    # reason, it is convenient to represent depths below surface as negative
    # numbers.
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## --- Define thermal subsidence model parameters

    # # # # # # # # # # Enter thermal subsidence parameter priors here! # # # # # # # # # #
    # Enter initial guesses for the beta factor and thermal subsidence onset age and their uncertainties
    Beta = 1.3
    Beta_sigma = 0.2
    T0 = 820
    T0_sigma = 20

    therm = NewThermalSubsidenceParameters()
    therm.Param = [Beta, T0] 
    therm.Sigma = [Beta_sigma, T0_sigma] 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## --- Run stratigraphic model

    # # # # # # # # # # Configure stratigraphic model here! # # # # # # # # # #
    # Configure the stratigraphic Monte Carlo model
    config = NewStratAgeModelConfiguration()
    # If you in doubt, you can probably leave these parameters as-is
    config.resolution = 1 # Same units as sample height. Smaller is slower!
    config.bounding = 0.1 # how far away do we place runaway bounds, as a fraction of total section height. Larger is slower.
    #(bottom, top) = extrema(smpl.Height)
    npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
    config.nsteps = 1500 # Number of steps to run in distribution MCMC #!!TRIAL RUN SETTING!! ACTUAL RUNS SHOULD BE 15000
    config.burnin = 3000*npoints_approx # Number to discard #!!TRIAL RUN SETTING!! ACTUAL RUNS SHOULD BE 10000
    config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Run the stratigraphic MCMC model
    @time (subsmdl, agedist, lldist, beta_t0dist) = SubsidenceStratMetropolis(smpl, config, therm, model_strat_heights)

    # Plot results (mean and 95% confidence interval for both model and data)
    hdl = plot([mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(round(Int,minimum(mdl.Height)),0.5,:blue), label="model")
    plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="", fg_color_legend=:white) # Center line
    t = smpl.Age_Sidedness .== 0 # Two-sided constraints (plot in black)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],smpl.Age_975CI[t]-smpl.Age[t]),label="data",seriestype=:scatter,color=:black)
    t = smpl.Age_Sidedness .== 1 # Minimum ages (plot in cyan)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],zeros(count(t))),label="",seriestype=:scatter,color=:cyan,msc=:cyan)
    any(t) && zip(smpl.Age[t], smpl.Age[t].+nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
    t = smpl.Age_Sidedness .== -1 # Maximum ages (plot in orange)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(zeros(count(t)),smpl.Age_975CI[t]-smpl.Age[t]),label="",seriestype=:scatter,color=:orange,msc=:orange)
    any(t) && zip(smpl.Age[t], smpl.Age[t].-nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
    plot!(hdl, xlabel="Age ($(smpl.Age_Unit))", ylabel="Height ($(smpl.Height_Unit))")
    savefig(hdl,"AgeDepthModel.pdf")
    display(hdl)
    
#=
    # More plots - NOT READY TO RUN - still need to readjust the format
    curve_ages = similar(Sμ)
    curve_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta*sin(pi/beta)))).+t0
    model_strat_heights_m = copy(-model_strat_heights[2:end]).*1000
    plot!(hdl, curve_ages, model_strat_heights_m, linecolor =:blue, label = "beta and t0", xlabel="Age (Ma) calculated based on beta and t0 from MCMC", ylabel="Thermal subsidence (m)")
    savefig(hdl, figpath*"Age_depth_curve_upper_section_comparison.png")

    Sᵣ = reverse(Sμ)
    target_ages = [815.29, 811.51, 810.7, 788.7, 752.7, 739.9]
    target_age_errors = [0.32, 0.25, 6.3, 0.24, 5.5, 6.1]
    target_heights = [-(2.232-0.25), -(2.232-0.425), -(2.232-0.435), -(2.232-1.14), -(2.232-1.93), -(2.232-2.125)]
    target_heights_m = copy(target_heights)*1000
    target_height_errors = [0.212, 0.1, 0.05, 0.2, 0.035, 0.035]
    target_height_errors_m = copy(target_height_errors)*1000
    target_subs_heights = [-Sᵣ[250], -Sᵣ[425], -Sᵣ[435], -Sᵣ[1140], -Sᵣ[1930], -Sᵣ[2125]]
    target_subs_height_errors = [(Sᵣ[250]-Sᵣ[38]), (Sᵣ[425]-Sᵣ[325]), (Sᵣ[435]-Sᵣ[385]), (Sᵣ[1140]-Sᵣ[940]), (Sᵣ[1930]-Sᵣ[1895]), (Sᵣ[2125]-Sᵣ[2090])]

    scatter!(age_depth_curve, target_ages, target_heights_m, xerr = target_age_errors, yerr = target_height_errors_m, label = "correlated ages")
    scatter!(thermal_subs_curve, target_ages, target_subs_heights, xerr = target_age_errors, yerr = target_subs_height_errors, label = "correlated ages")
    
    figpath = "C:/Users/zhtia/Documents/GitHub/Chron.jl/trials/"
    savefig(age_depth_curve, figpath*"Age_depth_curve_entire_section.png")
    savefig(thermal_subs_curve, figpath*"Thermal_subs_curve_entire_section.png")    
=#

    ## --- Optional: Stratigraphic model including hiata

    # Data about hiata
    nHiatuses = 2 # The number of hiata you have data for
    hiatus = NewHiatusData(nHiatuses) # Struct to hold data
    hiatus.Height         = [-371.5, -405.0 ]
    hiatus.Height_sigma   = [   0.0,    0.0 ]
    hiatus.Duration       = [ 100.0,   123.0]
    hiatus.Duration_sigma = [  30.5,    20.0]

    # Run the model. Note the additional `hiatus` arguments
    @time (subsmdl, agedist, hiatusdist, lldist, beta_t0dist) = StratMetropolis(smpl, hiatus, config, therm, model_strat_heights); sleep(0.5)

    # Plot results (mean and 95% confidence interval for both model and data)
    hdl = plot([mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(minimum(mdl.Height),0.5,:blue), label="model")
    plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="", fg_color_legend=:white)
    plot!(hdl, smpl.Age, smpl.Height, xerror=smpl.Age_sigma*2,label="data",seriestype=:scatter,color=:black)
    plot!(hdl, xlabel="Age ($(smpl.Age_Unit))", ylabel="Height ($(smpl.Height_Unit))")

#=
    # More plots - NOT READY TO RUN - still need to readjust the format
    hdl_hiatus = plot([mdl_age_025CI; reverse(mdl_age_975CI)],[mdl_height; reverse(mdl_height)], fill=(round(Int,minimum(mdl_height)),0.5,:blue), label="Model results - 2σ error envelope")
    plot!(hdl_hiatus, mdl_age, mdl_height, linecolor=:blue, label="Model results - mean", fg_color_legend=:white) # Center line
    t = Age_Sidedness .== 0 # Two-sided constraints (plot in black)
    any(t) && plot!(hdl_hiatus, Age[t], Height[t], xerror=Age_sigma[t],label="Input - age data",seriestype=:scatter,color=:black)
    age_at_height = (mdl_height,mdl_age,-2232)

    curve_ages = similar(Sμ)
    curve_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta*sin(pi/beta)))).+t0
    model_strat_heights_m = copy(-model_strat_heights[2:end]).*1000
        
    hiatus_thermal_subs_curve = plot(curve_ages, -Sμ, linecolor =:blue, xflip = true, label = "beta and t0", xlabel="Age (Ma) calculated based on beta and t0 from MCMC", ylabel="Thermal subsidence (m)")
    curve_1_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_025CI*sin(pi/beta_025CI)))).+t0_025CI
    curve_2_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_025CI*sin(pi/beta_025CI)))).+t0_975CI
    curve_3_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_975CI*sin(pi/beta_975CI)))).+t0_025CI
    curve_4_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_975CI*sin(pi/beta_975CI)))).+t0_975CI
    #curve_1_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_lowest_t025*sin(pi/beta_lowest_t025)))).+t0_025CI
    #curve_2_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_lowest_t975*sin(pi/beta_lowest_t975)))).+t0_975CI
    #curve_3_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_highest_t025*sin(pi/beta_highest_t025)))).+t0_025CI
    #curve_4_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_highest_t975*sin(pi/beta_highest_t975)))).+t0_975CI
    plot!(hiatus_thermal_subs_curve, curve_1_ages, -Sμ, linecolor =:purple, label = "lowest beta and t0 at 025CI")
    plot!(hiatus_thermal_subs_curve, curve_2_ages, -Sμ, linecolor =:red, label = "lowest beta and t0 at 975CI")
    plot!(hiatus_thermal_subs_curve, curve_3_ages, -Sμ, linecolor =:orange, label = "beta at 975CI and t0 at 025CI")
    plot!(hiatus_thermal_subs_curve, curve_4_ages, -Sμ, linecolor =:yellow, label = "beta at 975CI and t0 at 975CI")

    hiatus_age_depth_curve = plot([curve_1_ages; reverse(curve_4_ages)], [model_strat_heights_m; reverse(model_strat_heights_m)], fill = (round(Int,minimum(mdl_height)),0.5,:white), label = "Interpolated model results - 2σ error envelope")
    plot!(Age[2:4], Height[2:4], xerror=Age_sigma[2:4],label="Age data",seriestype=:scatter,color=:black, legend =:bottomleft, grid = false, dpi = 300)
    plot!(hiatus_age_depth_curve, curve_2_ages, model_strat_heights_m, linecolor =:red, label = "lowest beta and t0 at 975CI")
    plot!(hiatus_age_depth_curve, curve_3_ages, model_strat_heights_m, linecolor =:orange, label = "beta at 975CI and t0 at 025CI")
    plot!(hiatus_age_depth_curve, curve_4_ages, model_strat_heights_m, linecolor =:yellow, label = "beta at 975CI and t0 at 975CI")
    #plot!(age_depth_curve, mdl_age, mdl_height, linecolor=:blue, label="", fg_color_legend=:white)

    Sᵣ = reverse(Sμ)
    target_ages = [815.29, 811.51, 810.7, 788.7, 752.7, 739.9]
    target_age_errors = [0.32, 0.25, 6.3, 0.24, 5.5, 6.1]
    target_heights = [-(2.232-0.25), -(2.232-0.425), -(2.232-0.435), -(2.232-1.14), -(2.232-1.93), -(2.232-2.125)]
    target_heights_m = copy(target_heights)*1000
    target_height_errors = [0.212, 0.1, 0.05, 0.2, 0.035, 0.035]
    target_height_errors_m = copy(target_height_errors)*1000
    target_subs_heights = [-Sᵣ[250], -Sᵣ[425], -Sᵣ[435], -Sᵣ[1140], -Sᵣ[1930], -Sᵣ[2125]]
    target_subs_height_errors = [(Sᵣ[250]-Sᵣ[38]), (Sᵣ[425]-Sᵣ[325]), (Sᵣ[435]-Sᵣ[385]), (Sᵣ[1140]-Sᵣ[940]), (Sᵣ[1930]-Sᵣ[1895]), (Sᵣ[2125]-Sᵣ[2090])]

    scatter!(hiatus_age_depth_curve, target_ages, target_heights_m, xerr = target_age_errors, yerr = target_height_errors_m, label = "correlated ages")
    scatter!(hiatus_thermal_subs_curve, target_ages, target_subs_heights, xerr = target_age_errors, yerr = target_subs_height_errors, label = "correlated ages")
=#

## --- End of File
