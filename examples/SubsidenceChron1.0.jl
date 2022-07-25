## For subsidence modeling in extensional basins

## --- Load required pacages, install Chron if required

    try
        using Chron
    catch
        using Pkg
        Pkg.add(PackageSpec(url="https://github.com/brenhinkeller/Chron.jl"))
        using Chron
    end

    #import Pkg;
    #Pkg.add(["StatGeochem","Distributions","Plots","StatsBase","ProgressMeter","LsqFit","KernelDensity","Interpolations","SpecialFunctions"])
    using StatGeochem, Distributions, Plots, Statistics, StatsBase, SpecialFunctions
    #using StatsBase: fit, Histogram, percentile
    #using ProgressMeter: @showprogress, Progress, update!
    #using LsqFit: curve_fit
    #using KernelDensity: kde
    #using Interpolations 
    using Plots; gr();

## --- Define properties of the stratigraphy

    # # # # # # # # # # # Enter stratigraphic information here! # # # # # # # # # # # #
    data_csv = importdataset("examples/Test_DB_PerfectSubsidence.csv",',')
    nLayers = length(data_csv["Thickness"])

    strat = NewStratData(nLayers)
    strat.Lithology          = data_csv["Lithology"] 
    strat.Thickness         .= data_csv["Thickness"] 

    nsims = 1000
    res = 0.02

## --- Run the decompaction and backstripping MC model
    @time (Sₜ, Sμ, Sσ, model_strat_heights) = DecompactBackstrip(strat, nsims, res)

    # Plot results - tectonic subsidence in comparison with present day stratigraphic heights
    p1 = plot(1:541, Sμ, alpha = 1, yflip = true, xflip = true, label = "Tectonic subsidence", color = "blue")
    plot!(p1, 1:541, reverse((model_strat_heights)*1000), yflip = true, label = "Present-day thickness", color = "red")
    plot!(p1, 1:541, Sₜ[:,2:end], alpha = 0.01, label = "", yflip = true, color = "blue", fg_color_legend=:white)

    #Test Plot 1: how the distribution of Sₜ compare with target mean value?
    using DelimitedFiles
    Sₜ_test = readdlm("St_test.txt", ',', Float64)*1000
    testplot1_1 = histogram(Sₜ[1,:], label = "Distribution from MC", linecolor = "white")
    vline!(testplot1_1, [Sₜ_test[1]], label = "Target value", linecolor = "black", linewidth = 2)
    testplot1_2 = histogram(Sₜ[121,:], label = "Distribution from MC", linecolor = "white")
    vline!(testplot1_2, [Sₜ_test[121]], label = "Target value", linecolor = "black", linewidth = 2)
    testplot1_3 = histogram(Sₜ[361,:], label = "Distribution from MC", linecolor = "white")
    vline!(testplot1_3, [Sₜ_test[361]], label = "Target value", linecolor = "black", linewidth = 2)
    testplot1_4 = histogram(Sₜ[501,:], label = "Distribution from MC", linecolor = "white")
    vline!(testplot1_4, [Sₜ_test[501]], label = "Target value", linecolor = "black", linewidth = 2)
    testplot1 = plot(testplot1_1, testplot1_2, testplot1_3, testplot1_4, layout = 4, title = ["0 Ma" "80 Ma" "145 Ma" "245 Ma"], size = (1000, 1000))
    png(testplot1, "Test1_DecompactBackstripMC")

## --- Define properties of age constraints

    # # # # # # # # # # # Enter age constraint (sample) information here! # # # # # # # # # # # #
    # Input the number of samples we wish to model (must match below)
    nSamples = 4
    # Make an instance of a ChronSection object for nSamples
    smpl = NewChronAgeData(nSamples)
    smpl.Name          = ("Sample 1", "Sample 2", "Sample 3", "Sample 4") # Et cetera
    smpl.Age          .= [ 388.7,  362.8,  327.0,  106.8] # Measured ages
    smpl.Age_sigma    .= [  0.2,    0.2,    0.2,    0.2] # Measured 1-σ uncertainties
    smpl.Height       .= [ -1700,   -1100,   -700,   -100] # Depths below surface should be negative
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
    Beta = 1.4
    Beta_sigma = 0.2
    T0 = 400
    T0_sigma = 50

    therm = NewThermalSubsidenceParameters()
    therm.Param = [Beta, T0] 
    therm.Sigma = [Beta_sigma, T0_sigma] 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## --- Run stratigraphic model

    # # # # # # # # # # Configure stratigraphic model here! # # # # # # # # # #
    # Configure the stratigraphic Monte Carlo model
    config = NewStratAgeModelConfiguration()
    # If you in doubt, you can probably leave these parameters as-is
    config.resolution = res*1000 # Same units as sample height. Smaller is slower!
    config.bounding = 0.1 # how far away do we place runaway bounds, as a fraction of total section height. Larger is slower.
    (bottom, top) = extrema(smpl.Height)
    npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
    config.nsteps = 15000 # Number of steps to run in distribution MCMC #!!TRIAL RUN SETTING!! ACTUAL RUNS SHOULD BE 15000
    config.burnin = 30000*npoints_approx # Number to discard #!!TRIAL RUN SETTING!! ACTUAL RUNS SHOULD BE 10000
    config.sieve = round(Int,npoints_approx)*10 # Record one out of every nsieve steps
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    ## --- Option 1: Stratigraphic MCMC model (without hiatus)
        #Run the model
        @time (subsmdl, agedist, lldist, beta_t0dist, lldist_age, lldist_height, lldist_tsparam, lldist_subs) = SubsidenceStratMetropolis(smpl, config, therm, model_strat_heights, Sμ, Sσ, 0.05, -5)

        #Calculate modeled ages for target horizons
        target_height = [0.36, 0.98, 1.311, 1.84]
        target_index = findclosest(target_height, model_strat_heights)
        target_subs = copy(Sₜ[target_index,:])
        #target_subs_rand = Array{Float64,2}(undef, length(target_index), 20)
        #t_subset = Array{Float64,1}(undef, 20)
        beta_t0dist_filter = Array{Float64,2}(undef, 2, config.nsteps)
        idx = 1
        for i = 1:config.nsteps
            if beta_t0dist[1,i]<=subsmdl.Beta_975CI[1] && beta_t0dist[1,i]>=subsmdl.Beta_025CI[1] && beta_t0dist[2,i]<=subsmdl.T0_975CI[1] && beta_t0dist[2,i]>=subsmdl.T0_025CI[1]
                beta_t0dist_filter[1,idx] = beta_t0dist[1,i]
                beta_t0dist_filter[2,idx] = beta_t0dist[2,i]
                idx += 1
            end
        end
        beta_t0dist_95CI = beta_t0dist_filter[:,1:(idx-1)]
        
        beta_t0_sample_size = div((idx-1),30)
        Sₜ_sample_size = 30
        beta_t0_sampled = Array{Float64,2}(undef, 2, beta_t0_sample_size)
        predicted_ages = Array{Float64,2}(undef, length(target_index), beta_t0_sample_size*Sₜ_sample_size)
        τ = 50
        E₀ = 3165.6475782289444
        using ProgressMeter
        pgrs = Progress(beta_t0_sample_size, desc="Progress...")
        pgrs_interval = ceil(Int,sqrt(beta_t0_sample_size))
        for i = 1:beta_t0_sample_size
            beta_t0_sampled[:,i] = beta_t0dist_95CI[:,(rand(1:(idx-1)))] 
            print(beta_t0_sampled[:,i])
            for j = 1:length(target_index)
                k = 1
                while k <= Sₜ_sample_size
                    random_draw = rand(target_subs[j,:])
                    temp_calc = (random_draw*pi)/(E₀*beta_t0_sampled[1,i]*sin(pi/beta_t0_sampled[1,i]))
                    if temp_calc<1
                        predicted_ages[j,k+(i-1)*Sₜ_sample_size] = τ*log(1-temp_calc)+beta_t0_sampled[2,i]
                        k += 1
                    end
                end
            end
            mod(i,pgrs_interval)==0 && update!(pgrs, i)
        end
        update!(pgrs, beta_t0_sample_size)

        # Plot results (mean and 95% confidence interval for both model and data)
        hdl = plot([subsmdl.Age_025CI; reverse(subsmdl.Age_975CI)],[subsmdl.Height; reverse(subsmdl.Height)], fill=(round(Int,minimum(subsmdl.Height)),0.5,:blue), label="model")
        plot!(hdl, subsmdl.Age, subsmdl.Height, linecolor=:blue, label="", fg_color_legend=:white) # Center line
        t = smpl.Age_Sidedness .== 0 # Two-sided constraints (plot in black)
        any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],smpl.Age_975CI[t]-smpl.Age[t]),label="data",seriestype=:scatter,color=:black)
        t = smpl.Age_Sidedness .== 1 # Minimum ages (plot in cyan)
        any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],zeros(count(t))),label="",seriestype=:scatter,color=:cyan,msc=:cyan)
        any(t) && zip(smpl.Age[t], smpl.Age[t].+nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
        t = smpl.Age_Sidedness .== -1 # Maximum ages (plot in orange)
        any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(zeros(count(t)),smpl.Age_975CI[t]-smpl.Age[t]),label="",seriestype=:scatter,color=:orange,msc=:orange)
        any(t) && zip(smpl.Age[t], smpl.Age[t].-nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
        plot!(hdl, xlabel="Age ($(smpl.Age_Unit))", ylabel="Height ($(smpl.Height_Unit))")
        #plot!(hdl, [smpl.Age[1], smpl.Height[1]],[smpl.Age[3], smpl.Height[3]])
        savefig(hdl,"AgeDepthModel_test6-3_0721.pdf")
        #display(hdl)

        #Test Plot 2 (to see how well the model predicted beta and t0 matches the actual values):
        beta_range = 1:0.005:therm.Param[1]+therm.Sigma[1]*3 #three sigmas
        t0_range = therm.Param[2]-therm.Sigma[2]*3:1:therm.Param[2]+therm.Sigma[2]*3 #three sigmas
        beta_pdf = pdf.(Normal(Beta, Beta_sigma), beta_range)
        t0_pdf = pdf.(Normal(T0, T0_sigma), t0_range)
        testplot2_1 = plot(beta_range, beta_pdf, linecolor = "red", label = "prior", legend=:topright)
        vline!([Beta], label = "prior mean (= actual t0)", linecolor = "red", linewidth = 2) #1.4+/-0.2
        histogram!(twinx(), beta_t0dist[1,:], label = "posterior", color ="blue", alpha = 0.5, legend=:topleft)
        vline!([nanmean(beta_t0dist,dims=2)[1]], label = "posterior mean", linecolor = "blue", linewidth = 2)
        vline!([1.4], label = "actual beta", linecolor = "black", linewidth =2)
        savefig(testplot2_1,"beta_test6-3_0721.pdf")

        testplot2_2 = histogram(beta_t0dist[2,:], label = "posterior", color = "blue", alpha = 0.5, legend=:topright)
        vline!([nanmean(beta_t0dist,dims=2)[2]], label = "posterior mean", linecolor = "blue", linewidth = 2) #120+/-20
        plot!(twinx(), t0_range, t0_pdf, linecolor = "red", label = "prior", legend=:topleft)
        vline!([T0], label = "prior mean", linecolor = "red", linewidth = 2)
        vline!([400], label = "actual beta", linecolor = "black", linewidth =2)
        savefig(testplot2_2,"t0_test6-3_0721.pdf")

        #Test Plot 3 (to see how ll and the different ll components change throughout the whole run):
        testplot3_1 = plot(lldist_subs, label = "Subsidence" , legend=:bottomright)
        plot!(testplot3_1, lldist_age, label = "Age")
        plot!(testplot3_1, lldist_tsparam, label = "TS Parameters")
        plot!(testplot3_1, lldist_height, label = "Height")
        png(testplot3_1, "Test1_1_LL_components")

        testplot3_2 = plot(lldist, label = "ll of accepted proposals" , legend=:bottomright)
        png(testplot3_2, "Test1_1_LL")

        #Test Plot 4 (to see how well is the model doing matching the ideal subsidence curve)
        Sₜ_025CI = nanpctile(Sₜ, 2.5, dims = 2)[6:86]
        Sₜ_975CI = nanpctile(Sₜ, 97.5, dims = 2)[6:86]
        Sμ_crop = Sμ[6:86]
        Sμ_sample = [Sμ[6],Sμ[36],Sμ[56],Sμ[86]]
        testplot4 = plot([subsmdl.Age; reverse(subsmdl.Age)],[reverse(Sₜ_025CI); Sₜ_975CI], fill=(round(Int,minimum(Sₜ_025CI)),0.4,:blue), label="model")
        plot!(testplot4, subsmdl.Age, reverse(Sμ_crop), linecolor=:blue, label="", fg_color_legend=:white) # Center line
        t = smpl.Age_Sidedness .== 0 # Two-sided constraints (plot in black)
        any(t) && plot!(testplot4, smpl.Age[t], reverse(Sμ_sample)[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],smpl.Age_975CI[t]-smpl.Age[t]), label="data",seriestype=:scatter,color=:black)
        t = smpl.Age_Sidedness .== 1 # Minimum ages (plot in cyan)
        any(t) && plot!(testplot4, smpl.Age[t], reverse(Sμ_sample)[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],zeros(count(t))),label="",seriestype=:scatter,color=:cyan,msc=:cyan)
        any(t) && zip(smpl.Age[t], smpl.Age[t].+nanmean(smpl.Age_sigma[t])*4, reverse(Sμ_sample[t])) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
        t = smpl.Age_Sidedness .== -1 # Maximum ages (plot in orange)
        any(t) && plot!(testplot4, smpl.Age[t], reverse(Sμ_sample)[t], xerror=(zeros(count(t)),smpl.Age_975CI[t]-smpl.Age[t]),label="",seriestype=:scatter,color=:orange,msc=:orange)
        any(t) && zip(smpl.Age[t], smpl.Age[t].-nanmean(smpl.Age_sigma[t])*4, reverse(Sμ_sample)[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
        plot!(testplot4, xlabel="Age ($(smpl.Age_Unit))", ylabel="Tectonic Subsidence ($(smpl.Height_Unit))")
        E₀ = 3165.6475782289444
        beta = nanmean(beta_t0dist, dims=2)[1]
        t0 = nanmean(beta_t0dist, dims=2)[2]
        Sμ_calculated = (E₀*beta/pi)*sin(pi/beta).*(1 .-exp.(-(t0 .-subsmdl.Age)./50))
        plot!(testplot4, subsmdl.Age, Sμ_calculated, linecolor=:purple, label="TS curve based on posterior")
        Sμ_ideal = (E₀*1.4/pi)*sin(pi/1.4).*(1 .-exp.(-(400 .-subsmdl.Age)./50))
        plot!(testplot4, subsmdl.Age, Sμ_ideal, linecolor=:red, label="actual TS curve")
        savefig(testplot4, "SubsidenceCurveComparison_test6_1_0721.pdf")
    
        #Test Plot 5 (to see how the predicted ages for horizons-of-interest match the true values)
        testplot5 = histogram(predicted_ages[1,:], label = "predicted age for horizon 1", color ="blue", alpha = 0.5, legend=:topleft)
        vline!([282.60], label = "actual age for horizon 1", linecolor = "black", linewidth = 2) 
        vline!([nanmedian(predicted_ages,dims=2)[1]], label = "posterior median", linecolor = "blue", linewidth = 2)
        savefig(testplot5,"PredictedAge_h1_test6-1_0724.pdf")

#=
    # More plots - For Svalbard - still need to readjust the format

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

    ## --- Option 2: Stratigraphic MCMC model including a hiatus with unknown duration
        #Input strat height of hiatus
        hiatus_height = -500

        #Run the model with the additional term of "hiatus_height"
        @time (subsmdl, agedist, lldist, beta_t0dist) = SubsidenceStratMetropolis(smpl, config, therm, model_strat_heights, Sμ, Sσ, hiatus_height, 0.05, -0.5)

        # Plot results (mean and 95% confidence interval for both model and data)
        hdl = plot([subsmdl.Age_025CI; reverse(subsmdl.Age_975CI)],[subsmdl.Height; reverse(subsmdl.Height)], fill=(round(Int,minimum(subsmdl.Height)),0.5,:blue), label="model")
        plot!(hdl, subsmdl.Age, subsmdl.Height, linecolor=:blue, label="", fg_color_legend=:white) # Center line
        t = smpl.Age_Sidedness .== 0 # Two-sided constraints (plot in black)
        any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],smpl.Age_975CI[t]-smpl.Age[t]),label="data",seriestype=:scatter,color=:black)
        t = smpl.Age_Sidedness .== 1 # Minimum ages (plot in cyan)
        any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],zeros(count(t))),label="",seriestype=:scatter,color=:cyan,msc=:cyan)
        any(t) && zip(smpl.Age[t], smpl.Age[t].+nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
        t = smpl.Age_Sidedness .== -1 # Maximum ages (plot in orange)
        any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(zeros(count(t)),smpl.Age_975CI[t]-smpl.Age[t]),label="",seriestype=:scatter,color=:orange,msc=:orange)
        any(t) && zip(smpl.Age[t], smpl.Age[t].-nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
        plot!(hdl, xlabel="Age ($(smpl.Age_Unit))", ylabel="Height ($(smpl.Height_Unit))")
        savefig(hdl,"AgeDepthModel_hiatus_test7-1_0724.pdf")
        display(hdl)


    ## --- Option 3: Stratigraphic MCMC model including hiata with known durations
        # Data about hiata
        nHiatuses = 2 # The number of hiata you have data for
        hiatus = NewHiatusData(nHiatuses) # Struct to hold data
        hiatus.Height         = [-371.5, -405.0 ]
        hiatus.Height_sigma   = [   0.0,    0.0 ]
        hiatus.Duration       = [ 100.0,   123.0]
        hiatus.Duration_sigma = [  30.5,    20.0]

        # Run the model. Note the additional `hiatus` arguments
        @time (subsmdl, agedist, hiatusdist, lldist, beta_t0dist) = StratMetropolis(smpl, config, therm, model_strat_heights, Sμ, Sσ, hiatus); sleep(0.5)

        # Plot results (mean and 95% confidence interval for both model and data)
        hdl = plot([subsmdl.Age_025CI; reverse(subsmdl.Age_975CI)],[subsmdl.Height; reverse(subsmdl.Height)], fill=(minimum(subsmdl.Height),0.5,:blue), label="model")
        plot!(hdl, subsmdl.Age, subsmdl.Height, linecolor=:blue, label="", fg_color_legend=:white)
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
