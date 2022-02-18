## --- Stratigraphic MCMC model without hiatus # # # # # # # # # # # # # # # # #

import Pkg;
Pkg.add(["StatGeochem","Distributions","Plots"])
using StatGeochem, Distributions, Plots
using StatsBase: fit, Histogram, percentile
using ProgressMeter: @showprogress, Progress, update!
using LsqFit: curve_fit
using KernelDensity: kde
using Interpolations


# Part 1: decompaction and backstripping 

    # Define the decompaction function
        function decompact!(yₚ, y, ϕ₀, c, n, m; niterations = 10)
            # y[n+1] = y₂ = present-day depth for base of layer n 
            # y[n] = y₁ = present-day depth for top of layer n
            δy = y[n+1] - y[n]
            # yₚ[n-m+2] = y₂' = depth for base of layer n during decompaction (at t=tₘ)
            # yₚ[n-m+1] = y₁' = depth for top of layer n during decompaction (at t=tₘ)
            yₚ[n-m+2] = yₚ[n-m+1] + δy 
            # In this for loop, the initial approximatation for the thickness of a given layer is the present-day thickness of that layer (calculated above)
            # The true thickness is then reached by running the decompaction equation multiple times (niterations)
            @inbounds for idx=1:niterations
            yₚ[n-m+2] = yₚ[n-m+1] + δy + ϕ₀/c * ((exp(-c*yₚ[n-m+1]) - exp(-c*yₚ[n-m+2])) - (exp(-c*y[n]) - exp(-c*y[n+1])))
            end
        end

    # Import data from csv and assign parameters for each lithology
        # Import data
        #=
        @time data = importdataset("Svalbard.csv",',')
        @time lithology = data["Lithology"]
        thickness = data["Thickness"]::Array{Float64,1}
        strat_heights = cumsum([0; thickness])
        layer_count = length(lithology)
        =#
        @time model_data_csv = importdataset("Svalbard.csv",',')
        model_lithology_csv = model_data_csv["Lithology"]
        model_thickness_csv = model_data_csv["Thickness"]
        model_strat_heights_csv = cumsum([0; model_thickness_csv])
        layer_count_csv = length(model_strat_heights_csv)-1
        model_strat_heights_km = [0:0.001:maximum(model_strat_heights_csv);]
        model_thickness_max = maximum(model_strat_heights_km)
        model_layer_count = length(model_strat_heights_km)-1
        model_strat_heights_m = copy(model_strat_heights_km[2:end])*1000

    # Allocate parameters as distributions; each element/distribution represents a layer
        # porosity depth coefficient(c)
        c_dist = Array{Distribution,1}(undef, layer_count_csv)
        # surface porosity (ϕ₀)
        ϕ₀_dist = Array{Distribution,1}(undef, layer_count_csv)
        # sediment grain density (ρg)
        ρg = Array{Float64,1}(undef, layer_count_csv)

    # Find the correct c, ϕ₀, and ρg for each layer based on lithology
        for i = 1:layer_count_csv
            if model_lithology_csv[i] == "Shale"
                # c has a lower bound of 0 b/c porosity at depth should not be greater than porosity at surface
                c_dist[i] = truncated(LogNormal(log(0.51),0.15), 0, Inf)
                # ϕ₀ has a lower bound of 0 and an upper bound of 1 b/c of the definition of porosity
                ϕ₀_dist[i] = truncated(Normal(0.63, 0.15), 0, 1)
                ρg[i] = 2660
            elseif model_lithology_csv[i] == "Sandstone"
                c_dist[i] = truncated(LogNormal(log(0.27), 0.1), 0, Inf)
                ϕ₀_dist[i] = truncated(Normal(0.49, 0.1), 0, 1)
                ρg[i] = 2650
            #elseif model_lithology_csv[i] == "Chalk"
                #c_dist[i] = truncated(LogNormal(log(0.71), 0.15), 0, Inf)
                #ϕ₀_dist[i] = truncated(Normal(0.7, 0.15), 0, 1)
                #ρg[i] = 2710
            elseif model_lithology_csv[i] == "Limestone"
                c_dist[i] = truncated(LogNormal(log(0.6), 0.2), 0, Inf)
                ϕ₀_dist[i] = truncated(Normal(0.4, 0.17), 0, 1)
                ρg[i] = 2710
            elseif model_lithology_csv[i] == "Dolostone"
                c_dist[i] = truncated(LogNormal(log(0.6), 0.2), 0, Inf)
                ϕ₀_dist[i] = truncated(Normal(0.2, 0.1), 0, 1)
                ρg[i] = 2870
            #elseif model_lithology_csv[i] == "Anhydrite"
                #c_dist[i] = truncated(LogNormal(log(0.2), 0.1), 0, Inf)
                #ϕ₀_dist[i] = truncated(Normal(0.05, 0.05), 0, 1)
                #ρg[i] = 2960
            elseif model_lithology_csv[i] == "Quartzite"
                c_dist[i] = truncated(LogNormal(log(0.3), 0.1), 0, Inf)
                ϕ₀_dist[i] = truncated(Normal(0.2, 0.1), 0, 1)
                ρg[i] = 2650
            end
        end

    # Prep for decompaction and backstripping MC
        # Define parameters for decompaction (number of simulations; water and mantle densities)
        nsims = 1000
        ρw = 1030
        ρm = 3330
           
        # Allocate depth matricies (rows are strat horizons, columns are timesteps)
        # decompacted depth
        Y = fill(1E-18, (model_layer_count+1, model_layer_count))
        # decompacted depth after applying varies corrections
        Y_corr = fill(1E-18, (model_layer_count+1, model_layer_count))
        Y_max = fill(1E-18, (nsims, model_layer_count))

        c_highres = Array{Float64,1}(undef, model_layer_count)
        ϕ₀_highres = Array{Float64,1}(undef, model_layer_count)
        ρg_highres = Array{Float64,1}(undef, model_layer_count)
        
        #Allocate porosity, density and tectonic subsidence matricies   
        # porosity of a strat unit at any depth
        ϕ_avg = fill(1E-18, (model_layer_count, model_layer_count))
        # bulk density of a single layer
        ρ_bulk = fill(1E-18, (model_layer_count, model_layer_count))
        # intermediate step - bulk density*thickness of a single layer
        m_bulk = fill(1E-18, (model_layer_count, model_layer_count))
        # bulk density of the entire column
        ρ_bulk_column = Array{Float64,1}(undef, model_layer_count)
        # tectonic subsidence # this is the only one that will need to propagate outside of the loop
        Sₜ_km = Array{Float64,2}(undef, model_layer_count, nsims)
    
    # ***** Skip this block when testing MCMC - results are stored in txt files already *****
    # MC for decompaction and backstripping
        print("Decompaction and Backstripping: ", nsims, " steps\n")
        pgrs = Progress(nsims, desc="Decompaction and Backstripping...")
        pgrs_interval = ceil(Int,10)    
        @time for sim = 1:nsims
            # randomly select c and ϕ₀ vectors (for the whole column) from the distributions
            c = rand.(c_dist)
            ϕ₀ = rand.(ϕ₀_dist)

            for i = 1:layer_count_csv
                for j = 1:model_layer_count
                    if model_strat_heights_km[j+1]>model_strat_heights_csv[i]
                        c_highres[j]=c[i]
                        ϕ₀_highres[j]=ϕ₀[i]
                        ρg_highres[j]=ρg[i]
                    end
                end
            end

            # Fill the first column with modern observed values (present-day depths)
            Y[:,1] .= model_strat_heights_km
            # Fill the first row with zeros
            Y[1,:] .= 0

            # Decompact
            # i = time steps during decompaction, which runs from 2 to layer_count b/c column 1 is present day
            # j = layer number, which runs from i to layer_count b/c the ith column begins with y₁' of layer i
            for i = 2:model_layer_count
                for j = i:model_layer_count
                    decompact!(view(Y,:,i), model_strat_heights_km, ϕ₀_highres[j], c_highres[j], j, i)
                end
            end

            # Apply decompaction corrections
            # vector for paleo sea-level corrections (each element represents a time step)
            paleo_sl = zeros(1, model_layer_count) #using zeros for now but should turn this into a distribution
            # vector for paleo water-depth corrections (each element represents a time step)
            paleo_wd = zeros(1, model_layer_count) #using zeros for now but should correlate with lithology?
            Y_corr[:,:] .= Y[:,:].-paleo_sl.*(ρw/(ρm-ρw)).+paleo_wd.-paleo_sl
            Y_max[sim, :] = maximum(Y_corr, dims = 1)
            
            # Remove the effect of sediment load - same indexing logic as the decompaction loop
            # i = time steps (columns) = 1:layer_count b/c need to calculate these parameters for the present day column/layers
            # j = layer number = i:layer_count b/c the ith column begins with y₁' of layer i
            for i = 1:model_layer_count
                for j = i:model_layer_count
                    ϕ_avg[j-i+1,i] = (ϕ₀_highres[j]/c_highres[j])*(exp(-c_highres[j]*Y_corr[j-i+1,i])-exp(-c_highres[j]*Y_corr[j-i+2,i]))/(Y_corr[j-i+2,i]-Y_corr[j-i+1,i])
                    ρ_bulk[j-i+1,i] = ρw*ϕ_avg[j-i+1,i]+(1-ϕ_avg[j-i+1,i])*ρg_highres[j]
                    m_bulk[j-i+1,i] = (Y_corr[j-i+2,i]-Y_corr[j-i+1,i])*ρ_bulk[j-i+1,i]
                end
                    # maximum(Y_corr[:,i]) = total depth of column at time step i
                    ρ_bulk_column[i] = sum(m_bulk[:,i])/maximum(Y_corr[:,i])
                    Sₜ_km[i,sim] = maximum(Y_corr[:,i])*((ρm-ρ_bulk_column[i])/(ρm-ρw))        
            end
            mod(sim,pgrs_interval)==0 && update!(pgrs, sim)
        end
        update!(pgrs,nsims)

    # Calculate summary statistics (mean and standard deviation)
        Sₜ = copy(Sₜ_km).*1000
        Sμ = dropdims(nanmean(Sₜ, dims=2), dims=2)
        Sσ = dropdims(nanstd(Sₜ, dims=2), dims=2)

        using DelimitedFiles
        writedlm("St_120221.txt", Sₜ)
        writedlm("St_mu_120221.txt", Sμ)
        writedlm("St_sigma_120221.txt", Sσ)
    
    # ***** End of the block that will be skipping when testing MCMC

    # !!!MISSING STEP: convert txts to csvs!!!

    # Read the results from the decompaction-backstripping loop from these txts
        using DelimitedFiles
        Sₜ = readdlm("St_120221.csv", ',', Float64)
        Sμ = readdlm("St_mu_120221.csv", ',', Float64)
        Sσ = readdlm("St_sigma_120221.csv", ',', Float64)

        #Sσ_cumulative = cumsum(Sσ, dims = 1)

        p1 = plot(1:2232, Sμ, alpha = 1, yflip = true, xflip = true, label = "Tectonic subsidence", color = "blue")
        plot!(p1, 1:2232, reverse(model_strat_heights_m), yflip = true, label = "Present-day thickness", color = "red")
        plot!(p1, 1:2232, Sₜ[:,2:end], alpha = 0.01, label = "", yflip = true, color = "blue", fg_color_legend=:white)


# Part 2: Modified StratMetropolis

    # Define the function that calculates the ll of the thermal subsidence fit
        
        # Thermal subsidence model parameters
        y_litho= 125000
        ρ_mantle = 3330
        ρ_water = 1000
        αᵥ = 3.28*10^(-5)
        T_mantle = 1333
        τ = 50 #Myr
        E₀ = (4*y_litho*ρ_mantle*αᵥ*T_mantle)/(pi^2*(ρ_mantle-ρ_water))
    
        function subsidence_ll(model_St, model_St_sigma, model_t, beta_t0)
            # Calculate subsidence_model_heights given this unique smooth thermal subsidence model at each age in model_t
            subsidence_model_heights = (E₀*beta_t0[1]/pi)*sin(pi/beta_t0[1]).*(1 .-exp.(-(beta_t0[2] .-model_t)./τ))
            #-τ*log(1 .-((model_St*pi)/(E₀*beta_t0[1]*sin(pi/beta_t0[1])))) .+beta_t0[2]
            # Turn that into a log likelihood using some age uncertainty of the curve
            log_likelihood = normpdf_ll(model_St, model_St_sigma, subsidence_model_heights)
            return log_likelihood
        end

        #= 
        #From here to the end of Strat MCMC (w/o hiatus) used to be under:
        #function StratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration)

        # Run stratigraphic MCMC model
        #print("Generating stratigraphic age-depth model...\n")

        # Stratigraphic age constraints
        Age = copy(smpl.Age)::Array{Float64,1}
        Age_sigma = copy(smpl.Age_sigma)::Array{Float64,1}
        Height = copy(smpl.Height)::Array{Float64,1}
        Height_sigma = smpl.Height_sigma::Array{Float64,1} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Array{Float64,1} # Bottom is a maximum age and top is a minimum age
        (bottom, top) = extrema(Height)
        (youngest, oldest) = extrema(Age)
        dt_dH = (oldest-youngest)/(top-bottom)
        aveuncert = nanmean(Age_sigma)
        model_heights = bottom:resolution:top
        =#

        #= 
        Instead I'll do the following for now to let the rest of the model run
        Age = [791.1; 737.5; 717]
        # One sigma on age constraints
        Age_sigma = [2.45; 4.8; 0.4]
        Height = [-model_strat_heights_km[1284]; -model_strat_heights_km[121]; -model_strat_heights_km[1]].*1000
        Height_sigma = fill(0.01, 3)
        Age_Sidedness = zeros(3)
        (bottom, top) = extrema(Height)
        (youngest, oldest) = extrema(Age)
        dt_dH = (oldest-youngest)/(top-bottom)
        aveuncert = nanmean(Age_sigma)
        model_heights = copy(-model_strat_heights_km[2:end]).*1000
        =#

        #Now let's try something else - what if we just use two data points (excluding the 717 age)
        Age = [791.1; 737.5]
        Age_sigma = [2.45; 4.8]
        Height = [-model_strat_heights_km[1284]; -model_strat_heights_km[121]].*1000
        Height_sigma = fill(0.01, 2)
        Age_Sidedness = zeros(2)
        (bottom, top) = extrema(Height)
        (youngest, oldest) = extrema(Age)
        dt_dH = (oldest-youngest)/(top-bottom)
        aveuncert = nanmean(Age_sigma)
        model_heights = copy(-model_strat_heights_km[2:end]).*1000

        # Model configuration -- read from struct (Current version: Using direct inputs for now)
        resolution = 1
        bounding = 0.1
        npoints_approx = round(Int,length(bottom:resolution:top) * (1 + 2*bounding))
        nsteps = 1500 # Number of steps to run in distribution MCMC # !!!DIVIDED BY 10 FOR TRIAL RUNS
        burnin = 3000*npoints_approx# Number to discard # !!!DIVIDED BY 10 FOR TRIAL RUNS
        sieve = round(Int,npoints_approx) # Record one out of every nsieve steps
        
        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            offset = round((top-bottom)*bounding)
            Age = [oldest + offset*dt_dH; Age; youngest - offset*dt_dH]
            Age_sigma = [nanmean(Age_sigma)/10; Age_sigma; nanmean(Age_sigma)/10]
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
        end
        
        active_height_t = (model_heights .>= bottom) .& (model_heights .<= top)
        npoints = length(model_heights)

        # STEP 1: calculate log likelihood of the modeled ages (and heights) in the initial proposal
        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        model_ages = a .+ b .* collect(model_heights)

        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        # Don't actually need this step anymore!
        closest = findclosest(sample_height, model_heights)
        closest_model_ages = model_ages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_model_ages[i] - Age[i])
                closest_model_ages[i] = Age[i]
            end
        end

        #New addition - weight of each term 
        w_age = 1
        w_height = 1
        w_parameters = 1
        w_subsidence = 1

        ll = w_age*normpdf_ll(Age, Age_sigma, closest_model_ages)
        ll += w_height*normpdf_ll(Height, Height_sigma, sample_height)

        # STEP 2: calculate log likelihood of the subsidence model parameters in the initial proposal
        # Define thermal subsidence model parameters 
        ideal_subs_parameters = [1.3, 820]
        ideal_subs_parameters_sigma = [0.2,20]
                
        # Initial proposal for subsidence parameters - randomly pick a set of values from the distribution
        subs_parameters = [ideal_subs_parameters[1]+0.01, ideal_subs_parameters[2]-1]
        # Calculate log likelihood of this initial proposal for subsidence parameters
        ll += w_parameters*normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parameters)

        # STEP 3: calculate log likelihood for the fit of the thermal subsidence curve in the initial proposal
        ll += w_subsidence*subsidence_ll(Sμ, Sσ, model_ages, subs_parameters)/(length(model_strat_heights_km)-1)

        # Preallocate variables for MCMC proposals
        llₚ = ll
        model_agesₚ = copy(model_ages)
        closestₚ = copy(closest)
        sample_heightₚ = copy(sample_height)
        closest_model_agesₚ = copy(closest_model_ages)
        subs_parametersₚ = copy(subs_parameters)
        
        # Run burnin
        # acceptancedist = fill(false,burnin)
        print("Burn-in: ", burnin, " steps\n")
        pgrs = Progress(burnin, desc="Burn-in...")
        pgrs_interval = ceil(Int,sqrt(burnin))
        for n=1:burnin
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)
            copyto!(subs_parametersₚ, subs_parameters)

            # Propose adjustment to subsidence_parametersₚ
            subs_parametersₚ += randn() *ideal_subs_parameters_sigma

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points that are still stratigraphically below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points that are still stratigraphically above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = w_age*normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += w_height*normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += w_subsidence*subsidence_ll(Sμ, Sσ, model_agesₚ, subs_parametersₚ)/(length(model_strat_heights_km)-1)
            llₚ += w_parameters*normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)
            
            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
                # acceptancedist[i] = true
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs, burnin) # Finalize

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)
        beta_t0dist = Array{Float64}(undef,2,nsteps)

        # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        pgrs_interval = ceil(Int,sqrt(nsteps*sieve))
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)
            copyto!(subs_parametersₚ, subs_parameters)

            # Propose adjustment to subsidence_parametersₚ
            subs_parametersₚ += randn() *ideal_subs_parameters_sigma

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = w_age*normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += w_height*normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += w_subsidence*subsidence_ll(Sμ, Sσ, model_agesₚ, subs_parametersₚ)/(length(model_strat_heights_km)-1)
            llₚ += w_parameters*normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= model_ages
                beta_t0dist[:,n÷sieve] .= subs_parameters
                #predicted_ages[n÷sieve] = beta_t0dist[2]+τ*log(1-(Sμ[?]*pi)/(E₀*beta_t0dist[1]*sin(pi/beta_t0dist[2])))
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs,nsteps*sieve)

        # Crop the result
        agedist = agedist[active_height_t,:]

        using DelimitedFiles
        writedlm("agedist.txt", agedist)
        writedlm("beta_t0dist.txt", beta_t0dist)
        writedlm("lldist.txt", lldist)

        agedist = readdlm("agedist.csv", ',', Float64)
        beta_t0dist = readdlm("beta_t0dist.csv", ',', Float64)
        lldist = readdlm("lldist.csv", ',', Float64)

        #mdl_age_025CI = pctile(agedist,2.5, dim=2)
        #mdl_age_975CI = pctile(agedist,97.5, dim=2)
        mdl_height = model_heights[active_height_t]
        mdl_age = nanmean(agedist,dim=2)
        mdl_median = nanmedian(agedist, dim=2)
        mdl_age_025CI = nanpctile(agedist,2.5, dim=2)
        mdl_age_975CI = nanpctile(agedist, 97.5, dim=2)
        beta = nanmean(beta_t0dist, dims=2)[1]
        t0 = nanmean(beta_t0dist, dims=2)[2]
        beta_025CI = nanpctile(beta_t0dist[1,:],2.5)
        beta_975CI = nanpctile(beta_t0dist[1,:],97.5)
        t0_025CI = nanpctile(beta_t0dist[2,:],2.5)
        t0_975CI = nanpctile(beta_t0dist[2,:],97.5)
        beta_σ = nanstd(beta_t0dist, dims=2)[1]
        t0_σ = nanstd(beta_t0dist, dims=2)[2]

        #exclude out of bound betas (i.e. the beta values that won't be able to reach 871.99 thermal subsidence at top)
        beta_lowest = beta_025CI 
        beta_highest = beta
        for i = beta_025CI:0.001:beta
            while (maximum(Sμ)*pi)/(E₀*i*sin(pi/i))>1
                beta_lowest = i+0.001
            end
        end

        for i = beta_lowest:0.001:beta
            while τ*log(1-((maximum(Sμ)*pi)/(E₀*i*sin(pi/i))))+t0_975CI<(717-0.8)
                beta_lowest = i+0.001
            end
        end
        
        for i = beta_975CI:-0.001:beta_lowest
            while τ*log(1-((maximum(Sμ)*pi)/(E₀*i*sin(pi/i))))+t0_975CI>(737.5+9.6)
                beta_highest = i-0.001
            end
        end
        beta_lowest_t975 = 1.3403078
        beta_highest_t975 = 1.3831634
        beta_lowest_t025 = 1.36839353
        beta_highest_t025 = copy(beta_975CI)

        hdl = plot([mdl_age_025CI; reverse(mdl_age_975CI)],[mdl_height; reverse(mdl_height)], fill=(round(Int,minimum(mdl_height)),0.5,:blue), label="Model results - 2σ error envelope")
        plot!(hdl, mdl_age, mdl_height, linecolor=:blue, label="Model results - mean", fg_color_legend=:white) # Center line
        t = Age_Sidedness .== 0 # Two-sided constraints (plot in black)
        any(t) && plot!(hdl, Age[t], Height[t], xerror=Age_sigma[t],label="Input - age data",seriestype=:scatter,color=:black)
        age_at_height = (mdl_height,mdl_age,-2232)
        
        figpath = "C:/Users/zhtia/Documents/GitHub/Chron.jl/trials/"
        savefig(hdl, figpath*"Age_depth_curve_upper_section_T15.png")

        curve_ages = similar(Sμ)
        curve_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta*sin(pi/beta)))).+t0
        model_strat_heights_m = copy(-model_strat_heights_km[2:end]).*1000
        plot!(hdl, curve_ages, model_strat_heights_m, linecolor =:blue, label = "beta and t0", xlabel="Age (Ma) calculated based on beta and t0 from MCMC", ylabel="Thermal subsidence (m)")
        savefig(hdl, figpath*"Age_depth_curve_upper_section_comparison_T15.png")

        # Let's see with two age constraints, what would the MCMC predicted beta and t0 give us in terms of the age at the top of section
        curve_ages = τ.*log.(1 .-((Sμ[1].*pi)./(E₀*beta*sin(pi/beta)))).+t0

        #testing - what should beta and t0 be in order for the curves to match
        #beta = 1.341213987308592
        #t0 = 819.1213987308635
        y_litho= 125000
        ρ_mantle = 3330
        ρ_water = 1000
        αᵥ = 3.28*10^(-5)
        T_mantle = 1333
        τ = 50 #Myr
        E₀ = (4*y_litho*ρ_mantle*αᵥ*T_mantle)/(pi^2*(ρ_mantle-ρ_water))

        #curve_heights = 0:1:872
        #curve_heights_current = copy(model_strat_heights_km[2:end]).*1000
        curve_ages = similar(Sμ)
        curve_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta*sin(pi/beta)))).+t0
        model_strat_heights_m = copy(-model_strat_heights_km[2:end]).*1000
        
        thermal_subs_curve = plot(curve_ages, -Sμ, linecolor =:blue, xflip = true, label = "beta and t0", xlabel="Age (Ma) calculated based on beta and t0 from MCMC", ylabel="Thermal subsidence (m)")
        curve_1_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_025CI*sin(pi/beta_025CI)))).+t0_025CI
        curve_2_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_025CI*sin(pi/beta_025CI)))).+t0_975CI
        curve_3_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_975CI*sin(pi/beta_975CI)))).+t0_025CI
        curve_4_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_975CI*sin(pi/beta_975CI)))).+t0_975CI
        #curve_1_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_lowest_t025*sin(pi/beta_lowest_t025)))).+t0_025CI
        #curve_2_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_lowest_t975*sin(pi/beta_lowest_t975)))).+t0_975CI
        #curve_3_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_highest_t025*sin(pi/beta_highest_t025)))).+t0_025CI
        #curve_4_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta_highest_t975*sin(pi/beta_highest_t975)))).+t0_975CI
        plot!(thermal_subs_curve, curve_1_ages, -Sμ, linecolor =:purple, label = "lowest beta and t0 at 025CI")
        plot!(thermal_subs_curve, curve_2_ages, -Sμ, linecolor =:red, label = "lowest beta and t0 at 975CI")
        plot!(thermal_subs_curve, curve_3_ages, -Sμ, linecolor =:orange, label = "beta at 975CI and t0 at 025CI")
        plot!(thermal_subs_curve, curve_4_ages, -Sμ, linecolor =:yellow, label = "beta at 975CI and t0 at 975CI")

        age_depth_curve = plot([curve_1_ages; reverse(curve_4_ages)], [model_strat_heights_m; reverse(model_strat_heights_m)], fill = (round(Int,minimum(mdl_height)),0.5,:white), label = "Interpolated model results - 2σ error envelope")
        plot!(Age[2:4], Height[2:4], xerror=Age_sigma[2:4],label="Age data",seriestype=:scatter,color=:black, legend =:bottomleft, grid = false, dpi = 300)
        plot!(age_depth_curve, curve_2_ages, model_strat_heights_m, linecolor =:red, label = "lowest beta and t0 at 975CI")
        plot!(age_depth_curve, curve_3_ages, model_strat_heights_m, linecolor =:orange, label = "beta at 975CI and t0 at 025CI")
        plot!(age_depth_curve, curve_4_ages, model_strat_heights_m, linecolor =:yellow, label = "beta at 975CI and t0 at 975CI")
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

        scatter!(age_depth_curve, target_ages, target_heights_m, xerr = target_age_errors, yerr = target_height_errors_m, label = "correlated ages")
        scatter!(thermal_subs_curve, target_ages, target_subs_heights, xerr = target_age_errors, yerr = target_subs_height_errors, label = "correlated ages")
        
        figpath = "C:/Users/zhtia/Documents/GitHub/Chron.jl/trials/"
        savefig(age_depth_curve, figpath*"Age_depth_curve_entire_section_T15.png")
        savefig(thermal_subs_curve, figpath*"Thermal_subs_curve_entire_section_T15.png")
        

        #calc_ages = similar(Sμ)
        #calc_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta*sin(pi/beta)))).+t0
        #plot!(curve, calc_ages, Sμ, linecolor =:black)

        subsidence_model_heights_MCMC = (E₀*beta/pi)*sin(pi/beta).*(1 .-exp.(-(t0.-mdl_age)./τ))
        mdl_age_plot = plot(mdl_age, subsidence_model_heights_MCMC)
        plot!(mdl_age_plot, mdl_age, reverse(Sμ[1:1284]), linecolor =:black)
        plot!(mdl_age_plot, mdl_age, reverse(Sμ[1:1284]), yerror = reverse(Sσ[1:1284]), linecolor =:black, alpha = 0.1)

        calc_age_plot = plot(calculated_age, reverse(Sμ))

        plot(Sₜ, alpha = 0.01, yflip = true, xflip = true, leg = false, xlabel="Time step", ylabel="Depth(km)", color = "blue")

        interpolated_distribution = Array{Float64}(undef,size(agedist,2))
        for i=1:size(agedist,2)
            interpolated_distribution[i] = linterp1s(Height,agedist[:,i],-2232)
        end
        hdl = histogram(interpolated_distribution, nbins=50, label="")
        Age_Unit = "Ma"
        height = -2232
        plot!(hdl, xlabel="Age ($(Age_Unit)) at height=$height", ylabel="Likelihood (unnormalized)")
        savefig(hdl, "Interpolated age distribution.pdf")
        display(hdl)


        hdl = plot([mdl_age_025CI; reverse(mdl_age_975CI)],[mdl_height; reverse(mdl_height)], fill=(round(Int,minimum(mdl_height)),0.5,:blue), label="model")
        plot!(hdl, mdl_age, mdl_height, linecolor=:blue, label="", fg_color_legend=:white) # Center line
        t = smplge_Sidedness .== 0 # Two-sided constraints (plot in black)
        any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],smpl.Age_975CI[t]-smpl.Age[t]),label="data",seriestype=:scatter,color=:black)
        t = Age_Sidedness .== 1 # Minimum ages (plot in cyan)
        any(t) && plot!(hdl, Age[t], Height[t], xerror=(Age[t]-smpl.Age_025CI[t],zeros(count(t))),label="",seriestype=:scatter,color=:cyan,msc=:cyan)
        any(t) && zip(smpl.Age[t], smpl.Age[t].+nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
        t = Age_Sidedness .== -1 # Maximum ages (plot in orange)
        any(t) && plot!(hdl, Age[t], Height[t], xerror=(zeros(count(t)),smpl.Age_975CI[t]-smpl.Age[t]),label="",seriestype=:scatter,color=:orange,msc=:orange)
        any(t) && zip(Age[t], Age[t].-nanmean(Age_sigma[t])*4, Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
        plot!(hdl, xlabel="Age ($(smpl.Age_Unit))", ylabel="Height ($(smpl.Height_Unit))")
        savefig(hdl,"AgeDepthModel.pdf")
        display(hdl)

        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian(agedist,dim=2), # Median age
            pctile(agedist,2.5,dim=2), # 2.5th percentile
            pctile(agedist,97.5,dim=2) # 97.5th percentile
        )

        return mdl, agedist, lldist
    



## --- Stratigraphic MCMC model with hiata # # # # # # # # # # # # # # # # # # #
    #Same idea - trying to add the same thermal subsidence modifications on the model that incorporates hiatus
    #= 
    function StratMetropolis(smpl::ChronAgeData, hiatus::HiatusData, config::StratAgeModelConfiguration)
        # Run stratigraphic MCMC model, with hiata
        print("Generating stratigraphic age-depth model...\n")
    =#
        # Thermal subsidence model parameters
        y_litho= 125000
        ρ_mantle = 3330
        ρ_water = 1000
        αᵥ = 3.28*10^(-5)
        T_mantle = 1333
        τ = 50 #Myr
        E₀ = (4*y_litho*ρ_mantle*αᵥ*T_mantle)/(pi^2*(ρ_mantle-ρ_water))

        # Define the function that calculates the ll of the thermal subsidence fit
        function subsidence_ll(model_St, model_St_sigma, model_t, beta_t0)
            # Calculate subsidence_model_heights given this unique smooth thermal subsidence model at each age in model_t
            subsidence_model_heights = (E₀*beta_t0[1]/pi)*sin(pi/beta_t0[1]).*(1 .-exp.(-(beta_t0[2] .-model_t)./τ))
            #-τ*log(1 .-((model_St*pi)/(E₀*beta_t0[1]*sin(pi/beta_t0[1])))) .+beta_t0[2]
            # Turn that into a log likelihood using some age uncertainty of the curve
            log_likelihood = normpdf_ll(model_St, model_St_sigma, subsidence_model_heights)
            return log_likelihood
        end

        # Stratigraphic age constraints. Type assertions for stability
        Age = [791.1; 737.5; 717]
        Age_sigma = [2.45; 4.8; 0.4]
        Height = [-model_strat_heights_km[1284]; -model_strat_heights_km[121]; -model_strat_heights_km[1]].*1000
        Height_sigma = fill(0.01, 3)
        Age_Sidedness = zeros(3)
        (youngest, oldest) = extrema(Age)
        aveuncert = nanmean(Age_sigma)
        (bottom, top) = extrema(Height)
        dt_dH = (oldest-youngest)/(top-bottom)
        model_heights = copy(-model_strat_heights_km[2:end]).*1000

        # Model configuration -- read from struct
        resolution = 1
        bounding = 0.1
        npoints_approx = round(Int,length(bottom:resolution:top) * (1 + 2*bounding))
        burnin = 3000*npoints_approx
        nsteps = 1500
        sieve = round(Int,npoints_approx)

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            offset = round((top-bottom)*bounding)
            Age = [oldest + offset*dt_dH; Age; youngest - offset*dt_dH]
            Age_sigma = [nanmean(Age_sigma)/10; Age_sigma; nanmean(Age_sigma)/10]
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
        end
        active_height_t = (model_heights .>= bottom) .& (model_heights .<= top)
        npoints = length(model_heights)

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        model_ages = a .+ b .* collect(model_heights)


        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_model_ages = model_ages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_model_ages[i] - Age[i])
                closest_model_ages[i] = Age[i]
            end
        end
        
        #New addition - weight of each term    
        w_age = 1
        w_height = 1
        w_parameters = 1
        w_subsidence = 1
        w_hiatus = 1

        ll = w_age*normpdf_ll(Age, Age_sigma, closest_model_ages)
        ll += w_height*normpdf_ll(Height, Height_sigma, sample_height)
        
        ideal_subs_parameters = [1.3, 820]
        ideal_subs_parameters_sigma = [0.2,20]
        subs_parameters = [ideal_subs_parameters[1]+0.01, ideal_subs_parameters[2]-1]
        ll += w_parameters*normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parameters)
        ll += w_subsidence*subsidence_ll(Sμ, Sσ, model_ages, subs_parameters)/(length(model_strat_heights_km)-1)

        # Ensure there is only one effective hiatus at most for each height node
        Hiatus_data_height = -0.01
        Hiatus_data_duration = 5
        Hiatus_data_duration_sigma = 0.5
        
        closest_hiatus = findclosestabove((Hiatus_data_height),model_heights)
        closest_hiatus_unique = unique(closest_hiatus)
        Hiatus_height = Array{Float64}(undef,size(closest_hiatus_unique))
        Hiatus_duration = Array{Float64}(undef,size(closest_hiatus_unique))
        Hiatus_duration_sigma = Array{Float64}(undef,size(closest_hiatus_unique))
        for i=1:length(closest_hiatus_unique)
            Hiatus_height[i] = mean((Hiatus_data_height)[closest_hiatus.==closest_hiatus_unique[i]])
            Hiatus_duration[i] = sum((Hiatus_data_duration)[closest_hiatus.==closest_hiatus_unique[i]])
            Hiatus_duration_sigma[i] = sqrt(sum((Hiatus_data_duration_sigma)[closest_hiatus.==closest_hiatus_unique[i]].^2))
        end

        # Add log likelihood for hiatus duration
        duration = @. min(model_ages[closest_hiatus_unique - 1] - model_ages[closest_hiatus_unique], Hiatus_duration)
        ll += w_hiatus*normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, duration)

        # Preallocate variables for MCMC proposals
        llₚ=ll
        chosen_point=0
        model_agesₚ = copy(model_ages)
        closestₚ = copy(closest)
        durationₚ = copy(duration)
        sample_heightₚ = copy(sample_height)
        closest_model_agesₚ = copy(closest_model_ages)
        subs_parametersₚ = copy(subs_parameters)

        # Run burnin
        # acceptancedist = fill(false,burnin)
        print("Burn-in: ", burnin, " steps\n")
        pgrs = Progress(burnin, desc="Burn-in...")
        pgrs_interval = ceil(Int,sqrt(burnin))
        for n=1:burnin
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)
            copyto!(subs_parametersₚ, subs_parameters)

            subs_parametersₚ += randn() *ideal_subs_parameters_sigma

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points that are still stratigraphically below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points that are still stratigraphically above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end

                # If chosen_point is a hiatus point, let there be a 20 percent chance of
                # adjusting the point below the hiatus as well
                if rand() < 0.2
                    # if Hiatus_height_uncert>0
                    #     closest_hiatus = findclosestabove(h.Height+randn(size(h.Height)).*Hiatus_height_uncert,heights)
                    # end
                    if any(closest_hiatus_unique.==chosen_point)
                        chosen_point -= 1
                        model_agesₚ[chosen_point] = model_ages[chosen_point] + r
                        #Resolve conflicts
                        if r > 0 # If proposing increased age
                            @inbounds for i=1:chosen_point # younger points below
                                if model_agesₚ[i] < model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        else # if proposing decreased age
                            @inbounds for i=chosen_point:npoints # older points above
                                if model_agesₚ[i] > model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        end
                    end
                end
            end


            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = w_age*normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += w_height*normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += w_subsidence*subsidence_ll(Sμ, Sσ, model_agesₚ, subs_parametersₚ)/(length(model_strat_heights_km)-1)
            llₚ += w_parameters*normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)
            
            # Add log likelihood for hiatus duration
            @. durationₚ = min(model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique], Hiatus_duration)
            llₚ += w_hiatus*normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(duration, durationₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
                # acceptancedist[i] = true
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs, burnin) # Finalize

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)
        beta_t0dist = Array{Float64}(undef,2,nsteps)
        hiatusdist = Array{Float64}(undef,length(duration),nsteps)

        # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        pgrs_interval = ceil(Int,sqrt(nsteps*sieve))
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)
            copyto!(subs_parametersₚ, subs_parameters)

            subs_parametersₚ += randn() *ideal_subs_parameters_sigma

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points that are still stratigraphically below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points that are still stratigraphically above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end

                # If chosen_point is a hiatus point, let there be a 20 percent chance of
                # adjusting the point below the hiatus as well
                if rand() < 0.2
                    # if Hiatus_height_uncert>0
                    #     closest_hiatus = findclosestabove(h.Height+randn(size(h.Height)).*Hiatus_height_uncert,heights)
                    # end
                    if any(closest_hiatus_unique.==chosen_point)
                        chosen_point -= 1
                        model_agesₚ[chosen_point] = model_ages[chosen_point] + r
                        #Resolve conflicts
                        if r > 0 # If proposing increased age
                            @inbounds for i=1:chosen_point # younger points below
                                if model_agesₚ[i] < model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        else # if proposing decreased age
                            @inbounds for i=chosen_point:npoints # older points above
                                if model_agesₚ[i] > model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = w_age*normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += w_height*normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += w_subsidence*subsidence_ll(Sμ, Sσ, model_agesₚ, subs_parametersₚ)/(length(model_strat_heights_km)-1)
            llₚ += w_parameters*normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)

            # Add log likelihood for hiatus duration
            @. durationₚ = min(model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique], Hiatus_duration)
            llₚ += w_hiatus*normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(duration, durationₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= model_ages
                beta_t0dist[:,n÷sieve] .= subs_parameters
                hiatusdist[n÷sieve] = duration[1]
                # !!! Test in the previous model if this works !!!
                # predicted_ages[n÷sieve] = beta_t0dist[2]+τ*log(1-(Sμ[?]*pi)/(E₀*beta_t0dist[1]*sin(pi/beta_t0dist[2])))
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs,nsteps*sieve)

        # Crop the result
        agedist = agedist[active_height_t,:]

        mdl_height = model_heights[active_height_t]
        mdl_age = nanmean(agedist,dim=2)
        mdl_median = nanmedian(agedist, dim=2)
        mdl_age_025CI = nanpctile(agedist,2.5, dim=2)
        mdl_age_975CI = nanpctile(agedist, 97.5, dim=2)
        beta = nanmean(beta_t0dist, dims=2)[1]
        t0 = nanmean(beta_t0dist, dims=2)[2]
        beta_025CI = nanpctile(beta_t0dist[1,:],2.5)
        beta_975CI = nanpctile(beta_t0dist[1,:],97.5)
        t0_025CI = nanpctile(beta_t0dist[2,:],2.5)
        t0_975CI = nanpctile(beta_t0dist[2,:],97.5)
        beta_σ = nanstd(beta_t0dist, dims=2)[1]
        t0_σ = nanstd(beta_t0dist, dims=2)[2]
        mdl_hiatus_duration = nanmean(hiatusdist)
        mdl_hiatus_σ = nanstd(hiatusdist)
        mdl_hiatus_025CI = nanpctile(hiatusdist,2.5)
        mdl_hiatus_975CI = nanpctile(hiatusdist,97.5)

        hdl_hiatus = plot([mdl_age_025CI; reverse(mdl_age_975CI)],[mdl_height; reverse(mdl_height)], fill=(round(Int,minimum(mdl_height)),0.5,:blue), label="Model results - 2σ error envelope")
        plot!(hdl_hiatus, mdl_age, mdl_height, linecolor=:blue, label="Model results - mean", fg_color_legend=:white) # Center line
        t = Age_Sidedness .== 0 # Two-sided constraints (plot in black)
        any(t) && plot!(hdl_hiatus, Age[t], Height[t], xerror=Age_sigma[t],label="Input - age data",seriestype=:scatter,color=:black)
        age_at_height = (mdl_height,mdl_age,-2232)

        curve_ages = similar(Sμ)
        curve_ages = τ.*log.(1 .-((Sμ.*pi)./(E₀*beta*sin(pi/beta)))).+t0
        model_strat_heights_m = copy(-model_strat_heights_km[2:end]).*1000
        
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
        
        figpath = "C:/Users/zhtia/Documents/GitHub/Chron.jl/trials/"
        savefig(hiatus_age_depth_curve, figpath*"Age_depth_curve_entire_section_hiatusT20.png")
        savefig(hiatus_thermal_subs_curve, figpath*"Thermal_subs_curve_entire_section_hiatusT20.png")
        savefig(hdl_hiatus, figpath*"Age_depth_curve_upper_section_hiatusT20.png")


        #= Will need to use part of this section later but not for now
        using DelimitedFiles
        writedlm("agedist.txt", agedist)
        writedlm("beta_t0dist.txt", beta_t0dist)
        writedlm("lldist.txt", lldist)

        agedist = readdlm("agedist.csv", ',', Float64)
        beta_t0dist = readdlm("beta_t0dist.csv", ',', Float64)
        lldist = readdlm("lldist.csv", ',', Float64)
        =#

        # Still need to optimize the following part (result summary and plotting）and move them here
        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian(agedist,dim=2), # Median age
            pctile(agedist,2.5,dim=2), # 2.5th percentile
            pctile(agedist,97.5,dim=2) # 97.5th percentile
        )

        return mdl, agedist, lldist, hiatusdist
    #end

    



## --- Stratigraphic MCMC model without hiatus, with distribution LL # # # # # #

    function StratMetropolisDist(smpl::ChronAgeData, config::StratAgeModelConfiguration)
        # Run stratigraphic MCMC model
        print("Generating stratigraphic age-depth model...\n")

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints
        Age = copy(smpl.Age)::Array{Float64,1}
        Age_sigma = copy(smpl.Age_sigma)::Array{Float64,1}
        Height = copy(smpl.Height)::Array{Float64,1}
        Height_sigma = smpl.Height_sigma::Array{Float64,1} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Array{Float64,1} # Bottom is a maximum age and top is a minimum age
        p = copy(smpl.Params)::Array{Float64,2}
        (youngest, oldest) = extrema(Age)
        aveuncert = nanmean(Age_sigma)
        (bottom, top) = extrema(Height)
        dt_dH = (oldest-youngest)/(top-bottom)
        model_heights = bottom:resolution:top


        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            offset = (top-bottom)*bounding
            Age = [oldest + offset*dt_dH; Age; youngest - offset*dt_dH]
            Age_sigma = [nanmean(Age_sigma)/10; Age_sigma; nanmean(Age_sigma)/10]
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
            pl = ones(5); pl[2] = oldest + offset*dt_dH; pl[3] = nanmean(Age_sigma)/10
            pu = ones(5); pu[2] = youngest - offset*dt_dH; pu[3] = nanmean(Age_sigma)/10
            p = hcat(pl,p,pu) # Add parameters for upper and lower runaway bounds
        end
        active_height_t = (model_heights .>= bottom) .& (model_heights .<= top)
        npoints = length(model_heights)

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        model_ages = a .+ b .* collect(model_heights)


        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_model_ages = model_ages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_model_ages[i] - Age[i])
                closest_model_ages[i] = Age[i]
            end
        end
        ll = bilinear_exponential_ll(closest_model_ages ,p)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

        # Preallocate variables for MCMC proposals
        llₚ = ll
        model_agesₚ = copy(model_ages)
        closestₚ = copy(closest)
        sample_heightₚ = copy(sample_height)
        closest_model_agesₚ = copy(closest_model_ages)

        # Run burnin
        # acceptancedist = fill(false,burnin)
        print("Burn-in: ", burnin, " steps\n")
        pgrs = Progress(burnin, desc="Burn-in...")
        pgrs_interval = ceil(Int,sqrt(burnin))
        for n=1:burnin
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point
                        # younger points below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints
                        # older points above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end
            end


            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = bilinear_exponential_ll(closest_model_agesₚ, p)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                # acceptancedist[i] = true
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs, burnin) # Finalize

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)


        # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        pgrs_interval = ceil(Int,sqrt(nsteps*sieve))
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point
                        # younger points below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints
                        # older points above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end
            end


            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = bilinear_exponential_ll(closest_model_agesₚ, p)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= model_ages
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs,nsteps*sieve)

        # Crop the result
        agedist = agedist[active_height_t,:]

        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian(agedist,dim=2), # Median age
            pctile(agedist,2.5,dim=2), # 2.5th percentile
            pctile(agedist,97.5,dim=2) # 97.5th percentile
        )

        return mdl, agedist, lldist
    end

## --- Stratigraphic MCMC model with hiata, with distribution LL # # # # # # # #

    function StratMetropolisDist(smpl::ChronAgeData, hiatus::HiatusData, config::StratAgeModelConfiguration)
        # Run stratigraphic MCMC model, with hiata
        print("Generating stratigraphic age-depth model...\n")

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints
        Age = copy(smpl.Age)::Array{Float64,1}
        Age_sigma = copy(smpl.Age_sigma)::Array{Float64,1}
        Height = copy(smpl.Height)::Array{Float64,1}
        Height_sigma = smpl.Height_sigma::Array{Float64,1} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Array{Float64,1} # Bottom is a maximum age and top is a minimum age
        p = copy(smpl.Params)::Array{Float64,2}
        (bottom, top) = extrema(Height)
        (youngest, oldest) = extrema(Age)
        dt_dH = (oldest-youngest)/(top-bottom)
        aveuncert = nanmean(Age_sigma)
        model_heights = bottom:resolution:top

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            offset = (top-bottom)*bounding
            Age = [oldest + offset*dt_dH; Age; youngest - offset*dt_dH]
            Age_sigma = [nanmean(Age_sigma)/10; Age_sigma; nanmean(Age_sigma)/10]
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
            pl = ones(5); pl[2] = oldest + offset*dt_dH; pl[3] = nanmean(Age_sigma)/10
            pu = ones(5); pu[2] = youngest - offset*dt_dH; pu[3] = nanmean(Age_sigma)/10
            p = hcat(pl,p,pu) # Add parameters for upper and lower runaway bounds
        end
        active_height_t = (model_heights .>= bottom) .& (model_heights .<= top)
        npoints = length(model_heights)

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        model_ages = a .+ b .* collect(model_heights)


        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_model_ages = model_ages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_model_ages[i] - Age[i])
                closest_model_ages[i] = Age[i]
            end
        end
        ll = bilinear_exponential_ll(closest_model_ages ,p)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

        # Ensure there is only one effective hiatus at most for each height node
        closest_hiatus = findclosestabove((hiatus.Height::Array{Float64,1}),model_heights)
        closest_hiatus_unique = unique(closest_hiatus)
        Hiatus_height = Array{Float64}(undef,size(closest_hiatus_unique))
        Hiatus_duration = Array{Float64}(undef,size(closest_hiatus_unique))
        Hiatus_duration_sigma = Array{Float64}(undef,size(closest_hiatus_unique))
        for i=1:length(closest_hiatus_unique)
            Hiatus_height[i] = mean((hiatus.Height::Array{Float64,1})[closest_hiatus.==closest_hiatus_unique[i]])
            Hiatus_duration[i] = sum((hiatus.Duration::Array{Float64,1})[closest_hiatus.==closest_hiatus_unique[i]])
            Hiatus_duration_sigma[i] = sqrt(sum((hiatus.Duration_sigma::Array{Float64,1})[closest_hiatus.==closest_hiatus_unique[i]].^2))
        end

        # Add log likelihood for hiatus duration
        duration = @. min(model_ages[closest_hiatus_unique - 1] - model_ages[closest_hiatus_unique], Hiatus_duration)
        ll += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, duration)

        # Preallocate variables for MCMC proposals
        llₚ = ll
        chosen_point = 0
        model_agesₚ = copy(model_ages)
        closestₚ = copy(closest)
        durationₚ = copy(duration)
        sample_heightₚ = copy(sample_height)
        closest_model_agesₚ = copy(closest_model_ages)

        # Run burnin
        # acceptancedist = fill(false,burnin)
        print("Burn-in: ", burnin, " steps\n")
        pgrs = Progress(burnin, desc="Burn-in...")
        pgrs_interval = ceil(Int,sqrt(burnin))
        for n=1:burnin
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point
                        # younger points below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints
                        # older points above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end

                # If chosen_point is a hiatus point, let there be a 20 percent chance of
                # adjusting the point below the hiatus as well
                if rand() < 0.2
                    # if Hiatus_height_uncert>0
                    #     closest_hiatus = findclosestabove(h.Height+randn(size(h.Height)).*Hiatus_height_uncert,heights)
                    # end
                    if any(closest_hiatus_unique.==chosen_point)
                        chosen_point -= 1
                        model_agesₚ[chosen_point] = model_ages[chosen_point] + r
                        #Resolve conflicts
                        if r > 0 # If proposing increased age
                            @inbounds for i=1:chosen_point # younger points below
                                if model_agesₚ[i] < model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        else # if proposing decreased age
                            @inbounds for i=chosen_point:npoints # older points above
                                if model_agesₚ[i] > model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        end
                    end
                end
            end


            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = bilinear_exponential_ll(closest_model_agesₚ, p)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Add log likelihood for hiatus duration
            @. durationₚ = min(model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique], Hiatus_duration)
            llₚ += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                # acceptancedist[i] = true
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs, burnin) # Finalize

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)
        hiatusdist = Array{Float64}(undef,length(duration),nsteps)

        # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        pgrs_interval = ceil(Int,sqrt(nsteps*sieve))
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end

                # If chosen_point is a hiatus point, let there be a 20 percent chance of
                # adjusting the point below the hiatus as well
                if rand() < 0.2
                    # if Hiatus_height_uncert>0
                    #     closest_hiatus = findclosestabove(h.Height+randn(size(h.Height)).*Hiatus_height_uncert,heights)
                    # end
                    if any(closest_hiatus_unique.==chosen_point)
                        chosen_point -= 1
                        model_agesₚ[chosen_point] = model_ages[chosen_point] + r
                        #Resolve conflicts
                        if r > 0 # If proposing increased age
                            @inbounds for i=1:chosen_point # younger points below
                                if model_agesₚ[i] < model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        else # if proposing decreased age
                            @inbounds for i=chosen_point:npoints # older points above
                                if model_agesₚ[i] > model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = bilinear_exponential_ll(closest_model_agesₚ, p)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Add log likelihood for hiatus duration
            @. durationₚ = min(model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique], Hiatus_duration)
            llₚ += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(duration, durationₚ)
                copyto!(sample_height, sample_heightₚ)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= model_ages
                hiatusdist[:,n÷sieve] .= duration
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs,nsteps*sieve)

        # Crop the result
        agedist = agedist[active_height_t,:]

        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian(agedist,dim=2), # Median age
            pctile(agedist,2.5,dim=2), # 2.5th percentile
            pctile(agedist,97.5,dim=2) # 97.5th percentile
        )

        return mdl, agedist, lldist, hiatusdist
    end

## --- Stratigraphic MCMC model without hiatus, for radiocarbon ages # # # # # #

    function StratMetropolis14C(smpl::ChronAgeData, config::StratAgeModelConfiguration)
        # Run stratigraphic MCMC model
        print("Generating stratigraphic age-depth model...\n")

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints
        Age = copy(smpl.Age)::Array{Float64,1}
        Age_sigma = copy(smpl.Age_sigma)::Array{Float64,1}
        Height = copy(smpl.Height)::Array{Float64,1}
        Height_sigma = smpl.Height_sigma::Array{Float64,1} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Array{Float64,1} # Bottom is a maximum age and top is a minimum age
        p = copy(smpl.Params)::Array{Float64,2}
        (bottom, top) = extrema(Height)
        (youngest, oldest) = extrema(Age)
        dt_dH = (oldest-youngest)/(top-bottom)
        aveuncert = nanmean(Age_sigma)
        model_heights = bottom:resolution:top

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            offset = (top-bottom)*bounding
            Age = [oldest + offset*dt_dH; Age; youngest - offset*dt_dH]
            Age_sigma = [nanmean(Age_sigma)/10; Age_sigma; nanmean(Age_sigma)/10]
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
            boundsigma = nanmean(Age_sigma)/10
            pl = normpdf_ll.(oldest + offset*dt_dH, boundsigma, 1:50000)
            pu = normpdf_ll.(youngest - offset*dt_dH, boundsigma, 1:50000)
            p = hcat(pl,p,pu) # Add parameters for upper and lower runaway bounds
        end
        active_height_t = (model_heights .>= bottom) .& (model_heights .<= top)
        npoints = length(model_heights)

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        model_ages = a .+ b .* collect(model_heights)


        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_model_ages = model_ages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_model_ages[i] - Age[i])
                closest_model_ages[i] = Age[i]
            end
        end
        ll = interpolate_ll(closest_model_ages, p)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

        # Preallocate variables for MCMC proposals
        llₚ = ll
        model_agesₚ = copy(model_ages)
        closestₚ = copy(closest)
        sample_heightₚ = copy(sample_height)
        closest_model_agesₚ = copy(closest_model_ages)

        # Run burnin
        # acceptancedist = fill(false,burnin)
        print("Burn-in: ", burnin, " steps\n")
        pgrs = Progress(burnin, desc="Burn-in...")
        pgrs_interval = ceil(Int,sqrt(burnin))
        for n=1:burnin
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end
            end


            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = interpolate_ll(closest_model_agesₚ, p)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                # acceptancedist[i] = true
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs, burnin) # Finalize

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)


        # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        pgrs_interval = ceil(Int,sqrt(nsteps*sieve))
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = interpolate_ll(closest_model_agesₚ, p)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= model_ages
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs,nsteps*sieve)

        # Crop the result
        agedist = agedist[active_height_t,:]

        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian(agedist,dim=2), # Median age
            pctile(agedist,2.5,dim=2), # 2.5th percentile
            pctile(agedist,97.5,dim=2) # 97.5th percentile
        )

        return mdl, agedist, lldist
    end

## --- Stratigraphic MCMC model with hiatus, for radiocarbon ages # # # # # #

    function StratMetropolis14C(smpl::ChronAgeData, hiatus::HiatusData, config::StratAgeModelConfiguration)
        # Run stratigraphic MCMC model, with hiata
        print("Generating stratigraphic age-depth model...\n")

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints
        Age = copy(smpl.Age)::Array{Float64,1}
        Age_sigma = copy(smpl.Age_sigma)::Array{Float64,1}
        Height = copy(smpl.Height)::Array{Float64,1}
        Height_sigma = smpl.Height_sigma::Array{Float64,1} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Array{Float64,1} # Bottom is a maximum age and top is a minimum age
        p = copy(smpl.Params)::Array{Float64,2}
        (bottom, top) = extrema(Height)
        (youngest, oldest) = extrema(Age)
        dt_dH = (oldest-youngest)/(top-bottom)
        aveuncert = nanmean(Age_sigma)
        model_heights = bottom:resolution:top

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            offset = (top-bottom)*bounding
            Age = [oldest + offset*dt_dH; Age; youngest - offset*dt_dH]
            Age_sigma = [nanmean(Age_sigma)/10; Age_sigma; nanmean(Age_sigma)/10]
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
            boundsigma = nanmean(Age_sigma)/10
            pl = normpdf_ll.(oldest + offset*dt_dH, boundsigma, 1:50000)
            pu = normpdf_ll.(youngest - offset*dt_dH, boundsigma, 1:50000)
            p = hcat(pl,p,pu) # Add parameters for upper and lower runaway bounds
        end
        active_height_t = (model_heights .>= bottom) .& (model_heights .<= top)
        npoints = length(model_heights)

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        model_ages = a .+ b .* collect(model_heights)


        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_model_ages = model_ages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_model_ages[i] - Age[i])
                closest_model_ages[i] = Age[i]
            end
        end
        ll = interpolate_ll(closest_model_ages, p)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

        # Ensure there is only one effective hiatus at most for each height node
        closest_hiatus = findclosestabove((hiatus.Height::Array{Float64,1}),model_heights)
        closest_hiatus_unique = unique(closest_hiatus)
        Hiatus_height = Array{Float64}(undef,size(closest_hiatus_unique))
        Hiatus_duration = Array{Float64}(undef,size(closest_hiatus_unique))
        Hiatus_duration_sigma = Array{Float64}(undef,size(closest_hiatus_unique))
        for i=1:length(closest_hiatus_unique)
            Hiatus_height[i] = mean((hiatus.Height::Array{Float64,1})[closest_hiatus.==closest_hiatus_unique[i]])
            Hiatus_duration[i] = sum((hiatus.Duration::Array{Float64,1})[closest_hiatus.==closest_hiatus_unique[i]])
            Hiatus_duration_sigma[i] = sqrt(sum((hiatus.Duration_sigma::Array{Float64,1})[closest_hiatus.==closest_hiatus_unique[i]].^2))
        end

        # Add log likelihood for hiatus duration
        duration = @. min(model_ages[closest_hiatus_unique - 1] - model_ages[closest_hiatus_unique], Hiatus_duration)
        ll += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, duration)

        # Preallocate variables for MCMC proposals
        llₚ = ll
        chosen_point = 0
        model_agesₚ = copy(model_ages)
        closestₚ = copy(closest)
        durationₚ = copy(duration)
        sample_heightₚ = copy(sample_height)
        closest_model_agesₚ = copy(closest_model_ages)

        # Run burnin
        # acceptancedist = fill(false,burnin)
        print("Burn-in: ", burnin, " steps\n")
        pgrs = Progress(burnin, desc="Burn-in...")
        pgrs_interval = ceil(Int,sqrt(burnin))
        for n=1:burnin
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end

                # If chosen_point is a hiatus point, let there be a 20 percent chance of
                # adjusting the point below the hiatus as well
                if rand() < 0.2
                    # if Hiatus_height_uncert>0
                    #     closest_hiatus = findclosestabove(h.Height+randn(size(h.Height)).*Hiatus_height_uncert,heights)
                    # end
                    if any(closest_hiatus_unique.==chosen_point)
                        chosen_point -= 1
                        model_agesₚ[chosen_point] = model_ages[chosen_point] + r
                        #Resolve conflicts
                        if r > 0 # If proposing increased age
                            @inbounds for i=1:chosen_point # younger points below
                                if model_agesₚ[i] < model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        else # if proposing decreased age
                            @inbounds for i=chosen_point:npoints # older points above
                                if model_agesₚ[i] > model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        end
                    end
                end
            end


            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = interpolate_ll(closest_model_agesₚ, p)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Add log likelihood for hiatus duration
            @. durationₚ = min(model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique], Hiatus_duration)
            llₚ += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                # acceptancedist[i] = true
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs, burnin) # Finalize

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)
        hiatusdist = Array{Float64}(undef,length(duration),nsteps)

        # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        pgrs_interval = ceil(Int,sqrt(nsteps*sieve))
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_heightₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - model_heights[1])/resolution)+1
                    if closestₚ[i] < 1 # Check we're still within bounds
                        closestₚ[i] = 1
                    elseif closestₚ[i] > npoints
                        closestₚ[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                model_agesₚ[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if model_agesₚ[i] < model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if model_agesₚ[i] > model_agesₚ[chosen_point]
                            model_agesₚ[i] = model_agesₚ[chosen_point]
                        end
                    end
                end

                # If chosen_point is a hiatus point, let there be a 20 percent chance of
                # adjusting the point below the hiatus as well
                if rand() < 0.2
                    # if Hiatus_height_uncert>0
                    #     closest_hiatus = findclosestabove(h.Height+randn(size(h.Height)).*Hiatus_height_uncert,heights)
                    # end
                    if any(closest_hiatus_unique.==chosen_point)
                        chosen_point -= 1
                        model_agesₚ[chosen_point] = model_ages[chosen_point] + r
                        #Resolve conflicts
                        if r > 0 # If proposing increased age
                            @inbounds for i=1:chosen_point # younger points below
                                if model_agesₚ[i] < model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        else # if proposing decreased age
                            @inbounds for i=chosen_point:npoints # older points above
                                if model_agesₚ[i] > model_agesₚ[chosen_point]
                                    model_agesₚ[i] = model_agesₚ[chosen_point]
                                end
                            end
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = interpolate_ll(closest_model_agesₚ, p)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Add log likelihood for hiatus duration
            @. durationₚ = min(model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique], Hiatus_duration)
            llₚ += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(duration, durationₚ)
                copyto!(sample_height, sample_heightₚ)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= model_ages
                hiatusdist[:,n÷sieve] .= duration
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs,nsteps*sieve)

        # Crop the result
        agedist = agedist[active_height_t,:]

        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian(agedist,dim=2), # Median age
            pctile(agedist,2.5,dim=2), # 2.5th percentile
            pctile(agedist,97.5,dim=2) # 97.5th percentile
        )

        return mdl, agedist, lldist, hiatusdist
    end

## ---
