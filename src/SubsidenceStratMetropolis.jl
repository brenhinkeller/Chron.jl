## --- Stratigraphic MCMC model without hiatus # # # # # # # # # # # # # # # # #


# Part 1: Decompaction and Backstripping 

function DecompactBackstrip(strat::StratData, nsims, res)

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
        lithology_inputs = strat.Lithology
        height_inputs = cumsum([0; strat.Thickness])
        nlayer_input = length(strat.Thickness)
        model_strat_heights = [0:res:maximum(height_inputs);]
        model_nlayer = length(model_strat_heights)-1

    # Allocate parameters as distributions; each element/distribution represents a layer
        # porosity depth coefficient(c)
        c_dist = Array{Distribution,1}(undef, nlayer_input)
        # surface porosity (ϕ₀)
        ϕ₀_dist = Array{Distribution,1}(undef, nlayer_input)
        # sediment grain density (ρg)
        ρg = Array{Float64,1}(undef, nlayer_input)

    # Find the correct c, ϕ₀, and ρg for each layer based on lithology
    for i = 1:nlayer_input
        if lithology_inputs[i] == "Shale"
            # c has a lower bound of 0 b/c porosity at depth should not be greater than porosity at surface
            c_dist[i] = truncated(Normal(0.51, 0.15), 0, Inf)
            # ϕ₀ has a lower bound of 0 and an upper bound of 1 b/c of the definition of porosity
            ϕ₀_dist[i] = truncated(Normal(0.63, 0.15), 0, 1)
            ρg[i] = 2720
        elseif lithology_inputs[i] == "Sandstone"
            c_dist[i] = truncated(Normal(0.27, 0.1), 0, Inf)
            ϕ₀_dist[i] = truncated(Normal(0.49, 0.1), 0, 1)
            ρg[i] = 2650
        elseif lithology_inputs[i] == "Chalk"
            c_dist[i] = truncated(Normal(0.71, 0.15), 0, Inf)
            ϕ₀_dist[i] = truncated(Normal(0.7, 0.15), 0, 1)
            ρg[i] = 2710
        elseif lithology_inputs[i] == "Limestone"
            c_dist[i] = truncated(Normal(0.6, 0.2), 0, Inf)
            ϕ₀_dist[i] = truncated(Normal(0.4, 0.17), 0, 1)
            ρg[i] = 2710
        elseif lithology_inputs[i] == "Dolostone"
            c_dist[i] = truncated(Normal(0.6, 0.2), 0, Inf)
            ϕ₀_dist[i] = truncated(Normal(0.2, 0.1), 0, 1)
            ρg[i] = 2870
        elseif lithology_inputs[i] == "Anhydrite"
            c_dist[i] = truncated(Normal(0.2, 0.1), 0, Inf)
            ϕ₀_dist[i] = truncated(Normal(0.05, 0.05), 0, 1)
            ρg[i] = 2960
        elseif lithology_inputs[i] == "Quartzite"
            c_dist[i] = truncated(Normal(0.3, 0.1), 0, Inf)
            ϕ₀_dist[i] = truncated(Normal(0.2, 0.1), 0, 1)
            ρg[i] = 2650
        end
    end

    # Prep for decompaction and backstripping MC
        # Define parameters for decompaction (number of simulations; water and mantle densities)
        ρw = 1000
        ρm = 3330
           
        # Allocate depth matricies (rows are strat horizons, columns are timesteps)
        # decompacted depth
        Y = fill(1E-18, (model_nlayer+1, model_nlayer+1))
        # decompacted depth after applying varies corrections
        Y_corr = fill(1E-18, (model_nlayer+1, model_nlayer+1))
        #Y_max = fill(1E-18, (nsims, model_nlayer))

        c_highres = Array{Float64,1}(undef, model_nlayer)
        ϕ₀_highres = Array{Float64,1}(undef, model_nlayer)
        ρg_highres = Array{Float64,1}(undef, model_nlayer)
        
        #Allocate porosity, density and tectonic subsidence matricies   
        # porosity of a strat unit at any depth
        ϕ_avg = fill(1E-18, (model_nlayer, model_nlayer+1))
        # bulk density of a single layer
        ρ_bulk = fill(1E-18, (model_nlayer, model_nlayer+1))
        # intermediate step - bulk density*thickness of a single layer
        m_bulk = fill(1E-18, (model_nlayer, model_nlayer+1))
        # bulk density of the entire column
        ρ_bulk_column = Array{Float64,1}(undef, model_nlayer+1)
        # tectonic subsidence # this is the only one that will need to propagate outside of the loop
        Sₜ_km = Array{Float64,2}(undef, model_nlayer+1, nsims)
        
        paleo_wd_dist = Array{Float64,1}(undef, nsims)
    
    # MC for decompaction and backstripping
        print("Decompaction and Backstripping: ", nsims, " steps\n")
        pgrs = Progress(nsims, desc="Decompaction and Backstripping...")
        pgrs_interval = ceil(Int,10)    
        for sim = 1:nsims
            # randomly select c and ϕ₀ vectors (for the whole column) from the distributions
            c = rand.(c_dist)
            ϕ₀ = rand.(ϕ₀_dist)

            for i = 1:nlayer_input
                for j = 1:model_nlayer
                    if model_strat_heights[j+1]>height_inputs[i]
                        c_highres[j]=c[i]
                        ϕ₀_highres[j]=ϕ₀[i]
                        ρg_highres[j]=ρg[i]
                    end
                end
            end

            # Fill the first column with modern observed values (present-day depths)
            Y[:,1] .= model_strat_heights
            # Fill the first row with zeros
            Y[1,:] .= 0

            # Decompact
            # i = time steps during decompaction, which runs from 2 to layer_count b/c column 1 is present day
            # j = layer number, which runs from i to layer_count b/c the ith column begins with y₁' of layer i
            for i = 2:model_nlayer
                for j = i:model_nlayer
                    decompact!(view(Y,:,i), model_strat_heights, ϕ₀_highres[j], c_highres[j], j, i)
                end
            end

            # Apply decompaction corrections
            #=
            # vector for paleo sea-level corrections (each element represents a time step)
            paleo_sl = fill(0,nlayer_input+1) #using zeros for now but should turn this into a distribution
            # vector for paleo water-depth corrections (each element represents a time step)
            paleo_wd = fill(200*rand(),nlayer_input+1)  #using zeros for now but should correlate with lithology?
            =#
            #rand_wd = 0.2*rand()
            #paleo_wd_dist[sim]=rand_wd
            #paleo_sl_highres = fill(0, model_nlayer+1)
            paleo_wd_highres = fill(0, model_nlayer+1)
            #=
            paleo_sl_highres = Array{Float64,1}(undef, model_nlayer+1)
            paleo_wd_highres = Array{Float64,1}(undef, model_nlayer+1)
            paleo_sl_highres[1] = paleo_sl[1]
            paleo_wd_highres[1] = paleo_wd[1]


            for i = 1:nlayer_input+1
                for j = 1:model_nlayer+1
                    if model_strat_heights[j]>height_inputs[i]
                        paleo_sl_highres[j]=paleo_sl[i+1]
                        paleo_wd_highres[j]=paleo_wd[i+1]
                    end
                end
            end
            =#

            #Y_corr[end] = 0-paleo_sl_highres[end]*(ρw/(ρm-ρw))+paleo_wd_highres[end]-paleo_sl_highres[end]
            #Sₜ_km[end,:] .= copy(Y_corr[end])*(ρm/(ρm-ρw))
            #Y_max[sim, :] = maximum(Y_corr, dims = 1)
            
            # Remove the effect of sediment load - same indexing logic as the decompaction loop
            # i = time steps (columns) = 1:layer_count b/c need to calculate these parameters for the present day column/layers
            # j = layer number = i:layer_count b/c the ith column begins with y₁' of layer i
            for i = 1:model_nlayer+1
                for j = i:model_nlayer
                    ϕ_avg[j-i+1,i] = (ϕ₀_highres[j]/c_highres[j])*(exp(-c_highres[j]*Y[j-i+1,i])-exp(-c_highres[j]*Y[j-i+2,i]))/(Y[j-i+2,i]-Y[j-i+1,i])
                    ρ_bulk[j-i+1,i] = ρw*ϕ_avg[j-i+1,i]+(1-ϕ_avg[j-i+1,i])*ρg_highres[j]
                    m_bulk[j-i+1,i] = (Y[j-i+2,i]-Y[j-i+1,i])*ρ_bulk[j-i+1,i]
                end
            end
            
            m_bulk[:,end] .= 0
            for i = 1:model_nlayer+1        
                # maximum(Y_corr[:,i]) = total depth of column at time step i
                ρ_bulk_column[i] = sum(m_bulk[:,i])/maximum(Y[:,i])
                #Y_corr[:,i] = copy(Y[:,i]).-paleo_sl_highres[i]*(ρw/(ρm-ρw)).+paleo_wd_highres[i].-paleo_sl_highres[i]
                #Sₜ_km[i,sim] = copy(Y_corr[model_nlayer+2-i,i])*((ρm-ρ_bulk_column[i])/(ρm-ρw))
                Sₜ_km[i,sim] = copy(Y[model_nlayer+2-i,i])*((ρm-ρ_bulk_column[i])/(ρm-ρw)).+paleo_wd_highres[i]       
            end

            mod(sim,pgrs_interval)==0 && update!(pgrs, sim)
        end
        update!(pgrs,nsims)

        # Calculate summary statistics (mean and standard deviation)
        Sₜ = copy(Sₜ_km).*1000
        Sμ = dropdims(nanmean(Sₜ, dims=2), dims=2)
        Sσ = dropdims(nanstd(Sₜ, dims=2), dims=2)

    return Sₜ, Sμ, Sσ, model_strat_heights, paleo_wd_dist
end


# Part 2a: Modified StratMetropolis for extensional basins - without hiatus

function SubsidenceStratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, model_strat_heights, Sμ, Sσ, beta_ip, t0_ip)

    # Run stratigraphic MCMC model
    print("Generating stratigraphic age-depth model...\n")

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
        model_heights = copy(-model_strat_heights[2:end]).*1000

    # Model configuration -- read from struct
        resolution = config.resolution
        bounding = config.bounding
        nsteps = config.nsteps
        burnin = config.burnin
        sieve = config.sieve
        
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
        
        if any(i -> i == -2138, Height) == true #Will replace this 2138 with a more generic function to identify whether we have a age constraint at the base of section
            active_height_t = (model_heights .> bottom) .& (model_heights .<= top)
        elseif any(i -> i == -2138, Height) == false
            active_height_t = (model_heights .>= bottom) .& (model_heights .<= top)
        end
        npoints = length(model_heights)

    # STEP 1: calculate log likelihood of the modeled ages (and heights) in the initial proposal
        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        b = -0.04
        model_ages = a .+ b .* collect(model_heights)

        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_model_ages = model_ages[closest]

        #Age_twosided = Array{Float64}(undef, count(==(0), Age))
        #Age_sigma_twosided = Array{Float64}(undef, count(==(0), Age))

        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_model_ages[i] - Age[i])
                closest_model_ages[i] = Age[i]
            end
            #if Age_sidedness[i] == 0
            #    for j = 1:length(Age_twosided)
            #        Age_twosided[j] = Age[i]
            #        Age_sigma_twosided[j] = Age_sigma[i]
            #    end
            #end
        end

        ll = normpdf_ll(Age, Age_sigma, closest_model_ages)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

    # STEP 2: calculate log likelihood of the subsidence model parameters in the initial proposal
        # Define thermal subsidence model parameters 
        ideal_subs_parameters = therm.Param
        ideal_subs_parameters_sigma = therm.Sigma
                
        # Initial proposal for subsidence parameters - randomly pick a set of values from the distribution
        subs_parameters = [ideal_subs_parameters[1]+beta_ip, ideal_subs_parameters[2]+t0_ip]
        # Calculate log likelihood of this initial proposal for subsidence parameters
        ll += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parameters)

    # STEP 3: calculate log likelihood for the fit of the thermal subsidence curve in the initial proposal
        ts_model_ages = reverse(model_ages[active_height_t])
        model_heights_all = copy(-model_strat_heights).*1000
        if any(i -> i == -2138, Height) == true #Will replace this 2138 with a more generic function to identify whether we have a age constraint at the base of section
            ts_Sμ = Sμ[(model_heights_all .> bottom) .& (model_heights_all .<= top)]
            ts_Sσ = Sσ[(model_heights_all .> bottom) .& (model_heights_all .<= top)]
        elseif any(i -> i == -2138, Height) == false
            ts_Sμ = Sμ[(model_heights_all .>= bottom) .& (model_heights_all .<= top)]
            ts_Sσ = Sσ[(model_heights_all .>= bottom) .& (model_heights_all .<= top)]
        end

        ll += subsidence_ll(ts_Sμ, ts_Sσ, ts_model_ages, subs_parameters)

    # Preallocate variables for MCMC proposals
        llₚ = ll
        model_agesₚ = copy(model_ages)
        closestₚ = copy(closest)
        sample_heightₚ = copy(sample_height)
        closest_model_agesₚ = copy(closest_model_ages)
        subs_parametersₚ = copy(subs_parameters)
        
    # Run burnin
        # acceptancedist = fill(false,burnin)
        lldist_burnin = Array{Float64}(undef,burnin÷1000)

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
            subs_parametersₚ .+= randn.() .* ideal_subs_parameters_sigma

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
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)
            
            ts_model_agesₚ = reverse(model_agesₚ[active_height_t])
            llₚ += subsidence_ll(ts_Sμ, ts_Sσ, ts_model_agesₚ, subs_parametersₚ)
            
            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
                # acceptancedist[i] = true
            end
            if mod(n,1000) == 0
                lldist_burnin[n÷1000] = ll
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
            subs_parametersₚ .+= randn.() .*ideal_subs_parameters_sigma

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
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)

            ts_model_agesₚ = reverse(model_agesₚ[active_height_t])
            llₚ += subsidence_ll(ts_Sμ, ts_Sσ, ts_model_agesₚ, subs_parametersₚ)

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

        subsmdl = SubsidenceStratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian(agedist,dim=2), # Median age
            nanpctile(agedist,2.5,dim=2), # 2.5th percentile
            nanpctile(agedist,97.5,dim=2), # 97.5th percentile
            nanmean(beta_t0dist[1,:],dim=1), # Mean beta
            nanstd(beta_t0dist[1,:],dim=1), # Standard deviation
            nanmedian(beta_t0dist[1,:],dim=1), # Median beta
            nanpctile(beta_t0dist[1,:],2.5,dim=1), # 2.5th percentile
            nanpctile(beta_t0dist[1,:],97.5,dim=1), # 97.5th percentile
            nanmean(beta_t0dist[2,:],dim=1), # Mean T0
            nanstd(beta_t0dist[2,:],dim=1), # Standard deviation
            nanmedian(beta_t0dist[2,:],dim=1), # Median T0
            nanpctile(beta_t0dist[2,:],2.5,dim=1), # 2.5th percentile
            nanpctile(beta_t0dist[2,:],97.5,dim=1) # 97.5th percentile
        )

    return subsmdl, agedist, lldist, beta_t0dist, lldist_burnin
end



# Part 2b: Modified StratMetropolis for extensional basins - with hiatus (when we only know the strat height and do not have any previous knowledge about the duration)

function SubsidenceStratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, model_strat_heights, Sμ, Sσ, hiatus_height, beta_ip, t0_ip)

    # Run stratigraphic MCMC model
    print("Generating stratigraphic age-depth model...\n")

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
        model_heights = copy(-model_strat_heights[2:end]).*1000

    # Model configuration -- read from struct
        resolution = config.resolution
        bounding = config.bounding
        nsteps = config.nsteps
        burnin = config.burnin
        sieve = config.sieve
        
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
        closest = findclosest(sample_height, model_heights)
        closest_model_ages = model_ages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_model_ages[i] - Age[i])
                closest_model_ages[i] = Age[i]
            end
        end

        ll = normpdf_ll(Age, Age_sigma, closest_model_ages)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

    # STEP 2: calculate log likelihood of the subsidence model parameters in the initial proposal
        # Define thermal subsidence model parameters 
        ideal_subs_parameters = therm.Param
        ideal_subs_parameters_sigma = therm.Sigma
                
        # Initial proposal for subsidence parameters - randomly pick a set of values from the distribution
        subs_parameters = [ideal_subs_parameters[1]+beta_ip, ideal_subs_parameters[2]+t0_ip]
        # Calculate log likelihood of this initial proposal for subsidence parameters
        ll += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parameters)

    # STEP 3: calculate log likelihood for the fit of the thermal subsidence curve in the initial proposal
        active_height_t_hiatus = bottom .<= model_heights .< hiatus_height
        ts_model_ages_hiatus = reverse(model_ages[active_height_t_hiatus])
        model_heights_all = copy(-model_strat_heights).*1000
        ts_Sμ_hiatus = Sμ[bottom .<= model_heights_all .< hiatus_height]
        ts_Sσ_hiatus = Sσ[bottom .<= model_heights_all .< hiatus_height]
        ll += subsidence_ll(ts_Sμ_hiatus, ts_Sσ_hiatus, ts_model_ages_hiatus, subs_parameters)

    # Preallocate variables for MCMC proposals
        llₚ = ll
        model_agesₚ = copy(model_ages)
        closestₚ = copy(closest)
        sample_heightₚ = copy(sample_height)
        closest_model_agesₚ = copy(closest_model_ages)
        subs_parametersₚ = copy(subs_parameters)
        
    # Run burnin
        # acceptancedist = fill(false,burnin)
        lldist_burnin = Array{Float64}(undef,burnin÷10)
        
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
            subs_parametersₚ .+= randn.() .*ideal_subs_parameters_sigma

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
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)

            ts_model_ages_hiatusₚ = reverse(model_agesₚ[active_height_t_hiatus])
            llₚ += subsidence_ll(ts_Sμ_hiatus, ts_Sσ_hiatus, ts_model_ages_hiatusₚ, subs_parametersₚ)
            
            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                copyto!(subs_parameters, subs_parametersₚ)
                # acceptancedist[i] = true
            end
            if mod(n,100) == 0
                lldist_burnin[[n÷100]] = ll
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
        h = plot(xlabel="Age", ylabel="Height", framestyle=:box)
    
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        pgrs_interval = ceil(Int,sqrt(nsteps*sieve))
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)
            copyto!(subs_parametersₚ, subs_parameters)

            # Propose adjustment to subsidence_parametersₚ
            subs_parametersₚ .+= randn.() .*ideal_subs_parameters_sigma

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
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)
            
            ts_model_ages_hiatusₚ = reverse(model_agesₚ[active_height_t_hiatus])
            llₚ += subsidence_ll(ts_Sμ_hiatus, ts_Sσ_hiatus, ts_model_ages_hiatusₚ, subs_parametersₚ)

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
                if mod(n, sieve*4000) ==0
                    ts_model_heights = Sμ[6:86]
                    plot!(h, model_ages[active_height_t], reverse(ts_model_heights), label="", alpha=0.5, yflip = true, xflip = true)
                    yₛ = (3165.647578*subs_parameters[1]/pi)*sin(pi/subs_parameters[1]).*(1 .-exp.(-(subs_parameters[2] .-model_ages[active_height_t_hiatus])./50))
                    # y should be a smooth curve determined by the subsidence model with the current t0, beta, etc.
                    plot!(h, model_ages[active_height_t_hiatus], yₛ, label="", alpha=0.5)
                    scatter!(h, smpl.Age, [223.6125, 580.5430, 840.8018, 975.1254, 1030.9710, 1072.7511], label="data",seriestype=:scatter,color=:black)
                end
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs,nsteps*sieve)

        savefig(h, "SubsidenceModelAgeComparison_hiatus.pdf")

    # Crop the result
        agedist = agedist[active_height_t,:]

        subsmdl = SubsidenceStratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian(agedist,dim=2), # Median age
            nanpctile(agedist,2.5,dim=2), # 2.5th percentile
            nanpctile(agedist,97.5,dim=2), # 97.5th percentile
            nanmean(beta_t0dist[1,:],dim=1), # Mean beta
            nanstd(beta_t0dist[1,:],dim=1), # Standard deviation
            nanmedian(beta_t0dist[1,:],dim=1), # Median beta
            nanpctile(beta_t0dist[1,:],2.5,dim=1), # 2.5th percentile
            nanpctile(beta_t0dist[1,:],97.5,dim=1), # 97.5th percentile
            nanmean(beta_t0dist[2,:],dim=1), # Mean T0
            nanstd(beta_t0dist[2,:],dim=1), # Standard deviation
            nanmedian(beta_t0dist[2,:],dim=1), # Median T0
            nanpctile(beta_t0dist[2,:],2.5,dim=1), # 2.5th percentile
            nanpctile(beta_t0dist[2,:],97.5,dim=1) # 97.5th percentile
        )

    return subsmdl, agedist, lldist, beta_t0dist, lldist_burnin
end



# Part 2c: Modified StratMetropolis for extensional basins - with hiata (when we know the duration of hiata relatively precisely)
    
function SubsidenceStratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration, therm::ThermalSubsidenceParameters, model_strat_heights, Sμ, Sσ, hiatus::HiatusData)
        
    # Run stratigraphic MCMC model, with hiata
    print("Generating stratigraphic age-depth model...\n")

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
        Age = copy(smpl.Age)::Array{Float64,1}
        Age_sigma = copy(smpl.Age_sigma)::Array{Float64,1}
        Height = copy(smpl.Height)::Array{Float64,1}
        Height_sigma = smpl.Height_sigma::Array{Float64,1} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Array{Float64,1} # Bottom is a maximum age and top is a minimum age
        (bottom, top) = extrema(Height)
        (youngest, oldest) = extrema(Age)
        dt_dH = (oldest-youngest)/(top-bottom)
        aveuncert = nanmean(Age_sigma)
        model_heights = copy(-model_strat_heights[2:end]).*1000

    # Model configuration -- read from struct
        resolution = config.resolution
        bounding = config.bounding
        nsteps = config.nsteps
        burnin = config.burnin
        sieve = config.sieve

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
        closest = findclosest(sample_height, model_heights)
        closest_model_ages = model_ages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_model_ages[i] - Age[i])
                closest_model_ages[i] = Age[i]
            end
        end
        ll = normpdf_ll(Age, Age_sigma, closest_model_ages)
        ll += normpdf_ll(Height, Height_sigma, sample_height)
        
    # STEP 2: calculate log likelihood of the subsidence model parameters in the initial proposal
        # Define thermal subsidence model parameters - read from struct
        ideal_subs_parameters = therm.Param
        ideal_subs_parameters_sigma = therm.Sigma        
        # Initial proposal for subsidence parameters - randomly pick a set of values from the distribution
        subs_parameters = [ideal_subs_parameters[1]+0.01, ideal_subs_parameters[2]-1]
        # Calculate log likelihood of this initial proposal
        ll += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parameters)

    # STEP 3: calculate log likelihood for the fit of the thermal subsidence curve in the initial proposal
        ll += subsidence_ll(Sμ, Sσ, model_ages, subs_parameters)/(length(model_strat_heights)-1)

        # Ensure there is only one effective hiatus at most for each height node
        closest_hiatus = findclosestabove((hiatus.Height::Array{Float64,1}),model_heights)
        closest_hiatus_unique = unique(closest_hiatus)
        Hiatus_height = Array{Float64}(undef,size(closest_hiatus_unique))
        Hiatus_duration = Array{Float64}(undef,size(closest_hiatus_unique))
        Hiatus_duration_sigma = Array{Float64}(undef,size(closest_hiatus_unique))
        for i=1:length(closest_hiatus_unique)
            t = closest_hiatus.==closest_hiatus_unique[i]
            Hiatus_height[i] = mean((hiatus.Height::Array{Float64,1})[t])
            Hiatus_duration[i] = sum((hiatus.Duration::Array{Float64,1})[t])
            Hiatus_duration_sigma[i] = sqrt(sum((hiatus.Duration_sigma::Array{Float64,1})[t].^2))
        end

    # STEP 4: add log likelihood for hiatus duration
        duration = @. min(model_ages[closest_hiatus_unique - 1] - model_ages[closest_hiatus_unique], Hiatus_duration)
        ll += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, duration)

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

            # Propose adjustment to subsidence_parametersₚ
            subs_parametersₚ .+= randn.() .*ideal_subs_parameters_sigma

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
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += subsidence_ll(Sμ, Sσ, model_agesₚ, subs_parametersₚ)/(length(model_strat_heights)-1)
            llₚ += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)
            
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

            # Propose adjustment to subsidence_parametersₚ
            subs_parametersₚ .+= randn.() .*ideal_subs_parameters_sigma

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
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)
            llₚ += subsidence_ll(Sμ, Sσ, model_agesₚ, subs_parametersₚ)/(length(model_strat_heights)-1)
            llₚ += normpdf_ll(ideal_subs_parameters, ideal_subs_parameters_sigma, subs_parametersₚ)

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
                copyto!(subs_parameters, subs_parametersₚ)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= model_ages
                beta_t0dist[:,n÷sieve] .= subs_parameters
                hiatusdist[n÷sieve] = duration[1]
                # predicted_ages[n÷sieve] = beta_t0dist[2]+τ*log(1-(Sμ[?]*pi)/(E₀*beta_t0dist[1]*sin(pi/beta_t0dist[2])))
            end

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs,nsteps*sieve)

    # Crop the result
        agedist = agedist[active_height_t,:]

        subsmdl = SubsidenceStratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian(agedist,dim=2), # Median age
            nanpctile(agedist,2.5,dim=2), # 2.5th percentile
            nanpctile(agedist,97.5,dim=2), # 97.5th percentile
            nanmean(beta_t0dist[1,:],dim=1), # Mean beta
            nanstd(beta_t0dist[1,:],dim=1), # Standard deviation
            nanmedian(beta_t0dist[1,:],dim=1), # Median beta
            nanpctile(beta_t0dist[1,:],2.5,dim=1), # 2.5th percentile
            nanpctile(beta_t0dist[1,:],97.5,dim=1), # 97.5th percentile
            nanmean(beta_t0dist[2,:],dim=1), # Mean T0
            nanstd(beta_t0dist[2,:],dim=1), # Standard deviation
            nanmedian(beta_t0dist[2,:],dim=1), # Median T0
            nanpctile(beta_t0dist[2,:],2.5,dim=1), # 2.5th percentile
            nanpctile(beta_t0dist[2,:],97.5,dim=1) # 97.5th percentile
        )
    return subsmdl, agedist, lldist, hiatusdist, beta_t0dist
end

## --- End of File
