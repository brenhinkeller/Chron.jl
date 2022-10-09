## --- Stratigraphic MCMC model without hiatus # # # # # # # # # # # # # # # # #


    """
    ```julia
    StratMetropolis(smpl::ChronAgeData, [hiatus::HiatusData,] config::StratAgeModelConfiguration)
    ```
    Runs the main Chron.jl age-depth model routine for a stratigraphic set of
    samples defined by sample heights and simple Gaussian age constraints in the
    `smpl` struct, and an age-depth model configuration defined by the `config`
    struct.

    Optionally, if a `hiatus` struct is provided, the model will additionally
    incorporate information about the durations of known hiatuses at the
    specified model heights.

    ### Examples:
    ```julia
    (mdl, agedist, lldist) = StratMetropolis(smpl, config)
    ```
    ```julia
    (mdl, agedist, hiatusdist, lldist) = StratMetropolis(smpl, hiatus, config)
    ```
    """
    function StratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration)
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
        @inbounds for i ∈ eachindex(Age)
            if Age_Sidedness[i] == sign(closest_model_ages[i] - Age[i])
                closest_model_ages[i] = Age[i]
            end
        end
        ll = normpdf_ll(Age, Age_sigma, closest_model_ages)
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
                @inbounds for i ∈ eachindex(sample_heightₚ)
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
            @inbounds for i ∈ eachindex(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
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
                @inbounds for i ∈ eachindex(sample_heightₚ)
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
            @inbounds for i ∈ eachindex(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
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
        agedistₜ = copy(agedist)

        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian!(agedistₜ,dim=2), # Median age
            nanpctile!(agedistₜ,2.5,dim=2), # 2.5th percentile
            nanpctile!(agedistₜ,97.5,dim=2) # 97.5th percentile
        )

        return mdl, agedist, lldist
    end

## --- Stratigraphic MCMC model with hiata # # # # # # # # # # # # # # # # # # #

    function StratMetropolis(smpl::ChronAgeData, hiatus::HiatusData, config::StratAgeModelConfiguration)
        # Run stratigraphic MCMC model, with hiata
        print("Generating stratigraphic age-depth model...\n")

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints. Type assertions for stability
        Age = copy(smpl.Age)::Array{Float64,1}
        Age_sigma = copy(smpl.Age_sigma)::Array{Float64,1}
        Height = copy(smpl.Height)::Array{Float64,1}
        Height_sigma = smpl.Height_sigma::Array{Float64,1} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Array{Float64,1} # Bottom is a maximum age and top is a minimum age
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
        @inbounds for i ∈ eachindex(Age)
            if Age_Sidedness[i] == sign(closest_model_ages[i] - Age[i])
                closest_model_ages[i] = Age[i]
            end
        end
        ll = normpdf_ll(Age, Age_sigma, closest_model_ages)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

        # Ensure there is only one effective hiatus at most for each height node
        closest_hiatus = findclosestabove((hiatus.Height::Array{Float64,1}),model_heights)
        closest_hiatus_unique = unique(closest_hiatus)
        Hiatus_height = Array{Float64}(undef,size(closest_hiatus_unique))
        Hiatus_duration = Array{Float64}(undef,size(closest_hiatus_unique))
        Hiatus_duration_sigma = Array{Float64}(undef,size(closest_hiatus_unique))
        for i ∈ eachindex(closest_hiatus_unique)
            t = closest_hiatus.==closest_hiatus_unique[i]
            Hiatus_height[i] = mean((hiatus.Height::Array{Float64,1})[t])
            Hiatus_duration[i] = sum((hiatus.Duration::Array{Float64,1})[t])
            Hiatus_duration_sigma[i] = sqrt(sum((hiatus.Duration_sigma::Array{Float64,1})[t].^2))
        end

        # Add log likelihood for hiatus duration
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
                @inbounds for i ∈ eachindex(sample_heightₚ)
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
            @inbounds for i ∈ eachindex(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Add log likelihood for hiatus duration
            @. durationₚ = model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique]
            llₚ += normcdf_ll!(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(duration, durationₚ)
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
                @inbounds for i ∈ eachindex(sample_heightₚ)
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
            @inbounds for i ∈ eachindex(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = normpdf_ll(Age, Age_sigma, closest_model_agesₚ)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Add log likelihood for hiatus duration
            @. durationₚ = model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique]
            llₚ += normcdf_ll!(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

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
        agedistₜ = copy(agedist)

        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian!(agedistₜ,dim=2), # Median age
            nanpctile!(agedistₜ,2.5,dim=2), # 2.5th percentile
            nanpctile!(agedistₜ,97.5,dim=2) # 97.5th percentile
        )

        return mdl, agedist, lldist, hiatusdist
    end

## --- Stratigraphic MCMC model without hiatus, with distribution LL # # # # # #

    """
    ```julia
    StratMetropolisDist(smpl::ChronAgeData, [hiatus::HiatusData,] config::StratAgeModelConfiguration)
    ```
    Runs the main Chron.jl age-depth model routine for a stratigraphic set of
    samples defined by sample heights and fitted asymmetric age distributions
    (`bilinear_exponential`) in the `smpl` struct, and an age-depth model
    configuration defined by the `config` struct.

    Optionally, if a `hiatus` struct is provided, the model will additionally
    incorporate information about the durations of known hiatuses at the
    specified model heights.

    ### Examples:
    ```julia
    (mdl, agedist, lldist) = StratMetropolisDist(smpl, config)
    ```
    ```julia
    (mdl, agedist, hiatusdist, lldist) = StratMetropolisDist(smpl, hiatus, config)
    ```
    """
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
        @inbounds for i ∈ eachindex(Age)
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
                @inbounds for i ∈ eachindex(sample_heightₚ)
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
            @inbounds for i ∈ eachindex(Age)
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
                @inbounds for i ∈ eachindex(sample_heightₚ)
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
            @inbounds for i ∈ eachindex(Age)
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
        agedistₜ = copy(agedist)

        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian!(agedistₜ,dim=2), # Median age
            nanpctile!(agedistₜ,2.5,dim=2), # 2.5th percentile
            nanpctile!(agedistₜ,97.5,dim=2) # 97.5th percentile
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
        @inbounds for i ∈ eachindex(Age)
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
        for i ∈ eachindex(closest_hiatus_unique)
            t = closest_hiatus.==closest_hiatus_unique[i]
            Hiatus_height[i] = mean((hiatus.Height::Array{Float64,1})[t])
            Hiatus_duration[i] = sum((hiatus.Duration::Array{Float64,1})[t])
            Hiatus_duration_sigma[i] = sqrt(sum((hiatus.Duration_sigma::Array{Float64,1})[t].^2))
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
                @inbounds for i ∈ eachindex(sample_heightₚ)
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
            @inbounds for i ∈ eachindex(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = bilinear_exponential_ll(closest_model_agesₚ, p)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Add log likelihood for hiatus duration
            @. durationₚ = model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique]
            llₚ += normcdf_ll!(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

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
                @inbounds for i ∈ eachindex(sample_heightₚ)
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
            @inbounds for i ∈ eachindex(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = bilinear_exponential_ll(closest_model_agesₚ, p)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Add log likelihood for hiatus duration
            @. durationₚ = model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique]
            llₚ += normcdf_ll!(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

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
        agedistₜ = copy(agedist)

        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian!(agedistₜ,dim=2), # Median age
            nanpctile!(agedistₜ,2.5,dim=2), # 2.5th percentile
            nanpctile!(agedistₜ,97.5,dim=2) # 97.5th percentile
        )

        return mdl, agedist, lldist, hiatusdist
    end

## --- Stratigraphic MCMC model without hiatus, for radiocarbon ages # # # # # #

    """
    ```julia
    StratMetropolis14C(smpl::ChronAgeData, [hiatus::HiatusData,] config::StratAgeModelConfiguration)
    ```
    Runs the main Chron.jl age-depth model routine for a stratigraphic set of
    samples defined by sample heights and interpolated radiocarbon age
    constraints in the `smpl` struct, and an age-depth model configuration
    defined by the `config` struct.

    Optionally, if a `hiatus` struct is provided, the model will additionally
    incorporate information about the durations of known hiatuses at the
    specified model heights.

    ### Examples:
    ```julia
    (mdl, agedist, lldist) = StratMetropolis14C(smpl, config)
    ```
    ```julia
    (mdl, agedist, hiatusdist, lldist) = StratMetropolis14C(smpl, hiatus, config)
    ```
    """
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
            pl = normpdf_ll.(oldest + offset*dt_dH, boundsigma, 1:size(smpl.Params,1))
            pu = normpdf_ll.(youngest - offset*dt_dH, boundsigma, 1:size(smpl.Params,1))
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
        @inbounds for i ∈ eachindex(Age)
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
                @inbounds for i ∈ eachindex(sample_heightₚ)
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
            @inbounds for i ∈ eachindex(Age)
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
                @inbounds for i ∈ eachindex(sample_heightₚ)
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
            @inbounds for i ∈ eachindex(Age)
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
        agedistₜ = copy(agedist)

        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian!(agedistₜ,dim=2), # Median age
            nanpctile!(agedistₜ,2.5,dim=2), # 2.5th percentile
            nanpctile!(agedistₜ,97.5,dim=2) # 97.5th percentile
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
            pl = normpdf_ll.(oldest + offset*dt_dH, boundsigma, 1:size(smpl.Params,1))
            pu = normpdf_ll.(youngest - offset*dt_dH, boundsigma, 1:size(smpl.Params,1))
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
        @inbounds for i ∈ eachindex(Age)
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
        for i ∈ eachindex(closest_hiatus_unique)
            t = closest_hiatus.==closest_hiatus_unique[i]
            Hiatus_height[i] = mean((hiatus.Height::Array{Float64,1})[t])
            Hiatus_duration[i] = sum((hiatus.Duration::Array{Float64,1})[t])
            Hiatus_duration_sigma[i] = sqrt(sum((hiatus.Duration_sigma::Array{Float64,1})[t].^2))
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
                @inbounds for i ∈ eachindex(sample_heightₚ)
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
            @inbounds for i ∈ eachindex(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = interpolate_ll(closest_model_agesₚ, p)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Add log likelihood for hiatus duration
            @. durationₚ = model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique]
            llₚ += normcdf_ll!(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

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
                @inbounds for i ∈ eachindex(sample_heightₚ)
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
            @inbounds for i ∈ eachindex(Age)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
                if Age_Sidedness[i] == sign(closest_model_agesₚ[i] - Age[i])
                    closest_model_agesₚ[i] = Age[i]
                end
            end
            llₚ = interpolate_ll(closest_model_agesₚ, p)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Add log likelihood for hiatus duration
            @. durationₚ = model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique]
            llₚ += normcdf_ll!(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

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
        agedistₜ = copy(agedist)

        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist,dim=2), # Mean age
            nanstd(agedist,dim=2), # Standard deviation
            nanmedian!(agedistₜ,dim=2), # Median age
            nanpctile!(agedistₜ,2.5,dim=2), # 2.5th percentile
            nanpctile!(agedistₜ,97.5,dim=2) # 97.5th percentile
        )

        return mdl, agedist, lldist, hiatusdist
    end

## ---
