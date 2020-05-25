## --- Stratigraphic MCMC model without hiatus # # # # # # # # # # # # # # # # #

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
        mages = a .+ b .* collect(model_heights)


        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_mages = mages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_mages[i] - Age[i])
                closest_mages[i] = Age[i]
            end
        end
        ll = normpdf_ll(Age, Age_sigma, closest_mages)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

        # Preallocate variables for MCMC proposals
        ll_prop = ll
        mages_prop = copy(mages)
        closest_prop = copy(closest)
        sample_height_prop = copy(sample_height)
        closest_mages_prop = copy(closest_mages)

        # Run burnin
        # acceptancedist = fill(false,burnin)
        print("Burn-in: ", burnin, " steps\n")
        pgrs = Progress(burnin, desc="Burn-in...")
        for n=1:burnin
            # Prepare proposal
            copyto!(mages_prop, mages)
            copyto!(closest_prop, closest)
            copyto!(sample_height_prop, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_height_prop)
                    sample_height_prop[i] += randn() * Height_sigma[i]
                    closest_prop[i] = round(Int,(sample_height_prop[i] - model_heights[1])/resolution)+1
                    if closest_prop[i] < 1 # Check we're still within bounds
                        closest_prop[i] = 1
                    elseif closest_prop[i] > npoints
                        closest_prop[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if mages_prop[i] < mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if mages_prop[i] > mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                end
            end


            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_mages_prop[i] = mages_prop[closest_prop[i]]
                if Age_Sidedness[i] == sign(closest_mages_prop[i] - Age[i])
                    closest_mages_prop[i] = Age[i]
                end
            end
            ll_prop = normpdf_ll(Age, Age_sigma, closest_mages_prop)
            ll_prop += normpdf_ll(Height, Height_sigma, sample_height_prop)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                ll = ll_prop
                copyto!(mages, mages_prop)
                copyto!(closest, closest_prop)
                copyto!(sample_height, sample_height_prop)
                # acceptancedist[i] = true
            end

            # Update progress meter every 100 steps
            mod(n,100)==0 && update!(pgrs, n)
        end

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)


        # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(mages_prop, mages)
            copyto!(closest_prop, closest)
            copyto!(sample_height_prop, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_height_prop)
                    sample_height_prop[i] += randn() * Height_sigma[i]
                    closest_prop[i] = round(Int,(sample_height_prop[i] - model_heights[1])/resolution)+1
                    if closest_prop[i] < 1 # Check we're still within bounds
                        closest_prop[i] = 1
                    elseif closest_prop[i] > npoints
                        closest_prop[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if mages_prop[i] < mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if mages_prop[i] > mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_mages_prop[i] = mages_prop[closest_prop[i]]
                if Age_Sidedness[i] == sign(closest_mages_prop[i] - Age[i])
                    closest_mages_prop[i] = Age[i]
                end
            end
            ll_prop = normpdf_ll(Age, Age_sigma, closest_mages_prop)
            ll_prop += normpdf_ll(Height, Height_sigma, sample_height_prop)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                ll = ll_prop
                copyto!(mages, mages_prop)
                copyto!(closest, closest_prop)
                copyto!(sample_height, sample_height_prop)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= mages
            end

            # Update progress meter every 100 steps
            mod(n,100)==0 && update!(pgrs, n)
        end

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
        mages = a .+ b .* collect(model_heights)


        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_mages = mages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_mages[i] - Age[i])
                closest_mages[i] = Age[i]
            end
        end
        ll = normpdf_ll(Age, Age_sigma, closest_mages)
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
        duration = @. min(mages[closest_hiatus_unique - 1] - mages[closest_hiatus_unique], Hiatus_duration)
        ll += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, duration)

        # Preallocate variables for MCMC proposals
        ll_prop=ll
        chosen_point=0
        mages_prop = copy(mages)
        closest_prop = copy(closest)
        duration_prop = copy(duration)
        sample_height_prop = copy(sample_height)
        closest_mages_prop = copy(closest_mages)

        # Run burnin
        # acceptancedist = fill(false,burnin)
        print("Burn-in: ", burnin, " steps\n")
        pgrs = Progress(burnin, desc="Burn-in...")
        for n=1:burnin
            # Prepare proposal
            copyto!(mages_prop, mages)
            copyto!(closest_prop, closest)
            copyto!(sample_height_prop, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_height_prop)
                    sample_height_prop[i] += randn() * Height_sigma[i]
                    closest_prop[i] = round(Int,(sample_height_prop[i] - model_heights[1])/resolution)+1
                    if closest_prop[i] < 1 # Check we're still within bounds
                        closest_prop[i] = 1
                    elseif closest_prop[i] > npoints
                        closest_prop[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if mages_prop[i] < mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if mages_prop[i] > mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
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
                        mages_prop[chosen_point] = mages[chosen_point] + r
                        #Resolve conflicts
                        if r > 0 # If proposing increased age
                            @inbounds for i=1:chosen_point # younger points below
                                if mages_prop[i] < mages_prop[chosen_point]
                                    mages_prop[i] = mages_prop[chosen_point]
                                end
                            end
                        else # if proposing decreased age
                            @inbounds for i=chosen_point:npoints # older points above
                                if mages_prop[i] > mages_prop[chosen_point]
                                    mages_prop[i] = mages_prop[chosen_point]
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
                closest_mages_prop[i] = mages_prop[closest_prop[i]]
                if Age_Sidedness[i] == sign(closest_mages_prop[i] - Age[i])
                    closest_mages_prop[i] = Age[i]
                end
            end
            ll_prop = normpdf_ll(Age, Age_sigma, closest_mages_prop)
            ll_prop += normpdf_ll(Height, Height_sigma, sample_height_prop)

            # Add log likelihood for hiatus duration
            @. duration_prop = min(mages_prop[closest_hiatus_unique - 1] - mages_prop[closest_hiatus_unique], Hiatus_duration)
            ll_prop += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, duration_prop)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                ll = ll_prop
                copyto!(mages, mages_prop)
                copyto!(closest, closest_prop)
                copyto!(duration, duration_prop)
                copyto!(sample_height, sample_height_prop)
                # acceptancedist[i] = true
            end

            # Update progress meter every 100 steps
            mod(n,100)==0 && update!(pgrs, n)
        end

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)
        hiatusdist = Array{Float64}(undef,length(duration),nsteps)

        # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(mages_prop, mages)
            copyto!(closest_prop, closest)
            copyto!(sample_height_prop, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_height_prop)
                    sample_height_prop[i] += randn() * Height_sigma[i]
                    closest_prop[i] = round(Int,(sample_height_prop[i] - model_heights[1])/resolution)+1
                    if closest_prop[i] < 1 # Check we're still within bounds
                        closest_prop[i] = 1
                    elseif closest_prop[i] > npoints
                        closest_prop[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if mages_prop[i] < mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if mages_prop[i] > mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
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
                        mages_prop[chosen_point] = mages[chosen_point] + r
                        #Resolve conflicts
                        if r > 0 # If proposing increased age
                            @inbounds for i=1:chosen_point # younger points below
                                if mages_prop[i] < mages_prop[chosen_point]
                                    mages_prop[i] = mages_prop[chosen_point]
                                end
                            end
                        else # if proposing decreased age
                            @inbounds for i=chosen_point:npoints # older points above
                                if mages_prop[i] > mages_prop[chosen_point]
                                    mages_prop[i] = mages_prop[chosen_point]
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
                closest_mages_prop[i] = mages_prop[closest_prop[i]]
                if Age_Sidedness[i] == sign(closest_mages_prop[i] - Age[i])
                    closest_mages_prop[i] = Age[i]
                end
            end
            ll_prop = normpdf_ll(Age, Age_sigma, closest_mages_prop)
            ll_prop += normpdf_ll(Height, Height_sigma, sample_height_prop)

            # Add log likelihood for hiatus duration
            @. duration_prop = min(mages_prop[closest_hiatus_unique - 1] - mages_prop[closest_hiatus_unique], Hiatus_duration)
            ll_prop += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, duration_prop)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                ll = ll_prop
                copyto!(mages, mages_prop)
                copyto!(closest, closest_prop)
                copyto!(duration, duration_prop)
                copyto!(sample_height, sample_height_prop)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= mages
                hiatusdist[:,n÷sieve] .= duration
            end

            # Update progress meter every 100 steps
            mod(n,100)==0 && update!(pgrs, n)
        end

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
        mages = a .+ b .* collect(model_heights)


        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_mages = mages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_mages[i] - Age[i])
                closest_mages[i] = Age[i]
            end
        end
        ll = bilinear_exponential_ll(closest_mages ,p)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

        # Preallocate variables for MCMC proposals
        ll_prop = ll
        mages_prop = copy(mages)
        closest_prop = copy(closest)
        sample_height_prop = copy(sample_height)
        closest_mages_prop = copy(closest_mages)

        # Run burnin
        # acceptancedist = fill(false,burnin)
        print("Burn-in: ", burnin, " steps\n")
        pgrs = Progress(burnin, desc="Burn-in...")
        for n=1:burnin
            # Prepare proposal
            copyto!(mages_prop, mages)
            copyto!(closest_prop, closest)
            copyto!(sample_height_prop, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_height_prop)
                    sample_height_prop[i] += randn() * Height_sigma[i]
                    closest_prop[i] = round(Int,(sample_height_prop[i] - model_heights[1])/resolution)+1
                    if closest_prop[i] < 1 # Check we're still within bounds
                        closest_prop[i] = 1
                    elseif closest_prop[i] > npoints
                        closest_prop[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point
                        # younger points below
                        if mages_prop[i] < mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints
                        # older points above
                        if mages_prop[i] > mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                end
            end


            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_mages_prop[i] = mages_prop[closest_prop[i]]
                if Age_Sidedness[i] == sign(closest_mages_prop[i] - Age[i])
                    closest_mages_prop[i] = Age[i]
                end
            end
            ll_prop = bilinear_exponential_ll(closest_mages_prop, p)
            ll_prop += normpdf_ll(Height, Height_sigma, sample_height_prop)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                ll = ll_prop
                copyto!(mages, mages_prop)
                copyto!(closest, closest_prop)
                copyto!(sample_height, sample_height_prop)
                # acceptancedist[i] = true
            end

            # Update progress meter every 100 steps
            mod(n,100)==0 && update!(pgrs, n)
        end

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)


        # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(mages_prop, mages)
            copyto!(closest_prop, closest)
            copyto!(sample_height_prop, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_height_prop)
                    sample_height_prop[i] += randn() * Height_sigma[i]
                    closest_prop[i] = round(Int,(sample_height_prop[i] - model_heights[1])/resolution)+1
                    if closest_prop[i] < 1 # Check we're still within bounds
                        closest_prop[i] = 1
                    elseif closest_prop[i] > npoints
                        closest_prop[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point
                        # younger points below
                        if mages_prop[i] < mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints
                        # older points above
                        if mages_prop[i] > mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                end
            end


            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_mages_prop[i] = mages_prop[closest_prop[i]]
                if Age_Sidedness[i] == sign(closest_mages_prop[i] - Age[i])
                    closest_mages_prop[i] = Age[i]
                end
            end
            ll_prop = bilinear_exponential_ll(closest_mages_prop, p)
            ll_prop += normpdf_ll(Height, Height_sigma, sample_height_prop)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                ll = ll_prop
                copyto!(mages, mages_prop)
                copyto!(closest, closest_prop)
                copyto!(sample_height, sample_height_prop)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= mages
            end

            # Update progress meter every 100 steps
            mod(n,100)==0 && update!(pgrs, n)
        end

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
        mages = a .+ b .* collect(model_heights)


        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_mages = mages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_mages[i] - Age[i])
                closest_mages[i] = Age[i]
            end
        end
        ll = bilinear_exponential_ll(closest_mages ,p)
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
        duration = @. min(mages[closest_hiatus_unique - 1] - mages[closest_hiatus_unique], Hiatus_duration)
        ll += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, duration)

        # Preallocate variables for MCMC proposals
        ll_prop = ll
        chosen_point = 0
        mages_prop = copy(mages)
        closest_prop = copy(closest)
        duration_prop = copy(duration)
        sample_height_prop = copy(sample_height)
        closest_mages_prop = copy(closest_mages)

        # Run burnin
        # acceptancedist = fill(false,burnin)
        print("Burn-in: ", burnin, " steps\n")
        pgrs = Progress(burnin, desc="Burn-in...")
        for n=1:burnin
            # Prepare proposal
            copyto!(mages_prop, mages)
            copyto!(closest_prop, closest)
            copyto!(sample_height_prop, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_height_prop)
                    sample_height_prop[i] += randn() * Height_sigma[i]
                    closest_prop[i] = round(Int,(sample_height_prop[i] - model_heights[1])/resolution)+1
                    if closest_prop[i] < 1 # Check we're still within bounds
                        closest_prop[i] = 1
                    elseif closest_prop[i] > npoints
                        closest_prop[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point
                        # younger points below
                        if mages_prop[i] < mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints
                        # older points above
                        if mages_prop[i] > mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
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
                        mages_prop[chosen_point] = mages[chosen_point] + r
                        #Resolve conflicts
                        if r > 0 # If proposing increased age
                            @inbounds for i=1:chosen_point # younger points below
                                if mages_prop[i] < mages_prop[chosen_point]
                                    mages_prop[i] = mages_prop[chosen_point]
                                end
                            end
                        else # if proposing decreased age
                            @inbounds for i=chosen_point:npoints # older points above
                                if mages_prop[i] > mages_prop[chosen_point]
                                    mages_prop[i] = mages_prop[chosen_point]
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
                closest_mages_prop[i] = mages_prop[closest_prop[i]]
                if Age_Sidedness[i] == sign(closest_mages_prop[i] - Age[i])
                    closest_mages_prop[i] = Age[i]
                end
            end
            ll_prop = bilinear_exponential_ll(closest_mages_prop, p)
            ll_prop += normpdf_ll(Height, Height_sigma, sample_height_prop)

            # Add log likelihood for hiatus duration
            @. duration_prop = min(mages_prop[closest_hiatus_unique - 1] - mages_prop[closest_hiatus_unique], Hiatus_duration)
            ll_prop += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, duration_prop)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                ll = ll_prop
                copyto!(mages, mages_prop)
                copyto!(closest, closest_prop)
                copyto!(sample_height, sample_height_prop)
                # acceptancedist[i] = true
            end

            # Update progress meter every 100 steps
            mod(n,100)==0 && update!(pgrs, n)
        end

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)
        hiatusdist = Array{Float64}(undef,length(duration),nsteps)

        # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(mages_prop, mages)
            copyto!(closest_prop, closest)
            copyto!(sample_height_prop, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_height_prop)
                    sample_height_prop[i] += randn() * Height_sigma[i]
                    closest_prop[i] = round(Int,(sample_height_prop[i] - model_heights[1])/resolution)+1
                    if closest_prop[i] < 1 # Check we're still within bounds
                        closest_prop[i] = 1
                    elseif closest_prop[i] > npoints
                        closest_prop[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if mages_prop[i] < mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if mages_prop[i] > mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
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
                        mages_prop[chosen_point] = mages[chosen_point] + r
                        #Resolve conflicts
                        if r > 0 # If proposing increased age
                            @inbounds for i=1:chosen_point # younger points below
                                if mages_prop[i] < mages_prop[chosen_point]
                                    mages_prop[i] = mages_prop[chosen_point]
                                end
                            end
                        else # if proposing decreased age
                            @inbounds for i=chosen_point:npoints # older points above
                                if mages_prop[i] > mages_prop[chosen_point]
                                    mages_prop[i] = mages_prop[chosen_point]
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
                closest_mages_prop[i] = mages_prop[closest_prop[i]]
                if Age_Sidedness[i] == sign(closest_mages_prop[i] - Age[i])
                    closest_mages_prop[i] = Age[i]
                end
            end
            ll_prop = bilinear_exponential_ll(closest_mages_prop, p)
            ll_prop += normpdf_ll(Height, Height_sigma, sample_height_prop)

            # Add log likelihood for hiatus duration
            @. duration_prop = min(mages_prop[closest_hiatus_unique - 1] - mages_prop[closest_hiatus_unique], Hiatus_duration)
            ll_prop += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, duration_prop)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                ll = ll_prop
                copyto!(mages, mages_prop)
                copyto!(closest, closest_prop)
                copyto!(duration, duration_prop)
                copyto!(sample_height, sample_height_prop)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= mages
                hiatusdist[:,n÷sieve] .= duration
            end

            # Update progress meter every 100 steps
            mod(n,100)==0 && update!(pgrs, n)
        end

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
        mages = a .+ b .* collect(model_heights)


        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_mages = mages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_mages[i] - Age[i])
                closest_mages[i] = Age[i]
            end
        end
        ll = interpolate_ll(closest_mages, p)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

        # Preallocate variables for MCMC proposals
        ll_prop = ll
        mages_prop = copy(mages)
        closest_prop = copy(closest)
        sample_height_prop = copy(sample_height)
        closest_mages_prop = copy(closest_mages)

        # Run burnin
        # acceptancedist = fill(false,burnin)
        print("Burn-in: ", burnin, " steps\n")
        pgrs = Progress(burnin, desc="Burn-in...")
        for n=1:burnin
            copyto!(mages_prop, mages)
            copyto!(closest_prop, closest)
            copyto!(sample_height_prop, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_height_prop)
                    sample_height_prop[i] += randn() * Height_sigma[i]
                    closest_prop[i] = round(Int,(sample_height_prop[i] - model_heights[1])/resolution)+1
                    if closest_prop[i] < 1 # Check we're still within bounds
                        closest_prop[i] = 1
                    elseif closest_prop[i] > npoints
                        closest_prop[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if mages_prop[i] < mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if mages_prop[i] > mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                end
            end


            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_mages_prop[i] = mages_prop[closest_prop[i]]
                if Age_Sidedness[i] == sign(closest_mages_prop[i] - Age[i])
                    closest_mages_prop[i] = Age[i]
                end
            end
            ll_prop = interpolate_ll(closest_mages_prop, p)
            ll_prop += normpdf_ll(Height, Height_sigma, sample_height_prop)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                ll = ll_prop
                copyto!(mages, mages_prop)
                copyto!(closest, closest_prop)
                copyto!(sample_height, sample_height_prop)
                # acceptancedist[i] = true
            end

            # Update progress meter every 100 steps
            mod(n,100)==0 && update!(pgrs, n)
        end

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)


        # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(mages_prop, mages)
            copyto!(closest_prop, closest)
            copyto!(sample_height_prop, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_height_prop)
                    sample_height_prop[i] += randn() * Height_sigma[i]
                    closest_prop[i] = round(Int,(sample_height_prop[i] - model_heights[1])/resolution)+1
                    if closest_prop[i] < 1 # Check we're still within bounds
                        closest_prop[i] = 1
                    elseif closest_prop[i] > npoints
                        closest_prop[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if mages_prop[i] < mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if mages_prop[i] > mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            @inbounds for i=1:length(Age)
                closest_mages_prop[i] = mages_prop[closest_prop[i]]
                if Age_Sidedness[i] == sign(closest_mages_prop[i] - Age[i])
                    closest_mages_prop[i] = Age[i]
                end
            end
            ll_prop = interpolate_ll(closest_mages_prop, p)
            ll_prop += normpdf_ll(Height, Height_sigma, sample_height_prop)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                ll = ll_prop
                copyto!(mages, mages_prop)
                copyto!(closest, closest_prop)
                copyto!(sample_height, sample_height_prop)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= mages
            end

            # Update progress meter every 100 steps
            mod(n,100)==0 && update!(pgrs, n)
        end

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
        mages = a .+ b .* collect(model_heights)


        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_mages = mages[closest]
        @inbounds for i=1:length(Age)
            if Age_Sidedness[i] == sign(closest_mages[i] - Age[i])
                closest_mages[i] = Age[i]
            end
        end
        ll = interpolate_ll(closest_mages, p)
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
        duration = @. min(mages[closest_hiatus_unique - 1] - mages[closest_hiatus_unique], Hiatus_duration)
        ll += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, duration)

        # Preallocate variables for MCMC proposals
        ll_prop = ll
        chosen_point = 0
        mages_prop = copy(mages)
        closest_prop = copy(closest)
        duration_prop = copy(duration)
        sample_height_prop = copy(sample_height)
        closest_mages_prop = copy(closest_mages)

        # Run burnin
        # acceptancedist = fill(false,burnin)
        print("Burn-in: ", burnin, " steps\n")
        pgrs = Progress(burnin, desc="Burn-in...")
        for n=1:burnin
            # Prepare proposal
            copyto!(mages_prop, mages)
            copyto!(closest_prop, closest)
            copyto!(sample_height_prop, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_height_prop)
                    sample_height_prop[i] += randn() * Height_sigma[i]
                    closest_prop[i] = round(Int,(sample_height_prop[i] - model_heights[1])/resolution)+1
                    if closest_prop[i] < 1 # Check we're still within bounds
                        closest_prop[i] = 1
                    elseif closest_prop[i] > npoints
                        closest_prop[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if mages_prop[i] < mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if mages_prop[i] > mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
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
                        mages_prop[chosen_point] = mages[chosen_point] + r
                        #Resolve conflicts
                        if r > 0 # If proposing increased age
                            @inbounds for i=1:chosen_point # younger points below
                                if mages_prop[i] < mages_prop[chosen_point]
                                    mages_prop[i] = mages_prop[chosen_point]
                                end
                            end
                        else # if proposing decreased age
                            @inbounds for i=chosen_point:npoints # older points above
                                if mages_prop[i] > mages_prop[chosen_point]
                                    mages_prop[i] = mages_prop[chosen_point]
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
                closest_mages_prop[i] = mages_prop[closest_prop[i]]
                if Age_Sidedness[i] == sign(closest_mages_prop[i] - Age[i])
                    closest_mages_prop[i] = Age[i]
                end
            end
            ll_prop = interpolate_ll(closest_mages_prop, p)
            ll_prop += normpdf_ll(Height, Height_sigma, sample_height_prop)

            # Add log likelihood for hiatus duration
            @. duration_prop = min(mages_prop[closest_hiatus_unique - 1] - mages_prop[closest_hiatus_unique], Hiatus_duration)
            ll_prop += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, duration_prop)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                ll = ll_prop
                copyto!(mages, mages_prop)
                copyto!(closest, closest_prop)
                copyto!(sample_height, sample_height_prop)
                # acceptancedist[i] = true
            end

            # Update progress meter every 100 steps
            mod(n,100)==0 && update!(pgrs, n)
        end

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)
        hiatusdist = Array{Float64}(undef,length(duration),nsteps)

        # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(mages_prop, mages)
            copyto!(closest_prop, closest)
            copyto!(sample_height_prop, sample_height)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i=1:length(sample_height_prop)
                    sample_height_prop[i] += randn() * Height_sigma[i]
                    closest_prop[i] = round(Int,(sample_height_prop[i] - model_heights[1])/resolution)+1
                    if closest_prop[i] < 1 # Check we're still within bounds
                        closest_prop[i] = 1
                    elseif closest_prop[i] > npoints
                        closest_prop[i] = npoints
                    end
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] += r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    @inbounds for i=1:chosen_point # younger points below
                        if mages_prop[i] < mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
                        end
                    end
                else # if proposing decreased age
                    @inbounds for i=chosen_point:npoints # older points above
                        if mages_prop[i] > mages_prop[chosen_point]
                            mages_prop[i] = mages_prop[chosen_point]
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
                        mages_prop[chosen_point] = mages[chosen_point] + r
                        #Resolve conflicts
                        if r > 0 # If proposing increased age
                            @inbounds for i=1:chosen_point # younger points below
                                if mages_prop[i] < mages_prop[chosen_point]
                                    mages_prop[i] = mages_prop[chosen_point]
                                end
                            end
                        else # if proposing decreased age
                            @inbounds for i=chosen_point:npoints # older points above
                                if mages_prop[i] > mages_prop[chosen_point]
                                    mages_prop[i] = mages_prop[chosen_point]
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
                closest_mages_prop[i] = mages_prop[closest_prop[i]]
                if Age_Sidedness[i] == sign(closest_mages_prop[i] - Age[i])
                    closest_mages_prop[i] = Age[i]
                end
            end
            ll_prop = interpolate_ll(closest_mages_prop, p)
            ll_prop += normpdf_ll(Height, Height_sigma, sample_height_prop)

            # Add log likelihood for hiatus duration
            @. duration_prop = min(mages_prop[closest_hiatus_unique - 1] - mages_prop[closest_hiatus_unique], Hiatus_duration)
            ll_prop += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, duration_prop)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                ll = ll_prop
                copyto!(mages, mages_prop)
                copyto!(closest, closest_prop)
                copyto!(duration, duration_prop)
                copyto!(sample_height, sample_height_prop)
            end

            # Record sieved results
            if mod(n,sieve) == 0
                lldist[n÷sieve] = ll
                agedist[:,n÷sieve] .= mages
                hiatusdist[:,n÷sieve] .= duration
            end

            # Update progress meter every 100 steps
            mod(n,100)==0 && update!(pgrs, n)
        end

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
