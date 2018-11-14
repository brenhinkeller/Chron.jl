## --- Custom objects for holding stratigraphic age data

    # Define a type of struct object to hold data about geochron samples from
    # a stratigraphic section of interest
    mutable struct StratAgeData
        Name::Tuple
        Height::Array{Float64}
        Height_sigma::Array{Float64}
        Age::Array{Float64}
        Age_sigma::Array{Float64}
        Age_025CI::Array{Float64}
        Age_975CI::Array{Float64}
        Age_Sidedness::Array{Float64}
        Params::Array{Float64}
        Path::String
        inputSigmaLevel::Int
        Age_Unit::String
        Height_Unit::String
    end

    function NewStratAgeData(nSamples)
        smpl = StratAgeData(
            ("Sample Names",),
            fill(NaN,nSamples),  # Sample heights
            fill(NaN,nSamples),  # Height_sigma
            fill(NaN,nSamples),  # Sample ages
            fill(NaN,nSamples),  # Sample age uncertainty
            fill(NaN,nSamples),  # Sample age 2.5% CI
            fill(NaN,nSamples),  # Sample age 97.5% CI
            zeros(nSamples), # Sidedness (zeros by default, geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
            fill(NaN,5,nSamples), # Sample age distribution parameters
            "Relative path where we can find .csv data files",
            2, # i.e., are the data files 1-sigma or 2-sigma
            "Ma",
            "m",
        )
        return smpl
    end
    export NewStratAgeData

    # A type of object to hold data about hiatuses
    mutable struct HiatusData
        Height::Array{Float64}
        Height_sigma::Array{Float64}
        Duration::Array{Float64}
        Duration_sigma::Array{Float64}
        Age_Unit::String
        Height_Unit::String
    end

    function NewHiatusData(nHiatuses)
        hiatus = HiatusData(
            fill(NaN,nHiatuses),  # Height
            fill(NaN,nHiatuses),  # Height_sigma
            fill(NaN,nHiatuses),  # Duration
            fill(NaN,nHiatuses),  # Duration_sigma
            "Ma",
            "m",
        )
        return hiatus
    end
    export NewHiatusData

    # A type of object to specify the configuration of the stratigraphic model
    mutable struct StratAgeModelConfiguration
        resolution::Float64
        burnin::Int
        nsteps::Int
        sieve::Int
        bounding::Float64
    end

    function NewStratAgeModelConfiguration()
        return StratAgeModelConfiguration(NaN, 0, 0, 0, NaN)
    end
    export NewStratAgeModelConfiguration


    struct StratAgeModel
        Height::Array{Float64}
        Age::Array{Float64}
        Age_sigma::Array{Float64}
        Age_Median::Array{Float64}
        Age_025CI::Array{Float64}
        Age_975CI::Array{Float64}
    end

## --- Stratigraphic MCMC model without hiatus # # # # # # # # # # # # # # # # #

    function StratMetropolis(smpl::StratAgeData,config::StratAgeModelConfiguration)
        # Run stratigraphic MCMC model

        print("Generating stratigraphic age-depth model...\n")

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints
        (bottom, top) = extrema(smpl.Height)
        (youngest, oldest) = extrema(smpl.Age)
        dt_dH = (oldest-youngest)/(top-bottom)
        aveuncert = mean(smpl.Age_sigma)

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            offset = (top-bottom)*bounding
            Age = [oldest + offset*dt_dH; smpl.Age; youngest - offset*dt_dH]
            Age_sigma = [mean(smpl.Age_sigma)/10; smpl.Age_sigma; mean(smpl.Age_sigma)/10]
            Height = [bottom-offset; smpl.Height; top+offset]
            Height_sigma = [0; smpl.Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; smpl.Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
        else
            Age = smpl.Age
            Age_sigma = smpl.Age_sigma
            Height = smpl.Height
            Height_sigma = smpl.Height_sigma .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = smpl.Age_Sidedness # Bottom is a maximum age and top is a minimum age
            model_heights = bottom:resolution:top
        end
        active_height_t = (model_heights .>= bottom) .& (model_heights .<= top)
        npoints = length(model_heights)

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        mages = a .+ b .* collect(model_heights)


        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = Height
        closest = findclosest(sample_height,model_heights)
        agell = .- (mages[closest] .- Age).^2 ./ (2 .* Age_sigma.^2) .- log.(sqrt.(2*pi*Age_sigma))
        heightll = .-(sample_height .- Height).^2 ./ (2 .* Height_sigma.^2) .- log.(sqrt.(2*pi*Height_sigma))
        diff_sign = Age_Sidedness .!= sign.(mages[closest] .- Age)
        ll = sum(agell[diff_sign]) + sum(-log.(sqrt.(2*pi*Age_sigma[.~diff_sign]))) + sum(heightll)

        # Introduce variables so they will be accessible outside loop
        mages_prop = copy(mages)
        closest_prop = copy(closest)
        sample_height_prop = copy(sample_height)
        ll_prop = copy(ll)

        acceptancedist = fill(false,burnin)
        # Run burnin
        print("Burn-in: ", burnin, " steps\n")
        index = collect(1:npoints)
        @showprogress 5 "Burn-in..." for i=1:burnin
            mages_prop = copy(mages)
            closest_prop = copy(closest)
            sample_height_prop = copy(sample_height)

            if rand() < 0.1
                # Adjust heights
                sample_height_prop += randn(size(Height)) .* Height_sigma
                closest_prop = findclosest(sample_height_prop, model_heights)
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] = mages[chosen_point] + r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    younger_points_below = (mages_prop .< mages_prop[chosen_point]) .& (index .< chosen_point)
                    mages_prop[younger_points_below] .= mages_prop[chosen_point]
                else # if proposing decreased age
                    older_points_above = (mages_prop .> mages_prop[chosen_point]) .& (index .> chosen_point)
                    mages_prop[older_points_above] .= mages_prop[chosen_point]
                end
            end


            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            agell_prop = .- (mages_prop[closest_prop] .- Age).^2 ./ (2 .* Age_sigma.^2) -log.(sqrt.(2*pi*Age_sigma))
            heightll_prop = .- (sample_height_prop .- Height).^2 ./ (2 .* Height_sigma.^2) -log.(sqrt.(2*pi*Height_sigma))
            diff_sign = Age_Sidedness .!= sign.(mages_prop[closest_prop]-Age)
            ll_prop = sum(agell_prop[diff_sign]) + sum(-log.(sqrt.(2*pi*Age_sigma[.~diff_sign]))) + sum(heightll_prop)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                mages = copy(mages_prop)
                sample_height = copy(sample_height_prop)
                closest = copy(closest_prop)
                ll = copy(ll_prop)
                acceptancedist[i] = true
            end
        end

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)
        # acceptancedist = Array{Bool}(undef,nsteps)


        # Run the model
        @showprogress 5 "Collecting..." for i=1:(nsteps*sieve)
            mages_prop = copy(mages)
            closest_prop = copy(closest)
            sample_height_prop = copy(sample_height)

            if rand() < 0.1
                # Adjust heights
                sample_height_prop += randn(size(Height)) .* Height_sigma
                closest_prop = findclosest(sample_height_prop, model_heights)
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] = mages[chosen_point] + r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    younger_points_below = (mages_prop .< mages_prop[chosen_point]) .& (index .< chosen_point)
                    mages_prop[younger_points_below] .= mages_prop[chosen_point]
                else # if proposing decreased age
                    older_points_above = (mages_prop .> mages_prop[chosen_point]) .& (index .> chosen_point)
                    mages_prop[older_points_above] .= mages_prop[chosen_point]
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            agell_prop = .- (mages_prop[closest_prop] .- Age).^2 ./ (2 .* Age_sigma.^2) -log.(sqrt.(2*pi*Age_sigma))
            heightll_prop = .- (sample_height_prop .- Height).^2 ./ (2 .* Height_sigma.^2) -log.(sqrt.(2*pi*Height_sigma))
            diff_sign = Age_Sidedness .!= sign.(mages_prop[closest_prop]-Age)
            ll_prop = sum(agell_prop[diff_sign]) + sum(-log.(sqrt.(2*pi*Age_sigma[.~diff_sign]))) + sum(heightll_prop)


            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                mages = copy(mages_prop)
                sample_height = copy(sample_height_prop)
                closest = copy(closest_prop)
                ll = copy(ll_prop)
            end

            # Record sieved results
            if mod(i,sieve) == 0
                agedist[:,Int(i/sieve)] = mages
                lldist[Int(i/sieve)] = ll
            end
        end

        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist[active_height_t,:],dim=2), # Mean age
            nanstd(agedist[active_height_t,:],dim=2), # Standard deviation
            nanmedian(agedist[active_height_t,:],dim=2), # Median age
            pctile(agedist[active_height_t,:],2.5,dim=2), # 2.5th percentile
            pctile(agedist[active_height_t,:],97.5,dim=2) # 97.5th percentile
        )

        return mdl, agedist[active_height_t,:], lldist
    end

## --- Stratigraphic MCMC model with hiata # # # # # # # # # # # # # # # # # # #

    function StratMetropolisHiatus(smpl::StratAgeData,hiatus::HiatusData,config::StratAgeModelConfiguration)
        # Run stratigraphic MCMC model, with hiata

        print("Generating stratigraphic age-depth model...\n")

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints
        (bottom, top) = extrema(smpl.Height)
        (youngest, oldest) = extrema(smpl.Age)
        dt_dH = (oldest-youngest)/(top-bottom)
        aveuncert = mean(smpl.Age_sigma)

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            offset = (top-bottom)*bounding
            Age = [oldest + offset*dt_dH; smpl.Age; youngest - offset*dt_dH]
            Age_sigma = [mean(smpl.Age_sigma)/10; smpl.Age_sigma; mean(smpl.Age_sigma)/10]
            Height = [bottom-offset; smpl.Height; top+offset]
            Height_sigma = [0; smpl.Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; smpl.Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
        else
            Age = smpl.Age
            Age_sigma = smpl.Age_sigma
            Height = smpl.Height
            Height_sigma = smpl.Height_sigma .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = smpl.Age_Sidedness # Bottom is a maximum age and top is a minimum age
            model_heights = bottom:resolution:top
        end
        active_height_t = (model_heights .>= bottom) .& (model_heights .<= top)
        npoints = length(model_heights)

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        mages = a .+ b .* collect(model_heights)


        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = Height
        closest = findclosest(sample_height,model_heights)
        agell = .- (mages[closest] .- Age).^2 ./ (2 .* Age_sigma.^2) .- log.(sqrt.(2*pi*Age_sigma))
        heightll = .- (sample_height .- Height).^2 ./ (2 .* Height_sigma.^2) .- log.(sqrt.(2*pi*Height_sigma))
        diff_sign = Age_Sidedness .!= sign.(mages[closest] .- Age)
        ll = sum(agell[diff_sign]) + sum(-log.(sqrt.(2*pi*Age_sigma[.~diff_sign]))) + sum(heightll)

        # Ensure there is only one effective hiatus at most for each height node
        closest_h = findclosestabove(hiatus.Height,model_heights)
        closest_h_unique = unique(closest_h)
        hiatus_Height = Array{Float64}(undef,size(closest_h_unique))
        hiatus_Duration = Array{Float64}(undef,size(closest_h_unique))
        hiatus_Duration_sigma = Array{Float64}(undef,size(closest_h_unique))
        for i=1:length(closest_h_unique)
            hiatus_Height[i] = mean(hiatus.Height[closest_h.==closest_h_unique[i]])
            hiatus_Duration[i] = sum(hiatus.Duration[closest_h.==closest_h_unique[i]])
            hiatus_Duration_sigma[i] = sqrt(sum(hiatus.Duration_sigma[closest_h.==closest_h_unique[i]].^2))
        end

        # Add log likelihood for hiatuses
        duration_prop = abs.(mages[closest_h_unique .- 1] .- mages[closest_h_unique])
        ll += sum(-max.(hiatus_Duration .- duration_prop, 0).^2 ./ (2 .* hiatus_Duration_sigma.^2)) #-log.(sqrt.(2*pi*hiatus_Duration_sigma))
        duration = copy(duration_prop)

        # Introduce variables so they will be accessible outside loop
        mages_prop = copy(mages)
        closest_prop = copy(closest)
        sample_height_prop = copy(sample_height)
        ll_prop = copy(ll)
        chosen_point = 0

        acceptancedist = fill(false,burnin)
        # Run burnin
        print("Burn-in: ", burnin, " steps\n")
        index = collect(1:npoints)
        @showprogress 5 "Burn-in..." for i=1:burnin
            mages_prop = copy(mages)
            closest_prop = copy(closest)
            sample_height_prop = copy(sample_height)

            if rand() < 0.1
                # Adjust heights
                sample_height_prop += randn(size(Height)) .* Height_sigma
                closest_prop = findclosest(sample_height_prop, model_heights)
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] = mages[chosen_point] + r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    younger_points_below = (mages_prop .< mages_prop[chosen_point]) .& (index .< chosen_point)
                    mages_prop[younger_points_below] .= mages_prop[chosen_point]
                else # if proposing decreased age
                    older_points_above = (mages_prop .> mages_prop[chosen_point]) .& (index .> chosen_point)
                    mages_prop[older_points_above] .= mages_prop[chosen_point]
                end

                # If chosen_point is a hiatus point, let there be a 20 percent chance of
                # adjusting the point below the hiatus as well
                if rand() < 0.2
                    # if hiatus_height_uncert>0
                    #     closest_h = findclosestabove(h.Height+randn(size(h.Height)).*hiatus_height_uncert,heights)
                    # end
                    if any(closest_h_unique.==chosen_point)
                        mages_prop[chosen_point-1] = mages[chosen_point-1] + r
                        #Resolve conflicts
                        if r>0
                            mages_prop[(mages_prop .< mages_prop[chosen_point-1]) .& (index .< [chosen_point-1])] .= mages_prop[chosen_point-1]
                        else
                            mages_prop[(mages_prop .> mages_prop[chosen_point-1]) .& (index .> [chosen_point-1])] .= mages_prop[chosen_point-1]
                        end
                    end
                end
            end


            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            agell_prop = .- (mages_prop[closest_prop] .- Age).^2 ./ (2 .* Age_sigma.^2) -log.(sqrt.(2*pi*Age_sigma))
            heightll_prop = .- (sample_height_prop .- Height).^2 ./ (2 .* Height_sigma.^2) -log.(sqrt.(2*pi*Height_sigma))
            diff_sign = Age_Sidedness .!= sign.(mages_prop[closest_prop]-Age)
            ll_prop = sum(agell_prop[diff_sign]) + sum(-log.(sqrt.(2*pi*Age_sigma[.~diff_sign]))) + sum(heightll_prop)

            # if hiatus_height_uncert>0
            #     closest_h = findclosestabove(h.Height+randn(size(h.Height)).*hiatus_height_uncert,heights)
            # end
            duration_prop = abs.(mages_prop[closest_h_unique .- 1] .- mages_prop[closest_h_unique])
            ll_prop += sum(-max.(hiatus_Duration .- duration_prop, 0).^2 ./ (2 .* hiatus_Duration_sigma.^2)) #-log.(sqrt.(2*pi*hiatus_Duration_sigma))

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                mages = copy(mages_prop)
                sample_height = copy(sample_height_prop)
                closest = copy(closest_prop)
                ll = copy(ll_prop)
                acceptancedist[i] = true
            end
        end

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        hiatusdist = Array{Float64}(undef,length(duration),nsteps)
        lldist = Array{Float64}(undef,nsteps)
        # acceptancedist = Array{Bool}(undef,nsteps)


        # Run the model
        @showprogress 5 "Collecting..." for i=1:(nsteps*sieve)
            mages_prop = copy(mages)
            closest_prop = copy(closest)
            sample_height_prop = copy(sample_height)

            if rand() < 0.1
                # Adjust heights
                sample_height_prop += randn(size(Height)) .* Height_sigma
                closest_prop = findclosest(sample_height_prop, model_heights)
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] = mages[chosen_point] + r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    younger_points_below = (mages_prop .< mages_prop[chosen_point]) .& (index .< chosen_point)
                    mages_prop[younger_points_below] .= mages_prop[chosen_point]
                else # if proposing decreased age
                    older_points_above = (mages_prop .> mages_prop[chosen_point]) .& (index .> chosen_point)
                    mages_prop[older_points_above] .= mages_prop[chosen_point]
                end

                # If chosen_point is a hiatus point, let there be a 20 percent chance of
                # adjusting the point below the hiatus as well
                if rand() < 0.2
                    # if hiatus_height_uncert>0
                    #     closest_h = findclosestabove(h.Height+randn(size(h.Height)).*hiatus_height_uncert,heights)
                    # end
                    if any(closest_h_unique.==chosen_point)
                        mages_prop[chosen_point-1] = mages[chosen_point-1] + r
                        #Resolve conflicts
                        if r>0
                            younger_points_below = (mages_prop .< mages_prop[chosen_point-1]) .& (index .< (chosen_point-1))
                            mages_prop[younger_points_below] .= mages_prop[chosen_point-1]
                        else
                            older_points_above = (mages_prop .> mages_prop[chosen_point-1]) .& (index .> (chosen_point-1))
                            mages_prop[older_points_above] .= mages_prop[chosen_point-1]
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            agell_prop = .- (mages_prop[closest_prop] .- Age).^2 ./ (2 .* Age_sigma.^2) -log.(sqrt.(2*pi*Age_sigma))
            heightll_prop = .- (sample_height_prop .- Height).^2 ./ (2 .* Height_sigma.^2) -log.(sqrt.(2*pi*Height_sigma))
            diff_sign = Age_Sidedness .!= sign.(mages_prop[closest_prop]-Age)
            ll_prop = sum(agell_prop[diff_sign]) + sum(-log.(sqrt.(2*pi*Age_sigma[.~diff_sign]))) + sum(heightll_prop)


            # Add log likelihood for duration
            # if hiatus_height_uncert>0
            #     closest_h = findclosestabove(h.Height+randn(size(h.Height)).*hiatus_height_uncert,heights)
            # end
            duration_prop = abs.(mages_prop[closest_h_unique .- 1] .- mages_prop[closest_h_unique])
            ll_prop += sum(-max.(hiatus_Duration .- duration_prop, 0).^2 ./ (2 .* hiatus_Duration_sigma.^2)) #-log.(sqrt.(2*pi*hiatus_Duration_sigma))


            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                mages = copy(mages_prop)
                sample_height = copy(sample_height_prop)
                closest = copy(closest_prop)
                ll = copy(ll_prop)
                duration = copy(duration_prop)
            end

            # Record sieved results
            if mod(i,sieve) == 0
                agedist[:,Int(i/sieve)] = mages
                lldist[Int(i/sieve)] = ll
                hiatusdist[:,Int(i/sieve)] = duration
            end
        end

        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist[active_height_t,:],dim=2), # Mean age
            nanstd(agedist[active_height_t,:],dim=2), # Standard deviation
            nanmedian(agedist[active_height_t,:],dim=2), # Median age
            pctile(agedist[active_height_t,:],2.5,dim=2), # 2.5th percentile
            pctile(agedist[active_height_t,:],97.5,dim=2) # 97.5th percentile
        )

        return mdl, agedist[active_height_t,:], lldist, hiatusdist
    end

## --- Stratigraphic MCMC model without hiatus, with distribution LL # # # # # #

    function StratMetropolisDist(smpl::StratAgeData,config::StratAgeModelConfiguration)
        # Run stratigraphic MCMC model

        print("Generating stratigraphic age-depth model...\n")

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints
        (bottom, top) = extrema(smpl.Height)
        (youngest, oldest) = extrema(smpl.Age)
        dt_dH = (oldest-youngest)/(top-bottom)
        aveuncert = mean(smpl.Age_sigma)
        p = smpl.Params

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            offset = (top-bottom)*bounding
            Age = [oldest + offset*dt_dH; smpl.Age; youngest - offset*dt_dH]
            Age_sigma = [mean(smpl.Age_sigma)/10; smpl.Age_sigma; mean(smpl.Age_sigma)/10]
            Height = [bottom-offset; smpl.Height; top+offset]
            Height_sigma = [0; smpl.Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; smpl.Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
            pl = ones(5); pl[2] = oldest + offset*dt_dH; pl[3] = mean(smpl.Age_sigma)/10
            pu = ones(5); pu[2] = youngest - offset*dt_dH; pu[3] = mean(smpl.Age_sigma)/10
            p = hcat(pl,p,pu) # Add parameters for upper and lower runaway bounds
        else
            Age = smpl.Age
            Age_sigma = smpl.Age_sigma
            Height = smpl.Height
            Height_sigma = smpl.Height_sigma .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = smpl.Age_Sidedness # Bottom is a maximum age and top is a minimum age
            model_heights = bottom:resolution:top
        end
        active_height_t = (model_heights .>= bottom) .& (model_heights .<= top)
        npoints = length(model_heights)

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        mages = a .+ b .* collect(model_heights)


        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = Height
        closest = findclosest(sample_height,model_heights)
        agell = bilinear_exponential_LL(mages[closest],p)
        heightll = .- (sample_height .- Height).^2 ./ (2 .* Height_sigma.^2) .- log.(sqrt.(2*pi*Height_sigma))
        diff_sign = Age_Sidedness .!= sign.(mages[closest] .- Age)
        ll = sum(agell[diff_sign]) + sum(-log.(sqrt.(2*pi*Age_sigma[.~diff_sign]))) + sum(heightll)

        # Introduce variables so they will be accessible outside loop
        mages_prop = copy(mages)
        closest_prop = copy(closest)
        sample_height_prop = copy(sample_height)
        ll_prop = copy(ll)

        acceptancedist = fill(false,burnin)
        # Run burnin
        print("Burn-in: ", burnin, " steps\n")
        index = collect(1:npoints)
        @showprogress 5 "Burn-in..." for i=1:burnin
            mages_prop = copy(mages)
            closest_prop = copy(closest)
            sample_height_prop = copy(sample_height)

            if rand() < 0.1
                # Adjust heights
                sample_height_prop += randn(size(Height)) .* Height_sigma
                closest_prop = findclosest(sample_height_prop, model_heights)
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] = mages[chosen_point] + r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    younger_points_below = (mages_prop .< mages_prop[chosen_point]) .& (index .< chosen_point)
                    mages_prop[younger_points_below] .= mages_prop[chosen_point]
                else # if proposing decreased age
                    older_points_above = (mages_prop .> mages_prop[chosen_point]) .& (index .> chosen_point)
                    mages_prop[older_points_above] .= mages_prop[chosen_point]
                end
                # #Resolve conflicts with random linear slope
                # slope = rand() * (oldest-youngest)/npoints
                # offset = (index[chosen_point]-index) * slope
                # if r > 0 # If proposing increased age
                #     younger_points_below = (mages_prop .< mages_prop[chosen_point] + offset) .& (index .< chosen_point)
                #     mages_prop[younger_points_below] = mages_prop[chosen_point] + offset[younger_points_below]
                # else # if proposing decreased age
                #     older_points_above = (mages_prop .> mages_prop[chosen_point] + offset) .& (index .> chosen_point)
                #     mages_prop[older_points_above] = mages_prop[chosen_point] + offset[older_points_above]
                # end
            end


            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            agell_prop = bilinear_exponential_LL(mages_prop[closest_prop],p)
            heightll_prop = .- (sample_height_prop .- Height).^2 ./ (2 .* Height_sigma.^2) .- log.(sqrt.(2*pi*Height_sigma))
            diff_sign = Age_Sidedness .!= sign.(mages_prop[closest_prop] .- Age)
            ll_prop = sum(agell_prop[diff_sign]) + sum(-log.(sqrt.(2*pi*Age_sigma[.~diff_sign]))) + sum(heightll_prop)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                mages = copy(mages_prop)
                sample_height = copy(sample_height_prop)
                closest = copy(closest_prop)
                ll = copy(ll_prop)
                acceptancedist[i] = true
            end
        end

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        lldist = Array{Float64}(undef,nsteps)
        # acceptancedist = Array{Bool}(undef,nsteps)


        # Run the model
        @showprogress 5 "Collecting..." for i=1:(nsteps*sieve)
            mages_prop = copy(mages)
            closest_prop = copy(closest)
            sample_height_prop = copy(sample_height)

            if rand() < 0.1
                # Adjust heights
                sample_height_prop += randn(size(Height)) .* Height_sigma
                closest_prop = findclosest(sample_height_prop, model_heights)
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] = mages[chosen_point] + r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    younger_points_below = (mages_prop .< mages_prop[chosen_point]) .& (index .< chosen_point)
                    mages_prop[younger_points_below] .= mages_prop[chosen_point]
                else # if proposing decreased age
                    older_points_above = (mages_prop .> mages_prop[chosen_point]) .& (index .> chosen_point)
                    mages_prop[older_points_above] .= mages_prop[chosen_point]
                end
                # #Resolve conflicts with random linear slope
                # slope = rand() * (oldest-youngest)/npoints
                # offset = (index[chosen_point]-index) * slope
                # if r > 0 # If proposing increased age
                #     younger_points_below = (mages_prop .< mages_prop[chosen_point] + offset) .& (index .< chosen_point)
                #     mages_prop[younger_points_below] = mages_prop[chosen_point] + offset[younger_points_below]
                # else # if proposing decreased age
                #     older_points_above = (mages_prop .> mages_prop[chosen_point] + offset) .& (index .> chosen_point)
                #     mages_prop[older_points_above] = mages_prop[chosen_point] + offset[older_points_above]
                # end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            agell_prop = bilinear_exponential_LL(mages_prop[closest_prop],p)
            heightll_prop = .-(sample_height_prop .- Height).^2 ./ (2 .* Height_sigma.^2) .- log.(sqrt.(2*pi*Height_sigma))
            diff_sign = Age_Sidedness .!= sign.(mages_prop[closest_prop] .- Age)
            ll_prop = sum(agell_prop[diff_sign]) + sum(-log.(sqrt.(2*pi*Age_sigma[.~diff_sign]))) + sum(heightll_prop)


            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                mages = copy(mages_prop)
                sample_height = copy(sample_height_prop)
                closest = copy(closest_prop)
                ll = copy(ll_prop)
            end

            # Record sieved results
            if mod(i,sieve) == 0
                agedist[:,Int(i/sieve)] = mages
                lldist[Int(i/sieve)] = ll
            end
        end

        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist[active_height_t,:],dim=2), # Mean age
            nanstd(agedist[active_height_t,:],dim=2), # Standard deviation
            nanmedian(agedist[active_height_t,:],dim=2), # Median age
            pctile(agedist[active_height_t,:],2.5,dim=2), # 2.5th percentile
            pctile(agedist[active_height_t,:],97.5,dim=2) # 97.5th percentile
        )

        return mdl, agedist[active_height_t,:], lldist
    end

## --- Stratigraphic MCMC model with hiata, with distribution LL # # # # # # # #

    function StratMetropolisDistHiatus(smpl::StratAgeData,hiatus::HiatusData,config::StratAgeModelConfiguration)
        # Run stratigraphic MCMC model, with hiata

        print("Generating stratigraphic age-depth model...\n")

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints
        (bottom, top) = extrema(smpl.Height)
        (youngest, oldest) = extrema(smpl.Age)
        dt_dH = (oldest-youngest)/(top-bottom)
        aveuncert = mean(smpl.Age_sigma)
        p = smpl.Params

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            offset = (top-bottom)*bounding
            Age = [oldest + offset*dt_dH; smpl.Age; youngest - offset*dt_dH]
            Age_sigma = [mean(smpl.Age_sigma)/10; smpl.Age_sigma; mean(smpl.Age_sigma)/10]
            Height = [bottom-offset; smpl.Height; top+offset]
            Height_sigma = [0; smpl.Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; smpl.Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
            pl = ones(5); pl[2] = oldest + offset*dt_dH; pl[3] = mean(smpl.Age_sigma)/10
            pu = ones(5); pu[2] = youngest - offset*dt_dH; pu[3] = mean(smpl.Age_sigma)/10
            p = hcat(pl,p,pu) # Add parameters for upper and lower runaway bounds
        else
            Age = smpl.Age
            Age_sigma = smpl.Age_sigma
            Height = smpl.Height
            Height_sigma = smpl.Height_sigma .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = smpl.Age_Sidedness # Bottom is a maximum age and top is a minimum age
            model_heights = bottom:resolution:top
        end
        active_height_t = (model_heights .>= bottom) .& (model_heights .<= top)
        npoints = length(model_heights)

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        mages = a .+ b .* collect(model_heights)


        # Calculate log likelihood of initial proposal
        # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
        # proposals older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
        sample_height = Height
        closest = findclosest(sample_height,model_heights)
        agell = bilinear_exponential_LL(mages[closest],p)
        heightll = .- (sample_height .- Height).^2 ./ (2 .* Height_sigma.^2) .- log.(sqrt.(2*pi*Height_sigma))
        diff_sign = Age_Sidedness .!= sign.(mages[closest] .- Age)
        ll = sum(agell[diff_sign]) + sum(-log.(sqrt.(2*pi*Age_sigma[.~diff_sign]))) + sum(heightll)

        # Ensure there is only one effective hiatus at most for each height node
        closest_h = findclosestabove(hiatus.Height,model_heights)
        closest_h_unique = unique(closest_h)
        hiatus_Height = Array{Float64}(undef,size(closest_h_unique))
        hiatus_Duration = Array{Float64}(undef,size(closest_h_unique))
        hiatus_Duration_sigma = Array{Float64}(undef,size(closest_h_unique))
        for i=1:length(closest_h_unique)
            hiatus_Height[i] = mean(hiatus.Height[closest_h.==closest_h_unique[i]])
            hiatus_Duration[i] = sum(hiatus.Duration[closest_h.==closest_h_unique[i]])
            hiatus_Duration_sigma[i] = sqrt(sum(hiatus.Duration_sigma[closest_h.==closest_h_unique[i]].^2))
        end

        # Add log likelihood for hiatuses
        duration_prop = abs.(mages[closest_h_unique .- 1] .- mages[closest_h_unique])
        ll += sum(-max.(hiatus_Duration .- duration_prop, 0).^2 ./ (2 .* hiatus_Duration_sigma.^2)) #-log.(sqrt.(2*pi*hiatus_Duration_sigma))
        duration = copy(duration_prop)

        # Introduce variables so they will be accessible outside loop
        mages_prop = copy(mages)
        closest_prop = copy(closest)
        sample_height_prop = copy(sample_height)
        ll_prop = copy(ll)
        chosen_point = 0

        acceptancedist = fill(false,burnin)
        # Run burnin
        print("Burn-in: ", burnin, " steps\n")
        index = collect(1:npoints)
        @showprogress 5 "Burn-in..." for i=1:burnin
            mages_prop = copy(mages)
            closest_prop = copy(closest)
            sample_height_prop = copy(sample_height)

            if rand() < 0.1
                # Adjust heights
                sample_height_prop += randn(size(Height)) .* Height_sigma
                closest_prop = findclosest(sample_height_prop, model_heights)
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] = mages[chosen_point] + r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    younger_points_below = (mages_prop .< mages_prop[chosen_point]) .& (index .< chosen_point)
                    mages_prop[younger_points_below] .= mages_prop[chosen_point]
                else # if proposing decreased age
                    older_points_above = (mages_prop .> mages_prop[chosen_point]) .& (index .> chosen_point)
                    mages_prop[older_points_above] .= mages_prop[chosen_point]
                end

                # If chosen_point is a hiatus point, let there be a 20 percent chance of
                # adjusting the point below the hiatus as well
                if rand() < 0.2
                    # if hiatus_height_uncert>0
                    #     closest_h = findclosestabove(h.Height+randn(size(h.Height)).*hiatus_height_uncert,heights)
                    # end
                    if any(closest_h_unique.==chosen_point)
                        mages_prop[chosen_point-1] = mages[chosen_point-1] + r
                        #Resolve conflicts
                        if r>0
                            mages_prop[(mages_prop .< mages_prop[chosen_point-1]) .& (index .< [chosen_point-1])] .= mages_prop[chosen_point-1]
                        else
                            mages_prop[(mages_prop .> mages_prop[chosen_point-1]) .& (index .> [chosen_point-1])] .= mages_prop[chosen_point-1]
                        end
                    end
                end
            end


            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            agell_prop = bilinear_exponential_LL(mages_prop[closest_prop],p)
            heightll_prop = .- (sample_height_prop .- Height).^2 ./ (2 .* Height_sigma.^2) -log.(sqrt.(2*pi*Height_sigma))
            diff_sign = Age_Sidedness .!= sign.(mages_prop[closest_prop]-Age)
            ll_prop = sum(agell_prop[diff_sign]) + sum(-log.(sqrt.(2*pi*Age_sigma[.~diff_sign]))) + sum(heightll_prop)

            # if hiatus_height_uncert>0
            #     closest_h = findclosestabove(h.Height+randn(size(h.Height)).*hiatus_height_uncert,heights)
            # end
            duration_prop = abs.(mages_prop[closest_h_unique .- 1] .- mages_prop[closest_h_unique])
            ll_prop += sum(-max.(hiatus_Duration .- duration_prop, 0).^2 ./ (2 .* hiatus_Duration_sigma.^2)) #-log.(sqrt.(2*pi*hiatus_Duration_sigma))

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                mages = copy(mages_prop)
                sample_height = copy(sample_height_prop)
                closest = copy(closest_prop)
                ll = copy(ll_prop)
                acceptancedist[i] = true
            end
        end

        # Run Markov Chain Monte Carlo
        print("Collecting sieved stationary distribution: ", nsteps*sieve, " steps\n")
        agedist = Array{Float64}(undef,npoints,nsteps)
        hiatusdist = Array{Float64}(undef,length(duration),nsteps)
        lldist = Array{Float64}(undef,nsteps)
        # acceptancedist = Array{Bool}(undef,nsteps)


        # Run the model
        @showprogress 5 "Collecting..." for i=1:(nsteps*sieve)
            mages_prop = copy(mages)
            closest_prop = copy(closest)
            sample_height_prop = copy(sample_height)

            if rand() < 0.1
                # Adjust heights
                sample_height_prop += randn(size(Height)) .* Height_sigma
                closest_prop = findclosest(sample_height_prop, model_heights)
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * aveuncert # Generate a random adjustment
                chosen_point = ceil(Int, rand() * npoints) # Pick a point
                mages_prop[chosen_point] = mages[chosen_point] + r
                #Resolve conflicts
                if r > 0 # If proposing increased age
                    younger_points_below = (mages_prop .< mages_prop[chosen_point]) .& (index .< chosen_point)
                    mages_prop[younger_points_below] .= mages_prop[chosen_point]
                else # if proposing decreased age
                    older_points_above = (mages_prop .> mages_prop[chosen_point]) .& (index .> chosen_point)
                    mages_prop[older_points_above] .= mages_prop[chosen_point]
                end

                # If chosen_point is a hiatus point, let there be a 20 percent chance of
                # adjusting the point below the hiatus as well
                if rand() < 0.2
                    # if hiatus_height_uncert>0
                    #     closest_h = findclosestabove(h.Height+randn(size(h.Height)).*hiatus_height_uncert,heights)
                    # end
                    if any(closest_h_unique.==chosen_point)
                        mages_prop[chosen_point-1] = mages[chosen_point-1] + r
                        #Resolve conflicts
                        if r>0
                            younger_points_below = (mages_prop .< mages_prop[chosen_point-1]) .& (index .< (chosen_point-1))
                            mages_prop[younger_points_below] .= mages_prop[chosen_point-1]
                        else
                            older_points_above = (mages_prop .> mages_prop[chosen_point-1]) .& (index .> (chosen_point-1))
                            mages_prop[older_points_above] .= mages_prop[chosen_point-1]
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            # Proposals younger than age constraint are given a pass if Age_Sidedness is -1 (maximum age)
            # proposal older than age constraint are given a pass if Age_Sidedness is +1 (minimum age)
            agell_prop = bilinear_exponential_LL(mages_prop[closest_prop],p)
            heightll_prop = .- (sample_height_prop .- Height).^2 ./ (2 .* Height_sigma.^2) -log.(sqrt.(2*pi*Height_sigma))
            diff_sign = Age_Sidedness .!= sign.(mages_prop[closest_prop]-Age)
            ll_prop = sum(agell_prop[diff_sign]) + sum(-log.(sqrt.(2*pi*Age_sigma[.~diff_sign]))) + sum(heightll_prop)


            # Add log likelihood for duration
            # if hiatus_height_uncert>0
            #     closest_h = findclosestabove(h.Height+randn(size(h.Height)).*hiatus_height_uncert,heights)
            # end
            duration_prop = abs.(mages_prop[closest_h_unique .- 1] .- mages_prop[closest_h_unique])
            ll_prop += sum(-max.(hiatus_Duration .- duration_prop, 0).^2 ./ (2 .* hiatus_Duration_sigma.^2)) #-log.(sqrt.(2*pi*hiatus_Duration_sigma))


            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (ll_prop - ll)
                mages = copy(mages_prop)
                sample_height = copy(sample_height_prop)
                closest = copy(closest_prop)
                ll = copy(ll_prop)
                duration = copy(duration_prop)
            end

            # Record sieved results
            if mod(i,sieve) == 0
                agedist[:,Int(i/sieve)] = mages
                lldist[Int(i/sieve)] = ll
                hiatusdist[:,Int(i/sieve)] = duration
            end
        end

        mdl = StratAgeModel(
            model_heights[active_height_t], # Model heights
            nanmean(agedist[active_height_t,:],dim=2), # Mean age
            nanstd(agedist[active_height_t,:],dim=2), # Standard deviation
            nanmedian(agedist[active_height_t,:],dim=2), # Median age
            pctile(agedist[active_height_t,:],2.5,dim=2), # 2.5th percentile
            pctile(agedist[active_height_t,:],97.5,dim=2) # 97.5th percentile
        )

        return mdl, agedist[active_height_t,:], lldist, hiatusdist
    end
