## --- Stratigraphic MCMC model without hiatus # # # # # # # # # # # # # # # # #


    """
    ```julia
    StratMetropolis(smpl::ChronAgeData, [hiatus::HiatusData,] config::StratAgeModelConfiguration)
    StratMetropolis(smpl::GeneralAgeData, [hiatus::HiatusData,] config::StratAgeModelConfiguration)
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
    function StratMetropolis(smpl::ChronAgeData, config::StratAgeModelConfiguration, systematic=nothing)
        # Run stratigraphic MCMC model
        @info "Generating stratigraphic age-depth model..."

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints
        Age = copy(smpl.Age)::Vector{Float64}
        Age_sigma = copy(smpl.Age_sigma)::Vector{Float64}
        Height = copy(smpl.Height)::Vector{Float64}
        Height_sigma = smpl.Height_sigma::Vector{Float64} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Vector{Float64} # Bottom is a maximum age and top is a minimum age
        chronometer = smpl.Chronometer
        (bottom, top) = extrema(Height)
        model_heights = bottom:resolution:top

        aveuncert = nanmean(Age_sigma)
        absdiff = diff(sort!(Age[Age_Sidedness.==0]))
        maxdiff = isempty(absdiff) ? 0.0 : nanmaximum(absdiff)
        proposal_sigma = sqrt(aveuncert^2 + (maxdiff/10)^2)

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            (youngest, oldest) = extrema(Age)
            dt_dH = (oldest-youngest)/(top-bottom)
            offset = round((top-bottom)*bounding/resolution)*resolution
            Age = [oldest + offset*dt_dH; Age; youngest - offset*dt_dH]
            Age_sigma = [aveuncert/10; Age_sigma; aveuncert/10]
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
            chronometer = (:None, chronometer..., :None)
        end
        active_height_t = bottom .<= model_heights .<= top

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        model_ages = a .+ b .* collect(model_heights)

        # Select sidedness method
        sidedness = if smpl.Sidedness_Method === :fast || all(iszero, smpl.Age_Sidedness)
            FastSidedness(Age_Sidedness)
        else
            CDFSidedness(Age_Sidedness)
        end

        # Run the Markov chain
        ages = Normal.(Age, Age_sigma)
        agedist, lldist = stratmetropolis(Height, Height_sigma, model_heights, sidedness, ages, model_ages, proposal_sigma, burnin, nsteps, sieve, chronometer, systematic)

        # Crop the result
        agedist = agedist[active_height_t,:]
        model_heights = model_heights[active_height_t]
        mdl = StratAgeModel(model_heights, agedist)

        return mdl, agedist, lldist
    end
    function StratMetropolis(smpl::ChronAgeData, hiatus::HiatusData, config::StratAgeModelConfiguration, systematic=nothing)
        # Run stratigraphic MCMC model, with hiata
        @info "Generating stratigraphic age-depth model..."

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints. Type assertions for stability
        Age = copy(smpl.Age)::Vector{Float64}
        Age_sigma = copy(smpl.Age_sigma)::Vector{Float64}
        Height = copy(smpl.Height)::Vector{Float64}
        Height_sigma = smpl.Height_sigma::Vector{Float64} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Vector{Float64} # Bottom is a maximum age and top is a minimum age
        chronometer = smpl.Chronometer
        (bottom, top) = extrema(Height)
        model_heights = bottom:resolution:top

        aveuncert = nanmean(Age_sigma)
        absdiff = diff(sort!(Age[Age_Sidedness.==0]))
        maxdiff = isempty(absdiff) ? 0.0 : nanmaximum(absdiff)
        proposal_sigma = sqrt(aveuncert^2 + (maxdiff/10)^2)

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            (youngest, oldest) = extrema(Age)
            dt_dH = (oldest-youngest)/(top-bottom)
            offset = round((top-bottom)*bounding/resolution)*resolution
            Age = [oldest + offset*dt_dH; Age; youngest - offset*dt_dH]
            Age_sigma = [aveuncert/10; Age_sigma; aveuncert/10]
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
            chronometer = (:None, chronometer..., :None)
        end
        active_height_t = bottom .<= model_heights .<= top
        npoints = length(model_heights)

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        model_ages = a .+ b .* collect(model_heights)

        # Select sidedness method
        sidedness = if smpl.Sidedness_Method === :fast || all(iszero, smpl.Age_Sidedness)
            FastSidedness(Age_Sidedness)
        else
            CDFSidedness(Age_Sidedness)
        end
        
        # Run the Markov chain
        ages = Normal.(Age, Age_sigma)
        agedist, lldist, hiatusdist = stratmetropolis(hiatus, Height, Height_sigma, model_heights, sidedness, ages, model_ages, proposal_sigma, burnin, nsteps, sieve, chronometer, systematic)

        # Crop the result
        agedist = agedist[active_height_t,:]
        model_heights = model_heights[active_height_t]
        mdl = StratAgeModel(model_heights, agedist)

        return mdl, agedist, hiatusdist, lldist
    end
    function StratMetropolis(smpl::GeneralAgeData, config::StratAgeModelConfiguration, systematic=nothing)
        # Run stratigraphic MCMC model
        @info "Generating stratigraphic age-depth model..."

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints
        ages = unionize(smpl.Age)::Vector{<:Union{<:Distribution{Univariate, Continuous}}}
        Height = copy(smpl.Height)::Vector{Float64}
        Height_sigma = smpl.Height_sigma::Vector{Float64} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Vector{Float64} # Bottom is a maximum age and top is a minimum age
        chronometer = smpl.Chronometer
        (bottom, top) = extrema(Height)
        model_heights = bottom:resolution:top

        aveuncert = nanmean(std.(ages))
        absdiff = diff(sort!(mean.(ages[Age_Sidedness.==0])))
        maxdiff = isempty(absdiff) ? 0.0 : nanmaximum(absdiff)
        proposal_sigma = sqrt(aveuncert^2 + (maxdiff/10)^2)

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            (youngest, oldest) = extrema(mean.(ages))
            dt_dH = (oldest-youngest)/(top-bottom)
            offset = round((top-bottom)*bounding/resolution)*resolution
            ages = unionize([Normal(oldest+offset*dt_dH, aveuncert/10); 
                             ages; 
                             Normal(youngest-offset*dt_dH, aveuncert/10)])
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
            chronometer = (:None, chronometer..., :None)
        end
        active_height_t = bottom .<= model_heights .<= top

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ mean.(ages)
        model_ages = a .+ b .* collect(model_heights)

        # Select sidedness method
        sidedness = if smpl.Sidedness_Method === :fast || all(iszero, smpl.Age_Sidedness)
            FastSidedness(Age_Sidedness)
        else
            CDFSidedness(Age_Sidedness)
        end

        # Run the Markov chain
        agedist, lldist = stratmetropolis(Height, Height_sigma, model_heights, sidedness, ages, model_ages, proposal_sigma, burnin, nsteps, sieve, chronometer, systematic)

        # Crop the result
        agedist = agedist[active_height_t,:]
        model_heights = model_heights[active_height_t]
        mdl = StratAgeModel(model_heights, agedist)

        return mdl, agedist, lldist
    end
    function StratMetropolis(smpl::GeneralAgeData, hiatus::HiatusData, config::StratAgeModelConfiguration, systematic=nothing)
        # Run stratigraphic MCMC model
        @info "Generating stratigraphic age-depth model..."

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints
        ages = unionize(smpl.Age)::Vector{<:Union{<:Distribution{Univariate, Continuous}}}
        Height = copy(smpl.Height)::Vector{Float64}
        Height_sigma = smpl.Height_sigma::Vector{Float64} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Vector{Float64} # Bottom is a maximum age and top is a minimum age
        chronometer = smpl.Chronometer
        (bottom, top) = extrema(Height)
        model_heights = bottom:resolution:top

        aveuncert = nanmean(std.(ages))
        absdiff = diff(sort!(mean.(ages[Age_Sidedness.==0])))
        maxdiff = isempty(absdiff) ? 0.0 : nanmaximum(absdiff)
        proposal_sigma = sqrt(aveuncert^2 + (maxdiff/10)^2)

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            (youngest, oldest) = extrema(mean.(ages))
            dt_dH = (oldest-youngest)/(top-bottom)
            offset = round((top-bottom)*bounding/resolution)*resolution
            ages = unionize([Normal(oldest+offset*dt_dH, aveuncert/10); 
                             ages; 
                             Normal(youngest-offset*dt_dH, aveuncert/10)])
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
            chronometer = (:None, chronometer..., :None)
        end
        active_height_t = bottom .<= model_heights .<= top

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ mean.(ages)
        model_ages = a .+ b .* collect(model_heights)

        # Select sidedness method
        sidedness = if smpl.Sidedness_Method === :fast || all(iszero, smpl.Age_Sidedness)
            FastSidedness(Age_Sidedness)
        else
            CDFSidedness(Age_Sidedness)
        end

        # Run the Markov chain
        agedist, lldist, hiatusdist = stratmetropolis(hiatus, Height, Height_sigma, model_heights, sidedness, ages, model_ages, proposal_sigma, burnin, nsteps, sieve, chronometer, systematic)

        # Crop the result
        agedist = agedist[active_height_t,:]
        model_heights = model_heights[active_height_t]
        mdl = StratAgeModel(model_heights, agedist)

        return mdl, agedist, hiatusdist, lldist
    end

## --- Stratigraphic MCMC model without hiatus, with distribution LL # # # # # #

    """
    ```julia
    StratMetropolisDist(smpl::ChronAgeData, [hiatus::HiatusData,] config::StratAgeModelConfiguration)
    ```
    Runs the main Chron.jl age-depth model routine for a stratigraphic set of
    samples defined by sample heights and fitted asymmetric age distributions
    (`BilinearExponential`) in the `smpl` struct, and an age-depth model
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
    function StratMetropolisDist(smpl::ChronAgeData, config::StratAgeModelConfiguration, systematic=nothing)
        # Run stratigraphic MCMC model
        @info "Generating stratigraphic age-depth model..."

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints
        Age = copy(smpl.Age)::Vector{Float64}
        Age_sigma = copy(smpl.Age_sigma)::Vector{Float64}
        Height = copy(smpl.Height)::Vector{Float64}
        Height_sigma = smpl.Height_sigma::Vector{Float64} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Vector{Float64} # Bottom is a maximum age and top is a minimum age
        Chronometer = smpl.Chronometer
        p = copy(smpl.Params)::Matrix{Float64}
        (bottom, top) = extrema(Height)
        model_heights = bottom:resolution:top

        aveuncert = nanmean(Age_sigma)
        absdiff = diff(sort!(Age[Age_Sidedness.==0]))
        maxdiff = isempty(absdiff) ? 0.0 : nanmaximum(absdiff)
        proposal_sigma = sqrt(aveuncert^2 + (maxdiff/10)^2)

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            (youngest, oldest) = extrema(Age)
            dt_dH = (oldest-youngest)/(top-bottom)
            offset = round((top-bottom)*bounding/resolution)*resolution
            Age = [oldest + offset*dt_dH; Age; youngest - offset*dt_dH]
            Age_sigma = [aveuncert/10; Age_sigma; aveuncert/10]
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
            pl = ones(5); pl[2] = oldest + offset*dt_dH; pl[3] = aveuncert/10
            pu = ones(5); pu[2] = youngest - offset*dt_dH; pu[3] = aveuncert/10
            p = hcat(pl,p,pu) # Add parameters for upper and lower runaway bounds
            Chronometer = (:None, Chronometer..., :None)
        end
        active_height_t = bottom .<= model_heights .<= top
        npoints = length(model_heights)

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        model_ages = a .+ b .* collect(model_heights)

        # Select sidedness method
        sidedness = if smpl.Sidedness_Method === :fast || all(iszero, smpl.Age_Sidedness)
            FastSidedness(Age_Sidedness)
        else
            CDFSidedness(Age_Sidedness)
        end

        # Run the Markov chain
        ages = BilinearExponential.(eachcol(p))
        agedist, lldist = stratmetropolis(Height, Height_sigma, model_heights, sidedness, ages, model_ages, proposal_sigma, burnin, nsteps, sieve, Chronometer, systematic)

        # Crop the result
        agedist = agedist[active_height_t,:]
        model_heights = model_heights[active_height_t]
        mdl = StratAgeModel(model_heights, agedist)

        return mdl, agedist, lldist
    end
    function StratMetropolisDist(smpl::ChronAgeData, hiatus::HiatusData, config::StratAgeModelConfiguration)
        # Run stratigraphic MCMC model, with hiata
        @info "Generating stratigraphic age-depth model..."

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints
        Age = copy(smpl.Age)::Vector{Float64}
        Age_sigma = copy(smpl.Age_sigma)::Vector{Float64}
        Height = copy(smpl.Height)::Vector{Float64}
        Height_sigma = smpl.Height_sigma::Vector{Float64} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Vector{Float64} # Bottom is a maximum age and top is a minimum age
        p = copy(smpl.Params)::Matrix{Float64}
        (bottom, top) = extrema(Height)
        model_heights = bottom:resolution:top

        aveuncert = nanmean(Age_sigma)
        absdiff = diff(sort!(Age[Age_Sidedness.==0]))
        maxdiff = isempty(absdiff) ? 0.0 : nanmaximum(absdiff)
        proposal_sigma = sqrt(aveuncert^2 + (maxdiff/10)^2)

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            (youngest, oldest) = extrema(Age)
            dt_dH = (oldest-youngest)/(top-bottom)
            offset = round((top-bottom)*bounding/resolution)*resolution
            Age = [oldest + offset*dt_dH; Age; youngest - offset*dt_dH]
            Age_sigma = [aveuncert/10; Age_sigma; aveuncert/10]
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
            pl = ones(5); pl[2] = oldest + offset*dt_dH; pl[3] = aveuncert/10
            pu = ones(5); pu[2] = youngest - offset*dt_dH; pu[3] = aveuncert/10
            p = hcat(pl,p,pu) # Add parameters for upper and lower runaway bounds
        end
        active_height_t = bottom .<= model_heights .<= top
        npoints = length(model_heights)

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        model_ages = a .+ b .* collect(model_heights)

        # Select sidedness method
        sidedness = if smpl.Sidedness_Method === :fast || all(iszero, smpl.Age_Sidedness)
            FastSidedness(Age_Sidedness)
        else
            CDFSidedness(Age_Sidedness)
        end

        # Run the Markov chain
        ages = BilinearExponential.(eachcol(p))
        agedist, lldist, hiatusdist = stratmetropolis(hiatus, Height, Height_sigma, model_heights, sidedness, ages, model_ages, proposal_sigma, burnin, nsteps, sieve)

        # Crop the result
        agedist = agedist[active_height_t,:]
        model_heights = model_heights[active_height_t]
        mdl = StratAgeModel(model_heights, agedist)
        
        return mdl, agedist, hiatusdist, lldist
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
        @info "Generating stratigraphic age-depth model..."

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints
        Age = copy(smpl.Age)::Vector{Float64}
        Age_sigma = copy(smpl.Age_sigma)::Vector{Float64}
        Height = copy(smpl.Height)::Vector{Float64}
        Height_sigma = smpl.Height_sigma::Vector{Float64} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Vector{Float64} # Bottom is a maximum age and top is a minimum age
        p = copy(smpl.Params)::Matrix{Float64}
        (bottom, top) = extrema(Height)
        model_heights = bottom:resolution:top

        aveuncert = nanmean(Age_sigma)
        absdiff = diff(sort!(Age[Age_Sidedness.==0]))
        maxdiff = isempty(absdiff) ? 0.0 : nanmaximum(absdiff)
        proposal_sigma = sqrt(aveuncert^2 + (maxdiff/10)^2)

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            (youngest, oldest) = extrema(Age)
            dt_dH = (oldest-youngest)/(top-bottom)
            offset = round((top-bottom)*bounding/resolution)*resolution
            Age = [oldest + offset*dt_dH; Age; youngest - offset*dt_dH]
            Age_sigma = [aveuncert/10; Age_sigma; aveuncert/10]
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
            boundsigma = aveuncert/10
            pl = normpdf_ll.(oldest + offset*dt_dH, boundsigma, 1:size(smpl.Params,1))
            pu = normpdf_ll.(youngest - offset*dt_dH, boundsigma, 1:size(smpl.Params,1))
            p = hcat(pl,p,pu) # Add parameters for upper and lower runaway bounds
        end
        active_height_t = bottom .<= model_heights .<= top
        npoints = length(model_heights)

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        model_ages = a .+ b .* collect(model_heights)

        # Select sidedness method
        sidedness = if smpl.Sidedness_Method === :fast || all(iszero, smpl.Age_Sidedness)
            FastSidedness(Age_Sidedness)
        else
            CDFSidedness(Age_Sidedness)
        end

        # Run the Markov chain
        ages = Radiocarbon.(Age, Age_sigma, (collect(c) for c in eachcol(p)))
        agedist, lldist = stratmetropolis(Height, Height_sigma, model_heights, sidedness, ages, model_ages, proposal_sigma, burnin, nsteps, sieve)

        # Crop the result
        agedist = agedist[active_height_t,:]
        model_heights = model_heights[active_height_t]
        mdl = StratAgeModel(model_heights, agedist)

        return mdl, agedist, lldist
    end
    function StratMetropolis14C(smpl::ChronAgeData, hiatus::HiatusData, config::StratAgeModelConfiguration)
        # Run stratigraphic MCMC model, with hiata
        @info "Generating stratigraphic age-depth model..."

        # Model configuration -- read from struct
        resolution = config.resolution
        burnin = config.burnin
        nsteps = config.nsteps
        sieve = config.sieve
        bounding = config.bounding

        # Stratigraphic age constraints
        Age = copy(smpl.Age)::Vector{Float64}
        Age_sigma = copy(smpl.Age_sigma)::Vector{Float64}
        Height = copy(smpl.Height)::Vector{Float64}
        Height_sigma = smpl.Height_sigma::Vector{Float64} .+ 1E-9 # Avoid divide-by-zero issues
        Age_Sidedness = copy(smpl.Age_Sidedness)::Vector{Float64} # Bottom is a maximum age and top is a minimum age
        p = copy(smpl.Params)::Matrix{Float64}
        (bottom, top) = extrema(Height)
        model_heights = bottom:resolution:top

        aveuncert = nanmean(Age_sigma)
        absdiff = diff(sort!(Age[Age_Sidedness.==0]))
        maxdiff = isempty(absdiff) ? 0.0 : nanmaximum(absdiff)
        proposal_sigma = sqrt(aveuncert^2 + (maxdiff/10)^2)

        if bounding>0
            # If bounding is requested, add extrapolated top and bottom bounds to avoid
            # issues with the stratigraphic markov chain wandering off to +/- infinity
            (youngest, oldest) = extrema(Age)
            dt_dH = (oldest-youngest)/(top-bottom)
            offset = round((top-bottom)*bounding/resolution)*resolution
            Age = [oldest + offset*dt_dH; Age; youngest - offset*dt_dH]
            Age_sigma = [aveuncert/10; Age_sigma; aveuncert/10]
            Height = [bottom-offset; Height; top+offset]
            Height_sigma = [0; Height_sigma; 0] .+ 1E-9 # Avoid divide-by-zero issues
            Age_Sidedness = [-1.0; Age_Sidedness; 1.0;] # Bottom is a maximum age and top is a minimum age
            model_heights = (bottom-offset):resolution:(top+offset)
            boundsigma = aveuncert/10
            pl = normpdf_ll.(oldest + offset*dt_dH, boundsigma, 1:size(smpl.Params,1))
            pu = normpdf_ll.(youngest - offset*dt_dH, boundsigma, 1:size(smpl.Params,1))
            p = hcat(pl,p,pu) # Add parameters for upper and lower runaway bounds
        end
        active_height_t = bottom .<= model_heights .<= top
        npoints = length(model_heights)

        # Start with a linear fit as an initial proposal
        (a,b) = hcat(fill!(similar(Height), 1), Height) \ Age
        model_ages = a .+ b .* collect(model_heights)

        # Select sidedness method
        sidedness = if smpl.Sidedness_Method === :fast || all(iszero, smpl.Age_Sidedness)
            FastSidedness(Age_Sidedness)
        else
            CDFSidedness(Age_Sidedness)
        end

        # Run the Markov chain
        ages = Radiocarbon.(Age, Age_sigma, (collect(c) for c in eachcol(p)))
        agedist, lldist, hiatusdist = stratmetropolis(hiatus, Height, Height_sigma, model_heights, sidedness, ages, model_ages, proposal_sigma, burnin, nsteps, sieve)

        # Crop the result
        agedist = agedist[active_height_t,:]
        model_heights = model_heights[active_height_t]
        mdl = StratAgeModel(model_heights, agedist)

        return mdl, agedist, hiatusdist, lldist
    end

## --- # Internals of the Markov chain

    function stratmetropolis(Height, Height_sigma, model_heights::AbstractRange, Age_Sidedness, ages, model_ages, proposal_sigma, burnin::Integer, nsteps::Integer, sieve::Integer, Chronometer=nothing, systematic=nothing)
        resolution = step(model_heights)
        npoints = length(model_heights)

        # Calculate log likelihood of initial proposal

        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_model_ages = model_ages[closest]
        ll = strat_ll(closest_model_ages, ages, Age_Sidedness)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

        # Preallocate variables for MCMC proposals
        llₚ = ll
        agesₚ = copy(ages)
        model_agesₚ = copy(model_ages)
        closestₚ = copy(closest)
        sample_heightₚ = copy(sample_height)
        closest_model_agesₚ = copy(closest_model_ages)

        # Run burnin
        # acceptancedist = fill(false,burnin)
        @info "Burn-in: $burnin steps"
        pgrs = Progress(burnin, desc="Burn-in...")
        pgrs_interval = ceil(Int,sqrt(burnin))
        for n=1:burnin
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)
            isnothing(systematic) || copyto!(agesₚ, ages)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i ∈ eachindex(sample_heightₚ, closestₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - first(model_heights))/resolution)+1
                    closestₚ[i] = max(min(closestₚ[i], lastindex(model_agesₚ)), firstindex(model_agesₚ))
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * proposal_sigma # Generate a random adjustment
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
            adjust!(agesₚ, Chronometer, systematic)
            @inbounds for i ∈ eachindex(closest_model_agesₚ, closestₚ)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
            end
            llₚ = strat_ll(closest_model_agesₚ, agesₚ, Age_Sidedness)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(sample_height, sample_heightₚ)
                # acceptancedist[i] = true
            end
            isnan(ll) && error("Log likelihood is `NaN`")

            # Update progress meter every `pgrs_interval` steps
            mod(n,pgrs_interval)==0 && update!(pgrs, n)
        end
        update!(pgrs, burnin) # Finalize

        # Run Markov Chain Monte Carlo
        @info "Collecting sieved stationary distribution: $(nsteps*sieve) steps"
        agedist = fill(NaN, npoints,nsteps)
        lldist = fill(NaN, nsteps)

        # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        pgrs_interval = ceil(Int,sqrt(nsteps*sieve))
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)
            isnothing(systematic) || copyto!(agesₚ, ages)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i ∈ eachindex(sample_heightₚ, closestₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - first(model_heights))/resolution)+1
                    closestₚ[i] = max(min(closestₚ[i], lastindex(model_agesₚ)), firstindex(model_agesₚ))
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * proposal_sigma # Generate a random adjustment
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
            adjust!(agesₚ, Chronometer, systematic)
            @inbounds for i ∈ eachindex(closest_model_agesₚ, closestₚ)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
            end
            llₚ = strat_ll(closest_model_agesₚ, agesₚ, Age_Sidedness)
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

        return agedist, lldist
    end

    function stratmetropolis(hiatus::HiatusData, Height, Height_sigma, model_heights::AbstractRange, Age_Sidedness, ages, model_ages, proposal_sigma, burnin::Integer, nsteps::Integer, sieve::Integer, Chronometer=nothing, systematic=nothing)
        resolution = step(model_heights)
        npoints = length(model_heights)

        # Calculate log likelihood of initial proposal
        sample_height = copy(Height)
        closest = findclosest(sample_height, model_heights)
        closest_model_ages = model_ages[closest]
        ll = strat_ll(closest_model_ages, ages, Age_Sidedness)
        ll += normpdf_ll(Height, Height_sigma, sample_height)

        # Ensure there is only one effective hiatus at most for each height node
        closest_hiatus = findclosestabove((hiatus.Height::Vector{Float64}),model_heights)
        closest_hiatus_unique = unique(closest_hiatus)
        Hiatus_height = Array{Float64}(undef,size(closest_hiatus_unique))
        Hiatus_duration = Array{Float64}(undef,size(closest_hiatus_unique))
        Hiatus_duration_sigma = Array{Float64}(undef,size(closest_hiatus_unique))
        for i ∈ eachindex(closest_hiatus_unique)
            t = closest_hiatus.==closest_hiatus_unique[i]
            Hiatus_height[i] = mean((hiatus.Height::Vector{Float64})[t])
            Hiatus_duration[i] = sum((hiatus.Duration::Vector{Float64})[t])
            Hiatus_duration_sigma[i] = sqrt(sum((hiatus.Duration_sigma::Vector{Float64})[t].^2))
        end

        # Add log likelihood for hiatus duration
        duration = @. min(model_ages[closest_hiatus_unique - 1] - model_ages[closest_hiatus_unique], Hiatus_duration)
        ll += normpdf_ll(Hiatus_duration, Hiatus_duration_sigma, duration)

        # Preallocate variables for MCMC proposals
        llₚ=ll
        chosen_point=0
        agesₚ = copy(ages)
        model_agesₚ = copy(model_ages)
        closestₚ = copy(closest)
        durationₚ = copy(duration)
        sample_heightₚ = copy(sample_height)
        closest_model_agesₚ = copy(closest_model_ages)

        # Run burnin
        # acceptancedist = fill(false,burnin)
        @info "Burn-in: $burnin steps"
        pgrs = Progress(burnin, desc="Burn-in...")
        pgrs_interval = ceil(Int,sqrt(burnin))
        for n=1:burnin
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)
            isnothing(systematic) || copyto!(agesₚ, ages)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i ∈ eachindex(sample_heightₚ, closestₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - first(model_heights))/resolution)+1
                    closestₚ[i] = max(min(closestₚ[i], lastindex(model_agesₚ)), firstindex(model_agesₚ))
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * proposal_sigma # Generate a random adjustment
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

                # Move discontinuities to the nearest hiatus, if one not already present
                for _ in 1:(length(closest_hiatus_unique)÷2)
                    h = rand(closest_hiatus_unique)
                    if model_agesₚ[h-1] == model_agesₚ[h]
                        iₙ = findclosestunequal(model_agesₚ, h)
                        if iₙ < h
                            model_agesₚ[iₙ:h-1] .= model_agesₚ[iₙ]
                        elseif iₙ > h
                            model_agesₚ[h:iₙ] .= model_agesₚ[iₙ]
                        end
                    end
                end
            end


            # Calculate log likelihood of proposal
            adjust!(agesₚ, Chronometer, systematic)
            @inbounds for i ∈ eachindex(closest_model_agesₚ, closestₚ)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
            end
            llₚ = strat_ll(closest_model_agesₚ, agesₚ, Age_Sidedness)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Add log likelihood for hiatus duration
            @. durationₚ = model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique]
            llₚ += normcdf_ll(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

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
        @info "Collecting sieved stationary distribution: $(nsteps*sieve) steps"
        agedist = fill(NaN, npoints,nsteps)
        lldist = fill(NaN, nsteps)
        hiatusdist = fill(NaN, length(duration),nsteps)

        # Run the model
        pgrs = Progress(nsteps*sieve, desc="Collecting...")
        pgrs_interval = ceil(Int,sqrt(nsteps*sieve))
        for n=1:(nsteps*sieve)
            # Prepare proposal
            copyto!(model_agesₚ, model_ages)
            copyto!(closestₚ, closest)
            copyto!(sample_heightₚ, sample_height)
            isnothing(systematic) || copyto!(agesₚ, ages)

            if rand() < 0.1
                # Adjust heights
                @inbounds for i ∈ eachindex(sample_heightₚ, closestₚ)
                    sample_heightₚ[i] += randn() * Height_sigma[i]
                    closestₚ[i] = round(Int,(sample_heightₚ[i] - first(model_heights))/resolution)+1
                    closestₚ[i] = max(min(closestₚ[i], lastindex(model_agesₚ)), firstindex(model_agesₚ))
                end
            else
                # Adjust one point at a time then resolve conflicts
                r = randn() * proposal_sigma # Generate a random adjustment
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

                # Move discontinuities to the nearest hiatus, if one not already present
                for _ in 1:(length(closest_hiatus_unique)÷2)
                    h = rand(closest_hiatus_unique)
                    if model_agesₚ[h-1] == model_agesₚ[h]
                        iₙ = findclosestunequal(model_agesₚ, h)
                        if iₙ < h
                            model_agesₚ[iₙ:h-1] .= model_agesₚ[iₙ]
                        elseif iₙ > h
                            model_agesₚ[h:iₙ] .= model_agesₚ[iₙ]
                        end
                    end
                end
            end

            # Calculate log likelihood of proposal
            adjust!(agesₚ, Chronometer, systematic)
            @inbounds for i ∈ eachindex(closest_model_agesₚ, closestₚ)
                closest_model_agesₚ[i] = model_agesₚ[closestₚ[i]]
            end
            llₚ = strat_ll(closest_model_agesₚ, agesₚ, Age_Sidedness)
            llₚ += normpdf_ll(Height, Height_sigma, sample_heightₚ)

            # Add log likelihood for hiatus duration
            @. durationₚ = model_agesₚ[closest_hiatus_unique - 1] - model_agesₚ[closest_hiatus_unique]
            llₚ += normcdf_ll(Hiatus_duration, Hiatus_duration_sigma, durationₚ)

            # Accept or reject proposal based on likelihood
            if log(rand(Float64)) < (llₚ - ll)
                ll = llₚ
                copyto!(model_ages, model_agesₚ)
                copyto!(closest, closestₚ)
                copyto!(duration, durationₚ)
                copyto!(sample_height, sample_heightₚ)
            end
            isnan(ll) && error("Log likelihood is `NaN`")

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
        return agedist, lldist, hiatusdist
    end