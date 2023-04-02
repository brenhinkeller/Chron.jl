## --- Bootstrap prior distribution shape

    """
    ```julia
    BootstrapCrystDistributionKDE(smpl::ChronAgeData; cutoff=-0.05, [tpbloss=0])
    ```
    Bootstrap an estimate of the pre-eruptive (or pre-depositional) mineral
    crystallization distribution shape from a Chron.ChronAgeData object containing
    data for several samples, using a kernel density estimate of stacked sample data.

    If the samples provided as csv files in `smpl.Path` take the five-column form
    of U-Pb isotopic data files, the ages and uncertainties will be those of upper
    intercepts given a time of Pb-loss optionally specified as `tpbloss`.

    Uncertainties will be treated as one or two-sigma absolute based on the value
    of `smpl.inputSigmaLevel`.

    ### Examples
    ```julia
    BootstrappedDistribution = BootstrapCrystDistributionKDE(smpl)
    ```
    """
    function BootstrapCrystDistributionKDE(smpl::ChronAgeData; cutoff::Number=-0.05, tpbloss::Number=0)
        # Extact variables froms struct
        Name = collect(smpl.Name)::Vector{String}
        Path = smpl.Path::String
        DistType = smpl.Age_DistType::Vector{Float64}
        σstr = "$(smpl.inputSigmaLevel)-sigma 𝑎𝑏𝑠𝑜𝑙𝑢𝑡𝑒"

        # Load all data points and scale from 0 to 1
        allscaled = Float64[]
        for i ∈ eachindex(Name)
            if DistType[i]==0
                # Read data for each sample from file
                filepath = joinpath(Path, Name[i]*".csv")
                data = readclean(filepath, ',', Float64)::Matrix{Float64}

                if size(data, 2) == 5
                    @info "Interpreting the five columns of $(Name[i]).csv as:\n | ²⁰⁷Pb/²³⁵U | $σstr | ²⁰⁶Pb/²³⁸U | $σstr | correlation coefficient |"
                    data[:,2]./=smpl.inputSigmaLevel
                    data[:,4]./=smpl.inputSigmaLevel
                    analyses = UPbAnalysis.(eachcol(data)...,)
                    uis = upperintercept.(tpbloss, analyses)
                    μ, σ = Isoplot.val.(uis), Isoplot.err.(uis)
                else
                    @info "Interpreting first two columns of $(Name[i]).csv as \n | Age | Age $σstr |"
                    data[:,2]./=smpl.inputSigmaLevel
                    μ, σ = view(data, :, 1), view(data, :, 2)
                end

                # Maximum extent of expected analytical tail (beyond eruption/deposition)
                maxtaillength = nanmean(σ) .* norm_quantile(1 - 1/(1+countnotnans(μ)))
                included = (μ .- nanminimum(μ)) .>= maxtaillength
                included .|= μ .> nanmedian(μ) # Don't exclude more than half (could only happen in underdispersed datasets)
                included .&= .!isnan.(μ) # Exclude NaNs

                # Include and scale only those data not within the expected analytical tail
                if sum(included)>0
                    μₜ = data[included,1]
                    scaled = μₜ .- minimum(μₜ)
                    if maximum(scaled) > 0
                        scaled ./= maximum(scaled)
                    end
                    append!(allscaled, scaled)
                end
            end
        end

        # Calculate kernel density estimate, truncated at 0
        kd = kde(allscaled,npoints=2^7)
        t = kd.x .> cutoff # Ensure sharp cutoff at eruption / deposition
        return kd.density[t]
    end


    """
    ```julia
    BootstrapCrystDistributionKDE(data::AbstractArray, [sigma::AbstractArray]; cutoff=-0.05)
    ```
    Bootstrap an estimate of the pre-eruptive (or pre-depositional) mineral
    crystallization distribution shape from a 1- or 2-d array of sample ages
    (one row per sample, one column per datum, padded with NaNs as needed) and
    an equivalent-size array of one-sigma uncertainties, using a kernel density
    estimate of stacked sample data.

    ### Examples
    ```julia
    # Bootstrap crystallization distribution for a synthetic dataset with ages
    # [1,2,3,...10] Ma and uncertainties of 1 Ma each
    BootstrappedDistribution = BootstrapCrystDistributionKDE(1:10, ones(10))
    ```
    """
    function BootstrapCrystDistributionKDE(data::AbstractArray{T}; cutoff::Number=-0.05) where {T<:Number}
        # Load all data points and scale from 0 to 1
        allscaled = Vector{float(T)}()
        for i=1:size(data,2)
            scaled = data[:,i] .- minimum(data[:,i])
            if maximum(scaled) > 0
                scaled = scaled ./ maximum(scaled)
            end
            append!(allscaled, scaled)
        end

        # Calculate kernel density estimate, truncated at 0
        kd = kde(allscaled,npoints=2^7)
        t = kd.x .> cutoff
        return kd.density[t]
    end

    function BootstrapCrystDistributionKDE(data::AbstractArray{T}, sigma::AbstractArray{<:Number}; cutoff::Number=-0.05) where {T<:Number}
        # Array to hold stacked, scaled data
        allscaled = Vector{float(T)}()

        # For each row of data
        for i=1:size(data,2)
            μ, σ = data[:,i], sigma[:,i]

            # Maximum extent of expected analytical tail (beyond eruption/deposition/cutoff)
            maxtaillength = nanmean(σ) .* norm_quantile(1 - 1/(1+countnotnans(μ)))
            included = (μ .- nanminimum(μ)) .>= maxtaillength
            included .|= μ .> nanmedian(μ) # Don't exclude more than half (could only happen in underdispersed datasets)
            included .&= .!isnan.(μ) # Exclude NaNs

            # Include and scale only those data not within the expected analytical tail
            if sum(included)>0
                scaled = data[included,i] .- minimum(data[included,i])
                if maximum(scaled) > 0
                    scaled = scaled ./ maximum(scaled)
                end
                append!(allscaled, scaled)
            end
        end

        # Calculate kernel density estimate, truncated at 0
        kd = kde(allscaled,npoints=2^7)
        t = kd.x .> cutoff
        return kd.density[t]
    end

## --- Some useful distributions

    # Mostly outsourced to Isoplot.jl

    # Alias for distribution provided by Isoplot.jl
    const ArClosureDistribtuion = ExponentialDistribution

    # Superseded by ArClosureDistribution
    const EllisDistribution =
    [6.80942, 5.36792, 4.45867, 3.83457, 3.28267, 2.77244, 2.33403, 1.98717, 1.72219, 1.52427, 1.37856, 1.27023, 1.18445,1.10697, 1.03176, 0.958823, 0.888329, 0.820435, 0.755302, 0.693089, 0.633956, 0.578063, 0.525569, 0.476635, 0.431421,0.390065, 0.35249, 0.31849, 0.287855, 0.260378, 0.235851, 0.214064, 0.194811, 0.177882, 0.163069, 0.150163, 0.138957,0.129243, 0.120811, 0.113453, 0.106962, 0.101129, 0.095775, 0.0908504, 0.0863393, 0.0822255, 0.0784932, 0.0751262, 0.0721086, 0.0694244, 0.0670576]

    const TruncatedNormalDistribution =
    [1.15708, 1.2038, 1.25037, 1.29662, 1.34239, 1.38751, 1.43181, 1.47511, 1.51724, 1.55802, 1.5973, 1.63489, 1.67064, 1.70439, 1.73598, 1.76527, 1.79213, 1.81643, 1.83805, 1.8569, 1.87289, 1.88593, 1.89596, 1.90294, 1.90682, 1.9076, 1.90527, 1.89983, 1.89133, 1.87978, 1.86526, 1.84783, 1.82758, 1.80461, 1.77901, 1.75092, 1.72046, 1.68777, 1.65301, 1.61632, 1.57786, 1.53781, 1.49633, 1.45359, 1.40977, 1.36504, 1.31958, 1.27355, 1.22712, 1.18045, 1.13371, 1.08704, 1.04059, 0.9945, 0.948905, 0.903921, 0.859665, 0.816241, 0.773749, 0.732271, 0.691887, 0.652662, 0.614655, 0.577921, 0.542496, 0.508411, 0.475689, 0.444348, 0.414395, 0.38583, 0.358649, 0.332839, 0.308382, 0.285256, 0.263435, 0.242885, 0.223573, 0.205461, 0.188509, 0.172673, 0.157909, 0.144172, 0.131416, 0.119592, 0.108655, 0.0985572, 0.0892521, 0.0806935, 0.0728368, 0.0656377, 0.0590535, 0.0530436, 0.0475671, 0.0425867, 0.0380655, 0.0339688, 0.0302636, 0.0269185, 0.0239041, 0.0211926]


    ## --- End of File
