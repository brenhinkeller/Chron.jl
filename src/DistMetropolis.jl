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

    #[1.0, 1.0] -- maximum entropy:
    const UniformDistribution = ones(100)

    #[2.0,1.0,0.0]
    const TriangularDistribution = collect(range(2,0,length=100))

    # Best for survivorship processes, apparently including Ar-Ar partial degassing / entrainment / inheritance
    const ExponentialDistribution = [7.94734, 7.38851, 6.86897, 6.38597, 5.93693, 5.51946, 5.13135, 4.77053, 4.43508, 4.12322, 3.83329, 3.56375, 3.31315, 3.08018, 2.8636, 2.66224, 2.47504, 2.301, 2.1392, 1.98878, 1.84893, 1.71892, 1.59805, 1.48568, 1.38122, 1.28409, 1.1938, 1.10986, 1.03181, 0.959261, 0.891808, 0.829099, 0.7708, 0.7166, 0.666211, 0.619365, 0.575813, 0.535324, 0.497681, 0.462686, 0.430152, 0.399905, 0.371785, 0.345642, 0.321338, 0.298742, 0.277736, 0.258206, 0.24005, 0.22317, 0.207478, 0.192889, 0.179325, 0.166716, 0.154993, 0.144094, 0.133962, 0.124542, 0.115785, 0.107643, 0.100074, 0.0930371, 0.0864951, 0.080413, 0.0747586, 0.0695018, 0.0646147, 0.0600712, 0.0558472, 0.0519202, 0.0482693, 0.0448752, 0.0417197, 0.0387861, 0.0360588, 0.0335233, 0.031166, 0.0289745, 0.0269371, 0.025043, 0.023282, 0.0216449, 0.0201229, 0.0187079, 0.0173925, 0.0161695, 0.0150325, 0.0139755, 0.0129928, 0.0120791, 0.0112298, 0.0104401, 0.00970602, 0.00902352, 0.00838902, 0.00779913, 0.00725072, 0.00674087, 0.00626688, 0.00582621, 0.00541653, 0.00503566, 0.00468156, 0.00435237, 0.00404633, 0.0037618, 0.00349728, 0.00325137, 0.00302274, 0.00281019]
    const ArClosureDistribtuion = ExponentialDistribution

    const HalfNormalDistribution =
    [2.65307, 2.6516, 2.64718, 2.63984, 2.6296, 2.61648, 2.60054, 2.58182, 2.5604, 2.53633, 2.5097, 2.48059, 2.4491, 2.41532, 2.37936, 2.34133, 2.30135, 2.25954, 2.21603, 2.17095, 2.12442, 2.07657, 2.02755, 1.97749, 1.92653, 1.87479, 1.82242, 1.76954, 1.71629, 1.66279, 1.60917, 1.55555, 1.50205, 1.44878, 1.39584, 1.34335, 1.29139, 1.24006, 1.18946, 1.13965, 1.09071, 1.04272, 0.995728, 0.9498, 0.904985, 0.861327, 0.818864, 0.777631, 0.737653, 0.698954, 0.66155, 0.625452, 0.590667, 0.557197, 0.52504, 0.494189, 0.464635, 0.436363, 0.409356, 0.383594, 0.359054, 0.335711, 0.313537, 0.292503,0.272576, 0.253725, 0.235916, 0.219112, 0.20328, 0.188382, 0.174383, 0.161244, 0.14893, 0.137403, 0.126628, 0.116568,0.107188, 0.0984536, 0.0903304, 0.0827854, 0.0757864, 0.069302, 0.0633021, 0.0577574, 0.0526399, 0.0479225, 0.0435794, 0.0395859, 0.0359185, 0.0325546, 0.029473, 0.0266534, 0.0240769, 0.0217252, 0.0195815, 0.0176297, 0.0158548, 0.0142428, 0.0127805, 0.0114555, 0.0102566]

    const TruncatedNormalDistribution =
    [1.15708, 1.2038, 1.25037, 1.29662, 1.34239, 1.38751, 1.43181, 1.47511, 1.51724, 1.55802, 1.5973, 1.63489, 1.67064, 1.70439, 1.73598, 1.76527, 1.79213, 1.81643, 1.83805, 1.8569, 1.87289, 1.88593, 1.89596, 1.90294, 1.90682, 1.9076, 1.90527, 1.89983, 1.89133, 1.87978, 1.86526, 1.84783, 1.82758, 1.80461, 1.77901, 1.75092, 1.72046, 1.68777, 1.65301, 1.61632, 1.57786, 1.53781, 1.49633, 1.45359, 1.40977, 1.36504, 1.31958, 1.27355, 1.22712, 1.18045, 1.13371, 1.08704, 1.04059, 0.9945, 0.948905, 0.903921, 0.859665, 0.816241, 0.773749, 0.732271, 0.691887, 0.652662, 0.614655, 0.577921, 0.542496, 0.508411, 0.475689, 0.444348, 0.414395, 0.38583, 0.358649, 0.332839, 0.308382, 0.285256, 0.263435, 0.242885, 0.223573, 0.205461, 0.188509, 0.172673, 0.157909, 0.144172, 0.131416, 0.119592, 0.108655, 0.0985572, 0.0892521, 0.0806935, 0.0728368, 0.0656377, 0.0590535, 0.0530436, 0.0475671, 0.0425867, 0.0380655, 0.0339688, 0.0302636, 0.0269185, 0.0239041, 0.0211926]

    const MeltsZirconDistribution =
    [0.282361, 0.28919, 0.296019, 0.302849, 0.30968, 0.316567, 0.323614, 0.33064, 0.337727, 0.344848, 0.352146, 0.359642,0.367482, 0.375622, 0.384052, 0.392828, 0.401968, 0.411518, 0.421442, 0.43171, 0.44213, 0.45295, 0.464036, 0.47539, 0.486938, 0.498644, 0.51075, 0.523026, 0.535688, 0.548764, 0.562124, 0.575927, 0.590363, 0.604879, 0.620415, 0.636022, 0.652333, 0.669112, 0.686441, 0.704341, 0.72283, 0.742036, 0.761964, 0.782541, 0.803718, 0.825707, 0.848386, 0.871895,0.896139, 0.920462, 0.946071, 0.972964, 0.999905, 1.02776, 1.05664, 1.08637, 1.11731, 1.14919, 1.18202, 1.21582, 1.24956, 1.28342, 1.31828, 1.35427, 1.39153, 1.43006, 1.46879, 1.50812, 1.5477, 1.58888, 1.63149, 1.6748, 1.71724, 1.76126, 1.80668, 1.85101, 1.89546, 1.94144, 1.98379, 2.02785, 2.06738, 2.10669, 2.1377, 2.16306, 2.17843, 2.17924, 2.16073, 2.11744, 2.04444, 1.93323, 1.7923, 1.62527, 1.44425, 1.25401, 1.0528, 0.843628, 0.632687, 0.421876, 0.211064, 0.000252985]

    const MeltsVolcanicZirconDistribution =
    [0.54933, 0.556409, 0.563488, 0.570567, 0.577653, 0.584759, 0.591912, 0.599251, 0.606793, 0.614519, 0.622425, 0.630421, 0.63852, 0.646681, 0.654972, 0.663533, 0.672274, 0.681233, 0.690399, 0.699787, 0.709334, 0.719174, 0.729157, 0.739423, 0.749935, 0.760644, 0.771726, 0.782974, 0.794507, 0.806296, 0.818297, 0.830517, 0.842957, 0.855411, 0.866744, 0.878127, 0.889792, 0.901792, 0.914121, 0.926689, 0.939557, 0.952834, 0.966425, 0.980333, 0.994521, 1.00914, 1.02403, 1.03928, 1.05487, 1.0705, 1.08587, 1.10097, 1.11608, 1.13153, 1.1474, 1.16353, 1.18025, 1.19743, 1.21504, 1.23312, 1.25034, 1.26711, 1.28441, 1.30212, 1.32024, 1.33892, 1.35769, 1.37491, 1.3923, 1.41046, 1.42924, 1.44775, 1.46432, 1.48171, 1.49969, 1.51516, 1.53001, 1.54571, 1.5566, 1.56814, 1.57522, 1.58168, 1.58206, 1.57869, 1.56907, 1.55064, 1.51982, 1.4737, 1.40944, 1.32047, 1.21218, 1.09157, 0.965488, 0.834108, 0.697552, 0.558304, 0.418827, 0.279262, 0.139695, 0.000127237]

    # Superseded by ArClosureDistribution
    const EllisDistribution =
    [6.80942, 5.36792, 4.45867, 3.83457, 3.28267, 2.77244, 2.33403, 1.98717, 1.72219, 1.52427, 1.37856, 1.27023, 1.18445,1.10697, 1.03176, 0.958823, 0.888329, 0.820435, 0.755302, 0.693089, 0.633956, 0.578063, 0.525569, 0.476635, 0.431421,0.390065, 0.35249, 0.31849, 0.287855, 0.260378, 0.235851, 0.214064, 0.194811, 0.177882, 0.163069, 0.150163, 0.138957,0.129243, 0.120811, 0.113453, 0.106962, 0.101129, 0.095775, 0.0908504, 0.0863393, 0.0822255, 0.0784932, 0.0751262, 0.0721086, 0.0694244, 0.0670576]

    ## --- End of File
