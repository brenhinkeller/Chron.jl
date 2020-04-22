## -- Functions used for MCMC estimation of mineral eruption / deposition ages from observed crystallization age spectra


    # Return the log-likelihood of drawing the dataset 'data' from the distribution 'dist'
    function check_dist_ll(dist::Array{Float64}, data::Array{Float64}, uncert::Array{Float64}, tmin::Float64, tmax::Float64)
        # Define some frequently used variables
        loglikelihood = 0.0
        datarows = length(data)
        distrows = length(dist)
        dist_xscale = distrows-1.00000000001
        dist_yave = mean(dist)
        dt = abs(tmax-tmin)
        # Cycle through each datum in data array
        for j=1:datarows
            # If possible, prevent aliasing problems by interpolation
            if (uncert[j] < dt/dist_xscale) && data[j] > tmin && data[j] < tmax
                # Find (double) index
                ix = (data[j] - tmin) / dt * dist_xscale + 1
                # Interpolate corresponding distribution value
                f = floor(Int,ix)
                likelihood = (dist[f+1]*(ix-f) + dist[f]*(1-(ix-f))) / (dt * dist_yave)
                # Otherwise, sum contributions from Gaussians at each point in distribution
            else
                likelihood = 0
                for i=1:distrows
                    distx = tmin + dt*(i-1)/dist_xscale # time-position of distribution point
                    # Likelihood curve follows a Gaussian PDF. Note: dt cancels
                    likelihood += dist[i] / (dist_yave * distrows * uncert[j] * sqrt(2*pi)) *
                    exp( - (distx-data[j])*(distx-data[j]) / (2*uncert[j]*uncert[j]) )
                end
            end
            loglikelihood += log(likelihood)
        end
        return loglikelihood
    end

    # Return the log-likelihood of a proposed crystallization distribution, with adjustments to prevent runaway at low N
    function check_cryst_ll(dist::Array{Float64}, data::Array{Float64}, uncert::Array{Float64}, tmin::Float64, tmax::Float64)
        # Define some frequently used variables
        loglikelihood = 0.0
        datarows = length(data)
        distrows = length(dist)
        dist_xscale = distrows - 1.00000000001
        dist_yave = mean(dist)
        dt = abs(tmax-tmin)
        # Cycle through each datum in data array
        for j=1:datarows
            # If possible, prevent aliasing problems by interpolation
            if (uncert[j] < dt/dist_xscale) && data[j] > tmin && data[j] < tmax
                # Find (double) index
                ix = (data[j] - tmin) / dt * dist_xscale + 1
                # Interpolate corresponding distribution value
                f = floor(Int,ix)
                likelihood = (dist[f+1]*(ix-f) + dist[f]*(1-(ix-f))) / (dt * dist_yave)
                # Otherwise, sum contributions from Gaussians at each point in distribution
            else
                likelihood = 0
                for i=1:distrows
                    distx = tmin + dt*(i-1)/dist_xscale # time-position of distribution point
                    # Likelihood curve follows a Gaussian PDF. Note: dt cancels
                    likelihood += dist[i] / (dist_yave * distrows * uncert[j] * sqrt(2*pi)) *
                    exp( - (distx-data[j])*(distx-data[j]) / (2*uncert[j]*uncert[j]) )
                end
            end
            loglikelihood += log(likelihood)
        end
        # Calculate a weighted mean and examine our MSWD
        wm, wsigma, mswd = awmean(data, uncert)
        if datarows == 1 || mswd < 1
            Zf = 1
        elseif mswd*sqrt(datarows) > 1000
            Zf = 0
        else
            f = datarows-1.0
            # Height of MSWD distribution relative to height at MSWD = 1 (see Wendt and Carl, 1991, Chemical geology)
            Zf = exp((f/2-1)*log(mswd) - f/2*(mswd-1))
        end
        # To prevent instability / runaway of the MCMC for small datasets (low N),
        # favor the weighted mean interpretation at high Zf (MSWD close to 1) and
        # the youngest-zircon interpretation at low Zf (MSWD far from one). The
        # penalty factors used here are determined by training against synthetic datasets.
        return loglikelihood - (2/log(1+datarows)) * (                  # Scaling factor that decreases with log number of data points (i.e., no penalty at high N)
        log((abs(tmin - wm)+wsigma)/wsigma)*Zf +                        # Penalty for proposing tmin too far from the weighted mean at low MSWD (High Zf)
        log((abs(tmax - wm)+wsigma)/wsigma)*Zf +                        # Penalty for proposing tmax too far from the weighted mean at low MSWD (High Zf)
        log((abs(tmin - data[1])+uncert[1])/uncert[1])*(1-Zf) +         # Penalty for proposing tmin too far from youngest zircon at high MSWD (low Zf)
        log((abs(tmax - data[end])+uncert[end])/uncert[end])*(1-Zf) );  # Penalty for proposing tmax too far from oldest zircon at high MSWD (low Zf)
    end


    # Run a Metropolis sampler to estimate the extrema of a finite-range distribution from samples drawn
    # from that distribution -- e.g., estimate zircon saturation and eruption ages from a distribution of
    # zircon crystallization ages.
    function metropolis_minmax_cryst(nsteps::Int,dist::Array{Float64},data::Array{Float64},uncert::Array{Float64})
        # standard deviation of the proposal function is stepfactor * last step; this is tuned to optimize accetance probability at 50%
        stepfactor = 2.9
        # Sort the data array from youngest to oldest
        sI = sortperm(data)
        data = data[sI] # Sort data
        uncert = uncert[sI] # Sort uncertainty
        # These quantities will be used more than once
        datarows = length(data)
        tmin_obs = minimum(data)
        tmax_obs = maximum(data)
        # Step sigma for Gaussian proposal distributions
        dt = tmax_obs - tmin_obs + uncert[1] + uncert[end]
        tmin_step = dt / datarows
        tmax_step = dt / datarows
        # Use oldest and youngest zircons for initial proposal
        tmin = tmin_obs - uncert[1]
        tmax = tmax_obs + uncert[end]
        tmin_proposed = tmin
        tmax_proposed = tmax
        # Log likelihood of initial proposal
        ll =  check_cryst_ll(dist, data, uncert, tmin, tmax)
        ll_proposed = ll
        # Allocate ouput arrays
        tminDist = Array{Float64}(undef,nsteps)
        tmaxDist = Array{Float64}(undef,nsteps)
        llDist = Array{Float64}(undef,nsteps)
        acceptanceDist = falses(nsteps)
        # Step through each of the N steps in the Markov chain
        for i=1:nsteps
            tmin_proposed = copy(tmin)
            tmax_proposed = copy(tmax)
            # Adjust either upper or lower bound
            if rand()<0.5
                tmin_proposed += tmin_step*randn()
            else
                tmax_proposed += tmax_step*randn()
            end
            # Flip bounds if reversed
            if (tmin_proposed>tmax_proposed)
                r = tmin_proposed
                tmin_proposed = tmax_proposed
                tmax_proposed = r
            end
            # Calculate log likelihood for new proposal
            ll_proposed =  check_cryst_ll(dist, data, uncert, tmin_proposed, tmax_proposed)
            # Decide to accept or reject the proposal
            if log(rand(Float64)) < (ll_proposed-ll)
                if tmin_proposed != tmin
                    tmin_step = abs(tmin_proposed-tmin)*stepfactor
                end
                if tmax_proposed != tmax
                    tmax_step = abs(tmax_proposed-tmax)*stepfactor
                end

                tmin = copy(tmin_proposed)
                tmax = copy(tmax_proposed)
                ll = copy(ll_proposed)
                acceptanceDist[i]=true
            end
            tminDist[i] = tmin
            tmaxDist[i] = tmax
            llDist[i] = ll
        end
        return tminDist, tmaxDist, llDist, acceptanceDist
    end


## --- Some useful distributions

    UniformDistribution = ones(100) #[1.0, 1.0]
    TriangularDistribution = collect(linsp(2,0,100)) #[2.0,1.0,0.0]
    TruncatedNormalDistribution =
    [1.15708, 1.2038, 1.25037, 1.29662, 1.34239, 1.38751, 1.43181, 1.47511, 1.51724, 1.55802, 1.5973, 1.63489, 1.67064, 1.70439, 1.73598, 1.76527, 1.79213, 1.81643, 1.83805, 1.8569, 1.87289, 1.88593, 1.89596, 1.90294, 1.90682, 1.9076, 1.90527, 1.89983, 1.89133, 1.87978, 1.86526, 1.84783, 1.82758, 1.80461, 1.77901, 1.75092, 1.72046, 1.68777, 1.65301, 1.61632, 1.57786, 1.53781, 1.49633, 1.45359, 1.40977, 1.36504, 1.31958, 1.27355, 1.22712, 1.18045, 1.13371, 1.08704, 1.04059, 0.9945, 0.948905, 0.903921, 0.859665, 0.816241, 0.773749, 0.732271, 0.691887, 0.652662, 0.614655, 0.577921, 0.542496, 0.508411, 0.475689, 0.444348, 0.414395, 0.38583, 0.358649, 0.332839, 0.308382, 0.285256, 0.263435, 0.242885, 0.223573, 0.205461, 0.188509, 0.172673, 0.157909, 0.144172, 0.131416, 0.119592, 0.108655, 0.0985572, 0.0892521, 0.0806935, 0.0728368, 0.0656377, 0.0590535, 0.0530436, 0.0475671, 0.0425867, 0.0380655, 0.0339688, 0.0302636, 0.0269185, 0.0239041, 0.0211926]

    HalfNormalDistribution =
    [2.65307, 2.6516, 2.64718, 2.63984, 2.6296, 2.61648, 2.60054, 2.58182, 2.5604, 2.53633, 2.5097, 2.48059, 2.4491, 2.41532, 2.37936, 2.34133, 2.30135, 2.25954, 2.21603, 2.17095, 2.12442, 2.07657, 2.02755, 1.97749, 1.92653, 1.87479, 1.82242, 1.76954, 1.71629, 1.66279, 1.60917, 1.55555, 1.50205, 1.44878, 1.39584, 1.34335, 1.29139, 1.24006, 1.18946, 1.13965, 1.09071, 1.04272, 0.995728, 0.9498, 0.904985, 0.861327, 0.818864, 0.777631, 0.737653, 0.698954, 0.66155, 0.625452, 0.590667, 0.557197, 0.52504, 0.494189, 0.464635, 0.436363, 0.409356, 0.383594, 0.359054, 0.335711, 0.313537, 0.292503,0.272576, 0.253725, 0.235916, 0.219112, 0.20328, 0.188382, 0.174383, 0.161244, 0.14893, 0.137403, 0.126628, 0.116568,0.107188, 0.0984536, 0.0903304, 0.0827854, 0.0757864, 0.069302, 0.0633021, 0.0577574, 0.0526399, 0.0479225, 0.0435794, 0.0395859, 0.0359185, 0.0325546, 0.029473, 0.0266534, 0.0240769, 0.0217252, 0.0195815, 0.0176297, 0.0158548, 0.0142428, 0.0127805, 0.0114555, 0.0102566]

    EllisDistribution =
    [6.80942, 5.36792, 4.45867, 3.83457, 3.28267, 2.77244, 2.33403, 1.98717, 1.72219, 1.52427, 1.37856, 1.27023, 1.18445,1.10697, 1.03176, 0.958823, 0.888329, 0.820435, 0.755302, 0.693089, 0.633956, 0.578063, 0.525569, 0.476635, 0.431421,0.390065, 0.35249, 0.31849, 0.287855, 0.260378, 0.235851, 0.214064, 0.194811, 0.177882, 0.163069, 0.150163, 0.138957,0.129243, 0.120811, 0.113453, 0.106962, 0.101129, 0.095775, 0.0908504, 0.0863393, 0.0822255, 0.0784932, 0.0751262, 0.0721086, 0.0694244, 0.0670576]

    MeltsZirconDistribution =
    [0.282361, 0.28919, 0.296019, 0.302849, 0.30968, 0.316567, 0.323614, 0.33064, 0.337727, 0.344848, 0.352146, 0.359642,0.367482, 0.375622, 0.384052, 0.392828, 0.401968, 0.411518, 0.421442, 0.43171, 0.44213, 0.45295, 0.464036, 0.47539, 0.486938, 0.498644, 0.51075, 0.523026, 0.535688, 0.548764, 0.562124, 0.575927, 0.590363, 0.604879, 0.620415, 0.636022, 0.652333, 0.669112, 0.686441, 0.704341, 0.72283, 0.742036, 0.761964, 0.782541, 0.803718, 0.825707, 0.848386, 0.871895,0.896139, 0.920462, 0.946071, 0.972964, 0.999905, 1.02776, 1.05664, 1.08637, 1.11731, 1.14919, 1.18202, 1.21582, 1.24956, 1.28342, 1.31828, 1.35427, 1.39153, 1.43006, 1.46879, 1.50812, 1.5477, 1.58888, 1.63149, 1.6748, 1.71724, 1.76126, 1.80668, 1.85101, 1.89546, 1.94144, 1.98379, 2.02785, 2.06738, 2.10669, 2.1377, 2.16306, 2.17843, 2.17924, 2.16073, 2.11744, 2.04444, 1.93323, 1.7923, 1.62527, 1.44425, 1.25401, 1.0528, 0.843628, 0.632687, 0.421876, 0.211064, 0.000252985]

    MeltsVolcanicZirconDistribution =
    [0.54933, 0.556409, 0.563488, 0.570567, 0.577653, 0.584759, 0.591912, 0.599251, 0.606793, 0.614519, 0.622425, 0.630421, 0.63852, 0.646681, 0.654972, 0.663533, 0.672274, 0.681233, 0.690399, 0.699787, 0.709334, 0.719174, 0.729157, 0.739423, 0.749935, 0.760644, 0.771726, 0.782974, 0.794507, 0.806296, 0.818297, 0.830517, 0.842957, 0.855411, 0.866744, 0.878127, 0.889792, 0.901792, 0.914121, 0.926689, 0.939557, 0.952834, 0.966425, 0.980333, 0.994521, 1.00914, 1.02403, 1.03928, 1.05487, 1.0705, 1.08587, 1.10097, 1.11608, 1.13153, 1.1474, 1.16353, 1.18025, 1.19743, 1.21504, 1.23312, 1.25034, 1.26711, 1.28441, 1.30212, 1.32024, 1.33892, 1.35769, 1.37491, 1.3923, 1.41046, 1.42924, 1.44775, 1.46432, 1.48171, 1.49969, 1.51516, 1.53001, 1.54571, 1.5566, 1.56814, 1.57522, 1.58168, 1.58206, 1.57869, 1.56907, 1.55064, 1.51982, 1.4737, 1.40944, 1.32047, 1.21218, 1.09157, 0.965488, 0.834108, 0.697552, 0.558304, 0.418827, 0.279262, 0.139695, 0.000127237]

    ## --- End of File
