## --- Test construction of curves
@test intcal13.Age_14C[1:10] == [196.8, 194.6, 192.4, 190.2, 188.0, 184.8, 181.6, 178.4, 175.2, 172.0]
@test intcal20.Age_14C[1:10] == [197.0, 195.0, 193.0, 190.0, 188.0, 185.0, 181.0, 178.0, 174.0, 170.0]
@test shcal20.Age_14C[1:10] == [171.0, 168.0, 164.0, 161.0, 159.0, 157.0, 155.0, 154.0, 154.0, 154.0]
@test marine20.Age_14C[1:10] == [603.0, 603.0, 603.0, 603.0, 603.0, 603.0, 603.0, 603.0, 603.0, 603.0]

## --- Test age-depth model

# Input the number of samples we wish to model (must match below)
nsamples = 4
# Make an instance of a ChronSection object for nsamples
smpl = ChronAgeData(nsamples)
smpl.Name           = ("Sample 1", "Sample 2", "Sample 3", "Sample 4") # Et cetera
smpl.Age_14C       .= [ 6991,  7088,  7230,  7540,] # Measured ages
smpl.Age_14C_sigma .= [   30,    70,    50,    50,] # Measured 1-σ uncertainties
smpl.Height        .= [ -355,  -380,-397.0,-411.5,] # Depths below surface should be negative
smpl.Height_sigma  .= fill(0.01, nsamples) # Usually assume little or no sample height uncertainty
smpl.Age_Sidedness .= zeros(nsamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
smpl.Age_Unit = "Years BP" # Unit of measurement for ages
smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma

# Choice of radiocarbon calibration
# Options include: intcal13, intcal20, shcal20
calibration = intcal13

# Calculate calendar age PDFs for each sample
smpl.Params = fill(NaN, length(calibration.Age_Calendar), nsamples)
for i = 1:nsamples
    # The likelihood that a measured 14C age could result from a sample of
    # a given calendar age is proportional to the intergral of the product
    # of the two respective distributions
    likelihood = normproduct.(smpl.Age_14C[i], smpl.Age_14C_sigma[i], calibration.Age_14C, calibration.Age_sigma)
    likelihood ./= sum(likelihood) * calibration.dt # Normalize

    samples = draw_from_distribution(likelihood, 10^4) .* maximum(calibration.Age_Calendar)
    smpl.Age[i] = mean(samples)
    smpl.Age_sigma[i] = std(samples)

    # Populate smpl.Params with log likelihood for each sample
    smpl.Params[:,i] = normlogproduct.(smpl.Age_14C[i], smpl.Age_14C_sigma[i], calibration.Age_14C, calibration.Age_sigma)
end

# Run stratigraphic model
# Configure the stratigraphic Monte Carlo model
config = StratAgeModelConfiguration()
# If you in doubt, you can probably leave these parameters as-is
config.resolution = 5 # Same units as sample height. Smaller is slower!
config.bounding = 0.5 # how far away do we place runaway bounds, as a fraction of total section height
(bottom, top) = extrema(smpl.Height)
npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
config.nsteps = 100000 # Number of steps to run in distribution MCMC
config.burnin = 100000*npoints_approx # Number to discard
config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Run the stratigraphic MCMC model
@time (mdl, agedist, lldist) = StratMetropolis14C(smpl, config)


# Test that results match expectation, within some tolerance
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈ [8322.52, 8237.68, 8155.65, 8072.77, 8026.99, 7983.31, 7939.54, 7915.44, 7892.49, 7869.92, 7846.76, 7822.28] atol=20
@test mdl.Age_025CI ≈ [8205.9, 8045.4, 7999.0, 7975.6, 7894.3, 7860.6, 7839.4, 7802.5, 7780.8, 7762.1, 7747.9, 7736.0] atol=35
@test mdl.Age_975CI ≈ [8405.0, 8386.1, 8347.1, 8163.5, 8151.0, 8123.9, 8029.5, 8016.3, 8003.0, 7985.9, 7958.2, 7916.3] atol=35
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
@test all(!isnan, agedist)
@test -Inf < mean(lldist) < 0


# Data about hiatuses
nHiatuses = 2 # The number of hiatuses you have data for
hiatus = HiatusData(nHiatuses) # Struct to hold data
hiatus.Height         = [-371.5, -405.0 ]
hiatus.Height_sigma   = [   0.0,    0.0 ]
hiatus.Duration       = [ 100.0,  123.0 ]
hiatus.Duration_sigma = [  30.5,   20.0 ]

# Run the model. Note the additional `hiatus` arguments
@time (mdl, agedist, hiatusdist, lldist) = StratMetropolis14C(smpl, hiatus, config); sleep(0.5)

# Test that results match expectation, within some tolerance
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈ [8339.84, 8294.33, 8107.77, 8064.14, 8029.11, 7995.66, 7961.79, 7949.08, 7936.17, 7828.04, 7815.81, 7803.7] atol=20
@test mdl.Age_025CI ≈ [8217.14, 8161.73, 7991.82, 7974.87, 7921.84, 7891.41, 7874.27, 7858.89, 7844.64, 7739.61, 7730.51, 7721.97] atol=35
@test mdl.Age_975CI ≈ [8408.06, 8394.94, 8231.34, 8159.62, 8147.76, 8123.97, 8054.25, 8039.57, 8025.91, 7918.86, 7907.9, 7889.46] atol=35
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
@test all(!isnan, agedist)
@test size(hiatusdist) == (nHiatuses, config.nsteps)
@test mean(hiatusdist, dims=2) ≈ [108.4; 189.0;;] atol=20
@test -Inf < mean(lldist) < 0

## --- IntCal20

# Choice of radiocarbon calibration
# Options include: intcal13, intcal20, shcal20
calibration = intcal20

# Calculate calendar age PDFs for each sample
smpl.Params = fill(NaN, length(calibration.Age_Calendar), nsamples)
for i = 1:nsamples
    # The likelihood that a measured 14C age could result from a sample of
    # a given calendar age is proportional to the intergral of the product
    # of the two respective distributions
    likelihood = normproduct.(smpl.Age_14C[i], smpl.Age_14C_sigma[i], calibration.Age_14C, calibration.Age_sigma)
    likelihood ./= sum(likelihood) * calibration.dt # Normalize

    samples = draw_from_distribution(likelihood, 10^4) .* maximum(calibration.Age_Calendar)
    smpl.Age[i] = mean(samples)
    smpl.Age_sigma[i] = std(samples)

    # Populate smpl.Params with log likelihood for each sample
    smpl.Params[:,i] = normlogproduct.(smpl.Age_14C[i], smpl.Age_14C_sigma[i], calibration.Age_14C, calibration.Age_sigma)
end

# Run stratigraphic model
# Configure the stratigraphic Monte Carlo model
config = StratAgeModelConfiguration()
# If you in doubt, you can probably leave these parameters as-is
config.resolution = 5 # Same units as sample height. Smaller is slower!
config.bounding = 0.5 # how far away do we place runaway bounds, as a fraction of total section height
(bottom, top) = extrema(smpl.Height)
npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
config.nsteps = 100000 # Number of steps to run in distribution MCMC
config.burnin = 100000*npoints_approx # Number to discard
config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Run the stratigraphic MCMC model
@time (mdl, agedist, lldist) = StratMetropolis14C(smpl, config)


# Test that results match expectation, within some tolerance
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈ [8308.5, 8228.4, 8150.6, 8070.8, 8024.0, 7980.6, 7936.7, 7912.4, 7889.1, 7866.3, 7842.8, 7818.2] atol=20
@test mdl.Age_025CI ≈ [8200.7, 8034.8, 7991.0, 7970.6, 7888.0, 7854.5, 7832.5, 7795.1, 7771.7, 7754.0, 7740.8, 7724.9] atol=40
@test mdl.Age_975CI ≈ [8401.5, 8382.0, 8342.7, 8167.4, 8156.2, 8128.5, 8019.7, 8011.7, 8002.5, 7986.2, 7957.2, 7916.2] atol=40
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
@test all(!isnan, agedist)
@test -Inf < mean(lldist) < 0


# Data about hiatuses
nHiatuses = 2 # The number of hiatuses you have data for
hiatus = HiatusData(nHiatuses) # Struct to hold data
hiatus.Height         = [-371.5, -405.0 ]
hiatus.Height_sigma   = [   0.0,    0.0 ]
hiatus.Duration       = [ 100.0,  123.0 ]
hiatus.Duration_sigma = [  30.5,   20.0 ]

# Run the model. Note the additional `hiatus` arguments
@time (mdl, agedist, hiatusdist, lldist) = StratMetropolis14C(smpl, hiatus, config); sleep(0.5)

# Test that results match expectation, within some tolerance
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈ [8327.67, 8283.87, 8099.66, 8057.74, 8023.23, 7990.66, 7957.71, 7944.6, 7931.19, 7823.06, 7810.29, 7797.89] atol=20
@test mdl.Age_025CI ≈ [8210.55, 8152.11, 7983.06, 7968.49, 7913.95, 7885.02, 7868.98, 7852.67, 7838.42, 7730.2, 7719.48, 7710.82] atol=40
@test mdl.Age_975CI ≈ [8405.89, 8391.62, 8225.61, 8163.41, 8151.31, 8123.92, 8030.35, 8023.65, 8016.86, 7916.1, 7905.95, 7886.48] atol=40
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
@test all(!isnan, agedist)
@test size(hiatusdist) == (nHiatuses, config.nsteps)
@test mean(hiatusdist, dims=2) ≈ [110.5; 187.3;;]  atol=20
@test -Inf < mean(lldist) < 0