## --- Test construction of curves
@test intcal13.Age_14C[1:10] == [196.8, 194.6, 192.4, 190.2, 188.0, 184.8, 181.6, 178.4, 175.2, 172.0]
@test intcal20.Age_14C[1:10] == [197.0, 195.0, 193.0, 190.0, 188.0, 185.0, 181.0, 178.0, 174.0, 170.0]
@test shcal20.Age_14C[1:10] == [171.0, 168.0, 164.0, 161.0, 159.0, 157.0, 155.0, 154.0, 154.0, 154.0]


## --- Test age-depth model

# Input the number of samples we wish to model (must match below)
nSamples = 4
# Make an instance of a ChronSection object for nSamples
smpl = NewChronAgeData(nSamples)
smpl.Name           = ("Sample 1", "Sample 2", "Sample 3", "Sample 4") # Et cetera
smpl.Age_14C       .= [ 6991,  7088,  7230,  7540,] # Measured ages
smpl.Age_14C_sigma .= [   30,    70,    50,    50,] # Measured 1-σ uncertainties
smpl.Height        .= [ -355,  -380,-397.0,-411.5,] # Depths below surface should be negative
smpl.Height_sigma  .= fill(0.01, nSamples) # Usually assume little or no sample height uncertainty
smpl.Age_Sidedness .= zeros(nSamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
smpl.Age_Unit = "Years BP" # Unit of measurement for ages
smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma

# Choice of radiocarbon calibration
# Options include: intcal13, intcal20, shcal20
calibration = intcal20

# Calculate calendar age PDFs for each sample
smpl.Params = fill(NaN, length(calibration.Age_Calendar), nSamples)
for i = 1:nSamples
    # The likelihood that a measured 14C age could result from a sample of
    # a given calendar age is proportional to the intergral of the product
    # of the two respective distributions
    likelihood = normproduct.(smpl.Age_14C[i], smpl.Age_14C_sigma[i], calibration.Age_14C, calibration.Age_sigma)
    likelihood ./= sum(likelihood) * calibration.dt # Normalize

    samples = draw_from_distribution(likelihood, 10^4) .* maximum(calibration.Age_Calendar)
    smpl.Age[i] = mean(samples)
    smpl.Age_sigma[i] = std(samples)

    # Populate smpl.Params with log likelihood for each sample
    smpl.Params[:,i] = normproduct_ll.(smpl.Age_14C[i], smpl.Age_14C_sigma[i], calibration.Age_14C, calibration.Age_sigma)
end

# Run stratigraphic model
# Configure the stratigraphic Monte Carlo model
config = NewStratAgeModelConfiguration()
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
@test isa(mdl.Age, Array{Float64,1})
@test isapprox(mdl.Age, [8319.5, 8235.7, 8154.8, 8072.3, 8025.8, 7982.2, 7938.8, 7914.7, 7891.3, 7868.5, 7845.1], atol=30)
@test isapprox(mdl.Age_025CI, [8205.0, 8045.4, 7998.9, 7975.5, 7893.6, 7859.9, 7838.6, 7802.2, 7779.3, 7761.3, 7747.5], atol=35)
@test isapprox(mdl.Age_975CI, [8403.9, 8385.1, 8345.7, 8163.4, 8150.0, 8123.4, 8028.4, 8014.9, 8002.3, 7984.7, 7957.1], atol=35)


# Data about hiatuses
nHiatuses = 2 # The number of hiatuses you have data for
hiatus = NewHiatusData(nHiatuses) # Struct to hold data
hiatus.Height         = [-371.5, -405.0 ]
hiatus.Height_sigma   = [   0.0,    0.0 ]
hiatus.Duration       = [ 100.0,  123.0 ]
hiatus.Duration_sigma = [  30.5,   20.0 ]

# Run the model. Note the additional `hiatus` arguments
@time (mdl, agedist, hiatusdist, lldist) = StratMetropolis14C(smpl, hiatus, config); sleep(0.5)

# Test that results match expectation, within some tolerance
@test isa(mdl.Age, Array{Float64,1})
@test isapprox(mdl.Age, [8339.8, 8160.6, 8112.4, 8066.9, 8032.9, 7999.7, 7964.7, 7947.4, 7850.4, 7833.2, 7818.8], atol=30)
@test isapprox(mdl.Age_025CI, [8216.9, 8022.6, 7991.4, 7975.1, 7924.7, 7892.3, 7872.9, 7844.2, 7749.1, 7738.5, 7729.2], atol=35)
@test isapprox(mdl.Age_975CI, [8408.3, 8272.7, 8238.4, 8160.2, 8149.8, 8130.8, 8092.3, 8064.0, 7960.0, 7935.8, 7916.7], atol=35)
