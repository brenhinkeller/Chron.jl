nSamples = 5 # The number of samples you have data for
smpl = ChronAgeData(nSamples)
smpl.Name      =  ("KJ08-157", "KJ04-75", "KJ09-66",  "KJ04-72", "KJ04-70",)
smpl.Height   .=  [     -52.0,      44.0,      54.0,      82.0,      93.0,]
smpl.Height_sigma .= [    3.0,       1.0,       3.0,       3.0,       3.0,]
smpl.Age_Sidedness .= zeros(nSamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
smpl.Path = abspath("../examples/DenverUPbExampleData/") # Where are the data files?
smpl.inputSigmaLevel = 2 # i.e., are the data files 1-sigma or 2-sigma. Integer.
smpl.Age_Unit = "Ma" # Unit of measurement for ages and errors in the data files
smpl.Height_Unit = "cm" # Unit of measurement for Height and Height_sigma

# Remove outliers (if any)
smpl = screen_outliers(smpl, maxgap=50, make_plots=false)

# Distribution boostrapping from chron strat object
BootstrappedDistribution = BootstrapCrystDistributionKDE(smpl)
@test isa(BootstrappedDistribution, Array{Float64,1})

# Estimate the eruption age distributions for each sample  - - - - - - - -

# Configure distribution model here
distSteps = 10^5 # Number of steps to run in distribution MCMC
distBurnin = floor(Int,distSteps/2) # Number to discard

# Run MCMC to estimate saturation and eruption/deposition age distributions
@time tMinDistMetropolis(smpl,distSteps,distBurnin,BootstrappedDistribution; make_plots=false)

# Run stratigraphic model - - - - - - - - - - - - - - - - - - - - - - - - -
# Configure the stratigraphic Monte Carlo model
config = StratAgeModelConfiguration()
# If you in doubt, you can probably leave these parameters as-is
config.resolution = 10.0 # Same units as sample height. Smaller is slower!
config.bounding = 0.5 # how far away do we place runaway bounds, as a fraction of total section height
(bottom, top) = extrema(smpl.Height)
npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
config.nsteps = 100000 # Number of steps to run in distribution MCMC
config.burnin = 100000*npoints_approx # Number to discard
config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Run the stratigraphic MCMC model
println("StratMetropolisDist:")
@time (mdl, agedist, lldist) = StratMetropolisDist(smpl, config)

# Test that results match expectation, within some tolerance
@test isa(mdl.Age, Array{Float64,1})
@test isapprox(mdl.Age, [66.06, 66.05, 66.03, 66.02, 66.01, 66.0, 65.98, 65.97, 65.96, 65.94, 65.94, 65.93, 65.93, 65.9], atol=0.1)
@test isapprox(mdl.Age_025CI, [66.0, 65.98, 65.96, 65.95, 65.94, 65.93, 65.93, 65.92, 65.92, 65.91, 65.91, 65.89, 65.88, 65.82], atol=0.15)
@test isapprox(mdl.Age_975CI, [66.09, 66.09, 66.08, 66.08, 66.07, 66.07, 66.05, 66.04, 66.02, 65.97, 65.97, 65.96, 65.96, 65.95], atol=0.15)

# Try adding systematic uncertainties too
smpl.Chronometer = (:UPb, :UPb, :ArAr, :UPb, :UPb)
systematic=SystematicUncertainty()
systematic.ArAr = 0.005/2 # One-sigma
systematic.UPb = 0.005/2 # One-sigma

# Run the stratigraphic MCMC model
println("StratMetropolisDist with systematic uncertainties:")
@time (mdl, agedist, lldist) = StratMetropolisDist(smpl, config, systematic)

# Test that results match expectation, within some tolerance
@test isa(mdl.Age, Array{Float64,1})
@test isapprox(mdl.Age, [66.06, 66.05, 66.03, 66.02, 66.01, 66.0, 65.98, 65.97, 65.96, 65.94, 65.94, 65.93, 65.93, 65.9], atol=0.1)
@test isapprox(mdl.Age_025CI, [66.0, 65.98, 65.96, 65.95, 65.94, 65.93, 65.93, 65.92, 65.92, 65.91, 65.91, 65.89, 65.88, 65.82], atol=0.15)
@test isapprox(mdl.Age_975CI, [66.09, 66.09, 66.08, 66.08, 66.07, 66.07, 66.05, 66.04, 66.02, 65.97, 65.97, 65.96, 65.96, 65.95], atol=0.15)

# Data about hiatuses
nHiatuses = 2 # The number of hiatuses you have data for
hiatus = HiatusData(nHiatuses) # Struct to hold data
hiatus.Height         = [-7.0, 35.0 ]
hiatus.Height_sigma   = [ 0.0,  0.0 ]
hiatus.Duration       = [ 0.3,  0.3 ]
hiatus.Duration_sigma = [ 0.05, 0.05]

# Run the model. Note the additional `hiatus` arguments
println("StratMetropolisDist with hiata:")
@time (mdl, agedist, hiatusdist, lldist) = StratMetropolisDist(smpl, hiatus, config)

# Test that results match expectation, within some tolerance
@test isa(mdl.Age, Array{Float64,1})
@test isapprox(mdl.Age, [66.08, 66.07, 66.07, 66.07, 66.02, 66.01, 66.01, 66.01, 65.94, 65.94, 65.93, 65.93, 65.92, 65.9], atol=0.1)
@test isapprox(mdl.Age_025CI, [66.05, 66.04, 66.03, 66.02, 65.94, 65.94, 65.93, 65.93, 65.91, 65.9, 65.9, 65.89, 65.88, 65.82], atol=0.15)
@test isapprox(mdl.Age_975CI, [66.1, 66.1, 66.1, 66.1, 66.08, 66.08, 66.08, 66.08, 65.98, 65.96, 65.96, 65.96, 65.96, 65.95], atol=0.15)

## --- As above, but treat everything as a gaussian/weighted mean

# Tel tMinDistMetropolis to treat these as gaussians, using the first row of data file
smpl.Age_DistType.=1

# Run MCMC to estimate saturation and eruption/deposition age distributions
@time tMinDistMetropolis(smpl,distSteps,distBurnin,BootstrappedDistribution; make_plots=false)

# Run the stratigraphic MCMC model
println("StratMetropolisDist with fitted Gaussians:")
@time (mdl, agedist, lldist) = StratMetropolisDist(smpl, config)
@test isa(mdl.Age, Array{Float64,1})
@test isapprox(mdl.Age, [65.97, 65.97, 65.96, 65.95, 65.95, 65.94, 65.93, 65.92, 65.92, 65.91, 65.9, 65.89, 65.87, 65.85], atol=0.1)
@test isapprox(mdl.Age_025CI, [65.86, 65.86, 65.85, 65.85, 65.84, 65.84, 65.84, 65.83, 65.83, 65.83, 65.82, 65.77, 65.74, 65.72], atol=0.15)
@test isapprox(mdl.Age_975CI, [66.07, 66.06, 66.06, 66.06, 66.06, 66.05, 66.04, 66.03, 66.02, 66.0, 65.99, 65.98, 65.97, 65.96], atol=0.15)
