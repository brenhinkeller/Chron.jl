# Make an instance of a ChronAgeData object for nSamples
nSamples = 4
smpl = ChronAgeData(nSamples)
smpl.Name          = ("Sample 1", "Sample 2", "Sample 3", "Sample 4") # Et cetera
smpl.Age          .= [ 699.1,  708.8,  723.0,  754.0,] # Measured ages
smpl.Age_sigma    .= [   3.0,    7.0,    5.0,    5.0,] # Measured 1-σ uncertainties
smpl.Height       .= [ -355,  -380,-397.0,-411.5,] # Depths below surface should be negative
smpl.Height_sigma .= fill(0.01, nSamples) # Usually assume little or no sample height uncertainty
smpl.Age_Sidedness .= zeros(nSamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
smpl.Age_Unit = "Years BP" # Unit of measurement for ages
smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma

# Configure the stratigraphic Monte Carlo model
config = StratAgeModelConfiguration()
config.resolution = 5 # Same units as sample height. Smaller is slower!
config.bounding = 0.5 # how far away do we place runaway bounds, as a fraction of total section height. Larger is slower.
(bottom, top) = extrema(smpl.Height)
npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
config.nsteps = 100000 # Number of steps to run in distribution MCMC
config.burnin = 100000*npoints_approx # Number to discard
config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps

# Run the stratigraphic MCMC model
@time (mdl, agedist, lldist) = StratMetropolis(smpl, config)

# Test that results match expectation, within some tolerance
@test isa(mdl.Age, Array{Float64,1})
@test all(isapprox.(mdl.Age, [751.74, 742.70, 733.65, 724.51, 720.14, 715.96, 711.79, 709.14, 706.57, 703.70, 701.36], atol=1))
@test all(isapprox.(mdl.Age_025CI, [742.41, 723.83, 718.67, 715.69, 707.41, 703.72, 701.34, 698.16, 696.27, 694.94, 693.89], atol=1.5))
@test all(isapprox.(mdl.Age_975CI, [761.04, 757.86, 752.64, 733.58, 731.23, 728.11, 722.24, 720.58, 718.62, 716.13, 712.46], atol=1.5))
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])


# Data about hiatuses
nHiatuses = 2 # The number of hiatuses you have data for
hiatus = HiatusData(nHiatuses) # Struct to hold data
hiatus.Height         = [-371.5, -405.0 ]
hiatus.Height_sigma   = [   0.0,    0.0 ]
hiatus.Duration       = [  10.0,   12.3 ]
hiatus.Duration_sigma = [   3.1,    2.0 ]

# Run the model. Note the additional `hiatus` arguments
@time (mdl, agedist, hiatusdist, lldist) = StratMetropolis(smpl, hiatus, config); sleep(0.5)


# Test that results match expectation, within some tolerance
@test isa(mdl.Age, Array{Float64,1})
@test all(isapprox.(mdl.Age, [752.93, 734.7, 729.3, 724.28, 721.11, 718.03, 714.86, 713.05, 703.2, 701.25, 699.62], atol=1))
@test all(isapprox.(mdl.Age_025CI, [743.85, 721.93, 718.22, 715.98, 710.55, 707.7, 705.64, 703.03, 694.52, 693.66, 692.99], atol=2))
@test all(isapprox.(mdl.Age_975CI, [762.02, 746.97, 742.23, 732.72, 730.85, 728.43, 724.24, 723.16, 713.8, 711.07, 708.2], atol=2))
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
