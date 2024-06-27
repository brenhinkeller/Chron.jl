# Make an instance of a ChronAgeData object for nSamples
nSamples = 4
smpl = ChronAgeData(nSamples)
@test smpl isa ChronAgeData
smpl.Name          = ("Sample 1", "Sample 2", "Sample 3", "Sample 4") # Et cetera
smpl.Age          .= [ 699.1,  708.8,  723.0,  754.0,] # Measured ages
smpl.Age_sigma    .= [   3.0,    7.0,    5.0,    5.0,] # Measured 1-σ uncertainties
smpl.Height       .= [  -355,   -380, -397.0, -411.5,] # Depths below surface should be negative
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
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈ [751.87, 742.76, 733.65, 724.51, 720.15, 715.96, 711.8, 709.15, 706.6, 704.04, 701.42, 698.67] atol=1
@test mdl.Age_025CI ≈ [742.48, 723.81, 718.7, 715.73, 707.42, 703.73, 701.34, 698.2, 696.35, 695.04, 694.0, 693.05] atol=3
@test mdl.Age_975CI ≈ [761.17, 757.98, 752.64, 733.57, 731.21, 728.07, 722.24, 720.59, 718.6, 716.11, 712.48, 704.28] atol=3
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
@test all(!isnan, agedist)
@test -Inf < mean(lldist) < 0

@test NamedTuple(mdl) isa NamedTuple

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
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈ [752.76, 747.88, 728.92, 724.14, 720.8, 717.62, 714.44, 712.95, 711.41, 700.82, 699.4, 697.99] atol=1
@test mdl.Age_025CI ≈ [743.41, 734.96, 717.54, 715.41, 709.67, 707.05, 705.26, 703.72, 702.49, 693.27, 692.65, 692.06] atol=3
@test mdl.Age_975CI ≈ [761.77, 759.23, 741.37, 732.59, 730.64, 728.05, 723.52, 722.34, 720.9, 709.03, 706.78, 703.46] atol=3
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
@test all(!isnan, agedist)
@test size(hiatusdist) == (nHiatuses, config.nsteps)
@test mean(hiatusdist, dims=2) ≈ [10.580012942504894; 18.96167245288326;;] atol=2
@test -Inf < mean(lldist) < 0

## --- Strat only, general Distributions-based case

# Make an instance of a ChronAgeData object for nSamples
nSamples = 4
smpl = GeneralAgeData(nSamples)
@test smpl isa GeneralAgeData