# Make an instance of a ChronAgeData object for nSamples
nSamples = 4
smpl = ChronAgeData(nSamples)
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
@test mdl.Age_025CI ≈ [742.48, 723.81, 718.7, 715.73, 707.42, 703.73, 701.34, 698.2, 696.35, 695.04, 694.0, 693.05] atol=1.5
@test mdl.Age_975CI ≈ [761.17, 757.98, 752.64, 733.57, 731.21, 728.07, 722.24, 720.59, 718.6, 716.11, 712.48, 704.28] atol=1.5
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
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈ [753.23, 748.14, 729.73, 724.3, 721.09, 718.04, 714.96, 713.54, 711.82, 701.99, 699.95, 698.21] atol=1
@test mdl.Age_025CI ≈ [744.19, 734.06, 718.01, 715.91, 710.5, 707.73, 705.75, 703.97, 701.87, 694.01, 693.25, 692.56] atol=2
@test mdl.Age_975CI ≈ [762.31, 760.19, 743.82, 732.92, 730.96, 728.45, 724.28, 723.32, 722.16, 712.29, 708.49, 703.75] atol=2
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
