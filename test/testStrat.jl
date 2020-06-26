# Make an instance of a ChronSection object for nSamples
nSamples = 4
smpl = NewChronAgeData(nSamples)
smpl.Name          = ("Sample 1", "Sample 2", "Sample 3", "Sample 4") # Et cetera
smpl.Age          .= [ 699.1,  708.8,  723.0,  754.0,] # Measured ages
smpl.Age_sigma    .= [   3.0,    7.0,    5.0,    5.0,] # Measured 1-σ uncertainties
smpl.Height       .= [ -355,  -380,-397.0,-411.5,] # Depths below surface should be negative
smpl.Height_sigma .= fill(0.01, nSamples) # Usually assume little or no sample height uncertainty
smpl.Age_Sidedness .= zeros(nSamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
smpl.Age_Unit = "Years BP" # Unit of measurement for ages
smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma

# Configure the stratigraphic Monte Carlo model
config = NewStratAgeModelConfiguration()
config.resolution = 5 # Same units as sample height. Smaller is slower!
config.bounding = 0.5 # how far away do we place runaway bounds, as a fraction of total section height. Larger is slower.
(bottom, top) = extrema(smpl.Height)
npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
config.nsteps = 1000000 # Number of steps to run in distribution MCMC
config.burnin = 100000*npoints_approx # Number to discard
config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps

# Run the stratigraphic MCMC model
@time (mdl, agedist, lldist) = StratMetropolis(smpl, config)

# Test that results match expectation, within some tolerance
@test isapprox(mdl.Age, [751.74, 742.6997, 733.6532, 724.5082, 720.1385, 715.9552, 711.7907, 709.1437, 706.5679, 703.9969, 701.3627], atol=1)
@test isapprox(mdl.Age_025CI, [742.4111, 723.8341, 718.6988, 715.692, 707.4121, 703.7166, 701.3356, 698.159, 696.2699, 694.9445, 693.8905], atol=1)
@test isapprox(mdl.Age_975CI, [761.0398, 757.8618, 752.6434, 733.5778, 731.2341, 728.1138, 722.2433, 720.5815, 718.6199, 716.1298, 712.4586], atol=1)
