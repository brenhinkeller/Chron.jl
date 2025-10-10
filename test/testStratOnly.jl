# Make an instance of a ChronAgeData object for nsamples
nsamples = 6
smpl = ChronAgeData(nsamples)
@test smpl isa ChronAgeData
smpl.Name          = ("minimum age", "Sample 1", "Sample 2", "Sample 3", "Sample 4", "maximum age") # Et cetera
smpl.Age          .= [ 690.0,  699.1,  708.8,  723.0,  754.0,  812.0] # Measured ages
smpl.Age_sigma    .= [   7.0,    3.0,    7.0,    5.0,    5.0,    6.0] # Measured 1-σ uncertainties
smpl.Height       .= [-350.0, -355.0, -380.0, -397.0, -411.5, -420.0] # Depths below surface should be negative
smpl.Height_sigma .= fill(0.01, nsamples) # Usually assume little or no sample height uncertainty
smpl.Age_Sidedness .= zeros(nsamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
smpl.Age_Sidedness[1] = 1. # Minimum age
smpl.Age_Sidedness[end] = -1. # Maximum age
smpl.Age_Unit = "Years BP" # Unit of measurement for ages
smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma

# Configure the stratigraphic Monte Carlo model
config = StratAgeModelConfiguration()
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
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈  [775.2009717531292, 764.5613873063818, 753.7193603935707, 743.859942276044, 734.1903641532622, 724.4313056816793, 719.9612277736258, 715.7765014002005, 711.6562087664921, 709.0844186854755, 706.5904994086322, 704.1094515386092, 701.5311095414545, 698.7091903065589, 695.7689678220307]  atol=1
@test mdl.Age_025CI ≈ [750.0939493479998, 746.7621661695302, 744.0833309387968, 723.3669863534756, 718.3636351303046, 715.540624074875, 706.7137614573669, 703.1038143674268, 700.844732326315, 697.8358269342763, 696.1454906755058, 694.9648114134607, 694.0249776716082, 693.1542101854486, 684.3681541813636] atol=3
@test mdl.Age_975CI ≈ [807.2510236903488, 795.3902664676122, 763.3234305018414, 760.0376619032745, 754.7160594008493, 733.5082697327649, 731.2581768641392, 728.2692821103202, 722.3718766353179, 720.8356091969036, 718.9801564927342, 716.6713426324544, 713.1583260608256, 704.316125432296, 703.4023532440506] atol=3
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
@test all(!isnan, agedist)
@test -Inf < mean(lldist) < 0

@test NamedTuple(mdl) isa NamedTuple

# Data about hiatuses
nhiatuses = 2 # The number of hiatuses you have data for
hiatus = HiatusData(nhiatuses) # Struct to hold data
hiatus.Height         = [-371.5, -405.0 ]
hiatus.Height_sigma   = [   0.0,    0.0 ]
hiatus.Duration       = [  10.0,   12.3 ]
hiatus.Duration_sigma = [   3.1,    2.0 ]

# Run the model. Note the additional `hiatus` arguments
@time (mdl, agedist, hiatusdist, lldist) = StratMetropolis(smpl, hiatus, config); sleep(0.5)

# Test that results match expectation, within some tolerance
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈   [775.5494538893379, 765.0398661222231, 754.4553671556596, 749.1165989011665, 729.1802542341943, 724.1690788281638, 720.8671584198628, 717.7310672911206, 714.5176266505601, 713.1589898075201, 702.0799302600404, 700.7161043228054, 699.4056554863363, 698.0132290102006, 695.2079260036835] atol=1
@test mdl.Age_025CI ≈  [750.8736692043892, 747.7489306192067, 745.269061071627, 735.8291026862441, 717.8684206092925, 715.6903606224124, 710.0972278697697, 707.5386715468874, 705.8103447794734, 704.3928559634423, 694.2857415733332, 693.6245799703629, 693.0486461983832, 692.4919675548376, 684.2374067953986] atol=3
@test mdl.Age_975CI ≈  [807.3229235373376, 795.5005210793569, 763.8089856382381, 761.1489325004852, 742.73112549815, 732.7647859218115, 730.8975324858483, 728.4214438533553, 723.6903413737293, 722.624657136144, 711.3505494727046, 709.6442679173194, 707.3487791460136, 703.5758337472577, 702.7247254673973] atol=3
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
@test all(!isnan, agedist)
@test size(hiatusdist) == (nhiatuses, config.nsteps)
@test mean(hiatusdist, dims=2) ≈ [11.083103050893387; 19.955414496172235;;] atol=2
@test -Inf < mean(lldist) < 0
