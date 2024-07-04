# Make an instance of a ChronAgeData object for nSamples
nSamples = 6
smpl = ChronAgeData(nSamples)
@test smpl isa ChronAgeData
smpl.Name          = ("minimum age", "Sample 1", "Sample 2", "Sample 3", "Sample 4", "maximum age") # Et cetera
smpl.Age          .= [ 690.0,  699.1,  708.8,  723.0,  754.0,  812.0] # Measured ages
smpl.Age_sigma    .= [   7.0,    3.0,    7.0,    5.0,    5.0,    6.0] # Measured 1-σ uncertainties
smpl.Height       .= [-350.0, -355.0, -380.0, -397.0, -411.5, -420.0] # Depths below surface should be negative
smpl.Height_sigma .= fill(0.01, nSamples) # Usually assume little or no sample height uncertainty
smpl.Age_Sidedness .= zeros(nSamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
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
config.nsteps = 100000 # Number of steps to run in distribution MCMC
config.burnin = 100000*npoints_approx # Number to discard
config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps

# Run the stratigraphic MCMC model
@time (mdl, agedist, lldist) = StratMetropolis(smpl, config)

# Test that results match expectation, within some tolerance
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈  [774.446455402488, 764.138871634567, 753.6234690338521, 743.7736679540459, 734.1288211670727, 724.4222568319159, 719.9641790169545, 715.7853469234702, 711.659901307201, 709.0927902799843, 706.6073112860937, 704.1349331543528, 701.576806557881, 698.7594637725156, 695.9077192905107] atol=1
@test mdl.Age_025CI ≈ [749.7477463948217, 746.630273097855, 743.9690750302638, 723.4084389968576, 718.3588586689344, 715.5410003896017, 706.7023568283762, 703.1283100699234, 700.9213523438647, 697.904849247825, 696.2105445428034, 695.0433927602111, 694.1016935971145, 693.2625238063798, 684.8895558316834] atol=3
@test mdl.Age_975CI ≈ [806.7021880750214, 794.9197023904139, 763.2743319392429, 760.0094083209071, 754.7169964175238, 733.5482560401001, 731.2806438300145, 728.3050781269916, 722.3491038457139, 720.830748668773, 718.9846443130959, 716.6641608791327, 713.2397206761551, 704.3594827974073, 703.4763901911105] atol=3
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
@test mdl.Age ≈ [774.9621590393864, 764.7483478838244, 754.3928287806855, 749.1279876514806, 729.1819943547974, 724.1793043771663, 720.8913053718023, 717.7627759721249, 714.5509500060313, 713.2002472528285, 702.1296086101366, 700.7704707147774, 699.4643412956174, 698.0733121394825, 695.3696935693185] atol=1
@test mdl.Age_025CI ≈ [750.6328462444853, 747.6305055671014, 745.2219453018333, 735.8795895625685, 717.8606462193696, 715.6849749469166, 710.1425439567424, 707.5591665793866, 705.808337538727, 704.4299789235162, 694.3410849948168, 693.699614909472, 693.1108571564564, 692.5485716537205, 684.6340487749085] atol=3
@test mdl.Age_975CI ≈ [806.6715104282963, 794.8753648212958, 763.7844862197752, 761.1148678026058, 742.5015528417246, 732.771524797932, 730.9079222710706, 728.4137307378472, 723.6929471023955, 722.6480626953821, 711.3953411778708, 709.6936433474348, 707.3968615013218, 703.6205967015544, 702.7858695444451] atol=3
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