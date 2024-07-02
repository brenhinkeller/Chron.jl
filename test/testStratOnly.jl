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
@test mdl.Age ≈ [775.639803231674, 764.7522135364757, 753.6908761271258, 743.7955507171365, 734.1450626270426, 724.4181380914805, 719.9344128500521, 715.730084005377, 711.5837749049695, 708.9681700152063, 706.4443998136725, 703.9187181054767, 701.2964403515165, 698.4067159228537, 694.0564170150618] atol=1
@test mdl.Age_025CI ≈ [749.86998905391, 746.6946942599805, 744.0594596745302, 723.3563855942004, 718.369656399334, 715.5038496825758, 706.6306024739025, 703.0289240992544, 700.7652433271606, 697.635218301761, 695.8957733017851, 694.6899457591284, 693.7027382576433, 692.8123972703291, 681.3406483559285] atol=3
@test mdl.Age_975CI ≈ [811.2671899519365, 797.4625865888173, 763.3573868897032, 760.0445173042854, 754.7002448815252, 733.5275337761803, 731.2848593336894, 728.2781265900221, 722.32529630248, 720.7715897667372, 718.9216447983222, 716.5564071732707, 713.0417479766365, 704.0850806207267, 702.9561071444978] atol=3
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
@test mdl.Age ≈ [776.1645401715448, 765.3913415310423, 754.4578701374423, 749.1402317775008, 729.1734073022601, 724.1504941010227, 720.8064612101196, 717.6325237531024, 714.3911204374407, 712.9981743647901, 701.8475241735739, 700.4609977090519, 699.1253113675604, 697.6951927842167, 693.5283072170754] atol=1
@test mdl.Age_025CI ≈ [750.7546389575148, 747.721013702703, 745.2360260172313, 735.8600442968074, 717.8362687469073, 715.6571312160221, 709.9467580290137, 707.3352431871043, 705.6012864195188, 704.1631160813711, 693.9692554028434, 693.2871054139605, 692.6791996279173, 692.0920645898433, 681.1134640520288] atol=3
@test mdl.Age_975CI ≈ [811.4534097071294, 797.7943496043165, 763.8616858646202, 761.1731148399618, 742.6187531312476, 732.7757692120042, 730.8844293232036, 728.3748870708009, 723.607690594645, 722.5139094919781, 711.2193585524507, 709.48421725086, 707.154594077965, 703.3249584517113, 702.2469721879497] atol=3
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