# Input the number of samples we wish to model (must match below)
nsamples = 5
# Make an instance of a GeneralAgeData object for n samples
smpl = GeneralAgeData(nsamples)
smpl.Name             = (      "Sample 1",      "Sample 2",         "Sample 3",       "Sample 4",        "Sample 5",) # Et cetera
smpl.Age_Distribution = [Normal(39.5,0.1), Uniform(37, 38),  Normal(36.3, 0.1), Uniform(33.5,34), Normal(32.1, 0.1),] # Measured ages
smpl.Height           = [             100,             200,                300,              400,               500,] # Depths below surface should be negative
smpl.Age_Sidedness    = zeros(nsamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
smpl.Age_Unit = "Ma" # Unit of measurement for ages
smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma
@test smpl isa GeneralAgeData

# Configure the stratigraphic Monte Carlo model
config = StratAgeModelConfiguration()
config.resolution = 25 # Same units as sample height. Smaller is slower!
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
@test mdl.Age ≈ [39.49042255544478, 38.96722969630049, 38.48901032510139, 38.019348376156465, 37.547914861440425, 37.24037420316064, 36.9395916349458, 36.63320760328689, 36.28718539131596, 35.632664873531674, 35.01193213752288, 34.39264774277665, 33.76253985253119, 33.34359879740881, 32.94572880571547, 32.544987349727144, 32.10539570198641] atol=1
@test mdl.Age_025CI ≈ [39.29419983674923, 37.77309983527187, 37.36277515395373, 37.14051900815384, 37.040464130339146, 36.40376761147781, 36.233783270176595, 36.15245996064178, 36.092534551946756, 34.19714751528934, 33.75747057397605, 33.57285811598977, 33.514620066299244, 32.32781095001219, 32.07773907051565, 31.97549199229507, 31.90826940222938] atol=3
@test mdl.Age_975CI ≈ [39.68874305750388, 39.61422452075973, 39.495106359902124, 39.18726807289346, 37.97994610871827, 37.92642713833131, 37.79261699278054, 37.513276654055716, 36.48525650077192, 36.40149071129732, 36.2515851822035, 35.82327614288284, 33.98945821561441, 33.952571632565075, 33.84021405731478, 33.56418961260443, 32.30140448263967] atol=3
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
@test all(!isnan, agedist)
@test -Inf < mean(lldist) < 0

@test NamedTuple(mdl) isa NamedTuple

# Data about hiatuses
nhiatuses = 2 # The number of hiatuses you have data for
hiatus = HiatusData(nhiatuses) # Struct to hold data
hiatus.Height         = [ 150.0,  350.0 ]
hiatus.Height_sigma   = [   0.0,    0.0 ]
hiatus.Duration       = [   1.0,    2.3 ]
hiatus.Duration_sigma = [   1.0,    0.5 ]

# Run the model. Note the additional `hiatus` arguments
@time (mdl, agedist, hiatusdist, lldist) = StratMetropolis(smpl, hiatus, config); sleep(0.5)

# Test that results match expectation, within some tolerance
@test mdl.Age isa Vector{Float64}
@test mdl.Age ≈ [39.496472648071595, 39.091141330052544, 38.73438343576809, 37.76813820088203, 37.509275495475784, 37.21858934105048, 36.93167371762563, 36.63650162848444, 36.31150709915855, 36.0963194277336, 35.89241061256909, 33.940706599457855, 33.72854324734922, 33.32058089645568, 32.927941172657285, 32.533392897014124, 32.10500838284318] atol=1
@test mdl.Age_025CI ≈ [39.30060372750883, 38.06383585705967, 37.726993370228996, 37.0984623727761, 37.033368139697366, 36.410222472054656, 36.24633402828455, 36.170530411220994, 36.114657632203894, 35.38520367332072, 35.06911166768526, 33.53357860639156, 33.51023088836563, 32.318385132333965, 32.0764177648628, 31.975135830594798, 31.907829938419805] atol=3
@test mdl.Age_975CI ≈ [39.694290669563514, 39.63149256606031, 39.54227860914912, 38.61974398929089, 37.97314042095274, 37.9115793467882, 37.76835628918834, 37.49122332734656, 36.50783954925645, 36.46765270770563, 36.42101401410459, 34.65157396826377, 33.98338947575703, 33.935690559912885, 33.80903058704956, 33.53921927758842, 32.30330165071577] atol=3
# Test that all age-depth models are in stratigraphic order
@test all([issorted(x, rev=true) for x in eachcol(agedist)])
@test all(!isnan, agedist)
@test size(hiatusdist) == (nhiatuses, config.nsteps)
@test mean(hiatusdist, dims=2) ≈ [0.9662452348821302; 1.9517040131111487;;] atol=2
@test -Inf < mean(lldist) < 0
