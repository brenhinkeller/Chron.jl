## --- Test simple metropolis_min and metropolis_minmax functions

    mu, sigma = collect(100:0.1:101), 0.01*ones(11);
    dist = MeltsVolcanicZirconDistribution ./ vmean(MeltsVolcanicZirconDistribution)
    @test Chron.Isoplot.dist_ll(dist, mu, sigma, 100,101) ≈ -3.6933357793356576

    tmindist = metropolis_min(2*10^5, dist, mu, sigma, burnin=10^5)
    @test isapprox(nanmean(tmindist), 99.9228, atol=0.015)

    tmindist, tmaxdist, lldist, acceptancedist = metropolis_minmax(2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)
    @test isapprox(nanmean(tmindist), 99.92, atol=0.015)
    @test isapprox(nanmean(tmaxdist), 101.08, atol=0.015)

## -- Test distribution bootstrapping

    dist = BootstrapCrystDistributionKDE(1:10)
    @test all(isapprox.(dist, [0.465404, 0.503149, 0.540308, 0.576485, 0.611314, 0.644472, 0.675691, 0.704756, 0.731516, 0.75588, 0.777815, 0.797345, 0.81454, 0.82951, 0.842398, 0.85337, 0.862605, 0.87029, 0.876611, 0.881749, 0.885873, 0.889136, 0.891676, 0.893612, 0.895043, 0.896047, 0.896684, 0.896993, 0.896993, 0.896684, 0.896047, 0.895043, 0.893612, 0.891676, 0.889136, 0.885873, 0.881749, 0.876611, 0.87029, 0.862605, 0.85337, 0.842398, 0.82951, 0.81454, 0.797345, 0.777815, 0.75588, 0.731516, 0.704756, 0.675691, 0.644472, 0.611314, 0.576485, 0.540308, 0.503149, 0.465404, 0.427489, 0.389826, 0.352829, 0.316888, 0.282361, 0.249559, 0.218742, 0.19011, 0.163804, 0.139901, 0.118424, 0.099338, 0.082567, 0.067992, 0.055465, 0.044819, 0.03587, 0.028432, 0.022318, 0.017347, 0.01335, 0.010173, 0.007675, 0.005732, 0.004238, 0.003102, 0.002247, 0.001612, 0.001145, 0.000806, 0.000563, 0.000392, 0.000275, 0.000197, 0.00015, 0.000128], atol=0.00001))

    dist = BootstrapCrystDistributionKDE(1:10, ones(10))
    @test all(isapprox.(dist, [0.500228, 0.535143, 0.569056, 0.601647, 0.632641, 0.661805, 0.688961, 0.713982, 0.736794, 0.757375, 0.775748, 0.79198, 0.806169, 0.818443, 0.828947, 0.837842, 0.845291, 0.851459, 0.856502, 0.86057, 0.863796, 0.866299, 0.868179, 0.869517, 0.870375, 0.870794, 0.870794, 0.870375, 0.869517, 0.868179, 0.866299, 0.863796, 0.86057, 0.856502, 0.851459, 0.845291, 0.837842, 0.828947, 0.818443, 0.806169, 0.79198, 0.775748, 0.757375, 0.736794, 0.713982, 0.688961, 0.661805, 0.632641, 0.601647, 0.569056, 0.535143, 0.500228, 0.464656, 0.428798, 0.39303, 0.357728, 0.323251, 0.289935, 0.258079, 0.227938, 0.199721, 0.173583, 0.149624, 0.127894, 0.108393, 0.091075, 0.075858, 0.062627, 0.051244, 0.041553, 0.03339, 0.026585, 0.020972, 0.016391, 0.012691, 0.009734, 0.007396, 0.005566, 0.004149, 0.003064, 0.002241, 0.001623, 0.001166, 0.00083, 0.000587, 0.000415, 0.000295, 0.000215, 0.000166, 0.000143], atol=0.00001))

## -- Test systematic uncertainty propagation

    dist = 100 .+ randn(100000)

    dist_upb = add_systematic_uncert_UPb(dist)
    @test std(dist_upb) > std(dist)
    @test isapprox(std(dist_upb) - std(dist), 0.0015, atol=0.0006)

    dist_ar = add_systematic_uncert_ArAr(dist)
    @test std(dist_ar) > std(dist)
    @test isapprox(std(dist_ar) - std(dist), 0.005, atol=0.002)

    dist = (200 .+ randn(100000)) / 1000 # ~100 Ka

    dist_uth = add_systematic_uncert_UTh(dist)
    @test std(dist_uth) > std(dist)
    @test isapprox(std(dist_uth) - std(dist), 1e-4, atol=1e-4)


    dist2d = 100 .+ randn(1000,1000)

    dist_upb = add_systematic_uncert_UPb(dist2d)
    # Check that uncertainty is added in vertical bands
    diff = dist_upb .- dist2d
    @test all(std(diff,dims=2) .> std(diff,dims=1)')

    dist_ar = add_systematic_uncert_ArAr(dist2d)
    # Check that uncertainty is added in vertical bands
    diff = dist_ar .- dist2d
    @test all(std(diff,dims=2) .> std(diff,dims=1)')

    dist2d = (200 .+ randn(1000,1000)) ./ 1000 # ~200 Ka

    dist_uth = add_systematic_uncert_UTh(dist2d)
    # Check that uncertainty is added in vertical bands
    diff = dist_uth .- dist2d
    @test all(std(diff,dims=2) .> std(diff,dims=1)')


## ---
