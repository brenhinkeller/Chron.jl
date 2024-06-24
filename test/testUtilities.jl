## --- Test bilinear exponential function

    d = BilinearExponential([-5.9345103368858085, 66.00606672179812, 0.17739474265630253, 70.57882331291309, 0.6017142541555505])
    @test d isa BilinearExponential{Float64}
    @test pdf(d, 66.02) ≈ 0.00010243952455548608
    @test logpdf(d, 66.02) ≈ (-3.251727600957773 -5.9345103368858085)
    @test pdf(d, [66.02, 66.08, 66.15]) ≈ [0.00010243952455548608, 1.0372066408907184e-13, 1.0569878258022094e-28]

    @test eltype(d) === partype(d) === Float64
    @test params(d) == (-5.9345103368858085, 66.00606672179812, 0.17739474265630253, 70.57882331291309, 0.6017142541555505)
    @test d + 1 isa BilinearExponential{Float64}
    @test location(d + 1) == 66.00606672179812 + 1
    @test d * 2 isa BilinearExponential{Float64}
    @test location(d * 2) == 66.00606672179812 * 2
    @test scale(d * 2) == 0.17739474265630253 * 2


## --- Test Radiocarbon distribution

    d = Radiocarbon(1000, 10, intcal13)
    @test d isa Radiocarbon{Float64}
    @test pdf(d, 950) ≈ 0.0036210273644940918
    @test logpdf(d, 950) ≈ -3.4446721311475414 
    @test pdf(d, [900, 950, 1000]) ≈ [3.372067127114498e-7, 0.0036210273644940918, 1.5582187464154278e-12]

    @test eltype(d) === partype(d) === Float64
    @test location(d) ≈ 926.3094740428785 atol=0.1
    @test scale(d) ≈ 7.411299532288234 atol=0.1

## --- Other utility functions for log likelihoods

    @test strat_ll(66.02, BilinearExponential(-5.9345103368858085, 66.00606672179812, 0.17739474265630253, 70.57882331291309, 0.6017142541555505)) ≈  (-3.251727600957773 -5.9345103368858085)
    @test strat_ll(950, Radiocarbon(1000, 10, intcal13)) ≈ -3.4446721311475414 

    @test strat_ll([0.0, 0.0], [Normal(0,1), Normal(0,1)]) ≈ 0
    @test strat_ll([0.0, 0.5], [Normal(0,1), Uniform(0,1)]) ≈ -0.9189385332046728
    @test strat_ll([0.0, 0.5, 66.02], [Normal(0,1), Uniform(0,1), BilinearExponential(-5.9345103368858085, 66.00606672179812, 0.17739474265630253, 70.57882331291309, 0.6017142541555505)]) ≈ (-0.9189385332046728 -3.251727600957773 -5.9345103368858085)
    @test strat_ll([0.0, 0.5, 66.02, 900], [Normal(0,1), Uniform(0,1), BilinearExponential(-5.9345103368858085, 66.00606672179812, 0.17739474265630253, 70.57882331291309, 0.6017142541555505), Radiocarbon(1000, 10, intcal13)]) ≈ -22.831420814939655

## ---
