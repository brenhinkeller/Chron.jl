## --- Test bilinear exponential function

    @test bilinear_exponential(66.02, [-5.9345103368858085, 66.00606672179812, 0.17739474265630253, 70.57882331291309, 0.6017142541555505]) ≈ 0.00010243952455548608
    @test bilinear_exponential([66.02, 66.08, 66.15], [-5.9345103368858085, 66.00606672179812, 0.17739474265630253, 70.57882331291309, 0.6017142541555505]) ≈ [0.00010243952455548608, 1.0372066408907184e-13, 1.0569878258022094e-28]
    @test bilinear_exponential_ll(66.02, [-5.9345103368858085, 66.00606672179812, 0.17739474265630253, 70.57882331291309, 0.6017142541555505]) ≈ -3.251727600957773

    d = Chron.BilinearExponential([-5.9345103368858085, 66.00606672179812, 0.17739474265630253, 70.57882331291309, 0.6017142541555505])
    @test d isa Chron.BilinearExponential{Float64}
    @test pdf(d, 66.02) ≈ 0.00010243952455548608
    @test logpdf(d, 66.02) ≈ (-3.251727600957773 -5.9345103368858085)

    @test eltype(d) === partype(d) === Float64
    @test params(d) == (-5.9345103368858085, 66.00606672179812, 0.17739474265630253, 70.57882331291309, 0.6017142541555505)
    @test d + 1 isa Chron.BilinearExponential{Float64}
    @test location(d + 1) == 66.00606672179812 + 1
    @test d * 2 isa Chron.BilinearExponential{Float64}
    @test location(d * 2) == 66.00606672179812 * 2
    @test scale(d * 2) == 0.17739474265630253 * 2

## --- Other utility functions for log likelihoods

    @test Chron.strat_ll([0.0, 0.0], [Normal(0,1), Normal(0,1)]) ≈ 0
    @test Chron.strat_ll([0.0, 0.5], [Normal(0,1), Uniform(0,1)]) ≈ -0.9189385332046728

## ---
