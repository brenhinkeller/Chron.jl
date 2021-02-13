# Summary statistics: simple cases, Float64
    A = [1:10.0..., NaN]
    @test nanmean(A) == 5.5
    @test nanrange(A) == 9.0
    @test nanminimum(A) == 1.0
    @test nanmaximum(A) == 10.0
    @test nanextrema(A) == (1.0, 10.0)
    @test nanstd([1,2,3,NaN]) == 1.0
    @test nanmedian([1,2,3,NaN]) == 2.0
    @test pctile([0:100...,NaN],99) == 99.0

# Arrays containing only NaNs should yield NaN
    A = fill(NaN,10)
    @test isnan(nanmean(A))
    @test isnan(nanrange(A))
    @test isnan(nanminimum(A))
    @test isnan(nanmaximum(A))
    @test all(isnan.(nanextrema(A)))
    @test isnan(nanstd(A))

# Summary statistics: simple cases, Int64
    A = collect(1:10)
    @test nanmean(A) == 5.5
    @test nanrange(A) == 9.0
    @test nanminimum(A) == 1.0
    @test nanmaximum(A) == 10.0
    @test nanextrema(A) == (1.0, 10.0)
    @test nanstd([1,2,3]) == 1.0
    @test nanmedian([1,2,3]) == 2.0
    @test pctile([0:100...],99) == 99.0

# Summary statistics: simple cases, ranges
    A = 1:10
    @test nanmean(A) == 5.5
    @test nanrange(A) == 9.0
    @test nanminimum(A) == 1.0
    @test nanmaximum(A) == 10.0
    @test nanextrema(A) == (1.0, 10.0)
    @test nanstd(1:3) == 1.0
    @test nanmedian(1:3) == 2.0
    @test pctile(0:100,99) == 99.0

    A = 1:10.
    @test nanmean(A) == 5.5
    @test nanrange(A) == 9.0
    @test nanminimum(A) == 1.0
    @test nanmaximum(A) == 10.0
    @test nanextrema(A) == (1.0, 10.0)
    @test nanstd(1:3.) == 1.0
    @test nanmedian(1:3.) == 2.0
    @test pctile(0:100.,99) == 99.0

# Summary statistics: dimensional tests, Int64
    A = reshape(1:300,100,3)
    @test nanminimum(A, dims=1) == minimum(A, dims=1)
    @test nanminimum(A, dims=2) == minimum(A, dims=2)
    @test nanmaximum(A, dims=1) == maximum(A, dims=1)
    @test nanmaximum(A, dims=2) == maximum(A, dims=2)
    @test nanmean(A, dims=1) == mean(A, dims=1)
    @test nanmean(A, dims=2) == mean(A, dims=2)
    @test nanstd(A, dims=1) ≈ std(A, dims=1)
    @test nanstd(A, dims=2) ≈ std(A, dims=2)
    @test nanmedian(A, dims=1) == median(A, dims=1)
    @test nanmedian(A, dims=2) == median(A, dims=2)

# Summary statistics: dimensional tests, Float64
    A = reshape(1:300.,100,3)
    @test nanminimum(A, dims=1) == minimum(A, dims=1)
    @test nanminimum(A, dims=2) == minimum(A, dims=2)
    @test nanmaximum(A, dims=1) == maximum(A, dims=1)
    @test nanmaximum(A, dims=2) == maximum(A, dims=2)
    @test nanmean(A, dims=1) == mean(A, dims=1)
    @test nanmean(A, dims=2) == mean(A, dims=2)
    @test nanstd(A, dims=1) ≈ std(A, dims=1)
    @test nanstd(A, dims=2) ≈ std(A, dims=2)
    @test nanmedian(A, dims=1) == median(A, dims=1)
    @test nanmedian(A, dims=2) == median(A, dims=2)

# Interpolation
    @test linterp1(1:10,21:30,5:0.5:6) == [25.0, 25.5, 26.0]
    @test linterp1s(10:-1:1,21:30,5:0.5:6) == [26.0, 25.5, 25.0]

# Matching
    @test findclosest(3.3:5.3,1:10) == 3:5
    @test findclosestbelow(3.3:5.3,1:10) == 3:5
    @test findclosestabove(3.3:5.3,1:10) == 4:6

# Distributions
    @test normpdf.(0, 1,-1:1) ≈ [0.24197072451914337, 0.3989422804014327, 0.24197072451914337]
    @test normpdf.(1:10, 1:10, 1:10) ≈ normpdf(collect.((1:10, 1:10, 1:10))...)

    @test normpdf_ll.(0,1,-5:5) == -(-5:5).^2/2
    r = collect(-5:5)
    @test normpdf_ll(0,1,r) == normpdf_ll(0,ones(11),r) == normpdf_ll(zeros(11),ones(11),r) == sum(normpdf_ll.(0,1,r))
    @test normpdf_ll(ones(10),1,collect(1:10)) == normpdf_ll(collect(1:10),1,ones(10))

    @test normcdf(1,1,1) == 0.5
    @test normcdf(0,1,-2:2) ≈ normcdf.(0,1,-2:2) ≈ [0.02275013194817921, 0.15865525393145707, 0.5, 0.8413447460685429, 0.9772498680518208]
    @test normcdf.(1:10, 1:10, 1:10) == normcdf(collect.((1:10, 1:10, 1:10))...) == fill(0.5, 10)

    @test normproduct(0,1,0,1) === normpdf(0,sqrt(2),0) === 0.28209479177387814
    @test normproduct_ll(0,1,0,1) === normpdf_ll(0,1,0) === 0.0

    @test [-2,0,2] ≈ norm_quantile.([0.022750131948, 0.5, 0.977249868052])
    @test norm_quantile.(0:0.25:1) ≈ [-Inf, -0.6744897501960818, 0.0, 0.6744897501960818, Inf]

    @test isapprox(norm_width(390682215445)/2, 7, atol=1e-5)

    A = draw_from_distribution(ones(100), 10000)::AbstractArray
    @test length(A) == 10000
    @test isapprox(mean(A), 0.5, atol=0.08)
    @test isapprox(std(A), sqrt(1/12), atol=0.08)


# Weighted means

    # Plain vanilla weighted mean
    @test awmean([0,1,2,3],[1,1,1,1]) == (1.5, 0.5, 5/3)
    @test awmean(1:10,ones(10)) == (5.5, 0.31622776601683794, 9.166666666666666)

    # Weighted mean with geochronogists' MSWD-correction to uncertainty
    @test gwmean([0,1,2,3],[1,1,1,1]) == (1.5, 0.6454972243679028, 5/3)
    @test gwmean(1:10,ones(10)) == (5.5, 0.9574271077563381, 9.166666666666666)

    # Mean Square of Weighted Deviates (aka reduced chi-squared)
    @test Chron.MSWD(0:10, ones(11)) == 11
    @test Chron.MSWD([0,1,2],[1,1,1]) == 1.0
