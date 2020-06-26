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
