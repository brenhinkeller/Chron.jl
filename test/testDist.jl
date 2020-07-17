
mu = collect(100:0.1:101); sigma = 0.01*ones(11);
@test Chron.check_dist_ll(MeltsVolcanicZirconDistribution, mu, sigma, 100,101) ≈ -3.6933357793356576
@test isapprox.(nanmean.(metropolis_minmax(2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5))[1:2], (99.92,101.08), atol=0.015) |> all
@test isapprox(nanmean(metropolis_min(2*10^5, MeltsVolcanicZirconDistribution, mu, sigma, burnin=10^5)), 99.9228, atol=0.015)
