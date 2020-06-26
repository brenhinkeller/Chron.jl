@test all(isapprox.(nanmean.(metropolis_minmax(10^6, MeltsVolcanicZirconDistribution, collect(100:0.1:101), 0.01*ones(11), burnin=10^5))[1:2], (99.92,101.08), atol=0.01))
