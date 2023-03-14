## -- Generate dataset
mu = [8; 10:0.1:11; 12; 15]
nAnalyses = length(mu)
sigma = 0.05 * ones(nAnalyses)
plot_rankorder_errorbar(mu, sigma)

## Estimate eruption age
tmin_dist = metropolis_min(2*10^5, reverse(TriangularDistribution), mu, sigma, burnin=10^5)
tmin_sigma = nanstd(tmin_dist)
tmin_mu = nanmean(tmin_dist)

## --
using Plots
m = ones(nAnalyses).*tmin_mu
l = ones(nAnalyses).*(tmin_mu-tmin_sigma)
u = ones(nAnalyses).*(tmin_mu+tmin_sigma)
plot!(1:nAnalyses,l,fillto=u,fillalpha=0.6,linealpha=0, label="", framestyle=:box)
