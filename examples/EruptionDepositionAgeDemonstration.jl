## --- Load (and install if necessary) the Chron.jl package

    using Chron
    using Plots

## --- Test Bayesian eruption age estimation with a synthetic dataset

    dt_sigma = 30; # Timescale relative to analytical uncertainty
    N = 42; # Number of zircons

    # Draw set of pseudorandom ages from MELTS volcanic zircon distribution,
    # with true minimum == 0 and analytical sigma == 1
    ages = draw_from_distribution(MeltsVolcanicZirconDistribution,N).*dt_sigma + randn(N);
    uncert = ones(N)

    # Calculate the weighted mean age
    (wx, wsigma, mswd) = awmean(ages, uncert)

    h = plot(xlabel="N",
        ylabel="Time before eruption (sigma)",
        fg_color_legend=:white,
        legend=:topleft,
        framestyle=:box
    )
    plot!(h,1:length(ages),sort(ages),yerror=uncert*2,seriestype=:scatter, label="Synthetic zircon ages")
    plot!(h,collect(xlims()), [0,0], label="Eruption age")
    plot!(h,collect(xlims()), [wx,wx], label="Weighted mean, MSWD $(round(mswd,sigdigits=2))")

## --- Calculate bootstrapped $\ \mathcal{\vec{f}}_{xtal}(t_r)$

    # Bootstrap the crystallization distribution,
    # accounting for any expected analytical "tail" beyond eruption/deposition
    dist = BootstrapCrystDistributionKDE(ages, uncert)
    dist ./= nanmean(dist) # Normalize

    # Plot bootstrapped distribution
    x = range(0,1.3,length=length(dist))
    plot(x, dist,
        label="bootstrapped",
        xlabel="Time before eruption (scaled)",
        ylabel="Probability Density",
        fg_color_legend=:white,
        legend=:bottomleft,
        framestyle=:box
    )
    plot!(range(0,1,length=100),MeltsVolcanicZirconDistribution,label="original")

## ---  Run MCMC to estimate eruption/deposition age distribution of synthetic dataset
    # Configure model
    nsteps = 200000 # Length of Markov chain
    burnin = 100000 # Number of steps to discard at beginning of Markov chain

    # Run MCMC
    tminDist = metropolis_min(nsteps,dist,ages,uncert; burnin)

    # Print results
    AgeEst = nanmean(tminDist)
    AgeEst_sigma = nanstd(tminDist)
    print("\nEstimated eruption age of synthetic dataset:\n $AgeEst +/- $(2*AgeEst_sigma) Ma (2σ)\n (True synthetic age 0 Ma)")

    # Plot results
    h = histogram(tminDist[burnin:end],nbins=50,label="Posterior distribution",xlabel="Eruption Age (Ma)",ylabel="N",framestyle=:box)
    plot!(h,[0,0],collect(ylims()),line=(3),label="True (synthetic) age",fg_color_legend=:white)
    plot!(h,[wx,wx],collect(ylims()),line=(3),label="Weighted mean age",legend=:topright)
    display(h)

## --- Estimate eruption age for real zircon data
    # The example dataset here is from [Wotzlaw et al., 2013](https://doi.org/10.1130/G34366.1) FCT+MLX
    # Age and one-sigma uncertainty.
    ages = [28.196, 28.206, 28.215, 28.224, 28.232, 28.241, 28.246, 28.289, 28.308, 28.332, 28.341, 28.359, 28.379, 28.383, 28.395, 28.4, 28.405, 28.413, 28.415, 28.418, 28.42, 28.422, 28.428, 28.452, 28.454, 28.454, 28.458, 28.468, 28.471, 28.475, 28.482, 28.485, 28.502, 28.52, 28.551, 28.561, 28.565, 28.582, 28.584, 28.586, 28.611, 28.638, 28.655]
    uncert = [0.019, 0.0155, 0.019, 0.0215, 0.018, 0.023, 0.013, 0.029, 0.0175, 0.0315, 0.0095, 0.0245, 0.0255, 0.0175, 0.0235, 0.014, 0.021, 0.022, 0.0125, 0.0135, 0.016, 0.0195, 0.0175, 0.0125, 0.01, 0.014, 0.015, 0.0205, 0.0155, 0.011, 0.0115, 0.0185, 0.0255, 0.014, 0.0125, 0.013, 0.015, 0.014, 0.012, 0.016, 0.0215, 0.0125, 0.0215]

    # Sort by age (just to make rank-order plots prettier)
    t = sortperm(ages)
    ages = ages[t];
    uncert = uncert[t];

## --- Calculate bootstrapped $\ \mathcal{\vec{f}}_{xtal}(t_r)$
    # Bootstrap the crystallization distribution,
    # accounting for any expected analytical "tail" beyond eruption/deposition
    BootstrappedDistribution = BootstrapCrystDistributionKDE(ages, uncert)
    BootstrappedDistribution ./= nanmean(BootstrappedDistribution) # Normalize

    # Plot bootstrapped distribution
    plot(range(0,1,length=length(BootstrappedDistribution)),BootstrappedDistribution,
        label="Bootstrapped f_xtal",
        ylabel="Probability Density",
        xlabel="Time before eruption (unitless)",
        fg_color_legend=:white,
        framestyle=:box
    )

## --  Run MCMC to estimate eruption age
    # Configure model
    nsteps = 4000000; # Length of Markov chain
    burnin = 150000; # Number of steps to discard at beginning of Markov chain


    # Choose the form of the prior closure/crystallization distribution to use
    dist = BootstrappedDistribution
    ## You might alternatively consider:
    # dist = UniformDistribution              # A reasonable default
    # dist = MeltsVolcanicZirconDistribution  # A single magmatic pulse, truncated by eruption
    # dist = ExponentialDistribution          # Applicable for survivorship processes, potentially including inheritance/dispersion in Ar-Ar dates

    # Run MCMC
    tminDist = metropolis_min(nsteps,dist,ages,uncert; burnin);

    # Print results
    AgeEst = nanmean(tminDist);
    AgeEst_sigma = nanstd(tminDist);
    ci = CI(tminDist)
    flush(stdout)
    print("\nEstimated eruption age:\n $AgeEst +/- $(2*AgeEst_sigma) Ma (2σ)\n $ci Ma (95% CI)")

    # Plot results
    h = histogram(tminDist,nbins=100,
        label="Posterior distribution",
        xlabel="Eruption Age (Ma)",
        ylabel="N",
        legend=:topleft,
        fg_color_legend=:white,
        framestyle=:box
    )
    # plot!(h,[wx,wx],collect(ylims()),line=(3),label="Weighted mean age",legend=:topleft)
    display(h)

## ---  Plot eruption age estimate relative to rank-order plot of raw zircon ages

    h = plot(xlabel="N",
        ylabel="Age (Ma)",
        legend=:topleft,
        fg_color_legend=:white,
        framestyle=:box
    )
    plot!(h,1:length(ages),ages,yerror=uncert*2,seriestype=:scatter, markerstrokecolor=:auto, label="Observed ages")
    plot!(h,[length(ages)],[AgeEst],yerror=2*AgeEst_sigma, markerstrokecolor=:auto, label="Bayesian eruption age estimate",color=:red)
    plot!(h,collect(xlims()),[AgeEst,AgeEst],color=:red, label="")

## ---
