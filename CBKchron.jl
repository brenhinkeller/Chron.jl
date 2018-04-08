# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                 CBKchron.jl                                   #
#                                                                               #
#       A two-part framework for (1) estimating eruption/deposition age         #
#  distributions from complex mineral age spectra and (2) subsequently building #
#  a stratigraphic age model based on those distributions. Each step relies on  #
#  a Markov-Chain Monte Carlo model.                                            #
#                                                                               #
#    The first model uses an informative prior distribution to estimate the     #
#  times of first (i.e., saturation) and last  mineral crystallization (i.e.,   #
#  eruption/deposition).                                                        #
#                                                                               #
#    The second model uses the estimated (posterior) eruption/deposition ages   #
#  distributions along with the constraint of stratigraphic superposition to    #
#  produce an age-depth model                                                   #
#                                                                               #
#       This file contains ..., while ... contains the ...                      #
#                                                                               #
#   Last modified by C. Brenhin Keller 2018-03-19                               #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## --- Load external resources

    # Basic statistics and UI resources
    using StatsBase, ProgressMeter;
    # Weighted mean, etc
    include("Utilities.jl");
    # Functions for estimating extrema of a finite-range distribution
    include("DistMetropolis.jl");
    # Functions for stratigraphic modelling
    include("StratMetropolis.jl");
    # External packages for plotting, etc.
    using Plots, Colors;
    gr();

## --- Define sample properties

    # # # # # # # # # # # Enter sample information here!  # # # # # # # # # # #
    Name   =        ("KJ08-157", "KJ04-75", "KJ09-66", "KJ04-72", "KJ04-70");
    Height =        [     -52.0,      44.0,      54.0,      82.0,      93.0];
    Height_Sigma =  [       3.0,       1.0,       3.0,       3.0,       3.0];
    InputAgeSigmaLevel = 2;
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Count the names to see how many samples we have
    nSamples = length(Name);
    # Make an instance of a ChronSection object with
    smpl = StratAgeData(
        Name,  # Sample names
        Height,  # Sample heights
        Height_Sigma,  # Height_sigma
        fill(NaN,nSamples),  # Sample ages
        fill(NaN,nSamples),  # Sample age uncertainty
        fill(NaN,nSamples),  # Sample age 2.5% CI
        fill(NaN,nSamples),  # Sample age 97.5% CI
        zeros(nSamples), # Sidedness (zeros by default, geochron constraints are two-sided)
        fill(NaN,nSamples,nSamples), # Sample age distribution parameters
        );

## --- Calculate bootstrapped distribution

    # Load all data points and scale from 0 to 1
    allscaled = Array{Float64}([]);
    for i=1:length(smpl.Name)
        data = readcsv(string("examples/DenverUPbExample/", smpl.Name[i], ".csv"))
        scaled = data[:,1]-minimum(data[:,1]);
        scaled = scaled./maximum(scaled);
        allscaled = [allscaled; scaled]
    end

    # Calculate kernel density estimate, truncated at 0
    using KernelDensity
    kd = kde(allscaled,npoints=2^7);
    t = kd.x.>0;
    BootstrappedDistribution = kd.density[t];
    plot(BootstrappedDistribution)

## --- Estimate the eruption age distributions for each sample

    # Number of steps in distribution MCMC
    distSteps = 10^6;
    distBurnin = floor(Int,distSteps/2);

    # Choose a distribution to use.
    # A variety of potentially useful distributions are provided in
    # DistMetropolis.jl and include UniformDisribution, TriangularDistribution,
    # BootstrappedDistribution, and MeltsVolcanicZirconDistribution, etc.
    # or you can define your own.
    dist = BootstrappedDistribution;

    # Estimate the distribution for each sample
    print("Estimating eruption/deposition age distributions...\n");
    using LsqFit
    for i=1:nSamples
      # Load data for each sample
        data = readcsv(string("examples/DenverUPbExample/", smpl.Name[i], ".csv"))
        print(i, ": ", smpl.Name[i], "\n"); # Display progress

      # Run MCMC to estimate saturation and eruption/deposition age distributions
        (tminDist, tmaxDist, llDist, acceptanceDist) = crystMinMaxMetropolis(distSteps,dist,data[:,1],data[:,2]/InputAgeSigmaLevel);

        # Fill in the strat sample object with our new results
        smpl.Age[i] = mean(tminDist[distBurnin:end]);
        smpl.Age_Sigma[i] = std(tminDist[distBurnin:end]);
        smpl.Age_025CI[i] = percentile(tminDist[distBurnin:end],2.5);
        smpl.Age_975CI[i] = percentile(tminDist[distBurnin:end],97.5);

        # Rank-order plot of analyses and eruption/deposition age range
        nAnalyses = length(data[:,1]);
        h1 = plotRankOrderErrorbar(data[:,1],2*data[:,2],ylabel="Age (Ma)",label="Ages")
        m = ones(nAnalyses).*smpl.Age[i];
        l = m - 2*smpl.Age_Sigma[i];
        u = m + 2*smpl.Age_Sigma[i];
        plot!(h1,1:nAnalyses,l,fillto=u,fillalpha=0.6,linealpha=0, label="Model interpretation");
        plot!(h1,1:nAnalyses,m,linecolor=:black,linestyle=:dot,label="");
        savefig(h1,string(smpl.Name[i],"_rankorder.pdf"));

      # Fit custom many-parametric distribution function to histogram
        edges = linspace(minimum(tminDist[distBurnin:end]),maximum(tminDist[distBurnin:end]),101); # Vector of bin edges
        hobj = fit(Histogram,tminDist[distBurnin:end],edges,closed=:left) # Fit histogram object

        t = hobj.weights.>0; # Only look at bins with one or more results
        N = hobj.weights[t] ./ length(tminDist[distBurnin:end]) .* length(t); # Normalized number of MCMC steps per bin
        bincenters = cntr(hobj.edges[1])[t]; # Vector of bin centers

        # Initial guess for parameters
        p = ones(5);
        p[1] = maximum(N);
        p[2] = mean(tminDist[distBurnin:end]);
        p[3] = std(tminDist[distBurnin:end]);

        # Fit nonlinear model
        fobj = curve_fit(doubleLinearExponential,bincenters,N,p);

        # Plot model fit to histogram
        h2 = plot(bincenters,N,label="Histogram");
        plot!(h2,bincenters,doubleLinearExponential(bincenters,fobj.param),label="Curve fit");
        plot!(h2,legend=:topleft,xlabel="Age",ylabel="Likelihood");
        savefig(h2,string(smpl.Name[i],"_distribution.pdf"));

        smpl.p[:,i] = fobj.param;
    end

    using JLD
    @save "smpl.jld" smpl

## --- Run stratigraphic model

# # # # # # # # # # # Configure stratigraphic model here! # # # # # # # # # # #
  resolution = 2.0; # Same units as sample height. Smaller is slower!
  (bottom, top) = extrema(smpl.Height);
  bounding = 0.1;
  npoints_approx = round(Int,length(bottom:resolution:top) * (1+2*bounding))
  burnin = 20000*npoints_approx;
  nsteps = 25000;
  sieve = round(Int,npoints_approx);
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Pass configuration to a struct
    config = StratAgeModelConfiguration(
        resolution, # Resolution
        burnin, # Lenght of burnin
        nsteps, # Number of regular MC steps
        sieve, # Seive strength
    );

    # Run the model
    (mdl, agedist, lldist) = StratMetropolis(smpl, config);

    # Plot results
    hdl = plot([mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(minimum(mdl.Height),0.5,:blue), label="model")
    plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="")
    plot!(hdl, smpl.Age, smpl.Height, xerror=smpl.Age_Sigma*2,label="data",seriestype=:scatter,color=:black)
    plot!(hdl, xlabel="Age (Ma)", ylabel="Height (cm)")
    savefig(hdl,"AgeDepthModel.pdf");

## ---

# A type of object to hold data about hiatuses
hiatus = HiatusData(
    [20.0,], # Height
    [0.0,], # Height_Sigma
    [0.5,], # Duration
    [0.1,], # Duration_Sigma
);

(mdl, agedist, hiatusdist, lldist) = StratMetropolisHiatus(smpl, hiatus, config);
