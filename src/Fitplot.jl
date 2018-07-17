## --- Utility functions for plotting

    function plot_rankorder_errorbar(data,uncert; seriestype=:scatter,ylabel="",label="",xticks=[],xlabel="")
        sI = sortperm(data);
        h = plot(1:length(sI),data[sI],yerror=uncert[sI],seriestype=seriestype,
                 label=label,ylabel=ylabel,xlabel=xlabel,xticks=xticks)
        return h;
    end

## -- Bootstrap prior distribution shape

    # Bootstrap a KDE of the pre-eruptive (or pre-deposition) mineral crystallization
    # distribution shape from a 2-d array of sample ages using a KDE of stacked sample data
    function BootstrapCrystDistributionKDEfromStrat(smpl::StratAgeData);
        # Load all data points and scale from 0 to 1
        allscaled = Array{Float64}([]);
        for i=1:length(smpl.Name)
            data = readcsv(string(smpl.Path, smpl.Name[i], ".csv"))

            # Maximum extent of expected analytical tail (beyond eruption/deposition)
            maxTailLength = mean(data[:,2])/smpl.inputSigmaLevel * norm_quantile(1 - 1/(1+size(data,1)));
            included = (data[:,1]-minimum(data[:,1])) .>= maxTailLength;
            included .|= data[:,1] .> nanmedian(data[:,1]); # Don't exclude more than half (could only happen in underdispersed datasets)
            included .&= .~isnan.(data[:,1]); # Exclude NaNs

            # Include and scale only those data not within the expected analytical tail
            if sum(included)>0
                scaled = data[included,1]-minimum(data[included,1]);
                if maximum(scaled) > 0
                    scaled = scaled./maximum(scaled);
                end
                allscaled = [allscaled; scaled]
            end
        end

        # Calculate kernel density estimate, truncated at 0
        kd = kde(allscaled,npoints=2^7);
        t = kd.x .> -0.05; # Ensure sharp cutoff at eruption / deposition
        return kd.density[t];
    end

    # Bootstrap a KDE of the pre-eruptive (or pre-deposition) mineral crystallization
    # distribution shape from a 2-d array of sample ages using a KDE of stacked sample data
    function BootstrapCrystDistributionKDE(data::Array{Float64})
        # Load all data points and scale from 0 to 1
        allscaled = Array{Float64,1}();
        for i=1:size(data,2)
            scaled = data[:,i]-minimum(data[:,i]);
            if maximum(scaled) > 0
                scaled = scaled./maximum(scaled);
            end
            allscaled = [allscaled; scaled]
        end

        # Calculate kernel density estimate, truncated at 0
        kd = kde(allscaled,npoints=2^7);
        t = kd.x .> -0.05;
        return kd.density[t];
    end


## --- Fit and plot results from stationary distribution of depostion/eruption age distribution model

    # Process and fit stationary distribution for eruption age
    function tMinDistMetropolis(smpl::StratAgeData,nsteps::Int,burnin::Int,dist::Array{Float64});

        # Estimate the distribution for each sample
        print("Estimating eruption/deposition age distributions...\n");
        for i=1:length(smpl.Name)
            # Load data for each sample
            data = readcsv(string(smpl.Path, smpl.Name[i], ".csv"))
            print(i, ": ", smpl.Name[i], "\n"); # Display progress

            # Run MCMC to estimate saturation and eruption/deposition age distributions
            (tminDist, tmaxDist, llDist, acceptanceDist) = crystMinMaxMetropolis(nsteps,dist,data[:,1],data[:,2]/smpl.inputSigmaLevel);

            # Fill in the strat sample object with our new results
            smpl.Age[i] = mean(tminDist[burnin:end]);
            smpl.Age_Sigma[i] = std(tminDist[burnin:end]);
            smpl.Age_025CI[i] = percentile(tminDist[burnin:end],2.5);
            smpl.Age_975CI[i] = percentile(tminDist[burnin:end],97.5);

            # Fit custom many-parametric distribution function to histogram
            edges = linspace(minimum(tminDist[burnin:end]),maximum(tminDist[burnin:end]),101); # Vector of bin edges
            hobj = fit(Histogram,tminDist[burnin:end],edges,closed=:left) # Fit histogram object

            t = hobj.weights.>0; # Only look at bins with one or more results
            N = hobj.weights[t] ./ length(tminDist[burnin:end]) .* length(t); # Normalized number of MCMC steps per bin
            bincenters = cntr(hobj.edges[1])[t]; # Vector of bin centers

            # Initial guess for parameters
            p = ones(5);
            p[1] = maximum(N);
            p[2] = mean(tminDist[burnin:end]);
            p[3] = std(tminDist[burnin:end]);

            # Fit nonlinear model
            fobj = curve_fit(bilinear_exponential,bincenters,N,p);
            smpl.Params[:,i] = fobj.param;

            # Rank-order plot of analyses and eruption/deposition age range
            nAnalyses = length(data[:,1]);
            h1 = plot_rankorder_errorbar(data[:,1],2*data[:,2]/smpl.inputSigmaLevel,ylabel="Age (Ma)",label="Ages")
            m = ones(nAnalyses).*smpl.Age[i];
            l = ones(nAnalyses).*smpl.Age_025CI[i];
            u = ones(nAnalyses).*smpl.Age_975CI[i];
            plot!(h1,1:nAnalyses,l,fillto=u,fillalpha=0.6,linealpha=0, label="Model interpretation");
            plot!(h1,1:nAnalyses,m,linecolor=:black,linestyle=:dot,label="",legend=:topleft);
            savefig(h1,string(smpl.Path,smpl.Name[i],"_rankorder.pdf"));

            # Plot model fit to histogram
            h2 = plot(bincenters,N,label="Histogram");
            plot!(h2,bincenters,bilinear_exponential(bincenters,smpl.Params[:,i]),label="Curve fit");
            plot!(h2,legend=:topleft,xlabel="Age",ylabel="Likelihood");
            savefig(h2,string(smpl.Path,smpl.Name[i],"_distribution.pdf"));

        end

        return smpl
    end

## ---

    # Process and fit stationary distribution for eruption age
    function tMinDistMetropolisLA(smpl::StratAgeData,nsteps::Int,burnin::Int,dist::Array{Float64});

        # Estimate the distribution for each sample
        print("Estimating eruption/deposition age distributions...\n");
        for i=1:length(smpl.Name)
            # Load data for each sample
            data = readcsv(string(smpl.Path, smpl.Name[i], ".csv"))
            print(i, ": ", smpl.Name[i], "\n"); # Display progress

            # Run MCMC to estimate saturation and eruption/deposition age distributions
            (tminDist, tmaxDist, llDist, acceptanceDist) = crystMinMaxMetropolisLA(nsteps,dist,data[:,1],data[:,2]/smpl.inputSigmaLevel);

            # Fill in the strat sample object with our new results
            smpl.Age[i] = mean(tminDist[burnin:end]);
            smpl.Age_Sigma[i] = std(tminDist[burnin:end]);
            smpl.Age_025CI[i] = percentile(tminDist[burnin:end],2.5);
            smpl.Age_975CI[i] = percentile(tminDist[burnin:end],97.5);

            # Fit custom many-parametric distribution function to histogram
            edges = linspace(minimum(tminDist[burnin:end]),maximum(tminDist[burnin:end]),101); # Vector of bin edges
            hobj = fit(Histogram,tminDist[burnin:end],edges,closed=:left) # Fit histogram object

            t = hobj.weights.>0; # Only look at bins with one or more results
            N = hobj.weights[t] ./ length(tminDist[burnin:end]) .* length(t); # Normalized number of MCMC steps per bin
            bincenters = cntr(hobj.edges[1])[t]; # Vector of bin centers

            # Initial guess for parameters
            p = ones(5);
            p[1] = maximum(N);
            p[2] = mean(tminDist[burnin:end]);
            p[3] = std(tminDist[burnin:end]);

            # Fit nonlinear model
            fobj = curve_fit(bilinear_exponential,bincenters,N,p);
            smpl.Params[:,i] = fobj.param;

            # Rank-order plot of analyses and eruption/deposition age range
            nAnalyses = length(data[:,1]);
            h1 = plot_rankorder_errorbar(data[:,1],2*data[:,2]/smpl.inputSigmaLevel,ylabel="Age (Ma)",label="Ages")
            m = ones(nAnalyses).*smpl.Age[i];
            l = m - 2*smpl.Age_Sigma[i];
            u = m + 2*smpl.Age_Sigma[i];
            plot!(h1,1:nAnalyses,l,fillto=u,fillalpha=0.6,linealpha=0, label="Model interpretation");
            plot!(h1,1:nAnalyses,m,linecolor=:black,linestyle=:dot,label="");
            savefig(h1,string(smpl.Path,smpl.Name[i],"_rankorder.pdf"));

            # Plot model fit to histogram
            h2 = plot(bincenters,N,label="Histogram");
            plot!(h2,bincenters, bilinear_exponential(bincenters,smpl.Params[:,i]),label="Curve fit");
            plot!(h2,legend=:topleft,xlabel="Age",ylabel="Likelihood");
            savefig(h2,string(smpl.Path,smpl.Name[i],"_distribution.pdf"));

        end

        return smpl
    end
