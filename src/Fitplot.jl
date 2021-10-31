## --- Utility functions for plotting

    function plot_rankorder_errorbar(data,uncert;
            scale=1,
            seriestype=:scatter,
            ylabel="",
            label="",
            xticks=[],
            xlabel="",
            color=:auto,
            markersize=2
        )
        h = plot()
        i0 = 0
        plot_rankorder_errorbar!(h,data,uncert,i0;
                scale=scale,
                seriestype=seriestype,
                label=label,
                ylabel=ylabel,
                xlabel=xlabel,
                xticks=xticks,
                color=color,
                markersize=markersize
        )
        return h
    end

    function plot_rankorder_errorbar!(h::Plots.Plot,data,uncert,i0;
            scale=1,
            seriestype=:scatter,
            label="",
            ylabel="",
            xlabel="",
            xticks=[],
            color=:auto,
            markersize=2
        )
        sI = sortperm(data)
        plot!(h, i0 .+ (1:length(sI))*scale,data[sI],
            yerror=uncert[sI],
            seriestype=seriestype,
            markerstrokecolor=color,
            color=color,
            label=label,
            ylabel=ylabel,
            xlabel=xlabel,
            xticks=xticks,
            markersize=markersize
        )
    end


## --- Remove outliers

    function screen_outliers(smpl::ChronAgeData; maxgap=100, make_plots=true)
        # Variables from struct
        Name = collect(smpl.Name)::Array{String,1}
        Path = smpl.Path::String
        Age_Unit = smpl.Age_Unit::String

        # Create and populate directory of screened output
        screenedpath = joinpath(Path, "screened/")
        mkpath(screenedpath)
        for i=1:length(Name)
            # With screening
            # Maximum offset before cutoff
            data = readdlm(joinpath(Path, Name[i]*".csv"), ',', Float64)::Array{Float64,2}
            data = sortslices(data,dims=1) # Sort ages in ascending order
            nAnalyses = size(data,1)
            maxdt_sigma = maxgap*norm_width(nAnalyses)/nAnalyses

            # Filter data to exclude outliers
            sigma_mutual = nanmean(data[:,2]) / smpl.inputSigmaLevel * sqrt(2)
            for j=nAnalyses:-1:2
                dt_sigma = abs(data[j,1]-data[j-1,1]) / sigma_mutual # Time gap divided by relative sigma

                # If we exceed the maximum allowed dt/sigma, delete any points
                # below (older than) the gap
                if dt_sigma>maxdt_sigma && j>2
                    data=data[1:j-1,:]
                end
            end
            if make_plots
                # Rank-order plot of all individual ages for comparison
                hdl = plot(1:nAnalyses,data[:,1],yerror=data[:,2]*2/smpl.inputSigmaLevel, seriestype=:scatter, color=:red, markerstrokecolor=:red,label="rejected",legend=:topleft,framestyle=:box,fg_color_legend=:white)
                plot!(hdl, 1:size(data,1),data[:,1],yerror=data[:,2]*2/smpl.inputSigmaLevel, seriestype=:scatter, color=:blue,markerstrokecolor=:blue,label="included",xlabel="N",ylabel="Age ($(Age_Unit))")
                savefig(hdl, joinpath(screenedpath, Name[i]*"_screening.pdf"))
            end
            writedlm(joinpath(screenedpath, Name[i]*".csv"), data, ',')
        end
        smpl.Path = screenedpath
        return smpl
    end

## --- Bootstrap prior distribution shape

    """
    ```julia
    BootstrapCrystDistributionKDE(smpl::ChronAgeData)
    ```
    Bootstrap an estimate of the pre-eruptive (or pre-deposition) mineral
    crystallization distribution shape from a Chron.ChronAgeData object containing
    data for several samples, using a kernel density estimate of stacked sample data.
    """
    function BootstrapCrystDistributionKDE(smpl::ChronAgeData; cutoff::Number=-0.05)
        # Extact variables froms struct
        Name = collect(smpl.Name)::Array{String,1}
        Path = smpl.Path::String
        DistType = smpl.Age_DistType::Array{Float64,1}

        # Load all data points and scale from 0 to 1
        allscaled = Array{Float64}([])
        for i=1:length(Name)
            if DistType[i]==0
                data = readdlm(joinpath(Path, Name[i]*".csv"), ',', Float64)::Array{Float64,2}

                # Maximum extent of expected analytical tail (beyond eruption/deposition)
                maxTailLength = mean(data[:,2]) ./ smpl.inputSigmaLevel .* norm_quantile(1 - 1/(1+size(data,1)))
                included = (data[:,1] .- minimum(data[:,1])) .>= maxTailLength
                included .|= data[:,1] .> nanmedian(data[:,1]) # Don't exclude more than half (could only happen in underdispersed datasets)
                included .&= .~isnan.(data[:,1]) # Exclude NaNs

                # Include and scale only those data not within the expected analytical tail
                if sum(included)>0
                    scaled = data[included,1] .- minimum(data[included,1])
                    if maximum(scaled) > 0
                        scaled = scaled ./ maximum(scaled)
                    end
                    append!(allscaled, scaled)
                end
            end
        end

        # Calculate kernel density estimate, truncated at 0
        kd = kde(allscaled,npoints=2^7)
        t = kd.x .> cutoff # Ensure sharp cutoff at eruption / deposition
        return kd.density[t]
    end


    """
    ```julia
    BootstrapCrystDistributionKDE(data::AbstractArray, [sigma::AbstractArray])
    ```
    Bootstrap an estimate of thq;e pre-eruptive (or pre-deposition) mineral
    crystallization distribution shape from a 2-d array of sample ages (one row per
    sample, one column per datum, padded with NaNs as needed) and an equivalent-size
    array of one-sigma uncertainties,
    using a kernel density estimate of stacked sample data.
    """
    function BootstrapCrystDistributionKDE(data::AbstractArray{<:Number}; cutoff::Number=-0.05)
        # Load all data points and scale from 0 to 1
        allscaled = Array{float(eltype(data)),1}([])
        for i=1:size(data,2)
            scaled = data[:,i] .- minimum(data[:,i])
            if maximum(scaled) > 0
                scaled = scaled ./ maximum(scaled)
            end
            append!(allscaled, scaled)
        end

        # Calculate kernel density estimate, truncated at 0
        kd = kde(allscaled,npoints=2^7)
        t = kd.x .> cutoff
        return kd.density[t]
    end

    function BootstrapCrystDistributionKDE(data::AbstractArray{<:Number}, sigma::AbstractArray{<:Number}; cutoff::Number=-0.05)
        # Array to hold stacked, scaled data
        allscaled = Array{float(eltype(data)),1}([])

        # For each row of data
        for i=1:size(data,2)

            # Maximum extent of expected analytical tail (beyond eruption/deposition/cutoff)
            maxTailLength = nanmean(sigma[:,i]) .* norm_quantile(1 - 1/(1+count(.!isnan.(data[:,i]))))
            included = (data[:,i] .- nanminimum(data[:,i])) .>= maxTailLength
            included .|= data[:,i] .> nanmedian(data[:,i]) # Don't exclude more than half (could only happen in underdispersed datasets)
            included .&= .~isnan.(data[:,i]) # Exclude NaNs

            # Include and scale only those data not within the expected analytical tail
            if sum(included)>0
                scaled = data[included,i] .- minimum(data[included,i])
                if maximum(scaled) > 0
                    scaled = scaled ./ maximum(scaled)
                end
                append!(allscaled, scaled)
            end
        end

        # Calculate kernel density estimate, truncated at 0
        kd = kde(allscaled,npoints=2^7)
        t = kd.x .> cutoff
        return kd.density[t]
    end


## --- Fit and plot results from stationary distribution of depostion/eruption age distribution model

    """
    ```julia
    smpl = tMinDistMetropolis(smpl::ChronAgeData,nsteps::Int,burnin::Int,dist::Array{Float64})
    ```
    Calculate the minimum limiting (eruption/deposition) age of each sample in `smpl`
    using the `metropolis_min` function, assuming mineral ages for each sample are
    drawn from the source distribution `dist`. Fits a `bilinearexponential` function
    to the resulting stationary distribution for each sample.
    """
    function tMinDistMetropolis(smpl::ChronAgeData,nsteps::Int,burnin::Int,dist::Array{Float64}; make_plots=true)
        # Extract variables from struct
        Name = collect(smpl.Name)::Array{String,1}
        Path = smpl.Path::String
        Age_Unit = smpl.Age_Unit::String
        DistType = smpl.Age_DistType::Array{Float64,1}

        # Estimate the eruption/deposition distribution for each sample
        print("Estimating eruption/deposition age distributions...\n")
        for i=1:length(Name)
            if DistType[i] == 0 # A distribution to fit properly
                # Load data for each sample
                data = readdlm(joinpath(Path, Name[i]*".csv"), ',', Float64)::Array{Float64,2}
                print(i, ": ", Name[i], "\n") # Display progress

                # Run MCMC to estimate saturation and eruption/deposition age distributions
                tminDist = metropolis_min(nsteps,dist,data[:,1],data[:,2]/smpl.inputSigmaLevel; burnin=burnin) # Since we don't end up using any of the other distributions
                # (tminDist, tmaxDist, llDist, acceptanceDist) = metropolis_minmax(nsteps,dist,data[:,1],data[:,2]/smpl.inputSigmaLevel, burnin=burnin)

                # Fill in the strat sample object with our new results
                smpl.Age[i] = mean(tminDist)
                smpl.Age_sigma[i] = std(tminDist)
                smpl.Age_025CI[i] = percentile(tminDist,2.5)
                smpl.Age_975CI[i] = percentile(tminDist,97.5)
                smpl.Age_Distribution[i] = tminDist

                # Fit custom many-parametric distribution function to histogram
                edges = range(minimum(tminDist),maximum(tminDist),length=101) # Vector of bin edges
                hobj = fit(Histogram,tminDist,edges,closed=:left) # Fit histogram object

                t = hobj.weights.>0 # Only look at bins with one or more results
                N = hobj.weights[t] ./ length(tminDist) .* length(t) # Normalized number of MCMC steps per bin
                bincenters = cntr(hobj.edges[1])[t] # Vector of bin centers

                # Initial guess for parameters
                p = ones(5)
                p[1] = log(maximum(N))
                p[2] = mean(tminDist)
                p[3] = std(tminDist)

                # Fit nonlinear model
                fobj = curve_fit(bilinear_exponential,bincenters,N,p)
                smpl.Params[:,i] = fobj.param

                if make_plots
                    # Rank-order plot of analyses and eruption/deposition age range
                    nAnalyses = length(data[:,1])
                    h1 = plot_rankorder_errorbar(data[:,1],2*data[:,2]/smpl.inputSigmaLevel,ylabel="Age ($(Age_Unit))",label="Data (observed ages)")
                    m = ones(nAnalyses).*smpl.Age[i]
                    l = ones(nAnalyses).*smpl.Age_025CI[i]
                    u = ones(nAnalyses).*smpl.Age_975CI[i]
                    plot!(h1,1:nAnalyses,l,fillto=u,fillalpha=0.6,linealpha=0, label="Model ($(round(m[1],digits=3)) +$(round(u[1]-m[1],digits=3))/-$(round(m[1]-l[1],digits=3)) $(Age_Unit))")
                    plot!(h1,1:nAnalyses,m,linecolor=:black,linestyle=:dot,label="",legend=:topleft,fg_color_legend=:white,framestyle=:box)
                    savefig(h1,joinpath(Path, Name[i]*"_rankorder.pdf"))
                    savefig(h1,joinpath(Path, Name[i]*"_rankorder.svg"))

                    # Plot model fit to histogram
                    h2 = plot(bincenters,N,label="Histogram",fg_color_legend=:white,framestyle=:box)
                    plot!(h2,bincenters, bilinear_exponential(bincenters,smpl.Params[:,i]), label="Curve fit")
                    plot!(h2,legend=:topleft, xlabel="Age ($(Age_Unit))", ylabel="Probability density")
                    savefig(h2,joinpath(Path, Name[i]*"_distribution.pdf"))
                    savefig(h2,joinpath(Path, Name[i]*"_distribution.pdf"))
                end

            elseif DistType[i] == 1 # A single Gaussian
                # Load data for each sample
                data = readdlm(joinpath(Path, Name[i]*".csv"), ',', Float64)::Array{Float64,2}
                print(i, ": ", Name[i], "\n") # Display progress
                μ = data[1,1]
                σ = data[1,2]

                # Fill in the strat sample object with our new results
                smpl.Age[i] = μ
                smpl.Age_sigma[i] = σ
                smpl.Age_025CI[i] = μ - 1.95996398454*σ
                smpl.Age_975CI[i] = μ + 1.95996398454*σ

                # Initial guess for parameters
                p = ones(5)
                p[1] = log(normpdf(0,σ,0))
                p[2] = μ
                p[3] = σ

                # Fit nonlinear model
                x = μ .+ (-10σ:σ/10:10σ)
                fobj = curve_fit(bilinear_exponential, x, normpdf(μ, σ, x), p)
                fobj.param[5] = 1 # Must be symmetrical
                smpl.Params[:,i] = fobj.param
            end
        end

        # Save results as csv
        results = vcat(["Sample" "Age" "2.5% CI" "97.5% CI" "sigma"], hcat(Name,smpl.Age,smpl.Age_025CI,smpl.Age_975CI,smpl.Age_sigma))::Array{Any,2}
        writedlm(joinpath(Path, "distresults.csv"), results, ',')

        return smpl
    end

## ---
