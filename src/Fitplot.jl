## --- Utility function for file IO

    @inline function readclean(filepath, delim, T)
        io = open(filepath, "r")
        if read(io, Char) == '\ufeff'
            @warn """Skipping hidden \'\\ufeff\' (U+FEFF) character at start of input file.

            This character is often added to CSV files by Microsoft Excel (and some other
            Microsoft products) as what appears to be what we might call an "extension",
            which would would cause file parsing to fail if we didn't manually remove it.

            Try using open software like LibreOffice instead of Excel to make this warning go away.
            """
        else
            seekstart(io)
        end
        data = readdlm(io, delim, T)
        close(io)
        return data
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
        for i ∈ eachindex(Name)
            # With screening
            # Maximum offset before cutoff
            # Read data for each sample from file
            filepath = joinpath(Path, Name[i]*".csv")
            data_raw = readclean(filepath, ',', Float64)::Matrix{Float64}
            # Sort ages in ascending order
            data = sortslices(data_raw, dims=1)
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


## --- Fit and plot results from stationary distribution of depostion/eruption age distribution model

    """
    ```julia
    tMinDistMetropolis(smpl::ChronAgeData, nsteps::Int, burnin::Int, dist::Array{Float64})
    ```
    Calculate the minimum limiting (eruption/deposition) age of each sample defined
    in the `smpl` struct, using the `metropolis_min` function, assuming mineral
    ages for each sample are drawn from the source distribution `dist`. Fits a
    `bilinear_exponential` function to the resulting stationary distribution
    for each sample and stores the results in `smpl.Params` for use by the
    `StratMetropolisDist` function.

    ### Examples
    ```julia
    smpl = tMinDistMetropolis(smpl, 5*10^5, 2*10^5, TriangularDistribution)
    ```
    """
    function tMinDistMetropolis(smpl::ChronAgeData,nsteps::Int,burnin::Int,dist::DenseArray{Float64};
            include=ntuple(i->true, length(smpl.Name)),
            make_plots=true,
        )
        # Extract variables from struct
        Name = collect(smpl.Name)::Array{String,1}
        Path = smpl.Path::String
        Age_Unit = smpl.Age_Unit::String
        DistType = smpl.Age_DistType::Vector{Float64}

        # Estimate the eruption/deposition distribution for each sample
        print("Estimating eruption/deposition age distributions...\n")
        for i ∈ eachindex(Name)
            include[i] || continue # Only calculate if include is true
            if DistType[i] == 0 # A distribution to fit properly
                # Read data for each sample from file
                filepath = joinpath(Path, Name[i]*".csv")
                data = readclean(filepath, ',', Float64)::Matrix{Float64}
                print(i, ": ", Name[i], "\n") # Display progress

                # Run MCMC to estimate saturation and eruption/deposition age distributions
                μ̄ = view(data,:,1)
                σ̄ = view(data, :, 2)./=smpl.inputSigmaLevel
                tminDist = metropolis_min(nsteps, dist, μ̄, σ̄; burnin) # Since we don't end up using any of the other distributions
                # (tminDist, tmaxDist, llDist, acceptanceDist) = metropolis_minmax(nsteps,dist,data[:,1],data[:,2]/smpl.inputSigmaLevel, burnin=burnin)

                # Fill in the strat sample object with our new results
                tminDistₜ = copy(tminDist)
                smpl.Age[i] = nanmean(tminDist)
                smpl.Age_sigma[i] = nanstd(tminDist)
                smpl.Age_025CI[i] = nanpctile!(tminDistₜ,2.5)
                smpl.Age_975CI[i] = nanpctile!(tminDistₜ,97.5)
                smpl.Age_Distribution[i] = tminDist

                # Fit custom many-parametric distribution function to histogram
                edges = range(nanminimum(tminDist),nanmaximum(tminDist),length=101) # Vector of bin edges
                bincounts = histcounts(tminDist, edges)

                t = bincounts.>0 # Only look at bins with one or more results
                N = bincounts[t] ./ length(tminDist) .* count(t) # Normalized number of MCMC steps per bin
                bincenters = cntr(edges)[t] # Vector of bin centers

                # Initial guess for parameters
                p = ones(5)
                p[1] = log(maximum(N))
                p[2] = nanmean(tminDist)
                p[3] = nanstd(tminDist)

                # Fit nonlinear model
                fobj = curve_fit(bilinear_exponential,bincenters,N,p)
                smpl.Params[:,i] = fobj.param

                if make_plots
                    # Rank-order plot of analyses and eruption/deposition age range
                    nAnalyses = length(data[:,1])
                    h1 = rankorder(data[:,1],2*data[:,2]/smpl.inputSigmaLevel,ylabel="Age ($(Age_Unit))",label="Data (observed ages)")
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
                    savefig(h2,joinpath(Path, Name[i]*"_distribution.svg"))
                    savefig(h2,joinpath(Path, Name[i]*"_distribution.pdf"))
                end

            elseif DistType[i] == 1 # A single Gaussian
                # Read data for each sample from file
                filepath = joinpath(Path, Name[i]*".csv")
                data = readclean(filepath, ',', Float64)::Matrix{Float64}
                print(i, ": ", Name[i], "\n") # Display progress
                μ = data[1,1]
                σ = data[1,2]/smpl.inputSigmaLevel

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
    # Automatically do what people probably want if given improper types
    function tMinDistMetropolis(smpl::ChronAgeData, nsteps::Number, burnin::Number, dist::AbstractArray{Number}; make_plots=true)
        tMinDistMetropolis(smpl,Int(nsteps), Int(burnin), Float64.(dist); make_plots)
    end


## ---
