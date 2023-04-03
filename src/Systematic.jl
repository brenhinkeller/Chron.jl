## --- Dealing with systematic uncertainty

    function add_systematic_uncert_UPb(agedistmyr::Vector{<:AbstractFloat})
        λ238 = log(2)/44683 # Jaffey decay constant, 1/Myr
        λ238_ = log(2)./(44683 .+ 24 .* randn(size(agedistmyr))) # Jaffey decay with uncertainty 1/Myr
        σtracer = 1 .+ (0.03/2/100) .* randn(size(agedistmyr)) # Earthtime tracer uncertainty
        r = exp.(agedistmyr .* λ238) .- 1
        return log.(r .* σtracer .+ 1) ./ λ238_
    end
    function add_systematic_uncert_UPb(agedistmyr::Matrix{<:AbstractFloat})
        λ238 = log(2)/44683 # Jaffey decay constant, 1/Myr
        λ238_ = log(2)./(44683 .+ 24 .* randn(1, size(agedistmyr,2))) # Jaffey decay with uncertainty 1/Myr
        σtracer = 1 .+ (0.03/2/100) .* randn(1, size(agedistmyr,2)) # Earthtime tracer uncertainty
        r = exp.(agedistmyr .* λ238) .- 1
        return log.(r .* σtracer .+ 1) ./ λ238_
    end

    export add_systematic_uncert_UPb

    function add_systematic_uncert_ArAr(agedistmyr::Vector{<:AbstractFloat}; constants=:Renne)
        if constants==:Min
            # Min et al., 200 compilation
            κ = 1.6408E-03
            λϵ = 0.580E-10 # 1/years
            λβ = 4.884E-10 # 1/years
            # Uncertainties
            σκ = 0.0047E-03
            σλϵ = 0.007E-10 # 1/years
            σλβ = 0.049E-10 # 1/years
            # Covariances
            σκσλϵ = 0.0
            σκσλβ = 0.0
            σλϵσλβ = 0.0
        else
            # Optimization intercalibration values from Renne 2011 reply
            # Values
            κ = 1.6417E-03
            λϵ = 0.5757E-10 # 1/years
            λβ = 4.9548E-10 # 1/years
            # Uncertainties
            σκ = 0.0045E-03
            σλϵ = 0.0017E-10 # 1/years
            σλβ = 0.0134E-10 # 1/years
            # Covariances
            σκσλϵ =  7.1903E-19
            σκσλβ = -6.5839E-19
            σλϵσλβ = -3.4711E-26
        end

        # Create a multivariate normal distribution from the decay constant
        # means and covariance matrix
        μ = [κ, λϵ, λβ]
        Σ = [ σκ^2   σκσλϵ  σκσλβ;
              σκσλϵ  σλϵ^2 σλϵσλβ;
              σκσλβ σλϵσλβ  σλβ^2;]
        κλλ_ = rand(MvNormal(μ,Σ),length(agedistmyr)) # Draw from MvNormal
        # Exctract individual variables
        κ_ = κλλ_[1,:]
        λϵ_ = κλλ_[2,:]
        λβ_ = κλλ_[3,:]

        # Convert age distribution to R distribution. Assume Myr input
        λ = λβ + λϵ
        Rdist = (exp.(agedistmyr.*λ.*1000000) .- 1) .* λϵ./(κ.*λ)

        # Now convert back to age distribution
        λ_ = λβ_ .+ λϵ_
        return log.(Rdist.*κ_.*λ_./λϵ_ .+ 1)./(λ_.*1000000)
    end
    function add_systematic_uncert_ArAr(agedistmyr::Matrix{<:AbstractFloat}; constants=:Renne)
        if constants==:Min
            # Min et al., 200 compilation
            κ = 1.6408E-03
            λϵ = 0.580E-10 # 1/years
            λβ = 4.884E-10 # 1/years
            # Uncertainties
            σκ = 0.0047E-03
            σλϵ = 0.007E-10 # 1/years
            σλβ = 0.049E-10 # 1/years
            # Covariances
            σκσλϵ = 0.0
            σκσλβ = 0.0
            σλϵσλβ = 0.0
        else
            # Optimization intercalibration values from Renne 2011 reply
            # Values
            κ = 1.6417E-03
            λϵ = 0.5757E-10 # 1/years
            λβ = 4.9548E-10 # 1/years
            # Uncertainties
            σκ = 0.0045E-03
            σλϵ = 0.0017E-10 # 1/years
            σλβ = 0.0134E-10 # 1/years
            # Covariances
            σκσλϵ =  7.1903E-19
            σκσλβ = -6.5839E-19
            σλϵσλβ = -3.4711E-26
        end

        # Create a multivariate normal distribution from the decay constant
        # means and covariance matrix
        μ = [κ, λϵ, λβ]
        Σ = [ σκ^2   σκσλϵ  σκσλβ;
              σκσλϵ  σλϵ^2 σλϵσλβ;
              σκσλβ σλϵσλβ  σλβ^2;]
        κλλ_ = rand(MvNormal(μ,Σ),size(agedistmyr,2)) # Draw from MvNormal
        # Exctract individual variables
        κ_ = κλλ_[1:1,:]
        λϵ_ = κλλ_[2:2,:]
        λβ_ = κλλ_[3:3,:]

        # Convert age distribution to R distribution
        λ = λβ + λϵ
        Rdist = (exp.(agedistmyr.*λ.*1000000) .- 1) .* λϵ./(κ.*λ)

        # Now convert back to age distribution
        λ_ = λβ_ .+ λϵ_
        return log.(Rdist.*κ_.*λ_./λϵ_ .+ 1)./(λ_.*1000000)
    end

    export add_systematic_uncert_ArAr

## --- End of File

    function add_systematic_uncert_UTh(agedistmyr::Vector{<:AbstractFloat})
        # Age = -log(1-Th230ₐ/U238ₐ))/λ230Th
        # Th230ₐ/U238ₐ = 1-exp(-Age*λ230Th)
        # Th230/U238 = (1-exp(-Age*λ230Th)) * λ238U/λ230Th

        λ238 = val(λ238U)
        λ238_ = val(λ238U) .+ err(λ238U).*randn(size(agedistmyr))

        λ230 = val(λ230Th)
        λ230_ = val(λ230Th) .+ err(λ230Th).*randn(size(agedistmyr))

        # Convert age distribution to ratio distribution
        # (i.e., distribution of slope of two-point isochron)
        sTh230_U238_dist = @. (1-exp(-agedistmyr*λ230)) * λ238/λ230

        # Convert from ratio distribution back to age distribution, with uncertainty
        # Maximum possible activity ratio is 1.0, corresponding to secular equilibrium
        return @. -log(1 - min(sTh230_U238_dist*λ230_/λ238_, 1.0)) / λ230_
    end
    function add_systematic_uncert_UTh(agedistmyr::Matrix{<:AbstractFloat})
        # Age = -log(1-Th230ₐ/U238ₐ))/λ230Th
        # Th230ₐ/U238ₐ = 1-exp(-Age*λ230Th)
        # Th230/U238 = (1-exp(-Age*λ230Th)) * λ238U/λ230Th

        λ238 = val(λ238U)
        λ238_ = val(λ238U) .+ err(λ238U).*randn(1, size(agedistmyr,2))

        λ230 = val(λ230Th)
        λ230_ = val(λ230Th) .+ err(λ230Th).*randn(1, size(agedistmyr,2))

        # Convert age distribution to ratio distribution
        # (i.e., distribution of slope of two-point isochron)
        sTh230_U238_dist = @. (1-exp(-agedistmyr*λ230)) * λ238/λ230

        # Convert from ratio distribution back to age distribution, with uncertainty
        # Maximum possible activity ratio is 1.0, corresponding to secular equilibrium
        return @. -log(1 - min(sTh230_U238_dist*λ230_/λ238_, 1.0)) / λ230_
    end
    export add_systematic_uncert_UTh
