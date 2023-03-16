## --- Dealing with systematic uncertainty

    function add_systematic_uncert_UPb(agedistmyr::Vector{Float64})
        λ238 = log(2)/44683 # Jaffey decay constant, 1/Myr
        λ238_ = log(2)./(44683 .+ 24 .* randn(size(agedistmyr))) # Jaffey decay with uncertainty 1/Myr
        σtracer = 1 .+ (0.03/2/100) .* randn(size(agedistmyr)) # Earthtime tracer uncertainty
        r = exp.(agedistmyr .* λ238) .- 1
        return log.(r .* σtracer .+ 1) ./ λ238_
    end
    function add_systematic_uncert_UPb(agedistmyr::Matrix{Float64})
        λ238 = log(2)/44683 # Jaffey decay constant, 1/Myr
        λ238_ = log(2)./(44683 .+ 24 .* randn(1, size(agedistmyr,2))) # Jaffey decay with uncertainty 1/Myr
        σtracer = 1 .+ (0.03/2/100) .* randn(1, size(agedistmyr,2)) # Earthtime tracer uncertainty
        r = exp.(agedistmyr .* λ238) .- 1
        return log.(r .* σtracer .+ 1) ./ λ238_
    end

    export add_systematic_uncert_UPb

    function add_systematic_uncert_ArAr(agedistmyr::Vector{Float64})
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
    function add_systematic_uncert_ArAr(agedistmyr::Matrix{Float64})
        # Optimization intercalibration values from Renne 2011 reply
        # Values
        κ = 1.6417E-03
        λϵ = 0.5757E-10
        λβ = 4.9548E-10
        # Uncertainties
        σκ = 0.0045E-03
        σλϵ = 0.0017E-10
        σλβ = 0.0134E-10
        # Covariances
        σκσλϵ =  7.1903E-19
        σκσλβ = -6.5839E-19
        σλϵσλβ = -3.4711E-26

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
