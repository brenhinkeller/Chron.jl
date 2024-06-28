## --- Dealing with arbitrary distributions


## --- "bilinear exponential" distribution type

    """
    ```Julia
    struct BilinearExponential{T<:Real} <: ContinuousUnivariateDistribution
        A::T
        μ::T
        σ::T
        shp::T
        skw::T
    end
    BilinearExponential(p::AbstractVector)
    BilinearExponential(p::NTuple{5,T})
    ```
    A five-parameter pseudo-distribution which can be used to approximate various 
    asymmetric probability distributions found in geochronology (including as a result 
    of Bayesian eruption age estimation). 
    
    This "bilinear exponential" distribution, as the name might suggest, is defined as an 
    exponential function with two log-linear segments, joined by an arctangent sigmoid:
    ```math
    ℯ^{A} * ℯ^{v*xₛ*shp*skw - (1-v)*xₛ*shp/skw}
    ```
    where
    ```math
    v = 1/2 - atan(xₛ)/π
    ```
    is a sigmoid, positive on the left-hand side, and
    ```math
    xₛ = (x - μ)/σ
    ```
    is `x` scaled by the location parameter `μ` and scale parameter `σ`, 
    In addition to the scale parameter `σ`, the additional shape parameters `shp` and `skw` 
    (which control the sharpness and skew of the resulting distribution, respectively), are 
    all three required to be nonnegative.
    """
    struct BilinearExponential{T<:Real} <: ContinuousUnivariateDistribution
        A::T
        μ::T
        σ::T
        shp::T
        skw::T
    end
    function BilinearExponential(p::AbstractVector{T}) where {T}
        @assert length(p) == 5
        @assert all(x->!(x<0), Iterators.drop(p,1))
        BilinearExponential{T}(p...)
    end
    function BilinearExponential(p::NTuple{5,T}) where {T}
        @assert all(x->!(x<0), Iterators.drop(p,1))
        BilinearExponential{T}(p...)
    end

    ## Conversions
    Base.convert(::Type{BilinearExponential{T}}, d::Normal) where {T<:Real} = BilinearExponential{T}(T(d.A), T(d.μ), T(d.σ), T(d.shp), T(d.skw))
    Base.convert(::Type{BilinearExponential{T}}, d::Normal{T}) where {T<:Real} = d

    ## Parameters
    Distributions.params(d::BilinearExponential) = (d.A, d.μ, d.σ, d.shp, d.skw)
    @inline Distributions.partype(d::BilinearExponential{T}) where {T<:Real} = T

    Distributions.location(d::BilinearExponential) = d.μ
    Distributions.scale(d::BilinearExponential) = d.σ

    Base.eltype(::Type{BilinearExponential{T}}) where {T} = T

    ## Evaluation
    @inline function Distributions.pdf(d::BilinearExponential, x::Real)
        xs = (x - d.μ)/d.σ # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        return exp(d.A + d.shp*d.skw*xs*v - d.shp/d.skw*xs*(1-v))
    end

    @inline function Distributions.logpdf(d::BilinearExponential, x::Real)
        xs = (x - d.μ)/d.σ # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        return d.A + d.shp*d.skw*xs*v - d.shp/d.skw*xs*(1-v)
    end

    ## Affine transformations
    Base.:+(d::BilinearExponential{T}, c::Real) where {T} = BilinearExponential{T}(d.A, d.μ + c, d.σ, d.shp, d.skw)
    Base.:*(d::BilinearExponential{T}, c::Real) where {T} = BilinearExponential{T}(d.A, d.μ * c, d.σ * abs(c), d.shp, d.skw)

## --- Radiocarbon distribution type

    struct Radiocarbon{T<:Real} <: ContinuousUnivariateDistribution
        μ::T
        σ::T
        dist::Vector{T}
        ldist::Vector{T}
    end

    function Radiocarbon(μ::T, σ::T, ldist::Vector{T}) where {T<:Real}
        dist = exp.(ldist)
        dist ./= sum(dist) * 1 # Normalize
        return Radiocarbon{T}(μ, σ, dist, ldist)
    end

    function Radiocarbon(Age_14C::Real, Age_14C_sigma::Real, calibration::NamedTuple=intcal13)
        @assert calibration.Age_Calendar == 1:1:length(calibration.Age_14C)
        @assert step(calibration.Age_Calendar) == calibration.dt == 1

        ldist = normproduct_ll.(Age_14C, Age_14C_sigma, calibration.Age_14C, calibration.Age_sigma)

        dist = exp.(ldist)
        dist ./= sum(dist) * calibration.dt # Normalize
        μ = histmean(dist, calibration.Age_Calendar)
        σ = histstd(dist, calibration.Age_Calendar, corrected=false)

        return Radiocarbon(μ, σ, dist, ldist)
    end

    ## Conversions
    Base.convert(::Type{Radiocarbon{T}}, d::Normal) where {T<:Real} = Radiocarbon{T}(T(d.μ), T(d.σ), T.(d.ldist))
    Base.convert(::Type{Radiocarbon{T}}, d::Normal{T}) where {T<:Real} = d

    ## Parameters
    Distributions.params(d::Radiocarbon) = (d.dist, d.ldist)
    @inline Distributions.partype(d::Radiocarbon{T}) where {T<:Real} = T

    Distributions.location(d::Radiocarbon) = d.μ
    Distributions.scale(d::Radiocarbon) = d.σ

    Base.eltype(::Type{Radiocarbon{T}}) where {T} = T

    ## Evaluation
    @inline function Distributions.pdf(d::Radiocarbon{T}, x::Real) where {T}
        return linterp_at_index(d.dist, x, zero(T))
    end

    @inline function Distributions.logpdf(d::Radiocarbon{T}, x::Real) where {T}
        return linterp_at_index(d.ldist, x, -maxintfloat(T))
    end

    ## Statistics
    Distributions.mean(d::Radiocarbon) = d.μ
    Distributions.std(d::Radiocarbon) = d.σ

## --- 

    # Interpolate log likelihood from an array
    function interpolate_ll(x::AbstractVector,p::AbstractMatrix{T}) where {T<:Real}
        ll = zero(T)
        @inbounds for i ∈ eachindex(x)
            ll += linterp_at_index(view(p,:,i), x[i], -maxintfloat(T))
        end
        return ll
    end
    function interpolate_ll(x::AbstractVector,ages::AbstractVector{Radiocarbon{T}}) where T
        ll = zero(T)
        @inbounds for i ∈ eachindex(x,ages)
            ldist = ages[i].ldist
            ll += linterp_at_index(ldist, x[i], -maxintfloat(T))
        end
        return ll
    end

    function StatGeochemBase.normpdf_ll(x::AbstractVector, ages::AbstractVector{Normal{T}}) where T
        ll = zero(T)
        @inbounds for i ∈ eachindex(x,ages)
            age = ages[i]
            μᵢ, σᵢ = age.μ, age.σ
            ll -= (x[i]-μᵢ)^2 / (2*σᵢ*σᵢ)
        end
        return ll
    end

## --- End of File
