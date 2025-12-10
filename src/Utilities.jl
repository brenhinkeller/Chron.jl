## --- "Bilinear Exponential" distribution type

    """
    ```Julia
    BilinearExponential(loc, scl, shp, skw)
    BilinearExponential(p::AbstractVector)
    struct BilinearExponential{T<:Real} <: ContinuousUnivariateDistribution
        A::T
        loc::T
        scl::T
        shp::T
        skw::T
    end
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
    xₛ = (x - loc)/scl
    ```
    is `x` scaled by the location parameter `loc` and scale parameter `scl`, 
    In addition to the scale parameter `scl`, the additional shape parameters `shp` and `skw` 
    (which control the sharpness and skew of the resulting distribution, respectively), are 
    all three required to be nonnegative.

    If only four parameters `(loc, scl, shp, skw)` are specified, the normalization constant `A` 
    will be calculated such that the resulting distribution is normalized.
    """
    struct BilinearExponential{T<:Real} <: ContinuousUnivariateDistribution
        A::T
        loc::T
        scl::T
        shp::T
        skw::T
        μ::T # Not a parameter, but precalculated for later use
        σ::T # Not a parameter, but precalculated for later use
        function BilinearExponential{T}(A, loc, scl, shp, skw, μ=loc, σ=scl) where {T}
            new{T}(A, loc, abs(scl), abs(shp), abs(skw), μ, σ)
        end
        function BilinearExponential(A::T, loc::T, scl::T, shp::T, skw::T, μ::T=loc, σ::T=scl) where {T}
            new{T}(A, loc, abs(scl), abs(shp), abs(skw), μ, σ)
        end
    end
    function BilinearExponential(p::AbstractVector{T}) where {T}
        i₀ = firstindex(p)
        @assert ((eachindex(p) == i₀:i₀+3) || (eachindex(p) == i₀:i₀+4))
        if eachindex(p) == i₀:i₀+3
            BilinearExponential{T}(p[i₀], p[i₀+1], p[i₀+2], p[i₀+3])
        elseif eachindex(p) == i₀:i₀+4
            BilinearExponential{T}(p[i₀+1], p[i₀+2], p[i₀+3], p[i₀+4])
        end
    end
    BilinearExponential{T}(loc, scl, shp, skw) where {T<:Real} = BilinearExponential(T(loc), T(scl), T(shp), T(skw))
    function BilinearExponential(loc::T, scl::T, shp::T, skw::T) where {T<:Real}
        # Calculate A such that the resulting distribution is properly normalized
        function f(A, x)
            xs = (x - loc)/scl # X scaled by mean and variance
            v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
            exp(A + shp*skw*xs*v - shp/skw*xs*(1-v))
        end
        area = first(quadgk(x->f(zero(T), x), loc-200abs(scl/shp), loc+200abs(scl/shp)))
        A = -log(area)
        μ = first(quadgk(x->x*f(A,x), loc-100*scl/shp, loc+100*scl/shp, maxevals=1000))
        σ = sqrt(first(quadgk(x->(x-μ)^2*f(A,x), loc-100*scl/shp, loc+100*scl/shp, maxevals=1000)))
        BilinearExponential{float(T)}(A, loc, scl, shp, skw, μ, σ)
    end

    ## Conversions
    Base.convert(::Type{BilinearExponential{T}}, d::BilinearExponential) where {T<:Real} = BilinearExponential{T}(T(d.A), T(d.loc), T(d.scl), T(d.shp), T(d.skw))
    Base.convert(::Type{BilinearExponential{T}}, d::BilinearExponential{T}) where {T<:Real} = d

    ## Parameters
    Distributions.params(d::BilinearExponential) = (d.A, d.loc, d.scl, d.shp, d.skw)
    @inline Distributions.partype(d::BilinearExponential{T}) where {T<:Real} = T

    Distributions.location(d::BilinearExponential) = d.loc
    Distributions.scale(d::BilinearExponential) = d.scl

    Base.eltype(::Type{BilinearExponential{T}}) where {T} = T

    ## Evaluation
    @inline function Distributions.pdf(d::BilinearExponential, x::Real)
        xs = (x - d.loc)/d.scl # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        return exp(d.A + d.shp*d.skw*xs*v - d.shp/d.skw*xs*(1-v))
    end
    function Distributions.cdf(d::BilinearExponential{T}, x::Real) where {T}
        maxspan = 200*d.scl/d.shp
        if d.loc + maxspan < x
            one(float(T))
        elseif x < d.loc - maxspan
            zero(float(T))
        else
            l = min(d.loc, x)
            return first(quadgk(x->pdf(d,x), l-maxspan, x))
        end
    end
    function Distributions.ccdf(d::BilinearExponential{T}, x::Real) where {T}
        maxspan =  200*d.scl/d.shp
        if x < d.loc - maxspan
            one(float(T))
        elseif d.loc + maxspan < x
            zero(float(T))
        else
            u = max(d.loc, x)
            first(quadgk(x->pdf(d,x), x, u+maxspan))
        end
    end

    @inline function Distributions.logpdf(d::BilinearExponential, x::Real)
        xs = (x - d.loc)/d.scl # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        return d.A + d.shp*d.skw*xs*v - d.shp/d.skw*xs*(1-v)
    end
    function Distributions.logcdf(d::BilinearExponential{T}, x::Real) where {T}
        maxspan = 200*d.scl/d.shp
        if x > d.loc + maxspan
            zero(float(T))
        else
            l = min(d.loc, x)
            c = logpdf(d, l)
            area_shifted, e = quadgk(x->exp(logpdf(d,x)-c), l-maxspan, x)
            log(area_shifted)+c
        end
    end
    function Distributions.logccdf(d::BilinearExponential{T}, x::Real) where {T}
        maxspan = 200*d.scl/d.shp
        if x < d.loc - maxspan
            zero(float(T))
        else
            u = max(d.loc, x)
            c = logpdf(d, u)
            area_shifted, e = quadgk(x->exp(logpdf(d,x)-c), x, u+maxspan)
            log(area_shifted)+c
        end
    end


    ## Statistics
    Distributions.mean(d::BilinearExponential) = d.μ
    Distributions.var(d::BilinearExponential) = d.σ*d.σ
    Distributions.std(d::BilinearExponential) = d.σ
    Distributions.skewness(d::BilinearExponential) = first(quadgk(x->(x-d.μ)^3*pdf(d,x), d.loc-100*d.scl/d.shp, d.loc+100*d.scl/d.shp, maxevals=1000))/d.σ^3
    Distributions.kurtosis(d::BilinearExponential) = first(quadgk(x->(x-d.μ)^4*pdf(d,x), d.loc-100*d.scl/d.shp, d.loc+100*d.scl/d.shp, maxevals=1000))/d.σ^4 - 3

    ## Affine transformations
    Base.:+(d::BilinearExponential{T}, c::Real) where {T} = BilinearExponential{T}(d.A, d.loc + c, d.scl, d.shp, d.skw)
    Base.:*(d::BilinearExponential{T}, c::Real) where {T} = BilinearExponential{T}(d.A-log(abs(c)), d.loc * c, d.scl*abs(c), d.shp, d.skw)

    
## --- Extend `normpdf_ll` to deal with Distributions.Normal 

    function StatGeochemBase.normpdf_ll(x::AbstractVector, ages::AbstractVector{Normal{T}}) where T
        ll = zero(T)
        @inbounds for i ∈ eachindex(x,ages)
            age = ages[i]
            μᵢ, σᵢ = age.μ, age.σ
            ll -= (x[i]-μᵢ)^2 / (2*σᵢ*σᵢ)
        end
        return ll
    end

## --- log likelihood functions allowing for arbitrary Distributions

    # Use dispatch to let us reduce duplication
    function strat_ll(x, ages)
        ll = zero(float(eltype(x)))
        @assert eachindex(x) == eachindex(ages)
        @inbounds for i in eachindex(x)
            ll += fastlogpdf(ages[i], x[i])
        end
        return ll
    end
    function strat_ll(x, ages, sidedness::FastSidedness)
        ll = zero(float(eltype(x)))
        @assert eachindex(x) == eachindex(ages) == eachindex(sidedness)
        @inbounds for i in eachindex(x)
            μᵢ = mean(ages[i])
            ll += if sidedness[i] == sign(x[i] - μᵢ)
                fastlogpdf(ages[i], μᵢ)
            else
                fastlogpdf(ages[i], x[i])
            end
        end
        return ll
    end
    function strat_ll(x, ages, sidedness::CDFSidedness)
        ll = zero(float(eltype(x)))
        @assert eachindex(x) == eachindex(ages) == eachindex(sidedness)
        @inbounds for i in eachindex(x)
            ll += if sidedness[i] > 0 # Minimum age
                logcdf(ages[i], x[i])
            elseif sidedness[i] < 0  # Maximum age
                logccdf(ages[i], x[i])
            else
                fastlogpdf(ages[i], x[i])
            end
        end
        return ll
    end

    fastlogpdf(d, x::Real) = logpdf(d, x)
    function fastlogpdf(d::Normal, x::Real)
        δ, σ = (x - d.μ), d.σ
        - δ*δ/(2*σ*σ)
    end


## --- End of File
