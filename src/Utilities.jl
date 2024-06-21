## --- Fitting non-Gaussian distributions

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
    function Distributions.pdf(d::BilinearExponential, x::Real)
        xs = (x - d.μ)/d.σ # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        return exp(d.A + d.shp*d.skw*xs*v - d.shp/d.skw*xs*(1-v))
    end

    function Distributions.logpdf(d::BilinearExponential, x::Real)
        xs = (x - d.μ)/d.σ # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        return d.A + d.shp*d.skw*xs*v - d.shp/d.skw*xs*(1-v)
    end

    ## Affine transformations
    Base.:+(d::BilinearExponential{T}, c::Real) where {T} = BilinearExponential{T}(d.A, d.μ + c, d.σ, d.shp, d.skw)
    Base.:*(d::BilinearExponential{T}, c::Real) where {T} = BilinearExponential{T}(d.A, d.μ * c, d.σ * abs(c), d.shp, d.skw)

    """
    ```Julia
    bilinear_exponential(x::Number, p::AbstractVector)
    ```
    Evaluate the value of a "bilinear exponential" function defined by the parameters
    `p` at point `x`. This function, which can be used to approximate various asymmetric
    probability distributions found in geochronology, is defined by an exponential
    function with two log-linear segments joined by an `atan` sigmoid:
    ```math
    ℯ^{p_1} * ℯ^{v x_s p_4 p_5 - (1-v) x_s p_4/p_5}
    ```
    where
    ```math
    v = 1/2 - atan(x_s)/π
    ```
    is a sigmoid, positive on the left-hand side, and
    ```math
    x_s = (x - p_2)/p_3
    ```
    is `x` scaled by mean and standard deviation.

    The elements of the parameter array `p` may be considered to approximately represent\n
        p[1] # pre-exponential (normaliation constant)
        p[2] # mean (central moment)
        p[3] # standard deviation
        p[4] # sharpness
        p[5] # skew

    where all parameters `pᵢ` must be nonnegative.
    """
    function bilinear_exponential(x::Number, p::AbstractVector{<:Number})
        @assert length(p) == 5
        _, μ, σ, shp, skw = (abs(x) for x in p)
        xs = (x - μ)/σ # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        return exp(first(p) + shp*skw*xs*v - shp/skw*xs*(1-v))
    end
    function bilinear_exponential(x::AbstractVector, p::AbstractVector{<:Number})
        @assert length(p) == 5 && firstindex(p) == 1
        result = Array{float(eltype(x))}(undef,size(x))
        @inbounds for i ∈ eachindex(x)
            xs = (x[i] - p[2])/abs(p[3]) # X scaled by mean and variance
            v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
            shp = abs(p[4])
            skw = abs(p[5])
            result[i] = exp(p[1] + shp*skw*xs*v - shp/skw*xs*(1-v))
        end
        return result
    end

    """
    ```Julia
    bilinear_exponential_ll(x, p)
    ```
    Return the log likelihood corresponding to a `bilinear_exponential` distribution defined
    by the parameters `p` evaluate at point `x`.

    If `x` is provided as a number and `p` as a vector, a single evaluation will
    be returned; if `x` is provided as a vector and `p` as a matrix, the sum
    over all i for each distribution `p[:,i]` evaluated at `x[i]` will be returned
    instead.

    See also `bilinear_exponential`
    """
    function bilinear_exponential_ll(x::Number, p::AbstractVector{T}) where {T<:Number}
        @assert length(p) == 5
        isnan(p[2]) && return zero(T)
        _, μ, σ, shp, skw = (abs(x) for x in p)
        xs = (x - μ)/σ # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        return shp*skw*xs*v - shp/skw*xs*(1-v)
    end
    function bilinear_exponential_ll(x::AbstractVector, p::AbstractMatrix{T}) where {T<:Number}
        ll = zero(T)
        @inbounds for i ∈ eachindex(x)
            any(isnan, p[2,i]) && continue
            xs = (x[i]-p[2,i])/abs(p[3,i]) # X scaled by mean and variance
            v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
            shp = abs(p[4,i])
            skw = abs(p[5,i])
            ll += shp*skw*xs*v - shp/skw*xs*(1-v)
        end
        return lles
    end
    function bilinear_exponential_ll(x::AbstractVector, ages::AbstractVector{BilinearExponential{T}}) where T
        ll = zero(T)
        @inbounds for i ∈ eachindex(x,ages)
            pᵢ = ages[i]
            isnan(pᵢ.μ) && continue
            xs = (x[i]-pᵢ.μ)/(pᵢ.σ) # X scaled by mean and variance
            v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
            shp, skw = pᵢ.shp, pᵢ.skw
            ll += shp*skw*xs*v - shp/skw*xs*(1-v)
        end
        return ll
    end

    struct Radiocarbon{T}
        μ::T
        σ::T
        dist::Vector{T}
    end

    # Interpolate log likelihood from an array
    function interpolate_ll(x::AbstractVector,p::AbstractMatrix{T}) where {T<:Number}
        ll = zero(T)
        @inbounds for i ∈ eachindex(x)
            ll += linterp_at_index(view(p,:,i), x[i], -1e9)
        end
        return ll
    end
    function interpolate_ll(x::AbstractVector,ages::AbstractVector{Radiocarbon{T}}) where T
        ll = zero(T)
        @inbounds for i ∈ eachindex(x,ages)
            dist = ages[i].dist
            ll += linterp_at_index(dist, x[i], -1e9)
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

## --- Biquadratic distribution (like bilinear, with quadratic sides)


    struct BiquadraticExponential{T<:Real}
        A::T
        μ::T
        σ::T
        shp::T
        skw::T
    end
    function BiquadraticExponential(p::AbstractVector)
        @assert length(p) == 5
        @assert all(x->x>0, Iterators.drop(p,1))
        BiquadraticExponential(p...)
    end

    function biquadratic_exponential(x::Number, p::AbstractVector{<:Number})
        @assert length(p) == 5
        _, μ, σ, shp, skw = (abs(x) for x in p)
        xs = (x - μ)/σ # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        return exp(first(p) -0.5*(shp*skw*xs*v - shp/skw*xs*(1-v))^2)
    end
    function biquadratic_exponential(x::AbstractVector, p::AbstractVector{<:Number})
        @assert length(p) == 5
        _, μ, σ, shp, skw = (abs(x) for x in p)
        result = Array{float(eltype(x))}(undef,size(x))
        @inbounds for i ∈ eachindex(x)
            xs = (x[i] - μ)/σ # X scaled by mean and variance
            v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
            result[i] = exp(first(p) -0.5*(shp*skw*xs*v - shp/skw*xs*(1-v))^2)
        end
        return result
    end

    function biquadratic_exponential_ll(x::Number, p::AbstractVector{<:Number})
        xs = (x - p[2])/abs(p[3]) # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        shp = abs(p[4])
        skw = abs(p[5])
        return -(shp*skw*xs*v - shp/skw*xs*(1-v))^2
    end
    function biquadratic_exponential_ll(x::AbstractVector, ages::AbstractVector{BiquadraticExponential{T}}) where T
        ll = zero(T)
        @inbounds for i ∈ eachindex(x,ages)
            pᵢ = ages[i]
            xs = (x[i]-pᵢ.μ)/(pᵢ.σ) # X scaled by mean and variance
            v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
            shp = pᵢ.shp
            skw = pᵢ.skw
            ll -= (shp*skw*xs*v - shp/skw*xs*(1-v))^2
        end
        return ll
    end

## --- End of File
