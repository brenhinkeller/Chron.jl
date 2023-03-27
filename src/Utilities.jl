## --- Fitting non-Gaussian distributions

    struct BilinearExponential{T}
        A::T
        μ::T
        σ::T
        sharpness::T
        skew::T
    end
    function BilinearExponential(p::AbstractVector)
        @assert length(p) == 5
        @assert all(x->x>0, p)
        BilinearExponential(p...)
    end

    """
    ```Julia
    bilinear_exponential(x::Number, p::AbstractVector)
    ```
    Evaluate the value of a "bilinear exponential" function defined by the parameters
    `p` at point `x`. This function, which can be used to approximate various asymmetric
    probability distributions found in geochronology, is defined by an exponential
    function with two log-linear segments joined by an `atan` sigmoid:
    ```math
    ℯ^{p_1 + v x_s p_4^2 p_5^2 - (1-v) x_s p_4^2/p_5^2}
    ```
    where
    ```math
    v = 1/2 - atan(x_s)/π
    ```
    is a sigmoid, positive on the left-hand side, and
    ```math
    x_s = (x - p_2)/p_3^2
    ```
    is `x` scaled by mean and variance.

    The elements of the parameter array `p` may be considered to approximately represent\n
        p[1] # pre-exponential (log normaliation constant)
        p[2] # mean (central moment)
        p[3] # standard deviation
        p[4] # sharpness
        p[5] # skew
    """
    function bilinear_exponential(x::Number, p::AbstractVector{<:Number})
        @assert length(p) == 5
        A, μ, σ, shp, skw = (abs(x) for x in p)
        xs = (x - μ)/σ # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        return A * exp(shp*skw*xs*v - shp/skw*xs*(1-v))
    end
    function bilinear_exponential(x::AbstractVector, p::AbstractVector{<:Number})
        @assert length(p) == 5
        result = Array{float(eltype(x))}(undef,size(x))
        @inbounds for i ∈ eachindex(x)
            xs = (x[i] - p[2])/abs(p[3]) # X scaled by mean and variance
            v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
            shp = abs(p[4])
            skw = abs(p[5])
            result[i] = abs(p[1]) * exp(shp*skw*xs*v - shp/skw*xs*(1-v))
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
    function bilinear_exponential_ll(x::Number, p::AbstractVector{<:Number})
        @assert length(p) == 5
        A, μ, σ, shp, skw = (abs(x) for x in p)
        xs = (x - μ)/σ # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        return shp*skw*xs*v - shp/skw*xs*(1-v)
    end
    function bilinear_exponential_ll(x::AbstractVector, p::AbstractMatrix{<:Number})
        ll = 0.0
        @inbounds for i ∈ eachindex(x)
            xs = (x[i]-p[2,i])/abs(p[3,i]) # X scaled by mean and variance
            v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
            shp = abs(p[4,i])
            skw = abs(p[5,i])
            ll += shp*skw*xs*v - shp/skw*xs*(1-v)
        end
        return ll
    end
    function bilinear_exponential_ll(x::AbstractVector, ages::AbstractVector{BilinearExponential{T}}) where T
        ll = zero(T)
        @inbounds for i ∈ eachindex(x,ages)
            pᵢ = ages[i]
            xs = (x[i]-pᵢ.μ)/(pᵢ.σ) # X scaled by mean and variance
            v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
            shp, skw = pᵢ.sharpness, pᵢ.skew
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
            ll += linterp_at_index(view(p,:,i), x[i], -Inf)
        end
        return ll
    end
    function interpolate_ll(x::AbstractVector,ages::AbstractVector{Radiocarbon{T}}) where T
        ll = zero(T)
        @inbounds for i ∈ eachindex(x,ages)
            dist = ages[i].dist
            ll += linterp_at_index(dist, x[i], -Inf)
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


    struct Binormal{T}
        A::T
        μ::T
        σ::T
        sharpness::T
        skew::T
    end
    function Binormal(p::AbstractVector)
        @assert length(p) == 5
        @assert all(x->x>0, p)
        BilinearExponential(p...)
    end

    function binormal(x::Number, p::AbstractVector{<:Number})
        @assert length(p) == 5
        A, μ, σ, shp, skw = (abs(x) for x in p)
        xs = (x - μ)/σ # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        return A * exp(-0.5*(shp*skw*xs*v - shp/skw*xs*(1-v))^2)
    end
    function binormal(x::AbstractVector, p::AbstractVector{<:Number})
        @assert length(p) == 5
        A, μ, σ, shp, skw = (abs(x) for x in p)
        result = Array{float(eltype(x))}(undef,size(x))
        @inbounds for i ∈ eachindex(x)
            xs = (x[i] - μ)/σ # X scaled by mean and variance
            v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
            result[i] = A * exp(-0.5*(shp*skw*xs*v - shp/skw*xs*(1-v))^2)
        end
        return result
    end


    function binormal_ll(x::Number, p::AbstractVector{<:Number})
        xs = (x - p[2])/abs(p[3]) # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        shp = abs(p[4])
        skw = abs(p[5])
        return -(shp*skw*xs*v - shp/skw*xs*(1-v))^2
    end
    function binormal_ll(x::AbstractVector, ages::AbstractVector{Binormal{T}}) where T
        ll = zero(T)
        @inbounds for i ∈ eachindex(x,ages)
            pᵢ = ages[i]
            xs = (x[i]-pᵢ.μ)/(pᵢ.σ) # X scaled by mean and variance
            v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
            shp = pᵢ.sharpness
            skw = pᵢ.skew
            ll -= (shp*skw*xs*v - shp/skw*xs*(1-v))^2
        end
        return ll
    end
