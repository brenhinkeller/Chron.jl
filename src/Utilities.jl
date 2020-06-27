## --- Constructing arrays

    # Construct linearly spaced array with n points between l and u
    # (linspace replacement )
    if VERSION>=v"0.7"
        function linsp(l::Number,u::Number,n::Number)
            return range(l,stop=u,length=n)
        end
    else
        function linsp(l::Number,u::Number,n::Number)
            return linspace(l,u,n)
        end
    end

## --- Weighted means

    """
    ```julia
    (wx, wσ, mswd) = awmean(x, σ)
    ```
    Weighted mean, absent the geochonologist's MSWD correction to uncertainty.
    """
    function awmean(x, σ)
        n = length(x)

        sum_of_values = sum_of_weights = χ2 = 0.0
        @inbounds @simd for i=1:n
            sum_of_values += x[i] / (σ[i]*σ[i])
            sum_of_weights += 1 / (σ[i]*σ[i])
        end
        wx = sum_of_values / sum_of_weights

        @inbounds @simd for i=1:n
            χ2 += (x[i] - wx) * (x[i] - wx) / (σ[i] * σ[i])
        end
        mswd = χ2 / (n-1)
        wσ = sqrt(1.0 / sum_of_weights)
        return wx, wσ, mswd
    end
    function awmean(x::Array{<:Number}, σ::Array{<:Number})
        n = length(x)

        sum_of_values = sum_of_weights = χ2 = 0.0
        @avx for i=1:n
            sum_of_values += x[i] / (σ[i]*σ[i])
            sum_of_weights += 1 / (σ[i]*σ[i])
        end
        wx = sum_of_values / sum_of_weights

        @avx for i=1:n
            χ2 += (x[i] - wx) * (x[i] - wx) / (σ[i] * σ[i])
        end
        mswd = χ2 / (n-1)
        wσ = sqrt(1.0 / sum_of_weights)
        return wx, wσ, mswd
    end

    """
    ```julia
    (wx, wσ, mswd) = gwmean(x, σ)
    ```
    Geochronologist's weighted mean, with "MSWD correction" to uncertainty,
    i.e., wσ is increased by a factor of sqrt(mswd)
    """
    function gwmean(x, σ)
        n = length(x)

        sum_of_values = sum_of_weights = χ2 = 0.0
        @inbounds @simd for i=1:n
            sum_of_values += x[i] / (σ[i]*σ[i])
            sum_of_weights += 1 / (σ[i]*σ[i])
        end
        wx = sum_of_values / sum_of_weights

        @inbounds @simd for i=1:n
            χ2 += (x[i] - wx) * (x[i] - wx) / (σ[i] * σ[i])
        end
        mswd = χ2 / (n-1)
        wσ = sqrt(mswd / sum_of_weights)
        return wx, wσ, mswd
    end
    function gwmean(x::Array{<:Number}, σ::Array{<:Number})
        n = length(x)
        sum_of_values = sum_of_weights = χ2 = 0.0
        @avx for i=1:n
            sum_of_values += x[i] / (σ[i]*σ[i])
            sum_of_weights += 1 / (σ[i]*σ[i])
        end
        wx = sum_of_values / sum_of_weights

        @avx for i=1:n
            χ2 += (x[i] - wx) * (x[i] - wx) / (σ[i] * σ[i])
        end
        mswd = χ2 / (n-1)
        wσ = sqrt(mswd / sum_of_weights)
        return wx, wσ, mswd
    end

    """
    ```julia
    MSWD(x, σ)
    ```
    Return the Mean Square of Weighted Deviates (AKA the reduced chi-squared
    statistic) of a dataset with values `x` and one-sigma uncertainties `σ`
    """
    function MSWD(x, σ)
        sum_of_values = sum_of_weights = χ2 = 0.0
        n = length(x)

        @inbounds @simd for i=1:n
            w = 1 / (σ[i]*σ[i])
            sum_of_values += w * x[i]
            sum_of_weights += w
        end
        wx = sum_of_values / sum_of_weights

        @inbounds @simd for i=1:n
            χ2 += (x[i] - wx) * (x[i] - wx) / (σ[i] * σ[i])
        end

        return χ2 / (n-1)
    end
    function MSWD(x::Array{<:Number}, σ::Array{<:Number})
        sum_of_values = sum_of_weights = χ2 = 0.0
        n = length(x)

        @avx for i=1:n
            w = 1 / (σ[i]*σ[i])
            sum_of_values += w * x[i]
            sum_of_weights += w
        end
        wx = sum_of_values / sum_of_weights

        @avx for i=1:n
            χ2 += (x[i] - wx) * (x[i] - wx) / (σ[i] * σ[i])
        end

        return χ2 / (n-1)
    end

## --- Interpolation

    # Linearly interpolate vector y at index i, returning outboundsval if outside of bounds
    function linterp_at_index(y::AbstractArray, i::Number, outboundsval=float(eltype(y))(NaN))
        if i > 1 && i < length(y)
            i_below = floor(Int, i)
            i_above = i_below + 1
            f = i - i_below
            return @inbounds Float64(f*y[i_above] + (1-f)*y[i_below])
        else
            return Float64(outboundsval)
        end
    end

    # Interpolate y-value at xq
    # Linear interpolation, sorting inputs
    if VERSION>v"0.7"
        function linterp1(x,y,xq)
            itp = LinearInterpolation(x,y, extrapolation_bc = Line())
            yq = itp(xq) # Interpolate value of y at queried x values
            return yq
        end
    else
        function linterp1(x,y,xq)
            itp = interpolate((x,),y, Gridded(Linear()))
            yq = itp[xq] # Interpolate value of y at queried x values
            return yq
        end
    end

    # Sort x and interpolate y-value at xq
    if VERSION>v"0.7"
        function linterp1s(x,y,xq)
            sI = sortperm(x) # indices to construct sorted array
            itp = LinearInterpolation(x[sI], y[sI], extrapolation_bc = Line())
            yq = itp(xq) # Interpolate value of y at queried x values
            return yq
        end
    else
        function linterp1s(x,y,xq)
            sI = sortperm(x) # indices to construct sorted array
            itp = interpolate((x[sI],), y[sI], Gridded(Linear()))
            yq = itp[xq] # Interpolate value of y at queried x values
            return yq
        end
    end

## --- Working with Gaussian distributions

    """
    ```julia
    normpdf(mu,sigma,x)
    ```
    Probability density function of the Normal (Gaussian) distribution

    ``ℯ^{-(x-μ)^2 / (2σ^2)} / σ√2π``

    with mean `mu` and standard deviation `sigma`, evaluated at `x`
    """
    normpdf(mu,sigma,x) = @. exp(-(x-mu)*(x-mu) / (2*sigma*sigma)) / (sqrt(2*pi)*sigma)

    """
    ```julia
    normpdf_ll(mu, sigma, x)
    ```
    Fast log likelihood corresponding to a Normal (Gaussian) distribution
    with mean `mu` and standard deviation `sigma`, evaluated at `x`.

    If `x`, `mu`, and `sigma` are given as arrays, the sum of the log likelihood
    over all `x` will be returned.

    See also `normpdf`
    """
    normpdf_ll(mu::Number,sigma::Number,x::Number) = -(x-mu)*(x-mu) / (2*sigma*sigma)
    function normpdf_ll(mu::Number,sigma::Number,x::AbstractArray)
        ll = 0.0
        inv_s2 = 1/(2*sigma*sigma)
        @avx for i=1:length(x)
            ll -= (x[i]-mu)*(x[i]-mu) * inv_s2
        end
        return ll
    end
    function normpdf_ll(mu::AbstractArray,sigma::AbstractArray,x::AbstractArray)
        ll = 0.0
        @avx for i=1:length(x)
            ll -= (x[i]-mu[i])*(x[i]-mu[i]) / (2*sigma[i]*sigma[i])
        end
        return ll
    end

    """
    ```julia
    normcdf(mu,sigma,x)
    ```
    Cumulative density function of the Normal (Gaussian) distribution

    ``1/2 + erf(\frac{x-μ}{σ√2})/2``

    with mean `mu` and standard deviation `sigma`, evaluated at `x`.
    """
    function normcdf(mu,sigma,x)
        return @. 0.5 + erf((x-mu) / (sigma*sqrt(2))) / 2
    end
    function normcdf(mu::Number,sigma::Number,x::AbstractArray)
        result = Array{float(eltype(x))}(undef,length(x))
        sigma_sqrt = sigma*sqrt(2)
        @inbounds @simd for i = 1:length(x)
            result[i] = 0.5 + erf((x[i]-mu) / sigma_sqrt) / 2
        end
        return result
    end

    # How far away from the mean (in units of sigma) should we expect proportion
    # F of the samples to fall in a Normal (Gaussian) distribution
    function norm_quantile(F::Number)
        return sqrt(2)*erfinv(2*F-1)
    end

    # How dispersed (in units of sigma) should we expect a sample of N numbers
    # drawn from a Normal (Gaussian) distribution to be?
    function norm_width(N)
        F = 1 - 1/(N+1)
        return 2*norm_quantile(F)
    end

    # Integral of the product of two normal distributions N(μ1,σ1) * N(μ2,σ2)
    function normproduct(μ1::Number, σ1::Number, μ2::Number, σ2::Number)
        # The integral of the product of two normal distributions is itself just
        # another Normal distribution! Specifically, one with variance σ1^2 + σ2^2
        normpdf(μ1, sqrt(σ1^2 + σ2^2), μ2)
    end

    # Log likelihood corresponding to the integral of N(μ1,σ1) * N(μ2,σ2)
    function normproduct_ll(μ1::Number, σ1::Number, μ2::Number, σ2::Number)
        # As above, but using the fast log likelihood of a Normal distribution
        normpdf_ll(μ1, sqrt(σ1^2 + σ2^2), μ2)
    end


## --- Drawing from distributions

    """
    ```julia
    x = draw_from_distribution(dist::AbstractArray{<:AbstractFloat}, n::Integer)
    ```
    Draw `n` random floating point numbers from a continuous probability distribution
    specified by a vector `dist` defining the PDF curve thereof.
    """
    function draw_from_distribution(dist::AbstractArray{<:AbstractFloat}, n::Integer)
        # Draw n random floating-point numbers from the distribution 'dist'
        x = Array{eltype(dist)}(undef, n)
        dist_ymax = maximum(dist)
        dist_xmax = prevfloat(length(dist) - 1.0)

        @inbounds for i = 1:n
            while true
                # Pick random x value
                rx = rand(Float64) * dist_xmax
                # Interpolate corresponding distribution value
                f = floor(Int,rx)
                y = dist[f+2]*(rx-f) + dist[f+1]*(1-(rx-f))
                # See if x value is accepted
                ry = rand(Float64) * dist_ymax
                if (y > ry)
                    x[i] = rx / dist_xmax
                    break
                end
            end
        end
        return x
    end
    export draw_from_distribution

    """
    ```julia
    draw_from_distribution!(dist::AbstractArray{<:AbstractFloat}, x::Array{<:AbstractFloat})
    ```
    Fill an existing variable `x` with random floating point numbers drawn from
    a continuous probability distribution specified by a vector `dist`
    defining the PDF curve thereof.
    """
    function draw_from_distribution!(x::Array{<:AbstractFloat}, dist::AbstractArray{<:AbstractFloat})
        # Fill the array x with random numbers from the distribution 'dist'
        dist_ymax = maximum(dist)
        dist_xmax = prevfloat(length(dist) - 1.0)

        @inbounds for i=1:length(x)
            while true
                # Pick random x value
                rx = rand(eltype(x)) * dist_xmax
                # Interpolate corresponding distribution value
                f = floor(Int,rx)
                y = dist[f+2]*(rx-f) + dist[f+1]*(1-(rx-f))
                # See if x value is accepted
                ry = rand(Float64) * dist_ymax
                if (y > ry)
                    x[i] = rx / dist_xmax
                    break
                end
            end
        end
    end

## --- Fitting non-Gaussian distributions

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
    function bilinear_exponential(x::Number, p::AbstractVector)
        xs = (x - p[2])/p[3]^2 # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        return exp(p[1] + (p[4]^2)*(p[5]^2)*xs*v - (p[4]^2)/(p[5]^2)*xs*(1-v))
    end
    function bilinear_exponential(x::AbstractVector, p::AbstractVector)
        f = Array{float(eltype(x))}(undef,size(x))
        @avx for i = 1:length(x)
            xs = (x[i] - p[2])/p[3]^2 # X scaled by mean and variance
            v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
            f[i] = exp(p[1] + (p[4]^2)*(p[5]^2)*xs*v - (p[4]^2)/(p[5]^2)*xs*(1-v))
        end
        return f
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
    function bilinear_exponential_ll(x::Number, p::AbstractVector)
        xs = (x - p[2])/p[3]^2 # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        return p[1] + (p[4]^2)*(p[5]^2)*xs*v - (p[4]^2)/(p[5]^2)*xs*(1-v)
    end
    function bilinear_exponential_ll(x::AbstractVector, p::AbstractMatrix)
        ll = 0.0
        @avx for i=1:length(x)
            xs = (x[i]-p[2,i])/p[3,i]^2 # X scaled by mean and variance
            v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
            ll += p[1,i] + (p[4,i]^2)*(p[5,i]^2)*xs*v - (p[4,i]^2)/(p[5,i]^2)*xs*(1 - v)
        end
        return ll
    end

    # Interpolate log likelihood from an array
    function interpolate_ll(x::AbstractVector,p::AbstractMatrix)
        ll = 0
        @inbounds for i = 1:length(x)
            ll += linterp_at_index(view(p,:,i), x[i], -Inf)
        end
        return ll
    end

## --- Various other utility functions

    """
    ```julia
    cntr(edges::AbstractArray)
    ```
    Given an array of bin edges, return a corresponding vector of bin centers
    """
    function cntr(edges::AbstractArray)
        centers = (edges[1:end-1] + edges[2:end]) ./ 2
        return centers
    end

    """
    ```julia
    nanmask(A)
    ```
    Create a Boolean mask of dimensions `size(A)` that is false wherever `A` is `NaN`
    """
    nanmask(A) = nanmask!(Array{Bool}(undef,size(A)), A)
    export nanmask

    """
    ```julia
    nanmask!(mask, A)
    ```
    Fill a Boolean mask of dimensions `size(A)` that is false wherever `A` is `NaN`
    """
    function nanmask!(mask, A)
        @avx for i=1:length(A)
            mask[i] = !isnan(A[i])
        end
        return mask
    end
    # Special methods for arrays that cannot contain NaNs
    nanmask!(mask, A::AbstractArray{<:Integer}) = fill!(mask, true)
    nanmask!(mask, A::AbstractArray{<:Rational}) = fill!(mask, true)

    """
    ```julia
    nanmax(a,b)
    ```
    As `max(a,b)`, but if either argument is `NaN`, return the other one
    """
    nanmax(a, b) = ifelse(a > b, a, b)
    nanmax(a, b::AbstractFloat) = ifelse(a==a, ifelse(b > a, b, a), b)
    nanmax(a::SVec{N,<:Integer}, b::SVec{N,<:Integer}) where N = vifelse(a > b, a, b)
    nanmax(a::SVec{N,<:AbstractFloat}, b::SVec{N,<:AbstractFloat}) where N = vifelse(a==a, vifelse(b > a, b, a), b)

    """
    ```julia
    nanmin(a,b)
    ```
    As `min(a,b)`, but if either argument is `NaN`, return the other one
    """
    nanmin(a, b) = ifelse(a < b, a, b)
    nanmin(a, b::AbstractFloat) = ifelse(a==a, ifelse(b < a, b, a), b)
    nanmin(a::SVec{N,<:Integer}, b::SVec{N,<:Integer}) where N = vifelse(a < b, a, b)
    nanmin(a::SVec{N,<:AbstractFloat}, b::SVec{N,<:AbstractFloat}) where N = vifelse(a==a, vifelse(b < a, b, a), b)

    """
    ```julia
    pctile(A, p; dims)
    ```
    Find the `p`th percentile of an indexable collection `A`, ignoring NaNs,
    optionally along a dimension specified by `dims`.

    A valid percentile value must satisfy 0 <= `p` <= 100.
    """
    pctile(A, p; dims=:, dim=:) = __pctile(A, p, dims, dim)
    __pctile(A, p, dims, dim) = _pctile(A, p, dim) |> vec
    __pctile(A, p, dims, ::Colon) = _pctile(A, p, dims)
    function _pctile(A, p, ::Colon)
        t = nanmask(A)
        return any(t) ? percentile(A[t],p) : NaN
    end
    function _pctile(A, p, region)
        s = size(A)
        if region == 2
            t = Array{Bool}(undef, s[2])
            result = Array{float(eltype(A))}(undef, s[1], 1)
            for i=1:s[1]
                nanmask!(t, A[i,:])
                result[i] = any(t) ? percentile(A[i,t],p) : NaN
            end
        elseif region == 1
            t = Array{Bool}(undef, s[1])
            result = Array{float(eltype(A))}(undef, 1, s[2])
            for i=1:s[2]
                nanmask!(t, A[:,i])
                result[i] = any(t) ? percentile(A[t,i],p) : NaN
            end
        else
            result = _pctile(A, p, :)
        end
        return result
    end


    """
    ```julia
    nanminimum(A; dims)
    ```
    As `minimum` but ignoring `NaN`s: Find the smallest non-`NaN` value of an
    indexable collection `A`, optionally along a dimension specified by `dims`.
    """
    nanminimum(A; dims=:, dim=:) = __nanminimum(A, dims, dim)
    __nanminimum(A, dims, dim) = _nanminimum(A, dim) |> vec
    __nanminimum(A, dims, ::Colon) = _nanminimum(A, dims)
    _nanminimum(A, region) = reduce(nanmin, A, dims=region, init=float(eltype(A))(NaN))
    _nanminimum(A::AbstractArray{<:Number}, ::Colon) = vreduce(nanmin, A)

    """
    ```julia
    nanmaximum(A; dims)
    ```
    Find the largest non-NaN value of an indexable collection `A`, optionally
    along a dimension specified by `dims`.
    """
    nanmaximum(A; dims=:, dim=:) = __nanmaximum(A, dims, dim)
    __nanmaximum(A, dims, dim) = _nanmaximum(A, dim) |> vec
    __nanmaximum(A, dims, ::Colon) = _nanmaximum(A, dims)
    _nanmaximum(A, region) = reduce(nanmax, A, dims=region, init=float(eltype(A))(NaN))
    _nanmaximum(A::AbstractArray{<:Number}, ::Colon) = vreduce(nanmax, A)

    """
    ```julia
    nanextrema(A; dims)
    ```
    Find the extrema (maximum & minimum) of an indexable collection `A`,
    ignoring NaNs, optionally along a dimension specified by `dims`.
    """
    nanextrema(A; dims=:) = _nanextrema(A, dims)
    _nanextrema(A, region) = collect(zip(_nanminimum(A, region), _nanmaximum(A, region)))
    _nanextrema(A, ::Colon) = (_nanminimum(A, :), _nanmaximum(A, :))

    """
    ```julia
    nanrange(A; dims)
    ```
    Calculate the range (maximum - minimum) of an indexable collection `A`,
    ignoring NaNs, optionally along a dimension specified by `dims`.
    """
    nanrange(A; dims=:) = _nanmaximum(A, dims) - _nanminimum(A, dims)


    """
    ```julia
    nanmean(A, [W]; dims)
    ```
    Ignoring NaNs, calculate the mean (optionally weighted) of an indexable
    collection `A`, optionally along dimensions specified by `dims`.
    """
    nanmean(A; dims=:, dim=:) = __nanmean(A, dims, dim)
    __nanmean(A, dims, dim) = _nanmean(A, dim) |> vec
    __nanmean(A, dims, ::Colon) = _nanmean(A, dims)
    function _nanmean(A, region)
        mask = nanmask(A)
        return sum(A.*mask, dims=region) ./ sum(mask, dims=region)
    end
    # Fallback method for non-Arrays
    function _nanmean(A, ::Colon)
        n = 0
        m = zero(eltype(A))
        @inbounds @simd for i ∈ eachindex(A)
            Aᵢ = A[i]
            t = Aᵢ == Aᵢ
            n += t
            m += Aᵢ * t
        end
        return m / n
    end
    # Can't have NaNs if array is all Integers
    function _nanmean(A::Array{<:Integer}, ::Colon)
        m = zero(eltype(A))
        @avx for i ∈ eachindex(A)
            m += A[i]
        end
        return m / length(A)
    end
    # Optimized AVX version for floats
    function _nanmean(A::AbstractArray{<:AbstractFloat}, ::Colon)
        n = 0
        m = zero(eltype(A))
        @avx for i ∈ eachindex(A)
            Aᵢ = A[i]
            t = Aᵢ == Aᵢ
            n += t
            m += Aᵢ * t
        end
        return m / n
    end
    nanmean(A, W; dims=:, dim=:) = __nanmean(A, W, dims, dim)
    __nanmean(A, W, dims, dim) = _nanmean(A, W, dim) |> vec
    __nanmean(A, W, dims, ::Colon) = _nanmean(A, W, dims)
    function _nanmean(A, W, region)
        mask = nanmask(A)
        return sum(A.*W.*mask, dims=region) ./ sum(W.*mask, dims=region)
    end
    # Fallback method for non-Arrays
    function _nanmean(A, W, ::Colon)
        n = zero(eltype(W))
        m = zero(promote_type(eltype(W), eltype(A)))
        @inbounds @simd for i ∈ eachindex(A)
            Aᵢ = A[i]
            Wᵢ = W[i]
            t = Aᵢ == Aᵢ
            n += Wᵢ * t
            m += Wᵢ * Aᵢ * t
        end
        return m / n
    end
    # Can't have NaNs if array is all Integers
    function _nanmean(A::Array{<:Integer}, W, ::Colon)
        n = zero(eltype(W))
        m = zero(promote_type(eltype(W), eltype(A)))
        @avx for i ∈ eachindex(A)
            Wᵢ = W[i]
            n += Wᵢ
            m += Wᵢ * A[i]
        end
        return m / n
    end
    # Optimized AVX method for floats
    function _nanmean(A::AbstractArray{<:AbstractFloat}, W, ::Colon)
        n = zero(eltype(W))
        m = zero(promote_type(eltype(W), eltype(A)))
        @avx for i ∈ eachindex(A)
            Aᵢ = A[i]
            Wᵢ = W[i]
            t = Aᵢ == Aᵢ
            n += Wᵢ * t
            m += Wᵢ * Aᵢ * t
        end
        return m / n
    end


    """
    ```julia
    nanstd(A, [W]; dims)
    ```
    Calculate the standard deviation (optionaly weighted), ignoring NaNs, of an
    indexable collection `A`, optionally along a dimension specified by `dims`.
    """
    nanstd(A; dims=:, dim=:) = __nanstd(A, dims, dim)
    __nanstd(A, dims, dim) = _nanstd(A, dim) |> vec
    __nanstd(A, dims, ::Colon) = _nanstd(A, dims)
    function _nanstd(A, region)
        mask = nanmask(A)
        N = sum(mask, dims=region)
        s = sum(A.*mask, dims=region)./N
        d = A .- s # Subtract mean, using broadcasting
        @avx for i ∈ eachindex(d)
            dᵢ = d[i]
            d[i] = (dᵢ * dᵢ) * mask[i]
        end
        s .= sum(d, dims=region)
        @avx for i ∈ eachindex(s)
            s[i] = sqrt( s[i] / max((N[i] - 1), 0) )
        end
        return s
    end
    function _nanstd(A, ::Colon)
        n = 0
        m = zero(eltype(A))
        @inbounds @simd for i ∈ eachindex(A)
            Aᵢ = A[i]
            t = Aᵢ == Aᵢ # False for NaNs
            n += t
            m += Aᵢ * t
        end
        mu = m / n
        s = zero(typeof(mu))
        @inbounds @simd for i ∈ eachindex(A)
            Aᵢ = A[i]
            d = (Aᵢ - mu) * (Aᵢ == Aᵢ)# zero if Aᵢ is NaN
            s += d * d
        end
        return sqrt(s / max((n-1), 0))
    end
    function _nanstd(A::AbstractArray{<:AbstractFloat}, ::Colon)
        n = 0
        m = zero(eltype(A))
        @avx for i ∈ eachindex(A)
            Aᵢ = A[i]
            t = Aᵢ == Aᵢ # False for NaNs
            n += t
            m += Aᵢ * t
        end
        mu = m / n
        s = zero(typeof(mu))
        @avx for i ∈ eachindex(A)
            Aᵢ = A[i]
            d = (Aᵢ - mu) * (Aᵢ == Aᵢ)# zero if Aᵢ is NaN
            s += d * d
        end
        return sqrt(s / max((n-1), 0))
    end

    nanstd(A, W; dims=:, dim=:) = __nanstd(A, W, dims, dim)
    __nanstd(A, W, dims, dim) = _nanstd(A, W, dim) |> vec
    __nanstd(A, W, dims, ::Colon) = _nanstd(A, W, dims)
    function _nanstd(A, W, region)
        mask = nanmask(A)
        w = sum(W.*mask, dims=region)
        s = sum(A.*W.*mask, dims=region) ./ w
        d = A .- s # Subtract mean, using broadcasting
        @avx for i ∈ eachindex(d)
            dᵢ = d[i]
            d[i] = (dᵢ * dᵢ * W[i]) * mask[i]
        end
        s .= sum(d, dims=region)
        @avx for i ∈ eachindex(s)
            s[i] = sqrt( s[i] / w[i] )
        end
        return s
    end
    function _nanstd(A, W, ::Colon)
        w = zero(eltype(W))
        m = zero(promote_type(eltype(W), eltype(A)))
        @inbounds @simd for i ∈ eachindex(A)
            Aᵢ = A[i]
            Wᵢ = W[i]
            t = Aᵢ == Aᵢ
            w += Wᵢ * t
            m += Wᵢ * Aᵢ * t
        end
        mu = m / w
        s = zero(typeof(mu))
        @inbounds @simd for i ∈ eachindex(A)
            Aᵢ = A[i]
            d = Aᵢ - mu
            s += (d * d * W[i]) * (Aᵢ == Aᵢ) # Zero if Aᵢ is NaN
        end
        return sqrt(s / w)
    end
    function _nanstd(A::AbstractArray{<:AbstractFloat}, W, ::Colon)
        w = zero(eltype(W))
        m = zero(promote_type(eltype(W), eltype(A)))
        @avx for i ∈ eachindex(A)
            Aᵢ = A[i]
            Wᵢ = W[i]
            t = Aᵢ == Aᵢ
            w += Wᵢ * t
            m += Wᵢ * Aᵢ * t
        end
        mu = m / w
        s = zero(typeof(mu))
        @avx for i ∈ eachindex(A)
            Aᵢ = A[i]
            d = Aᵢ - mu
            s += (d * d * W[i]) * (Aᵢ == Aᵢ) # Zero if Aᵢ is NaN
        end
        return sqrt(s / w)
    end


    """
    ```julia
    nanmedian(A; dims)
    ```
    Calculate the median, ignoring NaNs, of an indexable collection `A`,
    optionally along a dimension specified by `dims`.
    """
    nanmedian(A; dims=:, dim=:) = _nanmedian(A, dims, dim)
    _nanmedian(A, region, ::Colon) = _nanmedian(A, region)
    _nanmedian(A, ::Colon, region) = _nanmedian(A, region) |> vec
    _nanmedian(A, ::Colon, ::Colon) = _nanmedian(A, :)
    function _nanmedian(A, ::Colon)
        t = nanmask(A)
        return any(t) ? median(A[t]) : float(eltype(A))(NaN)
    end
    function _nanmedian(A, region)
        s = size(A)
        if region == 2
            t = Array{Bool}(undef, s[2])
            result = Array{float(eltype(A))}(undef, s[1], 1)
            for i=1:s[1]
                nanmask!(t, A[i,:])
                result[i] = any(t) ? median(A[i,t]) : float(eltype(A))(NaN)
            end
        elseif region == 1
            t = Array{Bool}(undef, s[1])
            result = Array{float(eltype(A))}(undef, 1, s[2])
            for i=1:s[2]
                nanmask!(t, A[:,i])
                result[i] = any(t) ? median(A[t,i]) : float(eltype(A))(NaN)
            end
        else
            result = _nanmedian(A, :)
        end
        return result
    end

    """
    ```julia
    findclosest(source, target)
    ```
    Return the index of the numerically closest value in the indexable collection
    `target` for each value in `source`.
    If muliple values are equally close, the first one is used
    """
    function findclosest(source, target)
        index = Array{Int64}(undef, size(source))
        return findclosest!(index, source, target)
    end
    function findclosest!(index, source, target)
        # Find closest (numerical) match in target for each value in source
        @inbounds for i = 1:length(source)
            d = abs(target[1] - source[i])
            index[i] = 1
            for j = 2:length(target)
                d_prop = abs(target[j] - source[i])
                if d_prop < d
                    d = d_prop
                    index[i] = j
                end
            end
        end
        return index
    end

    """
    ```julia
    findclosestbelow(source, target)
    ```
    Return the index of the nearest value of the indexable collection `target`
    that is less than (i.e., "below") each value in `source`.
    If no such target values exist in `target`, returns an index of 0.
    """
    function findclosestbelow(source, target)
        index = Array{Int64}(undef, size(source))
        return findclosestbelow!(index, source, target)
    end
    function findclosestbelow!(index, source, target)
        @inbounds for i = 1:length(source)
            index[i] = d = j = 0
            while j < length(target)
                j += 1
                if target[j] < source[i]
                    d = source[i] - target[j]
                    index[i] = j
                    break
                end
            end
            while j < length(target)
                j += 1
                if target[j] < source[i]
                    d_prop = source[i] - target[j]
                    if d_prop < d
                        d = d_prop
                        index[i] = j
                    end
                end
            end
        end
        return index
    end

    """
    ```julia
    findclosestabove(source, target)
    ```
    Return the index of the nearest value of the indexable collection `target`
    that is greater than (i.e., "above") each value in `source`.
    If no such values exist in `target`, returns an index of 0.
    """
    function findclosestabove(source, target)
        index = Array{Int64}(undef, size(source))
        return findclosestabove!(index,source,target)
    end
    function findclosestabove!(index, source, target)
        @inbounds for i = 1:length(source)
            index[i] = d = j = 0
            while j < length(target)
                j += 1
                if target[j] > source[i]
                    d = target[j] - source[i]
                    index[i] = j
                    break
                end
            end
            while j < length(target)
                j += 1
                if target[j] > source[i]
                    d_prop = target[j] - source[i]
                    if d_prop < d
                        d = d_prop
                        index[i] = j
                    end
                end
            end
        end
        return index
    end


## --- Other

    # Direct system() access without stripping special characters like Julia's
    # rather protective "run()".
    function system(cmdstr)
        return ccall((:system,), Int, (Cstring,), cmdstr)
    end
    export system

    if VERSION>=v"1.0"
        function linspace(l::Number,u::Number,n::Number)
            return range(l,stop=u,length=n)
        end
        export linspace

        function repmat(A::AbstractArray, vert::Integer)
            return repeat(A, outer=vert)
        end
        function repmat(A::AbstractArray, vert::Integer, horiz::Integer)
            return repeat(A, outer=(vert, horiz))
        end
        export repmat

        function contains(haystack::AbstractString, needle::Union{AbstractString,Regex,AbstractChar})
            return occursin(needle::Union{AbstractString,Regex,AbstractChar}, haystack::AbstractString)
        end
        export contains
    end

## --- End of File
