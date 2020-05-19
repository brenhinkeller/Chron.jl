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
    function normpdf(mu,sigma,x)
        return @. exp(-(x-mu)*(x-mu) / (2*sigma*sigma)) / (sqrt(2*pi)*sigma)
    end

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
    function normpdf_ll(mu::Number,sigma::Number,x::Number)
        return -(x-mu)*(x-mu) / (2*sigma*sigma)
    end
    function normpdf_ll(mu::AbstractArray,sigma::AbstractArray,x::AbstractArray)
        ll = 0.0
        @avx for i=1:length(x)
            ll -= (x[i]-mu[i])*(x[i]-mu[i]) / (2*sigma[i]*sigma[i])
        end
        return ll
    end
    export normpdf_ll

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
    function nanmask(A)
        mask = Array{Bool}(undef,size(A))
        @avx for i=1:length(A)
            mask[i] = !isnan(A[i])
        end
        return mask
    end

    """
    ```julia
    nanmax(a,b)
    ```
    As `max(a,b)`, but if either argument is `NaN`, return the other one
    """
    function nanmax(a::AbstractFloat,b::AbstractFloat)
        ifelse(isnan(a), b, ifelse(a < b, b, a))
    end
    nanmax(a::AbstractFloat,b::Number) = nanmax(promote(a,b)...)
    nanmax(a::Number,b::AbstractFloat) = nanmax(promote(a,b)...)
    nanmax(a,b) = max(a,b) # Fallback method for non-Floats

    """
    ```julia
    nanmin(a,b)
    ```
    As `min(a,b)`, but if either argument is `NaN`, return the other one
    """
    function nanmin(a::AbstractFloat,b::AbstractFloat)
        ifelse(isnan(a), b, ifelse(a > b, b, a))
    end
    nanmin(a::AbstractFloat,b::Number) = nanmin(promote(a,b)...)
    nanmin(a::Number,b::AbstractFloat) = nanmin(promote(a,b)...)
    nanmin(a,b) = min(a,b) # Fallback method for non-Floats

    """
    ```julia
    pctile(A, p; dims)
    ```
    Find the `p`th percentile of an indexable collection `A`, ignoring NaNs,
    optionally along a dimension specified by `dims`.

    A valid percentile value must satisfy 0 <= `p` <= 100.
    """
    pctile(A, p; dims=:, dim=:) = _pctile(A, p, dims, dim)
    _pctile(A, p, region, ::Colon) = _pctile(A, p, region)
    _pctile(A, p, ::Colon, region) = _pctile(A, p, region) |> vec
    _pctile(A, p, ::Colon, ::Colon) = _pctile(A, p, :)
    function _pctile(A, p, ::Colon)
        t = .~ isnan.(A)
        return any(t) ? percentile(A[t],p) : NaN
    end
    function _pctile(A, p, region)
        s = size(A)
        if region == 2
            result = Array{float(eltype(A))}(undef, s[1], 1)
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? percentile(A[i,t],p) : NaN
            end
        elseif region == 1
            result = Array{float(eltype(A))}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                result[i] = any(t) ? percentile(A[t,i],p) : NaN
            end
        else
            result =  _pctile(A, p, :)
        end
        return result
    end

    """
    ```julia
    nansum(A; dims)
    ```
    Calculate the sum of an indexable collection `A`, ignoring NaNs, optionally
    along dimensions specified by `dims`.
    """
    nansum(A; dims=:, dim=:) = _nansum(A, dims, dim)
    _nansum(A, region, ::Colon) = _nansum(A, region)
    _nansum(A, ::Colon, region) = _nansum(A, region) |> vec
    _nansum(A, ::Colon, ::Colon) = _nansum(A, :)
    function _nansum(A,::Colon)
        m = zero(eltype(A))
        @simd for x in A
            m += x * !isnan(x)
        end
        return m
    end
    function _nansum(A, region)
        sum(A.*nanmask(A), dims=region)
    end

    """
    ```julia
    nanminimum(A; dims)
    ```
    As `minimum` but ignoring `NaN`s: Find the smallest non-`NaN` value of an
    indexable collection `A`, optionally along a dimension specified by `dims`.
    """
    nanminimum(A; dims=:, dim=:) = _nanminimum(A, dims, dim)
    _nanminimum(A, region, ::Colon) = _nanminimum(A, region)
    _nanminimum(A, ::Colon, region) = _nanminimum(A, region) |> vec
    _nanminimum(A, ::Colon, ::Colon) = _nanminimum(A, :)
    function _nanminimum(A, ::Colon)
        result = A[1]
        for x in A
            result = nanmin(x, result)
        end
        return result
    end
    _nanminimum(A, region) = reduce(nanmin, A, dims=region, init=float(eltype(A))(NaN))


    """
    ```julia
    nanmaximum(A; dims)
    ```
    Find the largest non-NaN value of an indexable collection `A`, optionally
    along a dimension specified by `dims`.
    """
    nanmaximum(A; dims=:, dim=:) = _nanmaximum(A, dims, dim)
    _nanmaximum(A, region, ::Colon) = _nanmaximum(A, region)
    _nanmaximum(A, ::Colon, region) = _nanmaximum(A, region) |> vec
    _nanmaximum(A, ::Colon, ::Colon) = _nanmaximum(A, :)
    function _nanmaximum(A, ::Colon)
        result = A[1]
        for x in A
            result = nanmax(x,result)
        end
        return result
    end
    _nanmaximum(A, region) = reduce(nanmax, A, dims=region, init=float(eltype(A))(NaN))

    """
    ```julia
    nanextrema(A; dims)
    ```
    Find the extrema (maximum & minimum) of an indexable collection `A`,
    ignoring NaNs, optionally along a dimension specified by `dims`.
    """
    nanextrema(A; dims=:) = _nanextrema(A, dims)
    _nanextrema(A, ::Colon) = (_nanminimum(A, :), _nanmaximum(A, :))
    _nanextrema(A, region) = collect(zip(_nanminimum(A, region), _nanmaximum(A, region)))

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
    nanmean(A; dims=:, dim=:) = _nanmean(A, dims, dim)
    nanmean(A, W; dims=:, dim=:) = _nanmean(A, W, dims, dim)
    function _nanmean(A, ::Colon, ::Colon)
        n = 0
        m = zero(eltype(A))
        @inbounds @simd for i=1:length(A)
            t = !isnan(A[i])
            n += t
            m += A[i] * t
        end
        return m / n
    end
    function _nanmean(A, region, ::Colon)
        mask = nanmask(A)
        return sum(A.*mask, dims=region) ./ sum(mask, dims=region)
    end
    _nanmean(A, ::Colon, region) = vec(_nanmean(A, region, :))
    function _nanmean(A, W, ::Colon, ::Colon)
        n = zero(eltype(W))
        m = zero(promote_type(eltype(W), eltype(A)))
        @inbounds @simd for i=1:length(A)
            t = !isnan(A[i])
            n += W[i] * t
            m += A[i] * W[i] * t
        end
        return m / n
    end
    function _nanmean(A, W, region, ::Colon)
        mask = nanmask(A)
        return sum(A.*W.*mask, dims=region) ./ sum(W.*mask, dims=region)
    end
    _nanmean(A, W, ::Colon, region) = vec(_nanmean(A, W, region, :))


    """
    ```julia
    nanstd(A, [W]; dims)
    ```
    Calculate the standard deviation (optionaly weighted), ignoring NaNs, of an
    indexable collection `A`, optionally along a dimension specified by `dims`.
    """
    nanstd(A; dims=:, dim=:) = _nanstd(A, dims, dim)
    nanstd(A, W; dims=:, dim=:) = _nanstd(A, W, dims, dim)
    function _nanstd(A, ::Colon, ::Colon)
        n = 0
        m = zero(eltype(A))
        @inbounds @simd for i=1:length(A)
            t = !isnan(A[i])
            n += t
            m += A[i] * t
        end
        mu = m / n
        s = zero(typeof(mu))
        @inbounds @simd for i=1:length(A)
            d = (A[i] - mu) * !isnan(A[i])
            s += d * d
        end
        return sqrt(s / (n-1))
    end
    function _nanstd(A, region, ::Colon)
        mask = nanmask(A)
        N = sum(mask, dims=region)
        s = sum(A.*mask, dims=region)./N
        d = A .- s # Subtract mean, using broadcasting
        @inbounds @simd for i = 1:length(d)
            d[i] = (d[i] * d[i]) * mask[i]
        end
        s .= sum(d, dims=region)
        return @avx @. sqrt( s / (N - 1) )
    end
    _nanstd(A, ::Colon, region) = vec(_nanstd(A, region, :))
    function _nanstd(A, W, ::Colon, ::Colon)
        w = zero(eltype(W))
        m = zero(promote_type(eltype(W), eltype(A)))
        @inbounds @simd for i=1:length(A)
            t = !isnan(A[i])
            w += W[i] * t
            m += A[i] * W[i] * t
        end
        mu = m / w
        s = zero(typeof(mu))
        @inbounds @simd for i=1:length(A)
            d = (A[i] - mu)
            s += d * d * W[i] * !isnan(A[i])
        end
        return sqrt(s / w)
    end
    function _nanstd(A, W, region, ::Colon)
        mask = nanmask(A)
        w = sum(W.*mask, dims=region)
        s = sum(A.*W.*mask, dims=region) ./ w
        d = A .- s # Subtract mean, using broadcasting
        @inbounds @simd for i = 1:length(d)
            d[i] = (d[i] * d[i] * W[i]) * mask[i]
        end
        s .= sum(d, dims=region)
        return @avx @. sqrt( s / w )
    end
    _nanstd(A, W, ::Colon, region) = vec(_nanstd(A, W, region, :))

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
        t = .~ isnan.(A)
        return any(t) ? median(A[t]) : float(eltype(A))(NaN)
    end
    function _nanmedian(A, region)
        s = size(A)
        if region == 2
            result = Array{float(eltype(A))}(undef, s[1], 1)
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                result[i] = any(t) ? median(A[i,t]) : float(eltype(A))(NaN)
            end
        elseif region == 1
            result = Array{float(eltype(A))}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
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
        # Allocate index and difference arrays
        index = Array{Int64}(undef, size(source))
        diff_type = promote_type(eltype(source), eltype(target))
        diff = Array{diff_type}(undef, length(target))
        # Find closest (numerical) match in target for each value in source
        @inbounds for i = 1:length(source)
            for j = 1:length(diff)
                diff[j] = abs(target[j] - source[i])
            end
            index[i] = argmin(diff)
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
        # Allocate output array
        index = Array{Int64}(undef, size(source))
        diff_type = promote_type(eltype(source), eltype(target))
        diff = Array{diff_type}(undef, length(target))
        t = Array{Bool}(undef,length(target))
        @inbounds for i = 1:length(source)
            j = 0
            closestbelow = 0
            while j < length(diff)
                j += 1
                if target[j] < source[i]
                    diff[j] = source[i] - target[j]
                    closestbelow = j
                    break
                end
            end
            while j < length(diff)
                j += 1
                if target[j] < source[i]
                    diff[j] = source[i] - target[j]
                    if diff[j] < diff[closestbelow]
                        closestbelow = j
                    end
                end
            end
            index[i] = closestbelow
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
        # Allocate output array
        index = Array{Int64}(undef, size(source))
        diff_type = promote_type(eltype(source), eltype(target))
        diff = Array{diff_type}(undef, length(target))
        t = Array{Bool}(undef,length(target))
        @inbounds for i = 1:length(source)
            j = 0
            closestabove = 0
            while j < length(diff)
                j += 1
                if target[j] > source[i]
                    diff[j] = target[j] - source[i]
                    closestabove = j
                    break
                end
            end
            while j < length(diff)
                j += 1
                if target[j] > source[i]
                    diff[j] = target[j] - source[i]
                    if diff[j] < diff[closestabove]
                        closestabove = j
                    end
                end
            end
            index[i] = closestabove
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
