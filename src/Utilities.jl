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

    # Calculate a weigted mean, including MSWD, with MSWD correction to uncertainty.
    function gwmean(x, sigma)
        n = length(x)
        s1 = 0.0; s2 = 0.0; s3 = 0.0

        if n==1
            wx = x[1]
            mswd = 0
            wsigma = sigma[1]
        else
            for i=1:n
                s1 += x[i] / (sigma[i]*sigma[i])
                s2 += 1 / (sigma[i]*sigma[i])
            end
            wx = s1/s2

            for i=1:n
                s3 += (x[i] - wx) * (x[i] - wx) / (sigma[i]*sigma[i])
            end
            mswd = s3 / (n-1)
            wsigma = sqrt(mswd/s2)
        end
        return wx, wsigma, mswd
    end

    # Calculate a weigted mean, including MSWD, but without MSWD correction to uncertainty
    function awmean(x, sigma)
        n = length(x)
        s1 = 0.0; s2 = 0.0; s3 = 0.0

        if n==1
            wx = x[1]
            mswd = 0
            wsigma = sigma[1]
        else
            for i=1:n
                s1 += x[i] / (sigma[i]*sigma[i])
                s2 += 1 / (sigma[i]*sigma[i])
            end
            wx = s1/s2

            for i=1:n
                s3 += (x[i] - wx) * (x[i] - wx) / (sigma[i]*sigma[i])
            end
            mswd = s3 / (n-1)
            wsigma = sqrt(1.0/s2)
        end
        return wx, wsigma, mswd
    end

## --- Interpolation

    # Linearly interpolate vector y at index i, returning outboundsval if out of bounds
    function linterp_at_index(y::AbstractArray, i::Number, outboundsval=NaN)
        if i > 1 && i < length(y)
            i_below = floor(Int, i)
            i_above = i_below + 1
            f = i - i_below
            return @inbounds Float64(f*y[i_above] + (1-f)*y[i_below])
        else
            return Float64(outboundsval)
        end
    end

    # Linear interpolation, sorting inputs
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

    # linear interpolation
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

    # Draw random numbers from a distribution specified by a vector of points
    # defining the PDF curve
    function draw_from_distribution(dist::Array{Float64}, n::Int)
        # Draw n random numbers from the distribution 'dist'
        x = Array{Float64}(undef,n)
        dist_ymax = maximum(dist)
        dist_xmax = length(dist)-1.0

        for i=1:n
            while true
                # Pick random x value
                rx = rand() * dist_xmax
                # Interpolate corresponding distribution value
                f = floor(Int,rx)
                y = dist[f+2]*(rx-f) + dist[f+1]*(1-(rx-f))
                # See if x value is accepted
                ry = rand() * dist_ymax
                if (y > ry)
                    x[i] = rx / dist_xmax
                    break
                end
            end
        end
        return x
    end

    # Fill an existing variable with random numbers from a distribution specified
    # by a vector of points defining the PDF curve
    function fill_from_distribution(dist::Array{Float64}, x::Array{Float64})
        # Fill the array x with random numbers from the distribution 'dist'
        dist_ymax = maximum(dist)
        dist_xmax = length(dist)-1.0

        for i=1:length(x)
            while true
                # Pick random x value
                rx = rand() * dist_xmax
                # Interpolate corresponding distribution value
                f = floor(Int,rx)
                y = dist[f+2]*(rx-f) + dist[f+1]*(1-(rx-f))
                # See if x value is accepted
                ry = rand() * dist_ymax
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

    # Returns a vector of bin centers if given a vector of bin edges
    function cntr(edges)
        centers = (edges[1:end-1]+edges[2:end])/2
        return centers
    end

    # Return percentile of an array along a specified dimension
    function pctile(A,p; dim=0)
        s = size(A)
        if dim == 2
            out = Array{eltype(A)}(undef,s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                out[i] = any(t) ? percentile(A[i,t],p) : NaN
            end
        elseif dim == 1
            out = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                out[i] = any(t) ? percentile(A[t,i],p) : NaN
            end
        else
            t = .~ isnan.(A)
            out = any(t) ? percentile(A[t],p) : eltype(A)(NaN)
        end
        return out
    end

    # Smallest non-NaN value of an array
    function nanminimum(A; dim=0)
        s = size(A)
        if dim == 2
            out = Array{eltype(A)}(undef,s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                out[i] = any(t) ? minimum(A[i,t]) : NaN
            end
        elseif dim == 1
            out = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                out[i] = any(t) ? minimum(A[t,i]) : NaN
            end
        else
            t = .~ isnan.(A)
            out = any(t) ? minimum(A[t]) : NaN
        end
        return out
    end

    # Largest non-NaN value of an array
    function nanmaximum(A; dim=0)
        s = size(A)
        if dim == 2
            out = Array{eltype(A)}(undef,s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                out[i] = any(t) ? maximum(A[i,t]) : NaN
            end
        elseif dim == 1
            out = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                out[i] = any(t) ? maximum(A[t,i]) : NaN
            end
        else
            t = .~ isnan.(A)
            out = any(t) ? maximum(A[t]) : NaN
        end
        return out
    end

    # Range (max-min) of an array, ignoring NaNs
    function nanrange(A; dim=0)
        s = size(A)
        if dim == 2
            out = Array{eltype(A)}(undef,s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                extr = extrema(A[i,t])
                out[i] = extr[2] - extr[1]
            end
        elseif dim == 1
            out = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                extr = extrema(A[t,i])
                out[i] = extr[2] - extr[1]
            end
        else
            t = .~ isnan.(A)
            extr = extrema(A[t])
            out = extr[2] - extr[1]
        end
        return out
    end

    # Mean, ignoring NaNs
    function nanmean(A; dim=0)
        s = size(A)
        if dim == 2
            out = Array{eltype(A)}(undef,s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                out[i] = any(t) ? mean(A[i,t]) : NaN
            end
        elseif dim == 1
            out = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                out[i] = any(t) ? mean(A[t,i]) : NaN
            end
        else
            t = .~ isnan.(A)
            out = any(t) ? mean(A[t]) : NaN
        end
        return out
    end

    # Standard deviation, ignoring NaNs
    function nanstd(A; dim=0)
        s = size(A)
        if dim == 2
            out = Array{eltype(A)}(undef,s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                out[i] = any(t) ? std(A[i,t]) : NaN
            end
        elseif dim == 1
            out = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                out[i] = any(t) ? std(A[t,i]) : NaN
            end
        else
            t = .~ isnan.(A)
            out = any(t) ? std(A[t]) : NaN
        end
        return out
    end

    # Median, ignoring NaNs
    function nanmedian(A; dim=0)
        s = size(A)
        if dim == 2
            out = Array{eltype(A)}(undef, s[1])
            for i=1:s[1]
                t = .~ isnan.(A[i,:])
                out[i] = any(t) ? median(A[i,t]) : NaN
            end
        elseif dim == 1
            out = Array{eltype(A)}(undef, 1, s[2])
            for i=1:s[2]
                t = .~ isnan.(A[:,i])
                out[i] = any(t) ? median(A[t,i]) : NaN
            end
        else
            t = .~ isnan.(A)
            out = any(t) ? median(A[t]) : NaN
        end
        return out
    end

    # Return the index of the closest value in 'target' for each value in 'source'
    # If muliple values are equally close, the first one is used
    function findclosest(source, target)
        # Allocate index and difference arrays
        index = Array{Int64}(undef, size(source))
        diff_type = promote_type(eltype(source), eltype(target))
        diff = Array{diff_type}(undef, length(target))
        # Find closest (numerical) match in target for each value in source
        for i = 1:length(source)
            for j = 1:length(diff)
                diff[j] = abs(target[j] - source[i])
            end
            index[i] = argmin(diff)
        end
        return index
    end

    # Return the index of the closest value of array 'target' below (less than)
    # each value in 'source'. Returns an index of 0 no such values exist
    function findclosestbelow(source, target)
        # Allocate output array
        index = Array{Int64}(undef, size(source))
        diff_type = promote_type(eltype(source), eltype(target))
        diff = Array{diff_type}(undef, length(target))
        t = Array{Bool}(undef,length(target))
        for i = 1:length(source)
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

    # Return the index of the closest value of the vector 'target' above (greater
    # than) each value in 'source'. Returns an index of 0 if no values exist
    function findclosestabove(source, target)
        # Allocate output array
        index = Array{Int64}(undef, size(source))
        diff_type = promote_type(eltype(source), eltype(target))
        diff = Array{diff_type}(undef, length(target))
        t = Array{Bool}(undef,length(target))
        for i = 1:length(source)
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
