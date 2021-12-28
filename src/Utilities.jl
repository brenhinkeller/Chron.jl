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
        xs = (x - p[2])/abs2(p[3]) # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        shp = abs2(p[4])
        skw = abs2(p[5])
        return exp(p[1] + shp*skw*xs*v - shp/skw*xs*(1-v))
    end
    function bilinear_exponential(x::AbstractVector, p::AbstractVector)
        result = Array{float(eltype(x))}(undef,size(x))
        @turbo for i ∈ eachindex(x)
            xs = (x[i] - p[2])/abs2(p[3]) # X scaled by mean and variance
            v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
            shp = abs2(p[4])
            skw = abs2(p[5])
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
    function bilinear_exponential_ll(x::Number, p::AbstractVector)
        xs = (x - p[2])/abs2(p[3]) # X scaled by mean and variance
        v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
        shp = abs2(p[4])
        skw = abs2(p[5])
        return p[1] + shp*skw*xs*v - shp/skw*xs*(1-v)
    end
    function bilinear_exponential_ll(x::AbstractVector, p::AbstractMatrix)
        ll = 0.0
        @turbo for i ∈ eachindex(x)
            xs = (x[i]-p[2,i])/abs2(p[3,i]) # X scaled by mean and variance
            v = 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
            shp = abs2(p[4,i])
            skw = abs2(p[5,i])
            ll += p[1,i] + shp*skw*xs*v - shp/skw*xs*(1-v)
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

## --- End of File
