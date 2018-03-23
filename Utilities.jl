## --- Various utility functions

    # Returns a vector of bin centers if given a vector of bin edges
    function cntr(edges)
        centers = (edges[1:end-1]+edges[2:end])/2;
        return centers;
    end

## --- Plotting

function plotRankOrderErrorbar(data,uncert; seriestype=:scatter,ylabel="",label="",xticks=[],xlabel="")
    sI = sortperm(data);
    h = plot(1:length(sI),data[sI],yerror=uncert[sI],seriestype=seriestype,
             label=label,ylabel=ylabel,xlabel=xlabel,xticks=xticks)
    return h;
end

## --- Weighted means

    # Calculate a weigted mean, including MSWD, with MSWD correction to uncertainty.
    function gwmean(x, sigma)
        n = length(x);
        s1 = 0.0; s2 = 0.0; s3 = 0.0;

        if n==1
            wx = x[1];
            mswd = 0;
            wsigma = sigma[1];
        else
            for i=1:n
                s1 += x[i] / (sigma[i]*sigma[i]);
                s2 += 1 / (sigma[i]*sigma[i]);
            end
            wx = s1/s2;

            for i=1:n
                s3 += (x[i] - wx) * (x[i] - wx) / (sigma[i]*sigma[i]);
            end
            mswd = s3 / (n-1);
            wsigma = sqrt(mswd/s2);
        end
        return wx, wsigma, mswd
    end

    # Calculate a weigted mean, including MSWD, but without MSWD correction to uncertainty
    function awmean(x, sigma)
        n = length(x);
        s1 = 0.0; s2 = 0.0; s3 = 0.0;

        if n==1
            wx = x[1];
            mswd = 0;
            wsigma = sigma[1];
        else
            for i=1:n
                s1 += x[i] / (sigma[i]*sigma[i]);
                s2 += 1 / (sigma[i]*sigma[i]);
            end
            wx = s1/s2;

            for i=1:n
                s3 += (x[i] - wx) * (x[i] - wx) / (sigma[i]*sigma[i]);
            end
            mswd = s3 / (n-1);
            wsigma = sqrt(1.0/s2);
        end
        return wx, wsigma, mswd
    end


## --- Drawing from distributions

    # Draw random numbers from a distribution specified by a vector of points
    # defining the PDF curve
    function drawFromDistribution(dist::Array{Float64}, n::Int)
        # Draw n random numbers from the distribution 'dist'
        x = Array{Float64}(n);
        dist_ymax = maximum(dist);
        dist_xmax = length(dist)-1.0;

        for i=1:n;
            while true
                # Pick random x value
                rx = rand() * dist_xmax;
                # Interpolate corresponding distribution value
                f = floor(Int,rx)
                y = dist[f+2]*(rx-f) + dist[f+1]*(1-(rx-f))
                # See if x value is accepted
                ry = rand() * dist_ymax;
                if (y > ry)
                    x[i] = rx / dist_xmax;
                    break
                end
            end
        end
        return x
    end

    # Fill an existing variable with random numbers from a distribution specified
    # by a vector of points defining the PDF curve
    function fillFromDistribution(dist::Array{Float64}, x::Array{Float64})
        # Fill the array x with random numbers from the distribution 'dist'
        dist_ymax = maximum(dist);
        dist_xmax = length(dist)-1.0;

        for i=1:length(x);
            while true
                # Pick random x value
                rx = rand() * dist_xmax;
                # Interpolate corresponding distribution value
                f = floor(Int,rx)
                y = dist[f+2]*(rx-f) + dist[f+1]*(1-(rx-f))
                # See if x value is accepted
                ry = rand() * dist_ymax;
                if (y > ry)
                    x[i] = rx / dist_xmax;
                    break
                end
            end
        end
    end


## --- End of File
