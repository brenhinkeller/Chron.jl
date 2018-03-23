## --- Custom objects

    # Define a type of struct object to hold data about geochron samples from
    # a stratigraphic section of interest
    struct StratAgeData
        Name::Tuple
        Height::Array{Float64}
        Height_Sigma::Array{Float64}
        Age::Array{Float64}
        Age_Sigma::Array{Float64}
        p::Array{Float64}
    end

    # Define another type of struct object to hold data about age models
    # from our stratigraphic section of interest
    struct StratAgeModel
        Height::Array{Float64}
        Age::Array{Float64}
        Age_Sigma::Array{Float64}
        Age_Median::Array{Float64}
        Age_025CI::Array{Float64}
        Age_975CI::Array{Float64}
    end
## ---
