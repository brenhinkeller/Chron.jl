## --- Custom objects for holding Chron age data

    # Define a type of struct object to hold data about geochron samples from
    # a region or stratigraphic section of interest
    mutable struct ChronAgeData
        Name::Tuple
        Height::Array{Float64,1}
        Height_sigma::Array{Float64,1}
        Age::Array{Float64,1}
        Age_sigma::Array{Float64,1}
        Age_025CI::Array{Float64,1}
        Age_975CI::Array{Float64,1}
        Age_14C::Array{Float64,1}
        Age_14C_sigma::Array{Float64,1}
        Age_Sidedness::Array{Float64,1}
        Age_Distribution::Array{<:Array,1}
        Params::Array{Float64,2}
        Path::String
        inputSigmaLevel::Int
        Age_Unit::String
        Height_Unit::String
    end

    function NewChronAgeData(nSamples)
        smpl = ChronAgeData(
            ("Sample Names",),
            fill(NaN,nSamples),  # Sample heights
            fill(NaN,nSamples),  # Height_sigma
            fill(NaN,nSamples),  # Sample ages
            fill(NaN,nSamples),  # Sample age uncertainty
            fill(NaN,nSamples),  # Sample age 2.5% CI
            fill(NaN,nSamples),  # Sample age 97.5% CI
            fill(NaN,nSamples),  # Sample 14C ages
            fill(NaN,nSamples),  # Sample 14C uncertainties
            zeros(nSamples), # Sidedness (zeros by default, geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
            Array{Array}(undef,nSamples), # Stationary distribution of eruption age
            fill(NaN,5,nSamples), # Sample age distribution parameters
            "./", # Relative path where we can find .csv data files
            2, # i.e., are the data files 1-sigma or 2-sigma
            "Ma",
            "m",
        )
        return smpl
    end
    export NewChronAgeData

    # Just for backwards compatibility (hopefully not neccessary)
    NewStratAgeData = NewChronAgeData
    export NewStratAgeData

    # A type of object to hold data about hiatuses
    mutable struct HiatusData
        Height::Array{Float64,1}
        Height_sigma::Array{Float64,1}
        Duration::Array{Float64,1}
        Duration_sigma::Array{Float64,1}
        Age_Unit::String
        Height_Unit::String
    end

    function NewHiatusData(nHiatuses)
        hiatus = HiatusData(
            fill(NaN,nHiatuses),  # Height
            fill(NaN,nHiatuses),  # Height_sigma
            fill(NaN,nHiatuses),  # Duration
            fill(NaN,nHiatuses),  # Duration_sigma
            "Ma",
            "m",
        )
        return hiatus
    end
    export NewHiatusData

    # A type of object to specify the configuration of the stratigraphic model
    mutable struct StratAgeModelConfiguration
        resolution::Float64
        burnin::Int
        nsteps::Int
        sieve::Int
        bounding::Float64
    end

    function NewStratAgeModelConfiguration()
        return StratAgeModelConfiguration(NaN, 0, 0, 0, NaN)
    end
    export NewStratAgeModelConfiguration


    struct StratAgeModel
        Height::Array{Float64,1}
        Age::Array{Float64,1}
        Age_sigma::Array{Float64,1}
        Age_Median::Array{Float64,1}
        Age_025CI::Array{Float64,1}
        Age_975CI::Array{Float64,1}
    end
