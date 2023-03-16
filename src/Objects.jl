## --- Custom objects for holding Chron age data

    # Define a type of struct object to hold data about geochron samples from
    # a region or stratigraphic section of interest
    mutable struct ChronAgeData{N}
        Name::NTuple{N, String}
        Height::Vector{Float64}
        Height_sigma::Vector{Float64}
        Age::Vector{Float64}
        Age_sigma::Vector{Float64}
        Age_025CI::Vector{Float64}
        Age_975CI::Vector{Float64}
        Age_14C::Vector{Float64}
        Age_14C_sigma::Vector{Float64}
        Age_Sidedness::Vector{Float64}
        Age_DistType::Vector{Float64}
        Age_Distribution::Vector{Vector{Float64}}
        Chronometer::NTuple{N, Symbol}
        Params::Matrix{Float64}
        Path::String
        inputSigmaLevel::Int
        Age_Unit::String
        Height_Unit::String
    end

    function ChronAgeData(nSamples::Integer)
        smpl = ChronAgeData{nSamples}(
            ntuple(i->"Sample name", nSamples),
            fill(NaN,nSamples),  # Sample heights
            fill(NaN,nSamples),  # Height_sigma
            fill(NaN,nSamples),  # Sample ages
            fill(NaN,nSamples),  # Sample age uncertainty
            fill(NaN,nSamples),  # Sample age 2.5% CI
            fill(NaN,nSamples),  # Sample age 97.5% CI
            fill(NaN,nSamples),  # Sample 14C ages
            fill(NaN,nSamples),  # Sample 14C uncertainties
            zeros(nSamples), # Sidedness (zeros by default, geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
            zeros(nSamples), # DistType (for Distribution-fitting only: 0=distribution to be fit, 1=single Gaussian)
            Vector{Vector{Float64}}(undef,nSamples), # Stationary distribution of eruption age
            ntuple(i->:Chronometer, nSamples), # Age Types (e.g., :UPb or :ArAr)
            fill(NaN,5,nSamples), # Sample age distribution parameters
            "./", # Relative path where we can find .csv data files
            2, # i.e., are the data files 1-sigma or 2-sigma
            "Ma",
            "m",
        )
        return smpl
    end
    # For backwards compatibility
    const NewChronAgeData = ChronAgeData
    export NewChronAgeData
    const NewStratAgeData = ChronAgeData
    export NewStratAgeData
    const StratAgeData = ChronAgeData
    export StratAgeData

    # One-sigma systematic uncertainty
    mutable struct SystematicUncertainty
        UPb::Float64
        ArAr::Float64
    end
    SystematicUncertainty() = SystematicUncertainty(NaN, NaN)
    export SystematicUncertainty

    # A type of object to hold data about hiatuses
    mutable struct HiatusData
        Height::Vector{Float64}
        Height_sigma::Vector{Float64}
        Duration::Vector{Float64}
        Duration_sigma::Vector{Float64}
        Age_Unit::String
        Height_Unit::String
    end

    function HiatusData(nHiatuses::Integer)
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
    const NewHiatusData = HiatusData
    export NewHiatusData

    # A type of object to specify the configuration of the stratigraphic model
    mutable struct StratAgeModelConfiguration
        resolution::Float64
        burnin::Int
        nsteps::Int
        sieve::Int
        bounding::Float64
    end

    function StratAgeModelConfiguration()
        return StratAgeModelConfiguration(NaN, 0, 0, 0, NaN)
    end
    const NewStratAgeModelConfiguration = StratAgeModelConfiguration
    export NewStratAgeModelConfiguration


    struct StratAgeModel
        Height::Vector{Float64}
        Age::Vector{Float64}
        Age_sigma::Vector{Float64}
        Age_Median::Vector{Float64}
        Age_025CI::Vector{Float64}
        Age_975CI::Vector{Float64}
    end
