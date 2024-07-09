## --- Types to hold age constraints for age-depth modelling, before or after eruption age modelling

    # Define a type of struct object to hold data about geochron samples from
    # a region or stratigraphic section of interest
    @kwdef mutable struct ChronAgeData{N}
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
        Sidedness_Method::Symbol
    end

    function ChronAgeData(nsamples::Integer)
        ChronAgeData{nsamples}(
            ntuple(i->"Sample name", nsamples),
            collect(1.0:nsamples),  # Sample Height
            zeros(nsamples),     # Sample Height_sigma
            fill(NaN,nsamples),  # Sample ages
            fill(NaN,nsamples),  # Sample age uncertainty
            fill(NaN,nsamples),  # Sample age 2.5% CI
            fill(NaN,nsamples),  # Sample age 97.5% CI
            fill(NaN,nsamples),  # Sample 14C ages
            fill(NaN,nsamples),  # Sample 14C uncertainties
            zeros(nsamples), # Sidedness (zeros by default, geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
            zeros(nsamples), # DistType (for Distribution-fitting only: 0=distribution to be fit, 1=single Gaussian)
            Vector{Vector{Float64}}(undef,nsamples), # Stationary distribution of eruption age
            ntuple(i->:Chronometer, nsamples), # Age Types (e.g., :UPb or :ArAr)
            fill(NaN,5,nsamples), # Sample age distribution parameters
            "./", # Relative path where we can find .csv data files
            2, # i.e., are the data files 1-sigma or 2-sigma
            "Ma",
            "m",
            :cdf,
        )
    end

    @kwdef mutable struct GeneralAgeData{N}
        Name::NTuple{N, String}
        Height::Vector{Float64}
        Height_sigma::Vector{Float64}
        Age::Vector{<:Union{<:Distribution{Univariate, Continuous}}}
        Age_Sidedness::Vector{Float64}
        Chronometer::NTuple{N, Symbol}
        Age_Unit::String
        Height_Unit::String
        Sidedness_Method::Symbol
    end

    function GeneralAgeData(nsamples::Integer)
        GeneralAgeData{nsamples}(
            ntuple(i->"Sample name", nsamples),
            collect(1.0:nsamples),  # Sample Height
            zeros(nsamples),     # Sample Height_sigma
            fill(Uniform(0,4567), nsamples),  # Sample ages
            zeros(nsamples), # Sidedness (zeros by default, geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
            ntuple(i->:Chronometer, nsamples), # Age Types (e.g., :UPb or :ArAr)
            "Ma",
            "m",
            :cdf,
        )
    end

## --- `HiatusData` type to hold data about hiatuses

    mutable struct HiatusData
        Height::Vector{Float64}
        Height_sigma::Vector{Float64}
        Duration::Vector{Float64}
        Duration_sigma::Vector{Float64}
        Age_Unit::String
        Height_Unit::String
    end

    function HiatusData(nhiatuses::Integer)
        HiatusData(
            fill(NaN, nhiatuses),  # Height
            fill(NaN, nhiatuses),  # Height_sigma
            fill(NaN, nhiatuses),  # Duration
            fill(NaN, nhiatuses),  # Duration_sigma
            "Ma",
            "m",
        )
    end

## --- Types type to specify the configuration of and hold results from the stratigraphic age-depth model

    mutable struct StratAgeModelConfiguration
        resolution::Float64
        burnin::Int
        nsteps::Int
        sieve::Int
        bounding::Float64
    end
    StratAgeModelConfiguration() = StratAgeModelConfiguration(NaN, 0, 0, 0, NaN)

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

    function StratAgeModel(heights, agedist)
        agedistₜ = copy(agedist)
        StratAgeModel(
            heights, # Model heights
            nanmean(agedistₜ,dim=2), # Mean age
            nanstd(agedistₜ,dim=2), # Standard deviation
            nanmedian!(agedistₜ,dim=2), # Median age
            nanpctile!(agedistₜ,2.5,dim=2), # 2.5th percentile
            nanpctile!(agedistₜ,97.5,dim=2) # 97.5th percentile
        )
    end

    Base.NamedTuple(mdl::StratAgeModel) = NamedTuple{fieldnames(typeof(mdl))}(ntuple(i->getfield(mdl, i), nfields(mdl)))

## --- Types to allow dispatch to different age sidedness methods
    abstract type Sidedness end
    struct CDFSidedness{T<:Real} <: Sidedness
        s::Vector{T}
    end
    struct FastSidedness{T<:Real} <: Sidedness
        s::Vector{T}
    end
    # Forward all array-like methods to the wrapped vector
    Base.getindex(A::Sidedness, args...) = getindex(A.s, args...)
    Base.setindex!(A::Sidedness, args...) = setindex!(A.s, args...)
    Base.eachindex(A::Sidedness) = eachindex(A.s)
    Base.size(A::Sidedness, args...) = size(A.s, args...)
    Base.length(A::Sidedness) = length(A.s)
    Base.ndims(A::Sidedness) = ndims(A.s)
    Base.ndims(::Type{<:Sidedness}) = 1
    Base.copyto!(A::Sidedness, args...) = copyto!(A.s, args...)

## --- Types to allow specification of systematic uncertainties

    mutable struct SystematicUncertainty
        UPb::Float64
        ArAr::Float64
    end
    SystematicUncertainty() = SystematicUncertainty(NaN, NaN)