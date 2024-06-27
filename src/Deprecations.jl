    # For backwards compatibility with old scripts and publications
    function NewChronAgeData(nSamples::Integer)
        @warn "`NewChronAgeData` is deprecated, use `ChronAgeData` instead"
        ChronAgeData(nSamples)
    end
    export NewChronAgeData
    function StratAgeData(nSamples::Integer)
        @warn "`StratAgeData` is deprecated, use `ChronAgeData` instead"
        ChronAgeData(nSamples)
    end
    export StratAgeData
    function NewStratAgeData(nSamples::Integer)
        @warn "`NewStratAgeData` is deprecated, use `ChronAgeData` instead"
        ChronAgeData(nSamples)
    end
    export NewStratAgeData


    function NewHiatusData(nHiatuses::Integer)
        @warn "`NewHiatusData` is deprecated, use `HiatusData` instead"
        HiatusData(nHiatuses)
    end