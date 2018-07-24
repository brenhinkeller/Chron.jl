# Chron.jl

A two-part framework for (1) estimating eruption/deposition age distributions from complex mineral age spectra and (2) subsequently building a stratigraphic age model based on those distributions. Each step relies on a Markov-Chain Monte Carlo model.

The first (distribution) MCMC model uses information about the possible shape of the true mineral crystallization (or closure) age distribution (e.g., no crystallization possible after eruption or deposition). In this first model, the true eruption or deposition age is a parameter of this scaled crystallization distribution. The stationary distribution of this first MCMC model then gives an estimate of the eruption/deposition age.

The second (stratigraphic) MCMC  model uses the estimated (posterior) eruption/deposition age distributions along with the constraint of stratigraphic superposition to produce an age-depth model

## Installation

In the Julia REPL
```Julia
julia> Pkg.clone("https://github.com/brenhinkeller/Chron.jl")
```

## Usage

Run [examples.jl](examples/examples.jl) to see how the code works.


[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/brenhinkeller/Chron.jl/master?filepath=demo.ipynb)

