# Chron.jl

[![DOI](https://github.com/brenhinkeller/Chron.jl/blob/master/osf.io:TQX3F.svg)](https://doi.org/10.17605/OSF.IO/TQX3F)

A two-part framework for (1) estimating eruption/deposition age distributions from complex mineral age spectra and (2) subsequently building a stratigraphic age model based on those distributions. Each step relies on a Markov-Chain Monte Carlo model.

The first (distribution) MCMC model is based on the work of [Keller, Schoene, and Samperton (2018)]( https://doi.org/10.7185/geochemlet.1826) and uses information about the possible shape of the true mineral crystallization (or closure) age distribution (e.g., no crystallization possible after eruption or deposition). In this first model, the true eruption or deposition age is a parameter of this scaled crystallization distribution. The stationary distribution of this first MCMC model then gives an estimate of the eruption/deposition age.

The second (stratigraphic) MCMC  model uses the estimated (posterior) eruption/deposition age distributions along with the constraint of stratigraphic superposition to produce an age-depth model

## Installation

Chron.jl is written in the [Julia](https://julialang.org/) language.

In the Julia package manager (type `]` in the REPL)
```Julia
(v1.0) pkg> add "https://github.com/brenhinkeller/Chron.jl"
```
or for previous versions of Julia, in the REPL
```Julia
julia> Pkg.clone("https://github.com/brenhinkeller/Chron.jl")
```

## Usage
### Online / notebook usage
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/brenhinkeller/Chron.jl/master?filepath=examples%2Fdemo.ipynb)

For a quick test (without having to install anything), try the [interactive online Jupyter notebook](https://mybinder.org/v2/gh/brenhinkeller/Chron.jl/master?filepath=examples%2Fdemo.ipynb) (note: it'll take a few minutes for the notebook to launch).

This runs [examples/demo.ipynb](examples/demo.ipynb) on a [JupyterHub](https://github.com/jupyterhub/jupyterhub) server hosted by the [Binder](https://mybinder.org) project. If you make changes to the interactive online notebook, you can save them with `File` > `Download as` > `Notebook (.ipynb)` To run a downloaded notebook locally, use [IJulia](https://github.com/JuliaLang/IJulia.jl)

```Julia
julia> using IJulia
julia> notebook()
```

### Standard usage

After installing [Julia](https://julialang.org/downloads/) with or without [Juno](http://junolab.org/), and Chron.jl (above), run [examples/demo.jl](examples/demo.jl) to see how the code works.
