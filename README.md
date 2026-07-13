# HyperCalibration.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://miguelmaso.github.io/HyperCalibration.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://miguelmaso.github.io/HyperCalibration.jl/dev/)
[![Build Status](https://github.com/miguelmaso/HyperCalibration.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/miguelmaso/HyperCalibration.jl/actions/workflows/ci.yml?branch=main)
[![Aqua QA](https://juliatesting.github.io/Aqua.jl/dev/assets/badge.svg)](https://juliatesting.github.io/Aqua.jl/stable/)

Tools for the calibration of hyperelastic constitutive models, built on top of [HyperFEM.jl](https://github.com/MultiSimOLab/HyperFEM.jl).

The main components of the package are:
- **ExperimentData**: a data structure holding the definition of an experimental test.
- **constitutive modelling**: a set of functions that evaluate a *constitutive model* from HyperFEM with the data from an experiment and predicts the dependent variable.
- **objective functions**: predefined *loss* function that evaluates the difference from observed data and predicted data with a given *constitutive model*, a set of parameters, and a set of experiments.

## How to use

The package is a registered package, and can be installed with `Pkg.add`.
```julia
julia> using Pkg; Pkg.add("Optim")
```
or through the pkg REPL mode by typing
```
] add HyperCalibration
```

## Example

```julia
using HyperCalibration, HyperFEM
using Optim
using CSV
using Plots

experiments = CSV.read(joinpath(@__DIR__, "quasi_static_tests.csv"), UniaxialQuasiStaticTest, decimal='.')

build_model(μ, N) = EightChain(μ=μ, N=N)
pn = [  "μ",   "N"]  # Parameter names
p0 = [  1e4,   1.0]  # Initial seed

f(p) = loss(build_model, p, experiments)

result = optimize(f, p0, NelderMead())
model = build_model(result.minimizer...)
plot(model, experiments[1], xlabel="Stretch [-]", ylabel="Stress [KPa]", units_scale=1e-3)
```

## Project funded by

- Grants PID2022-141957OA-C22/PID2022-141957OB-C22 funded by MCIN/AEI/ 10.13039/501100011033
