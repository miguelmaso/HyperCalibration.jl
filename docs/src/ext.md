
HyperCalibration.jl uses Julia's package extension mechanism (`[weakdeps]` /
`[extensions]` in `Project.toml`) to add optional functionality that only
loads when the corresponding package is also loaded. Neither extension
changes the core API described in [ExperimentsData](@ref) /
[ConstitutiveModelling](@ref) / [ObjectiveFunctions](@ref) — they add
convenience methods on top of it.

## CSV loading — `ExperimentsCSVExt`

Loaded automatically when [`CSV`](https://csv.juliadata.org/stable/) is loaded alongside `HyperCalibration`. Adds
a `CSV.read` method that parses a CSV file directly into a vector of
[`ExperimentData`](@ref):

```julia
using HyperCalibration, CSV

experiments = CSV.read("quasi_static_tests.csv", UniaxialQuasiStaticTest)
```

Rows are grouped by an `id` column, and one `ExperimentData` is built per
group. The expected columns depend on the target test type — e.g.
`UniaxialQuasiStaticTest` expects `id`, `stretch`, `stress`; a thermal
variant such as `UniaxialThermalCyclicLoadingTest` additionally expects
`rate` and `temp`. A `thickness` keyword can be passed for tests that need a
specimen thickness (`UniaxialThermoElectricCyclicLoadingTest`); it defaults
to `0.0` otherwise.

In the case the header in the CSV data does not coincide with the expected keywords,
the `header` option can be passed to the reader, e.g.:

```julia
data = IOBuffer("#, λ, σ
1, 1.0, 0.0
1, 1.091769158, 5823.352789
1, 1.19820246, 7726.97553
1, 1.297540208, 9757.927377
")
experiments = CSV.read(data, UniaxialQuasiStaticTest, header=[:id, :stretch, :stress], skipto=2)
```

Supported target types are the [test presets](@ref "Test presets") in
`ExperimentsInterface.jl`. Requesting an unsupported type raises an error
prompting you to add a `build_experiment` method for it.


## Plotting — `HyperCalibrationPlots`

Loaded automatically when [`Plots`](https://docs.juliaplots.org/stable/) is loaded alongside `HyperCalibration`.
Adds `Plots.jl` recipes for visually comparing a constitutive model's
predictions against experimental data:

```julia
using HyperCalibration, HyperFEM, Plots

plot(model, experiments[1], xlabel="Stretch [-]", ylabel="Stress [KPa]", units_scale=1e-3)
```

- `plot(model::PhysicalModel, data::ExperimentData)` draws the model
  prediction as a line and the measured data as a scatter series over the
  data's [`independent_variable`](@ref).
- `plot(model::PhysicalModel, datasets::AbstractVector{<:ExperimentData})`
  overlays one series pair per experiment; pass `label` and/or `seriescolor`
  as vectors (one entry per experiment) to control each series individually.
- `units_scale` rescales both predicted and measured values before plotting
  (e.g. `1e-3` to display stress in kPa instead of Pa).
