
## ExperimentsData

Data structures for defining an experimental test: the measurement, the
loading/environmental **protocol**, the ambient **conditions**, and the
specimen **geometry**. These four pieces are combined into an
[`ExperimentData`](@ref) instance, which is the unit consumed by
[`ConstitutiveModelling`](@ref) and [`ObjectiveFunctions`](@ref).

```@autodocs
Modules = [HyperCalibration]
Pages = ["ExperimentsData.jl"]
Order = [:type, :function]
```

### Test presets

Ready-made [`ExperimentData`](@ref) specializations (e.g.
[`UniaxialQuasiStaticTest`](@ref), [`DielectricSpectroscopyTest`](@ref)) and
their convenience constructors, for the common combinations of
measurement/protocol/condition/geometry.

```@autodocs
Modules = [HyperCalibration]
Pages = ["ExperimentsInterface.jl"]
Order = [:type, :function]
```


## ConstitutiveModelling

Functions that evaluate a `HyperFEM` constitutive model against the protocol,
condition, and geometry of an [`ExperimentData`](@ref), predicting the
dependent variable (stress, specific heat capacity, permittivity...) that
[`ObjectiveFunctions`](@ref) compares against the measured data.

```@autodocs
Modules = [HyperCalibration]
Pages = ["ConstitutiveModelling.jl"]
Order = [:type, :function]
```

### Reduced kinematics

Internal kinematic maps (deformation gradient construction, thermal Jacobian
solve) used by [`ConstitutiveModelling`](@ref) to reduce a 3D constitutive model to
the 1D/2D loading states of a given [`Kinematics`](@ref) (uniaxial, biaxial).

```@autodocs
Modules = [HyperCalibration]
Pages = ["ReducedKinematics.jl"]
Order = [:type, :function]
```


## ObjectiveFunctions

Predefined loss functions that evaluate the mismatch between observed data
and the data predicted by [`ConstitutiveModelling`](@ref), for a given
constitutive model, a set of parameters, and one or more
[`ExperimentData`](@ref). These are the functions you minimize (e.g. with
`Optim.jl`) during calibration.

```@autodocs
Modules = [HyperCalibration]
Pages = ["ObjectiveFunctions.jl"]
Order = [:type, :function]
```
