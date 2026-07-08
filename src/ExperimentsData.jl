
abstract type AbstractMeasurement end
abstract type AbstractProtocol end
abstract type AbstractCondition end
abstract type AbstractGeometry end

# --- Kinematics ---

abstract type Kinematics end

"""
Uniaxial kinematics allow to evaluate the stress response at one Gauss point.
"""
struct Uniaxial <: Kinematics end

"""
Biaxial kinematics allow to evaluate the stress response at one Gauss point.
"""
struct Biaxial <: Kinematics end

# --- Measurements ---

"""
The nominal stress recoded during a tensile test.
"""
struct TensileMeasurement <: AbstractMeasurement
  σ::Vector{Float64}
end

"""
The specific heat capacity recoded during a thermal test.
"""
struct ThermalMeasurement <: AbstractMeasurement
  cv::Vector{Float64}
end

"""
The dielectric parmittivity recoded during a dielectric test.
"""
struct DielectricMeasurement <: AbstractMeasurement
  ε::Vector{Float64}
end

# --- Protocols ---

abstract type MechanicalProtocol{K<:Kinematics} <: AbstractProtocol end

"""
A sequential protocol defines a series of sub-protocols to be executed in order.
E.g. a relaxation test is defined by a fast loading stage (CyclicLoadingProtocol) followed by a long holding stage (CreepProtocol).
"""
struct SequentialProtocol{K<:Kinematics} <: MechanicalProtocol{K}
  stages::Vector{MechanicalProtocol{K}}
end

"""
A sequence of stretches.
"""
struct QuasiStaticProtocol{K<:Kinematics} <: MechanicalProtocol{K}
  λ::Vector{Float64}
end

"""
A sequence of stretches at a constant rate and time step.
"""
struct CyclicLoadingProtocol{K<:Kinematics} <: MechanicalProtocol{K}
  λ::Vector{Float64}
  v::Float64
  Δt::Float64
end

"""
A sequence of times at a constant stretch.
"""
struct CreepProtocol{K<:Kinematics} <: MechanicalProtocol{K}
  t::Vector{Float64}
  λ::Float64
end

"""
A sequence of temperatures at constant rate.
"""
struct TemperatureSweepProtocol <: AbstractProtocol
  θ::Vector{Float64}
  v::Float64
end

"""
A sequence of frequencies.
"""
struct FrequencySweepProtocol <: AbstractProtocol
  f::Vector{Float64}
end

"""
Return the sequence of stretches for a given protocol.
"""
stretches(::AbstractProtocol) = throw(ArgumentError("stretches not defined for this protocol"))

stretches(p::SequentialProtocol) = vcat(map(stretches, p.stages)...)

stretches(p::QuasiStaticProtocol) = p.λ

stretches(p::CyclicLoadingProtocol) = p.λ

"""
Return the constant time step for a given protocol.
"""
time_step(::AbstractProtocol) = throw(ArgumentError("time_step not defined for this protocol"))

time_step(p::CyclicLoadingProtocol) = p.Δt

time_step(p::CreepProtocol) = diff(p.t)[1]

"""
Return the sequence of temperatures for a given protocol.
"""
temperatures(::AbstractProtocol) = throw(ArgumentError("temperatures not defined for this protocol"))

temperatures(p::TemperatureSweepProtocol) = p.θ

# --- Conditions ---

"""
A standard condition with no additional parameters.
"""
struct StandardCondition <: AbstractCondition
end

"""
An experiment at a constant temperature.
"""
struct IsothermalCondition <: AbstractCondition
  θ::Float64
end

"""
An experiment at a constant voltage.
"""
struct ElectricalCondition <: AbstractCondition
  V::Float64
end

"""
An experiment at a constant temperature and voltage.
"""
struct ThermoElectricalCondition <: AbstractCondition
  θ::Float64
  V::Float64
end

"""
Return the room temperature for a given condition.
"""
temperature(::AbstractCondition) = 293.15

temperature(c::IsothermalCondition) = c.θ

temperature(c::ThermoElectricalCondition) = c.θ

"""
Return the applied voltage for a given condition.
"""
voltage(::AbstractCondition) = 0.0

voltage(c::ElectricalCondition) = c.V

voltage(c::ThermoElectricalCondition) = c.V

# --- Geometries ---

"""
A specimen geometry with a reference thickness and cross-sectional area.
"""
struct PlateGeometry <: AbstractGeometry
  t0::Float64
  A0::Float64
end

PlateGeometry() = PlateGeometry(0.0, 0.0)

PlateGeometry(t0) = PlateGeometry(t0, 0.0)

"""
Return the thickness for a given geometry.
"""
thickness(::AbstractGeometry) = throw(ArgumentError("thickness not defined for this geometry"))

thickness(g::PlateGeometry) = g.t0

"""
Return the electric field for a given condition and geometry.
"""
electric_field(::AbstractCondition, ::AbstractGeometry) = voltage(c) / thickness(g)

# --- Experiments ---

"""
    ExperimentData{M, P, C, G}

An experiment holds a measurement M that is a function of a Protocol P -defining the independent variables-
and conditionc C -environtmental static variables- and geometry G -for the reference specimen configuration-.
Additionally, a weight can be assigned to the experiment to define its importance in the calibration process.
"""
mutable struct ExperimentData{M<:AbstractMeasurement, P<:AbstractProtocol, C<:AbstractCondition, G<:AbstractGeometry}
  const measurement::M
  const protocol::P
  const condition::C
  const geometry::G
  const id::Int
  weight::Float64
end

"""
Experiment data with mechanical measurement (stress) and mechanical protocol (stretch).
"""
MechanicalTest{C<:AbstractCondition, G<:AbstractGeometry} = ExperimentData{TensileMeasurement, MechanicalProtocol, C, G}

"""
An experiment data with thermal measurement (specific heat capacity) and thermal protocol (temperature).
"""
ThermalTest{C<:AbstractCondition, G<:AbstractGeometry} = ExperimentData{ThermalMeasurement, TemperatureSweepProtocol, C, G}

"""
An experiment data with dielectric measurement (permittivity) and thermal protocol (temperature).
"""
ThermoDielectricTest{C<:AbstractCondition, G<:AbstractGeometry} = ExperimentData{DielectricMeasurement, TemperatureSweepProtocol, C, G}

"""
Experiment data with uniaxial quasi-static stretch-stress and standard conditions.
"""
const UniaxialQuasiStaticTest = ExperimentData{TensileMeasurement, QuasiStaticProtocol{Uniaxial}, StandardCondition, PlateGeometry}

"""
Experiment data with constant rate uniaxial stretch-stress and standard conditions.
"""
const UniaxialCyclicLoadingTest = ExperimentData{TensileMeasurement, CyclicLoadingProtocol{Uniaxial}, StandardCondition, PlateGeometry}

"""
Experiment data with uniaxial relaxation stressess and standard conditions.
"""
const UniaxialRelaxationTest = ExperimentData{TensileMeasurement, SequentialProtocol{Uniaxial}, StandardCondition, PlateGeometry}

"""
Experiment data with uniaxial quasi-static stretch-stress and thermal conditions.
"""
const UniaxialThermalQuasiStaticTest = ExperimentData{TensileMeasurement, QuasiStaticProtocol{Uniaxial}, IsothermalCondition, PlateGeometry}

"""
Experiment data with constant rate uniaxial stretch-stress and thermal conditions.
"""
const UniaxialThermalCyclicLoadingTest = ExperimentData{TensileMeasurement, CyclicLoadingProtocol{Uniaxial}, IsothermalCondition, PlateGeometry}

"""
Experiment data with uniaxial relaxation stressess and thermal conditions.
"""
const UniaxialThermalRelaxationTest = ExperimentData{TensileMeasurement, SequentialProtocol{Uniaxial}, IsothermalCondition, PlateGeometry}

"""
Experiment data with constant rate uniaxial stretch-stress and coupled thermo-electrical conditions.
"""
const UniaxialThermoElectricCyclicLoadingTest = ExperimentData{TensileMeasurement, CyclicLoadingProtocol{Uniaxial}, ThermoElectricalCondition, PlateGeometry}

"""
Experiment data with temperature-specific heat capacity (DSC).
"""
const DifferentialScanningCalorimetryTest = ExperimentData{ThermalMeasurement, TemperatureSweepProtocol, StandardCondition, PlateGeometry}

"""
Experiment data with frequency-dielectric permittivity (BDS).
"""
const DielectricSpectroscopyTest = ExperimentData{DielectricMeasurement, FrequencySweepProtocol, StandardCondition, PlateGeometry}
