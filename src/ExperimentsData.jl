
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
  TensileMeasurement(σ::AbstractVector) = new(collect(Float64, σ))
end

"""
The specific heat capacity recoded during a thermal test.
"""
struct ThermalMeasurement <: AbstractMeasurement
  cv::Vector{Float64}
  ThermalMeasurement(cv::AbstractVector) = new(collect(Float64, cv))
end

"""
The dielectric parmittivity recoded during a dielectric test.
"""
struct DielectricMeasurement <: AbstractMeasurement
  ε::Vector{Float64}
  DielectricMeasurement(ε::AbstractVector) = new(collect(Float64, ε))
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
  QuasiStaticProtocol{K}(λ::AbstractVector) where {K<:Kinematics} = new{K}(collect(Float64, λ))
end

"""
A sequence of stretches at a constant rate and time step.
"""
struct CyclicLoadingProtocol{K<:Kinematics} <: MechanicalProtocol{K}
  λ::Vector{Float64}
  v::Float64
  Δt::Float64
  CyclicLoadingProtocol{K}(λ::AbstractVector, v::Real, Δt::Real) where {K<:Kinematics} = new{K}(collect(Float64, λ), Float64(v), Float64(Δt))
end

"""
A sequence of times at a constant stretch.
"""
struct CreepProtocol{K<:Kinematics} <: MechanicalProtocol{K}
  t::Vector{Float64}
  λ::Float64
  CreepProtocol{K}(t::AbstractVector, λ::Real) where {K<:Kinematics} = new{K}(collect(Float64, t), Float64(λ))
end

"""
A sequence of temperatures at constant rate.
"""
struct TemperatureSweepProtocol <: AbstractProtocol
  θ::Vector{Float64}
  v::Float64
  TemperatureSweepProtocol(θ::AbstractVector, v::Real) = new(collect(Float64, θ), Float64(v))
end

"""
A sequence of frequencies.
"""
struct FrequencySweepProtocol <: AbstractProtocol
  f::Vector{Float64}
  FrequencySweepProtocol(f::AbstractVector) = new(collect(Float64, f))
end

"""
Return the sequence of stretches for a given protocol.
"""
stretches(::AbstractProtocol) = throw(ArgumentError("stretches not defined for this protocol."))

stretches(p::SequentialProtocol) = vcat(map(stretches, p.stages)...)

stretches(p::QuasiStaticProtocol) = p.λ

stretches(p::CyclicLoadingProtocol) = p.λ

"""
Return the sequence of temperatures for a given protocol.
"""
temperatures(::AbstractProtocol) = throw(ArgumentError("temperatures not defined for this protocol."))

temperatures(p::TemperatureSweepProtocol) = p.θ

"""
Return the sequence of frequencies for a given protocol.
"""
frequencies(::AbstractProtocol) = throw(ArgumentError("frequencies not defined for this protocol."))

frequencies(p::FrequencySweepProtocol) = p.f

"""
Return the sequence of the independent variable for a given protocol.
"""
independent_variable(::AbstractProtocol) = throw(ArgumentError("independent_variable not implemented for this protocol."))

independent_variable(p::MechanicalProtocol) = stretches(p)

independent_variable(p::TemperatureSweepProtocol) = temperatures(p)

independent_variable(p::FrequencySweepProtocol) = frequencies(p)

"""
Return the constant time step for a given protocol.
"""
time_step(::AbstractProtocol) = throw(ArgumentError("time_step not defined for this protocol."))

time_step(p::CyclicLoadingProtocol) = p.Δt

time_step(p::CreepProtocol) = diff(p.t)[1]

"""
Return the rate of independent variable increment for a given protocol.
"""
rate(::AbstractProtocol) = throw(ArgumentError("rate not defined for this protocol."))

rate(p::CyclicLoadingProtocol) = p.v

rate(p::TemperatureSweepProtocol) = p.v

"""
Return the maximum stretch for a givben protocol.
"""
max_stretch(::AbstractProtocol) = throw(ArgumentError("max_stretch not defined for this protocol."))

max_stretch(p::MechanicalProtocol) = maximum(stretches(p))


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
  IsothermalCondition(θ::Real) = new(Float64(θ))
end

"""
An experiment at a constant voltage.
"""
struct ElectricalCondition <: AbstractCondition
  V::Float64
  ElectricalCondition(V::Real) = new(Float64(V))
end

"""
An experiment at a constant temperature and voltage.
"""
struct ThermoElectricalCondition <: AbstractCondition
  θ::Float64
  V::Float64
  ThermoElectricalCondition(θ::Real, V::Real) = new(Float64(θ), Float64(V))
end

"""
Return the room temperature for a given condition.
"""
temperature(::AbstractCondition) = throw(ArgumentError("temperature not defined for this condition."))

temperature(c::IsothermalCondition) = c.θ

temperature(c::ThermoElectricalCondition) = c.θ

"""
Return the applied voltage for a given condition.
"""
voltage(::AbstractCondition) = throw(ArgumentError("voltage not defined for this condition."))

voltage(c::ElectricalCondition) = c.V

voltage(c::ThermoElectricalCondition) = c.V


# --- Geometries ---

"""
A specimen geometry with a reference thickness and cross-sectional area.
"""
struct PlateGeometry <: AbstractGeometry
  t0::Float64
  A0::Float64
  PlateGeometry(t0::Real=0.0, A0::Real=0.0) = new(Float64(t0), Float64(A0))
end


"""
Return the thickness for a given geometry.
"""
thickness(::AbstractGeometry) = throw(ArgumentError("thickness not defined for this geometry"))

thickness(g::PlateGeometry) = g.t0

"""
Return the electric field for a given condition and geometry.
"""
electric_field(c::AbstractCondition, g::AbstractGeometry) = voltage(c) / thickness(g)


# --- Experiments ---

"""
    ExperimentData{M, P, C, G}

An experiment holds a measurement M that is a function of a Protocol P -defining the independent variables-
and conditionc C -environtmental static variables- and geometry G -for the reference specimen configuration-.
Additionally, a weight can be assigned to the experiment to define its importance in the calibration process.
"""
mutable struct ExperimentData{M<:AbstractMeasurement, P<:AbstractProtocol, C<:AbstractCondition, G<:AbstractGeometry}
  const id::Int
  const measurement::M
  const protocol::P
  const condition::C
  const geometry::G
  weight::Float64
end

"""
Super type for ExperimentData with mechanical measurement (stress) and mechanical protocol (stretch).
"""
MechanicalTest{P<:MechanicalProtocol, C<:AbstractCondition, G<:AbstractGeometry} = ExperimentData{TensileMeasurement, P, C, G}

"""
Super type fo ExperimentData with thermal measurement (specific heat capacity) and thermal protocol (temperature).
"""
ThermalTest{C<:AbstractCondition, G<:AbstractGeometry} = ExperimentData{ThermalMeasurement, TemperatureSweepProtocol, C, G}

"""
Super type for ExperimentData with dielectric measurement (permittivity) and thermal protocol (temperature).
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


# --- Experiment getters ---

stretches(d::ExperimentData)            = stretches(d.protocol)
temperatures(d::ExperimentData)         = temperatures(d.protocol)
frequencies(d::ExperimentData)          = frequencies(d.protocol)
independent_variable(d::ExperimentData) = independent_variable(d.protocol)
time_step(d::ExperimentData)            = time_step(d.protocol)
rate(d::ExperimentData)                 = rate(d.protocol)
max_stretch(d::ExperimentData)          = max_stretch(d.protocol)

temperature(d::ExperimentData)          = temperature(d.condition)
voltage(d::ExperimentData)              = voltage(d.condition)

thickness(d::ExperimentData)            = thickness(d.geometry)
electric_field(d::ExperimentData)       = electric_field(d.condition, d.geometry)


# --- Experiment labels ---

pretty_label(f::Function, d)            = string(f(d))
pretty_label(f::typeof(rate), d)        = @sprintf("%.2f ", f(d)) * rate_units(d)
pretty_label(f::typeof(max_stretch), d) = @sprintf("%3.0f %%", 100*(f(d)-1))
pretty_label(f::typeof(temperature), d) = @sprintf("%2.0f ºC", f(d)-273.15)
pretty_label(f::typeof(voltage), d)     = @sprintf("%4d V", f(d))
pretty_label(fs::Tuple{Vararg{Function}}, d) = join(map(f -> pretty_label(f, d), fs), ", ")

rate_units(d::ExperimentData) = rate_units(d.protocol)
rate_units(::TemperatureSweepProtocol) = "K/min"
rate_units(::MechanicalProtocol) = "/s"
