
abstract type AbstractMeasurement end
abstract type AbstractProtocol end
abstract type AbstractCondition end
abstract type AbstractGeometry end

# --- Kinematics ---

abstract type Kinematics end
struct Uniaxial <: Kinematics end
struct Biaxial <: Kinematics end

# --- Measurements ---

struct TensileMeasurement <: AbstractMeasurement
  σ::Vector{Float64}
end

struct ThermalMeasurement <: AbstractMeasurement
  cv::Vector{Float64}
end

struct DielectricMeasurement <: AbstractMeasurement
  ε::Vector{Float64}
end

# --- Protocols ---

abstract type MechanicalProtocol{K<:Kinematics} <: AbstractProtocol end

struct SequentialProtocol <: AbstractProtocol
  stages::Vector{AbstractProtocol}
end

struct QuasiStaticProtocol{K<:Kinematics} <: MechanicalProtocol{K}
  λ::Vector{Float64}
end

struct CyclicLoadingProtocol{K<:Kinematics} <: MechanicalProtocol{K}
  v::Float64
  Δt::Float64
  λ::Vector{Float64}
end

struct CreepProtocol{K<:Kinematics} <: MechanicalProtocol{K}
  λ::Float64
  t::Vector{Float64}
end

struct TemperatureSweepProtocol <: AbstractProtocol
  θ::Vector{Float64}
end

struct FrequencySweepProtocol <: AbstractProtocol
  f::Vector{Float64}
end

stretches(::AbstractProtocol) = throw(ArgumentError("stretches not defined for this protocol"))

stretches(p::SequentialProtocol) = vcat(map(stretches, p.stages)...)

stretches(p::QuasiStaticProtocol) = p.λ

stretches(p::CyclicLoadingProtocol) = p.λ

time_step(::AbstractProtocol) = throw(ArgumentError("time_step not defined for this protocol"))

time_step(p::CyclicLoadingProtocol) = p.Δt

time_step(p::CreepProtocol) = diff(p.t)[1]

temperatures(::AbstractProtocol) = throw(ArgumentError("temperatures not defined for this protocol"))

temperatures(p::TemperatureSweepProtocol) = p.θ

# --- Conditions ---

struct StandardCondition <: AbstractCondition
end

struct IsothermalCondition <: AbstractCondition
  θ::Float64
end

struct ElectricalCondition <: AbstractCondition
  V::Float64
end

struct ThermoElectricalCondition <: AbstractCondition
  θ::Float64
  V::Float64
end

temperature(::AbstractCondition) = 293.15

temperature(c::IsothermalCondition) = c.θ

temperature(c::ThermoElectricalCondition) = c.θ

voltage(::AbstractCondition) = 0.0

voltage(c::ElectricalCondition) = c.V

voltage(c::ThermoElectricalCondition) = c.V

# --- Geometries ---

struct PlateGeometry <: AbstractGeometry
  t0::Float64  # reference thickness
  A0::Float64  # Reference cross sectional area
end

thickness(::AbstractGeometry) = throw(ArgumentError("thickness not defined for this geometry"))

thickness(g::PlateGeometry) = g.t0

electric_field(::AbstractCondition, ::AbstractGeometry) = voltage(c) / thickness(g)

# --- Experiments ---

mutable struct ExperimentData{M<:AbstractMeasurement, P<:AbstractProtocol, C<:AbstractCondition, G<:AbstractGeometry}
  const measurement::M
  const protocol::P
  const condition::C
  const geometry::G
  weight::Float64
end
