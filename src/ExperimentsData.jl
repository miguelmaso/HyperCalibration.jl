
abstract type AbstractMeasurement end
abstract type AbstractProtocol end
abstract type AbstractCondition end

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

struct SequentialProtocol <: AbstractProtocol
  stages::Vector{AbstractProtocol}
end

struct QuasiStaticProtocol{K<:Kinematics} <: AbstractProtocol
  λ::Vector{Float64}
end

struct CyclicLoadingProtocol{K<:Kinematics} <: AbstractProtocol
  v::Float64
  Δt::Float64
  λ::Vector{Float64}
end

struct CreepProtocol{K<:Kinematics} <: AbstractProtocol
  t::Vector{Float64}
end

struct TemperatureSweepProtocol <: AbstractProtocol
  θ::Vector{Float64}
end

struct FrequencySweepProtocol <: AbstractProtocol
  f::Vector{Float64}
end

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

# --- Experiments ---

mutable struct ExperimentData{M<:AbstractMeasurement, P<:AbstractProtocol, C<:AbstractCondition}
  const measurement::M
  const protocol::P
  const condition::C
  weight::Float64
end
