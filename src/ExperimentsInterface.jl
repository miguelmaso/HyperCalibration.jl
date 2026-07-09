
# --- Constructors ---

function UniaxialQuasiStaticTest(id::Int, λ::AbstractArray, σ::AbstractArray, weight=1.0)
  measurement = TensileMeasurement(σ)
  protocol = QuasiStaticProtocol{Uniaxial}(λ)
  condition = StandardCondition()
  geometry = PlateGeometry()
  ExperimentData(measurement, protocol, condition, geometry, id, weight)
end

function UniaxialCyclicLoadingTest(id::Int, λ::AbstractArray, v::Real, σ::AbstractArray, weight=1.0)
  i_load = argmax(λ)-1  # This is the loading branch
  Δt = (λ[i_load]-λ[1]) / (i_load-1) / v
  measurement = TensileMeasurement(σ)
  protocol = CyclicLoadingProtocol{Uniaxial}(λ, v, Δt)
  condition = StandardCondition()
  geometry = PlateGeometry()
  ExperimentData(measurement, protocol, condition, geometry, id, weight)
end

function UniaxialThermalQuasiStaticTest(id::Int, λ::AbstractArray, σ::AbstractArray, θ::Real, weight=1.0)
  measurement = TensileMeasurement(σ)
  protocol = QuasiStaticProtocol{Uniaxial}(λ)
  condition = IsothermalCondition(θ)
  geometry = PlateGeometry()
  ExperimentData(measurement, protocol, condition, geometry, id, weight)
end

function UniaxialThermalCyclicLoadingTest(id::Int, λ::AbstractArray, v::Real, σ::AbstractArray, θ::Real, weight=1.0)
  i_load = argmax(λ)-1  # This is the loading branch
  Δt = (λ[i_load]-λ[1]) / (i_load-1) / v
  measurement = TensileMeasurement(σ)
  protocol = CyclicLoadingProtocol{Uniaxial}(λ, v, Δt)
  condition = IsothermalCondition(θ)
  geometry = PlateGeometry()
  ExperimentData(measurement, protocol, condition, geometry, id, weight)
end

function UniaxialThermoElectricCyclicLoadingTest(id::Int, λ::AbstractArray, v::Real, σ::AbstractArray, θ::Real, V::Real, thickness::Real, weight=1.0)
  i_load = argmax(λ)-1  # This is the loading branch
  Δt = (λ[i_load]-λ[1]) / (i_load-1) / v
  measurement = TensileMeasurement(σ)
  protocol = CyclicLoadingProtocol{Uniaxial}(λ, v, Δt)
  condition = ThermoElectricalCondition(θ, V)
  geometry = PlateGeometry(thickness)
  ExperimentData(measurement, protocol, condition, geometry, id, weight)
end

function DifferentialScanningCalorimetryTest(id::Int, θ::AbstractArray, v::Real, cv::AbstractArray, weight=1.0)
  measurement = ThermalMeasurement(cv)
  protocol = TemperatureSweepProtocol(θ, v)
  condition = StandardCondition()
  geometry = PlateGeometry()
  ExperimentData(measurement, protocol, condition, geometry, id, weight)
end

function DielectricSpectroscopyTest(id::Int, f::AbstractArray, ε::AbstractArray, θ::Real, weight=1.0)
  measurement = DielectricMeasurement(ε)
  protocol = FrequencySweepProtocol(f)
  condition = IsothermalCondition(θ)
  geometry = PlateGeometry()
  ExperimentData(measurement, protocol, condition, geometry, id, weight)
end


# --- Print ---

function Base.show(io::IO, d::ExperimentData)
  details = experiment_details(d)
  print(io, "[$(d.id)]: $(details) w=$(d.weight)")
end

function Base.summary(io::IO, data::AbstractVector{T}) where {T<:ExperimentData}
  print(io, "$(length(data))-element Vector{$(experiment_type_name(T))}")
end

experiment_details(::ExperimentData) = ""

function experiment_details(d::UniaxialQuasiStaticTest)
  λ_max = pretty_label(max_stretch, d)
  "$(λ_max), "
end

function experiment_details(d::UniaxialCyclicLoadingTest)
  λ_max = pretty_label(max_stretch, d)
  r     = pretty_label(rate, d)
  "$(λ_max), $(r), "
end

function experiment_details(d::UniaxialThermalQuasiStaticTest)
  λ_max = pretty_label(max_stretch, d)
  θ     = pretty_label(temperature, d)
  "$(λ_max), $(θ), "
end

function experiment_details(d::UniaxialThermalCyclicLoadingTest)
  λ_max = pretty_label(max_stretch, d)
  r     = pretty_label(rate, d)
  θ     = pretty_label(temperature, d)
  "$(λ_max), $(r), $(θ), "
end

function experiment_details(d::UniaxialThermoElectricCyclicLoadingTest)
  λ_max = pretty_label(max_stretch, d)
  r     = pretty_label(rate, d)
  θ     = pretty_label(temperature, d)
  V     = pretty_label(voltage, d)
  "$(λ_max), $(r), $(θ), $(V), "
end

function experiment_details(d::DifferentialScanningCalorimetryTest)
  r = pretty_label(rate, d)
  "$(r), "
end

function experiment_details(d::DielectricSpectroscopyTest)
  θ = pretty_label(temperature, d)
  "$(θ), "
end


experiment_type_name(::Type{T}) where {T<:ExperimentData} = string(T)

for test_type in (
  :UniaxialQuasiStaticTest,
  :UniaxialCyclicLoadingTest,
  :UniaxialRelaxationTest,
  :UniaxialThermalQuasiStaticTest,
  :UniaxialThermalCyclicLoadingTest,
  :UniaxialThermalRelaxationTest,
  :UniaxialThermoElectricCyclicLoadingTest,
  :DifferentialScanningCalorimetryTest,
  :DielectricSpectroscopyTest
)
  @eval experiment_type_name(::Type{$test_type}) = $(string(test_type))
end
