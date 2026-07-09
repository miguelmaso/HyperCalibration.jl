
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
