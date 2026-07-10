
# --- Specialisation of ExperimentData ---

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


# --- Constructors ---

function UniaxialQuasiStaticTest(id::Integer, λ::AbstractArray, σ::AbstractArray, weight=1.0)
  measurement = TensileMeasurement(σ)
  protocol = QuasiStaticProtocol{Uniaxial}(λ)
  condition = StandardCondition()
  geometry = PlateGeometry()
  ExperimentData(id, measurement, protocol, condition, geometry, weight)
end
UniaxialQuasiStaticTest(λ::AbstractArray, σ::AbstractArray, weight=1.0) = UniaxialQuasiStaticTest(0, λ, σ, weight)

function UniaxialCyclicLoadingTest(id::Integer, λ::AbstractArray, v::Real, σ::AbstractArray, weight=1.0)
  i_load = argmax(λ)-1  # This is the loading branch
  Δt = (λ[i_load]-λ[1]) / (i_load-1) / v
  measurement = TensileMeasurement(σ)
  protocol = CyclicLoadingProtocol{Uniaxial}(λ, v, Δt)
  condition = StandardCondition()
  geometry = PlateGeometry()
  ExperimentData(id, measurement, protocol, condition, geometry, weight)
end
UniaxialCyclicLoadingTest(λ::AbstractArray, v::Real, σ::AbstractArray, weight=1.0) = UniaxialCyclicLoadingTest(0, λ, v, σ, weight)

function UniaxialThermalQuasiStaticTest(id::Integer, λ::AbstractArray, σ::AbstractArray, θ::Real, weight=1.0)
  measurement = TensileMeasurement(σ)
  protocol = QuasiStaticProtocol{Uniaxial}(λ)
  condition = IsothermalCondition(θ)
  geometry = PlateGeometry()
  ExperimentData(id, measurement, protocol, condition, geometry, weight)
end
UniaxialThermalQuasiStaticTest(λ::AbstractArray, σ::AbstractArray, θ::Real, weight=1.0) = UniaxialThermalQuasiStaticTest(0, λ, σ, θ, weight)

function UniaxialThermalCyclicLoadingTest(id::Integer, λ::AbstractArray, v::Real, σ::AbstractArray, θ::Real, weight=1.0)
  i_load = argmax(λ)-1  # This is the loading branch
  Δt = (λ[i_load]-λ[1]) / (i_load-1) / v
  measurement = TensileMeasurement(σ)
  protocol = CyclicLoadingProtocol{Uniaxial}(λ, v, Δt)
  condition = IsothermalCondition(θ)
  geometry = PlateGeometry()
  ExperimentData(id, measurement, protocol, condition, geometry, weight)
end
UniaxialThermalCyclicLoadingTest(λ::AbstractArray, v::Real, σ::AbstractArray, θ::Real, weight=1.0) = UniaxialThermalCyclicLoadingTest(0, λ, v, σ, θ, weight)

function UniaxialThermoElectricCyclicLoadingTest(id::Integer, λ::AbstractArray, v::Real, σ::AbstractArray, θ::Real, V::Real, thickness::Real, weight=1.0)
  i_load = argmax(λ)-1  # This is the loading branch
  Δt = (λ[i_load]-λ[1]) / (i_load-1) / v
  measurement = TensileMeasurement(σ)
  protocol = CyclicLoadingProtocol{Uniaxial}(λ, v, Δt)
  condition = ThermoElectricalCondition(θ, V)
  geometry = PlateGeometry(thickness)
  ExperimentData(id, measurement, protocol, condition, geometry, weight)
end
UniaxialThermoElectricCyclicLoadingTest(λ::AbstractArray, v::Real, σ::AbstractArray, θ::Real, V::Real, thickness::Real, weight=1.0) = UniaxialThermoElectricCyclicLoadingTest(0, λ, v, σ, θ, V, thickness, weight)

function DifferentialScanningCalorimetryTest(id::Integer, θ::AbstractArray, v::Real, cv::AbstractArray, weight=1.0)
  measurement = ThermalMeasurement(cv)
  protocol = TemperatureSweepProtocol(θ, v)
  condition = StandardCondition()
  geometry = PlateGeometry()
  ExperimentData(id, measurement, protocol, condition, geometry, weight)
end
DifferentialScanningCalorimetryTest(θ::AbstractArray, v::Real, cv::AbstractArray, weight=1.0) = DifferentialScanningCalorimetryTest(0, θ, v, cv, weight)

function DielectricSpectroscopyTest(id::Integer, f::AbstractArray, ε::AbstractArray, θ::Real, weight=1.0)
  measurement = DielectricMeasurement(ε)
  protocol = FrequencySweepProtocol(f)
  condition = IsothermalCondition(θ)
  geometry = PlateGeometry()
  ExperimentData(id, measurement, protocol, condition, geometry, weight)
end
DielectricSpectroscopyTest(f::AbstractArray, ε::AbstractArray, θ::Real, weight=1.0) = DielectricSpectroscopyTest(0, f, ε, θ, weight)


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
