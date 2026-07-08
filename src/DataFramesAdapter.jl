using CSV, DataFrames

function UniaxialQuasiStaticTest(df; weight=1.0)
  id = df.id[1]
  λ  = df.stretch
  σ  = df.stress
  measurement = TensileMeasurement(σ)
  protocol = QuasiStaticProtocol{Uniaxial}(λ)
  condition = StandardCondition()
  geometry = PlateGeometry()
  ExperimentData(measurement, protocol, condition, geometry, id, weight)
end

function UniaxialCyclicLoadingTest(df; weight=1.0)
  id = df.id[1]
  v  = df.vel[1]
  λ  = df.stretch
  σ  = df.stress
  i_load = argmax(λ)-1  # This is the loading branch
  Δt = (λ[i_load]-λ[1]) / (i_load-1) / v
  measurement = TensileMeasurement(σ)
  protocol = CyclicLoadingProtocol{Uniaxial}(λ, v, Δt)
  condition = StandardCondition()
  geometry = PlateGeometry()
  ExperimentData(measurement, protocol, condition, geometry, id, weight)
end

function UniaxialThermalQuasiStaticTest(df; weight=1.0)
  id = df.id[1]
  λ  = df.stretch
  σ  = df.stress
  θ  = df.temp[1]
  measurement = TensileMeasurement(σ)
  protocol = QuasiStaticProtocol{Uniaxial}(λ)
  condition = IsothermalCondition(θ)
  geometry = PlateGeometry()
  ExperimentData(measurement, protocol, condition, geometry, id, weight)
end

function UniaxialThermalCyclicLoadingTest(df; weight=1.0)
  id = df.id[1]
  v  = df.vel[1]
  λ  = df.stretch
  σ  = df.stress
  θ  = df.temp[1]
  i_load = argmax(λ)-1  # This is the loading branch
  Δt = (λ[i_load]-λ[1]) / (i_load-1) / v
  measurement = TensileMeasurement(σ)
  protocol = CyclicLoadingProtocol{Uniaxial}(λ, v, Δt)
  condition = IsothermalCondition(θ)
  geometry = PlateGeometry()
  ExperimentData(measurement, protocol, condition, geometry, id, weight)
end

function UniaxialThermoElectricCyclicLoadingTest(df; thickness, weight=1.0)
  id = df.id[1]
  v  = df.vel[1]
  λ  = df.stretch
  σ  = df.stress
  θ  = df.temp[1]
  V  = df.voltage[1]
  i_load = argmax(λ)-1  # This is the loading branch
  Δt = (λ[i_load]-λ[1]) / (i_load-1) / v
  measurement = TensileMeasurement(σ)
  protocol = CyclicLoadingProtocol{Uniaxial}(λ, v, Δt)
  condition = ThermoElectricalCondition(θ, V)
  geometry = PlateGeometry(thickness)
  ExperimentData(measurement, protocol, condition, geometry, id, weight)
end

function DifferentialScanningCalorimetryTest(df, weight=1.0)
  id = df.id[1]
  r  = df.rate[1]
  θ  = df.temp
  cv = df.cv
  measurement = ThermalMeasurement(cv)
  protocol = TemperatureSweepProtocol(θ, r)
  condition = StandardCondition()
  geometry = PlateGeometry()
  ExperimentData(measurement, protocol, condition, geometry, id, weight)
end

function DielectricSpectroscopyTest(df, weight=1.0)
  id = df.id[1]
  θ = df.temp[1]
  f = df.freq
  ε = df.dielec
  measurement = DielectricMeasurement(ε)
  protocol = FrequencySweepProtocol(f)
  condition = StandardCondition()
  geometry = PlateGeometry()
  ExperimentData(measurement, protocol, condition, geometry, id, weight)
end

function load_data(filepath::String, experiment_type::Type; decimal='.', kwargs...)
  df = CSV.read(filepath, DataFrame; decimal=decimal)
  grouped = groupby(df, :id)
  [experiment_type(g; kwargs...) for g in grouped]
end
