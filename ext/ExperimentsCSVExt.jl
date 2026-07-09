module ExperimentsCSVExt

using HyperCalibration
using CSV


function CSV.read(filepath::String, ::Type{T}; thickness=0.0, kwargs...) where {T <: ExperimentData}
  data = CSV.File(filepath; kwargs...)

  # Group row indices by id. We use get! to initialize an empty Int array for new IDs, then push! the row index
  grouped_indices = Dict{eltype(data.id), Vector{Int}}()
  for (i, id) in enumerate(data.id)
      push!(get!(grouped_indices, id, Int[]), i)
  end

  return [build_experiment(T, data, id, indices, thickness) for (id, indices) in pairs(grouped_indices)]
end


function build_experiment(::Type{T}, data, id, indices, thickness) where {T}
  error("No build_experiment method defined for type $T")
end

function build_experiment(::Type{UniaxialQuasiStaticTest}, data, id, indices, thickness)
  λ = @view data.stretch[indices]
  σ = @view data.stress[indices]
  UniaxialQuasiStaticTest(id, λ, σ)
end

function build_experiment(::Type{UniaxialCyclicLoadingTest}, data, id, indices, thickness)
  vel = data.vel[indices[1]]
  λ = @view data.stretch[indices]
  σ = @view data.stress[indices]
  UniaxialCyclicLoadingTest(id, λ, vel, σ)
end

function build_experiment(::Type{UniaxialThermalQuasiStaticTest}, data, id, indices, thickness)
  θ = data.temp[indices[1]]
  λ = @view data.stretch[indices]
  σ = @view data.stress[indices]
  UniaxialThermalQuasiStaticTest(id, λ, σ, θ)
end

function build_experiment(::Type{UniaxialThermalCyclicLoadingTest}, data, id, indices, thickness)
  vel = data.vel[indices[1]]
  θ = data.temp[indices[1]]
  λ = @view data.stretch[indices]
  σ = @view data.stress[indices]
  UniaxialThermalCyclicLoadingTest(id, λ, vel, σ, θ)
end

function build_experiment(::Type{UniaxialThermoElectricCyclicLoadingTest}, data, id, indices, thickness)
  vel = data.vel[indices[1]]
  θ = data.temp[indices[1]]
  V = data.voltage[indices[1]]
  λ = @view data.stretch[indices]
  σ = @view data.stress[indices]
  UniaxialThermoElectricCyclicLoadingTest(id, λ, vel, σ, θ, V, thickness)
end

function build_experiment(::Type{DifferentialScanningCalorimetryTest}, data, id, indices, thickness)
  rate = data.rate[indices[1]]
  θ  = @view data.temp[indices]
  cv = @view data.cv[indices]
  DifferentialScanningCalorimetryTest(id, θ, rate, cv)
end

function build_experiment(::Type{DielectricSpectroscopyTest}, data, id, indices, thickness)
  θ = data.temp[indices[1]]
  f = @view data.freq[indices]
  ε = @view data.dielec[indices]
  DielectricSpectroscopyTest(id, f, ε, θ)
end

end