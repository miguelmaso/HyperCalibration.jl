module ExperimentsCSVExt

using HyperCalibration
using CSV


function CSV.read(filepath::String, experiment_type::Type{<:ExperimentData}; thickness=0.0, kwargs...)

  data = CSV.File(filepath; kwargs...)

  # We use get! to initialize an empty Int array for new IDs, then push! the row index
  grouped_indices = Dict{eltype(data.id), Vector{Int}}()
  for (i, id) in enumerate(data.id)
      push!(get!(grouped_indices, id, Int[]), i)
  end

  experiments = []
  for (id, idxs) in pairs(grouped_indices)

    if experiment_type === UniaxialQuasiStaticTest
      λ = @view data.stretch[idxs]
      σ = @view data.stress[idxs]

      push!(experiments, UniaxialQuasiStaticTest(id, λ, σ))

    elseif experiment_type === UniaxialCyclicLoadingTest
      vel = data.vel[idxs[1]]
      λ = @view data.stretch[idxs]
      σ = @view data.stress[idxs]

      push!(experiments, UniaxialCyclicLoadingTest(id, λ, vel, σ))

    elseif experiment_type === UniaxialThermalQuasiStaticTest
      θ = data.temp[idxs[1]]
      λ = @view data.stretch[idxs]
      σ = @view data.stress[idxs]

      push!(experiments, UniaxialThermalQuasiStaticTest(id, λ, σ, θ))

    elseif experiment_type === UniaxialThermalCyclicLoadingTest
      vel = data.vel[idxs[1]]
      θ = data.temp[idxs[1]]
      λ = @view data.stretch[idxs]
      σ = @view data.stress[idxs]

      push!(experiments, UniaxialThermalCyclicLoadingTest(id, λ, vel, σ, θ))

    elseif experiment_type === UniaxialThermoElectricCyclicLoadingTest
      vel = data.vel[idxs[1]]
      θ = data.temp[idxs[1]]
      V = data.voltage[idxs[1]]
      λ = @view data.stretch[idxs]
      σ = @view data.stress[idxs]

      push!(experiments, UniaxialThermoElectricCyclicLoadingTest(id, λ, vel, σ, θ, V, thickness))

    elseif experiment_type === DifferentialScanningCalorimetryTest
      rate = data.rate[idxs[1]]
      θ  = @view data.temp[idxs]
      cv = @view data.cv[idxs]

      push!(experiments, DifferentialScanningCalorimetryTest(id, θ, rate, cv))

    elseif experiment_type === DifferentialScanningCalorimetryTest
      θ = data.temp[idxs[1]]
      f = @view data.freq[idxs]
      ε = @view data.dielec[idxs]

      push!(experiments, DifferentialScanningCalorimetryTest(id, f, ε, θ))

    end
  end
  
  return experiments
end

end