
"""
Return the `y_true` and `y_pred` sampling points for a given constitutive model and experimental test.
"""
function experiment_prediction(model::PhysicalModel, data::MechanicalTest)
  y_true = data.measurement.σ
  y_pred = evaluate_stress(model, data.protocol, data.condition, data.geometry)
  return y_true, y_pred
end

function experiment_prediction(model::PhysicalModel, data::ThermalTest)
  y_true = data.measurement.cv
  y_pred = evaluate_cv(model, data.protocol, data.condition, data.geometry)
  return y_true, y_pred
end

function experiment_prediction(model::PhysicalModel, data::ThermoDielectricTest)
  y_true = data.measurement.ε
  y_pred = evaluate_epsilon(model, data.protocol, data.condition, data.geometry)
  return y_true, y_pred
end

"""
Compute the loss function for a given constitutive model and one/multiple experimental tests
as the weighted normalized squared difference between the predicted and the measured values.
"""
function loss(model::PhysicalModel, data::ExperimentData)
  y_true, y_pred = experiment_prediction(model, data)
  y_max = maximum(abs.(y_true))
  s2 = zero(eltype(y_true))
  @inbounds for i in eachindex(y_true)
    s2 += abs2( (y_true[i] - y_pred[i]) / y_max )
  end
  return data.weight * s2 / length(y_true)
end

function loss(model::PhysicalModel, data::Vector{<:ExperimentData})
  sum(loss(model, d) for d in data) / sum(d.weight for d in data)
end

function loss(model_builder, params, data)
  model = model_builder(params...)
  loss(model, data)
end
