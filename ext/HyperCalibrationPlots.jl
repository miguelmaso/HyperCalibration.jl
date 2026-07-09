module HyperCalibrationPlots

using HyperCalibration
using HyperFEM
using Plots.RecipesBase

import HyperFEM.PhysicalModels: PhysicalModel  # FIXME: This import must be removed after the HyperFEM 0.0.4 release

@recipe function f(model::PhysicalModel, data::ExperimentData)

  x_data = independent_variable(data.protocol)
  y_true, y_pred = experiment_prediction(model, data)

  user_label = get(plotattributes, :label, nothing)
  user_color = get(plotattributes, :seriescolor, nothing)

  split_series = (user_label isa AbstractArray && length(user_label) > 1) || 
                 (user_color isa AbstractArray && length(user_color) > 1)

  # 1. Prediction line
  @series begin
    seriestype := :line
    x_data, y_pred
  end

  # 2. Experimental data scatter.
  @series begin
    seriestype := :scatter
    primary    := split_series
    x_data, y_true
  end
end

end
