module HyperCalibrationPlots

using HyperCalibration
using HyperFEM
using Plots.RecipesBase

import HyperFEM.PhysicalModels: PhysicalModel  # FIXME: This import must be removed after the HyperFEM 0.0.4 release


@recipe function f(model::PhysicalModel, datasets::AbstractVector{<:ExperimentData})

  user_labels = get(plotattributes, :label, nothing)
  user_colors = get(plotattributes, :seriescolor, nothing)

  for (i, data) in enumerate(datasets)
    @series begin

      if user_labels isa AbstractArray && length(user_labels) == length(datasets)
        label := user_labels[mod1(i, length(user_labels))]
      end

      if user_colors isa AbstractArray && length(user_colors) == length(datasets)
        seriescolor := user_colors[mod1(i, length(user_colors))]
      end

      model, data  # Re-routes each item back to the single-experiment recipe
    end
  end

end


@recipe function f(model::PhysicalModel, data::ExperimentData)

  scale = pop!(plotattributes, :units_scale, 1.0)

  x_data = independent_variable(data.protocol)
  y_true, y_pred = experiment_prediction(model, data) .* scale

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
