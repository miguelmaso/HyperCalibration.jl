module HyperCalibrationPlots

using HyperCalibration
using HyperFEM
using RecipesBase  # Provided by Plots

import HyperFEM.PhysicalModels: PhysicalModel

# global COLOR_INDEX = 1

# function Plots.plot(model::PhysicalModel, data::ExperimentData; kwargs...)
#   p = plot()
#   plot!(model, data; kwargs...)
#   p
# end

# function Plots.plot!(model::PhysicalModel, data::DifferentialScanningCalorimetryTest; kwargs...)
#   get!(kwargs, :color, COLOR_INDEX)
#   x_data = temperatures(data.protocol)
#   y_true, y_pred = experiment_prediction(model, data)
#   p = plot!(x_data, y_pred; kwargs...)
#   p = scatter!(x_data, y_true; kwargs...)
#   global COLOR_INDEX += 1
#   p
# end

@recipe function f(model::PhysicalModel, data::DifferentialScanningCalorimetryTest)
    
  x_data = temperatures(data.protocol)
  y_true, y_pred = experiment_prediction(model, data)
  
  # 1. Plot the Prediction Line
  @series begin
    seriestype := :line
    label --> "Prediction"
    x_data, y_pred
  end

  # 2. Plot the Experimental Data Scatter. Plots naturally groups the two series under the same color cycle index.
  @series begin
    seriestype := :scatter
    label --> "Experimental True"
    x_data, y_true
  end
end

end
