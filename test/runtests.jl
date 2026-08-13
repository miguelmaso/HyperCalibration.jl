using Aqua
using Test
using HyperCalibration
using HyperFEM

@testset "HyperCalibration" begin

  @testset "code quality (Aqua.jl)" begin
    Aqua.test_all(HyperCalibration)
  end

  include("test_experiments_data.jl")
  include("test_constitutive_modelling.jl")
  include("test_objective_functions.jl")
end
