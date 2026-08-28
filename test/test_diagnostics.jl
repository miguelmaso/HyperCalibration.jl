using Test

@testset "Diagnostics" begin

  struct DummyModel <: PhysicalModel
    a::Float64
    b::Float64
  end

  builder = (a, b) -> DummyModel(a, b)
  params = [2.0, 1.0]

  data = [[1.0, 2.0, 3.0, 4.0]]
  HyperCalibration.npoints(d::Vector{Vector{Float64}}) = length(first(d))
  HyperCalibration.experiment_prediction(m::DummyModel, d::Vector{Vector{Float64}}) = experiment_prediction(m, first(d))
  HyperCalibration.experiment_prediction(m::DummyModel, d::Vector{Float64}) = (
    [2.1, 4.1, 5.9, 8.2],  # y_true
    m.a * d .+ m.b         # y_pred
  )

  res = CalibrationResult(builder, params, data)

  # Covariance & R²
  cov_mat, JtJ, sse = covariance_matrix(res)
  @test size(cov_mat) == (2, 2)
  @test size(JtJ) == (2, 2)
  @test sse > 0.0

  r2 = r2_score(res)
  @test 0.0 <= r2 <= 1.0

  # Parameter Stats & Pretty Print
  stats = parameter_stats(res; names=["a", "b"])

  io_plain = IOBuffer()
  show(io_plain, MIME("text/plain"), stats)
  @test !isempty(String(take!(io_plain)))

  io_latex = IOBuffer()
  show(io_latex, MIME("text/latex"), stats)
  @test !isempty(String(take!(io_latex)))

  # Parameter Sampling
  n_samples = 10
  samples = sample_parameters(res, n_samples)
  @test size(samples) == (length(params), n_samples)

end
