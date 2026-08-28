using Test

@testset "Diagnostics" begin

  struct DummyModel <: PhysicalModel
    a::Float64
    b::Float64
  end

  builder = (a, b) -> DummyModel(a, b)
  params = [2.0, 1.0]

  data = [1, 2, 3, 4]
  HyperCalibration.npoints(d::Vector) = length(d)
  HyperCalibration.experiment_prediction(m::DummyModel, d::Vector) = (
    Float64[2.1, 4.1, 5.9, 8.2], # y_true
    m.a .* d .+ m.b              # y_pred
  )

  res = CalibrationResult(builder, params, [data])

  @testset "Residuals & Jacobian" begin
    r = residuals(res.model, data)
    @test length(r) == 4
    @test r isa Vector{Float64}

    J = finite_difference_jacobian(builder, params, data)
    @test size(J) == (4, 2)
  end

  @testset "Covariance & R²" begin
    cov_mat, JtJ, sse = covariance_matrix(res)
    @test size(cov_mat) == (2, 2)
    @test size(JtJ) == (2, 2)
    @test sse > 0.0

    r2 = r2_score(res)
    @test 0.0 <= r2 <= 1.0
  end

  @testset "Parameter Stats & Pretty Print" begin
    stats = parameter_stats(res; names=["a", "b"])
    
    @test length(stats.std_errs) == 2
    @test length(stats.ci_bounds) == 2
    @test length(stats.sensitivities) == 2

    io_plain = IOBuffer()
    show(io_plain, MIME("text/plain"), stats)
    @test !isempty(String(take!(io_plain)))

    io_latex = IOBuffer()
    show(io_latex, MIME("text/latex"), stats)
    @test !isempty(String(take!(io_latex)))
  end

  @testset "Parameter Sampling" begin
    n_samples = 50
    samples = sample_parameters(res, n_samples)
    @test size(samples) == (length(params), n_samples)
  end

end
