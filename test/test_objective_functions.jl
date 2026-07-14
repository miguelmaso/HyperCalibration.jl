
@testset "ObjectiveFunctions" begin
  model = NeoHookean3D(λ=100.0, μ=10.0)
  protocol = QuasiStaticProtocol{Uniaxial}([1.0, 1.1, 1.2])
  predicted = evaluate_stress(model, protocol, StandardCondition(), PlateGeometry())

  data = UniaxialQuasiStaticTest(12, stretches(protocol), [0.0, 1.0, 2.0])
  observed, replayed = experiment_prediction(model, data)
  @test observed == data.measurement.σ
  @test replayed == predicted
  @test loss(model, data) ≥ 0
  @test loss(model, [data]) == loss(model, data)
end
