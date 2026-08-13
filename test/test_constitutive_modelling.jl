
@testset "ConstitutiveModelling" begin
  model = NeoHookean3D(λ=100.0, μ=10.0)
  protocol = QuasiStaticProtocol{Uniaxial}([1.0, 1.1, 1.2])
  predicted = evaluate_stress(model, protocol, StandardCondition(), PlateGeometry())

  @test length(predicted) == length(stretches(protocol))
  @test all(isfinite, predicted)
  @test predicted[1] ≈ 0.0 atol=1e-12
end
