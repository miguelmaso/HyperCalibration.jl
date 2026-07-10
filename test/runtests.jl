using Test
using HyperCalibration
using HyperFEM

@testset "HyperCalibration" begin
  @testset "experiment data contracts" begin
    @test TensileMeasurement([1, 2]).σ == [1.0, 2.0]
    @test QuasiStaticProtocol{Uniaxial}([1, 2]).λ == [1.0, 2.0]
    @test IsothermalCondition(293).θ == 293.0
    @test PlateGeometry(1, 2) == PlateGeometry(1.0, 2.0)

    λ = [1, 11 // 10, 6 // 5]
    σ = [0, 12, 25]
    data = UniaxialQuasiStaticTest(7, view(λ, :), view(σ, :), 0.5)

    @test data isa UniaxialQuasiStaticTest
    @test data.id == 7
    @test data.weight == 0.5
    @test stretches(data) == [1.0, 1.1, 1.2]
    @test independent_variable(data) == stretches(data)
    @test max_stretch(data) == 1.2
    @test data.measurement.σ == [0.0, 12.0, 25.0]
    @test pretty_label(max_stretch, data) == " 20 %"
    @test occursin("UniaxialQuasiStaticTest", sprint(summary, [data]))
    @test occursin("[7]:  20 %,  w=0.5", sprint(show, data))
    @test UniaxialQuasiStaticTest([1.0, 1.1], [0.0, 1.0]).id == 0

    cyclic = UniaxialCyclicLoadingTest(8, [1, 1.1, 1.2, 1.1, 1], 0.2, zeros(5))
    @test rate(cyclic) == 0.2
    @test time_step(cyclic) ≈ 0.5
    @test pretty_label(rate, cyclic) == "0.20 /s"

    thermal = DifferentialScanningCalorimetryTest(9, [293, 303], 2, [1, 2])
    @test temperatures(thermal) == [293.0, 303.0]
    @test rate(thermal) == 2.0
    @test pretty_label(rate, thermal) == "2.00 K/min"

    dielectric = DielectricSpectroscopyTest(10, [1, 10], [2, 3], 293)
    @test frequencies(dielectric) == [1.0, 10.0]
    @test temperature(dielectric) == 293.0

    coupled = UniaxialThermoElectricCyclicLoadingTest(11, [1, 1.1, 1.2, 1.1, 1], 0.2, zeros(5), 300, 50, 0.5)
    @test HyperCalibration.electric_field(coupled.condition, coupled.geometry) == 100.0
    @test_throws ArgumentError temperature(StandardCondition())
    @test_throws ArgumentError voltage(StandardCondition())
  end

  @testset "reduced-kinematics constitutive integration" begin
    model = NeoHookean3D(λ=100.0, μ=10.0)
    protocol = QuasiStaticProtocol{Uniaxial}([1.0, 1.1, 1.2])
    predicted = evaluate_stress(model, protocol, StandardCondition(), PlateGeometry())

    @test length(predicted) == length(stretches(protocol))
    @test all(isfinite, predicted)
    @test predicted[1] ≈ 0.0 atol=1e-12

    data = UniaxialQuasiStaticTest(12, stretches(protocol), [0.0, 1.0, 2.0])
    observed, replayed = experiment_prediction(model, data)
    @test observed == data.measurement.σ
    @test replayed == predicted
    @test loss(model, data) ≥ 0
    @test loss(model, [data]) == loss(model, data)
  end
end
