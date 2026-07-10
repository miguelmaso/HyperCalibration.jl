using HyperCalibration, HyperFEM
using Optim
using CSV
using Plots

experiments = CSV.read(joinpath(@__DIR__, "data/set 2 quasi-static.csv"), UniaxialQuasiStaticTest, decimal=',')

build_model(μ, N) = EightChain(μ=μ, N=N)
pn = [  "μ",   "N"]  # Parameter names
p0 = [  1e4,   1.0]  # Initial seed

f(p) = loss(build_model, p, experiments)

results = optimize(f, p0, NelderMead())
model = build_model(results.minimizer...)
plot(model, experiments[1], xlabel="Stretch [-]", ylabel="Stress [KPa]", units_scale=1e-3)
