using Revise
using HyperCalibration
using HyperFEM
using Optimization, OptimizationOptimJL, OptimizationMetaheuristics
using Plots


## Laboratory constants
θr = 293.15    # Reference temperature, ºK
t0 = 0.0005    # Specimen thickness, m (0.5mm)
αr = 1.8e-4    # Thermal expansion, /ºK (extracted from 3M VHB technical data sheet)
ε0 = 8.85e-12  # Air permittivity


## Load data

set_1_cal   = load_data(abspath(dirname(@__FILE__), "data/set 1 calorimetry.csv"), DifferentialScanningCalorimetryTest, decimal=',')
set_2_quasi = load_data(abspath(dirname(@__FILE__), "data/set 2 quasi-static.csv"), UniaxialQuasiStaticTest, decimal=',')
# set_3_creep = load_data(abspath(dirname(@__FILE__), "data/set 3 creep.csv"), UniaxialRelaxationTest, decimal=',')
set_4_load  = load_data(abspath(dirname(@__FILE__), "data/set 4 loading.csv"), UniaxialThermalCyclicLoadingTest, decimal=',')
set_5_elec  = load_data(abspath(dirname(@__FILE__), "data/set 5 dielectric.csv"), DielectricSpectroscopyTest, decimal=',')
set_6_coupl = load_data(abspath(dirname(@__FILE__), "data/set 6 coupled.csv"), UniaxialThermoElectricCyclicLoadingTest, decimal=',', thickness=t0)


## Step 1: Thermal characterization

build_heat(cv0, γv, κr) = ThermalVolumetric(cv0=cv0, θr=θr, α=αr, κr=κr, κ=1.0, γ=γv)

pn = ["cv0", "γv",  "κr" ]  # Parameter names
p0 = [1.0e6,  0.5, 1.0e9 ]  # Initial seed
lb = [ 10.0,  0.0, 1.0e8 ]  # Minimum search limits
ub = [1.0e8,  1.0, 1.0e10]  # Maximum search limits

opt_func = OptimizationFunction((p, d) -> loss(build_heat, p, d))
opt_prob = OptimizationProblem(opt_func, p0, set_1_cal, lb=lb, ub=ub)
opt_heat = solve(opt_prob, ParticleSwarm(lower=lb, upper=ub, n_particles=100), maxiters=1000, maxtime=30)
sol_heat = opt_heat.u

model = build_heat(sol_heat...)
plot(model, set_1_cal[1], label=["Prediction" "Experiment"], xlabel="T [ºC]", ylabel="cv [J/m³·ºK]")


# r2 = stats(build_heat, sol_heat, set_1_cal, pn)


## Step 2: Hyperelastic characterization

build_longterm(μ, N) = EightChain(μ=μ, N=N)
pn = [  "μ",   "N"]  # Parameter names
p0 = [  1e4,   1.0]  # Initial seed
lb = [  1e3,   0.0]  # Lower search limits
ub = [  1e5,  50.0]  # Upper search limits

build_longterm(μ1, μ2, α1, α2) = NonlinearMooneyRivlin3D(λ=0.0, μ1=μ1, μ2=μ2, α1=α1, α2=α2)
pn = [ "μ1", "μ2", "α1", "α2"]  # Parameter names
p0 = [  1e4,  1e4,  0.8,  0.8]  # Initial seed
lb = [  1e2,  1e3,  0.5,  0.5]  # Lower search limits
ub = [  1e5,  1e5,  3.0,  2.0]  # Upper search limits

build_longterm(C1, C2, C3) = Yeoh3D(λ=0.0, C10=C1, C20=C2, C30=C3)
pn = ["C10",  "C20",  "C30"]  # Parameter names
p0 = [  3e4,   -2e2,    3e0]  # Initial seed
lb = [1.0e3, -2.0e3,  0.0e0]  # Minimum search limits
ub = [2.0e5,  2.0e3,  2.0e2]  # Maximum search limits

opt_func = OptimizationFunction((p,d) -> loss(build_longterm, p, d))
opt_prob_ps = OptimizationProblem(opt_func, p0, set_2_quasi, lb=lb, ub=ub)
opt_long_ps = solve(opt_prob_ps, ParticleSwarm(lower=lb, upper=ub, n_particles=1000), maxiters=1000, maxtime=30)
opt_prob_nm = OptimizationProblem(opt_func, opt_long_ps.u, set_2_quasi)
opt_long_nm = solve(opt_prob_nm, NelderMead(), maxiters=100, maxtime=30)
sol_long = opt_long_nm.u

model = build_longterm(sol_long...)
plot(model, set_2_quasi[1])

# r2 = stats(build_longterm, sol_long, set_2_quasi, pn)
# text_par = text(join(map((n,v) -> @sprintf("%s=%.2g",n,v), pn, sol_long), "\n"), 12, :left)
# text_r2 = text(@sprintf("R² = %.1f %%", 100*r2), 12, :left)

# p = plot(xlabel="Stretch [-]", ylabel="Stress [KPa]")
# plot_experiment!(model, getfirst(r -> r.θ ≈ θr, set_2_quasi))
# annotate!((0.05, 0.8), text_r2, relative=true)
# display(p);

