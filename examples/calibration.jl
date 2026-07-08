using Revise
using HyperCalibration
using HyperFEM
using Optimization, OptimizationOptimJL, OptimizationMetaheuristics


## Laboratory constants
θr = 293.15    # Reference temperature, ºK
t0 = 0.0005    # Specimen thickness, m (0.5mm)
αr = 1.8e-4    # Thermal expansion, /ºK (extracted from 3M VHB technical data sheet)
ε0 = 8.85e-12  # Air permittivity


## Load data

set_1_cal   = load_data(abspath(dirname(@__FILE__), "data/set 1 calorimetry.csv"), DifferentialScanningCalorimetryTest, decimal=',')
set_2_load  = load_data(abspath(dirname(@__FILE__), "data/set 2 loading.csv"), UniaxialCyclicLoadingTest, decimal=',')
# set_3_creep = load_data(abspath(dirname(@__FILE__), "data/set 3 creep.csv"), UniaxialRelaxationTest, decimal=',')
set_4_quasi = load_data(abspath(dirname(@__FILE__), "data/set 4 quasi-static.csv"), UniaxialQuasiStaticTest, decimal=',')
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
opt_heat = solve(opt_prob, ParticleSwarm(lower=lb, upper=ub, n_particles=100), maxiters=1000, maxtime=60)
sol_heat = opt_heat.u

# model = build_heat(sol_heat...)
# r2 = stats(build_heat, sol_heat, set_1_cal, pn)
# text_r2 = text(@sprintf("R² = %.0f %%", 100*r2), 12, :left)

# p = plot(xlabel="T [ºC]", ylabel="cv [J/m³·ºK]")
# plot_experiment!(model, set_1_cal[1])
# annotate!((0.05, 0.8), text_r2, relative=true)
# display(p);

