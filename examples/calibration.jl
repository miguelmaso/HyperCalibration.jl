using Revise
using HyperCalibration


set_1_cal   = load_data(abspath(dirname(@__FILE__), "data/set 1 calorimetry.csv"), DifferentialScanningCalorimetryTest, decimal=',')
set_2_load  = load_data(abspath(dirname(@__FILE__), "data/set 2 loading.csv"), UniaxialCyclicLoadingTest, decimal=',')
# set_3_creep = load_data(abspath(dirname(@__FILE__), "data/set 3 creep.csv"), UniaxialRelaxationTest, decimal=',')
set_4_quasi = load_data(abspath(dirname(@__FILE__), "data/set 4 quasi-static.csv"), UniaxialQuasiStaticTest, decimal=',')
set_5_elec  = load_data(abspath(dirname(@__FILE__), "data/set 7 dielectric.csv"), DielectricSpectroscopyTest, decimal=',')
set_6_coupl = load_data(abspath(dirname(@__FILE__), "data/set 8 coupled.csv"), UniaxialThermoElectricCyclicLoadingTest, decimal=',', thickness=0.0005)
