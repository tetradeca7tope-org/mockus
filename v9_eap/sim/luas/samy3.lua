load_platform("input/S5 - 4ft square plate w box and sides fixed.nec")
add_antenna("input/A3 - 150dipole fixed.nec")
add_point(  0.101600,   0.101600,   0.203200)
add_antenna("input/A1 - 75mono fixed.nec")
add_point(  0.304800,  -0.304800,   0.000000)
add_antenna("input/A5 - 300dipole fixed.nec")
add_point(  0.609600,   0.304800,   0.000000)
params = {
    mutation = 0.0,
    exp_weight = 2,
    algorithm = "EX",
    run_simulator = 0,
    max_gain = 1,
    max_coup = 1,
    min_coup = 1,
    auto_seed = 1, 
}
