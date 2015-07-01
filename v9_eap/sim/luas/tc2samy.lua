load_platform("input/S6 - 4ft square plate w box and sides sloped front fixed.nec")
add_antenna("input/A1 - 75mono fixed.nec")
add_point( -0.508000,   0.508000,   0.000000)
add_antenna("input/A2 - 150mono fixed.nec")
add_point(  0.304800,   0.508000,   0.000000)
add_antenna("input/A5 - 300dipole fixed.nec")
add_point( -0.101600,   0.101600,   0.203200)
params = {
    mutation = 0.0,
    exp_weight = 2,
    algorithm = "EX",
    run_simulator = 1,
    max_gain = 1,
    max_coup = 1,
    min_coup = 1,
    auto_seed = 1, 
}
