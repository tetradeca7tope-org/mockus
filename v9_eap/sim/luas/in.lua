load_platform("S6 - 4ft square plate w box and sides sloped front fixed.nec")
add_antenna("input/A1 - 75mono fixed.nec")
add_point( 0.5844147, 0.9189849, 0.3114814 )
add_antenna("input/A2 - 150mono fixed.nec")
add_point( -0.9285766, 0.6982586, 0.8679865 )
add_antenna("input/A3 - 150dipole fixed.nec")
add_point( 0.3574703, 0.5154803, 0.4862649 )
params = {mutation = 0.0,exp_weight = 2,algorithm = "EX",run_simulator = 1,max_gain = 1,max_coup = 1,min_coup = 1,auto_seed = 1}