DEM Code in Rust

Hertz contact model for single spheres, no friction

No Boundary conditions 

Setup to test Bonded sphere model with 3 spheres and 3 bonds in a triangle

Example Input File

START 0.005 0.005 0.005 1 1 1
DAMPING 0.40
GRAVITY 0.0 0.0 0.0
MATERIAL 1 0.025E-04 7.295e3 4.3e10 0.33
FORCE 1
CYC 100000000 100 0


What each command inputs are
START x_domain y_domain z_domain x_axis_collision_box y_axis_collision_box z_axis_collision_box
DAMPING restitution_coefficient
GRAVITY g_x g_y g_z
MATERIAL material_id radius density younge_mod poissions_ratio
FORCE material_id (sets up the tringle simulation)
CYC number_of_cycles vtp_print_rate(update_rate variable in code if you need to print out every so many time steps) stress_averaging_and_print_rate


