DEM Code in Rust

Hertz contact model for single spheres, no friction

Boundary conditions are period in x and z, and a Lees-Edwards boundary condition in y

Example Input File:

START 0.002 0.002 0.001 1 1 1 
DAMPING 0.95
LEB 100.0
GRAVITY 0.0 0.0 0.0
MATERIAL 1 2500 8.7e9 0.30
RIG rigid_body_data.txt
RGB 222 1
RELAX
CYC 10000000 400 5000

Example Rigid Body Data File: (i.e. rigid_body_data.txt)
0.0 0.0 0.0 6e-5
12e-5 0.0 0.0 6e-5



What each command inputs are
START x_domain y_domain z_domain x_axis_collision_box y_axis_collision_box z_axis_collision_box
DAMPING restitution_coefficient
GRAVITY g_x g_y g_z
MATERIAL material_id density younge_mod poissions_ratio
RIG filename
RGB number_of_particles_to_generate material_id
RELAX (needed after RGP to remove overlaps)
CYC number_of_cycles vtp_print_rate stress_averaging_and_print_rate



TODO!
fix collision boxes (they aren't being used right now)



