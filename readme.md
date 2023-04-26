# DEM Code in Rust

A discrete element method code written in Rust. This code only runs Shear cells with a Lees-Edwards boundary condition. Results are matching granular kinetic theory. 

Hertz contact model for single spheres, no friction

Boundary conditions are period in x and z, and a Lees-Edwards boundary condition in y

Example Input File
```
START 0.002 0.002 0.001 2 2 2
DAMPING 0.95
LEB 100.0
GRAVITY 0.0 0.0 0.0
MATERIAL 1 6e-5 2500 8.7e9 0.30
RGP 444 1
RELAX
CYC 10000000 400 5000
```

```
What each command inputs are
START x_domain y_domain z_domain x_axis_collision_box y_axis_collision_box z_axis_collision_box
DAMPING restitution_coefficient
GRAVITY g_x g_y g_z
MATERIAL material_id radius density younge_mod poissions_ratio
RGP number_of_particles_to_generate material_id
RELAX (needed after RGP to remove overlaps)
CYC number_of_cycles vtp_print_rate stress_averaging_and_print_rate
```

