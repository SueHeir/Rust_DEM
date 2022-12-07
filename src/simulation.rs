mod calculations;
pub(crate) mod grid;
mod print;
use core::f64::consts::PI;
use std::process::exit;

use nalgebra::{Matrix3, Vector3};
use rand::Rng;

use crate::domain;
use crate::rigid_body;
use crate::simulation::print::print_vtp;

pub fn handle_commands(
    command_stack: Vec<String>,
    mut d_data: domain::DomainData,
    mut p_data: rigid_body::RigidBodiesData,
) {
    // generate_material_map(&mut p_data);

    for command in command_stack {
        let results: Vec<&str> = command.split_whitespace().collect();

        match &results[0][0..3] {
            "REL" => relax(&mut d_data, &mut p_data),
            "CYC" => {
                //If updateRate and clear rate are not set, the default is used
                let mut update_rate = 2500;
                let mut clear_rate = 75000;

                if !results[2].is_empty() {
                    update_rate = results[2].parse::<i32>().unwrap();
                }
                if !results[3].is_empty() {
                    clear_rate = results[3].parse::<i32>().unwrap();
                }

                cycle(
                    &mut d_data,
                    &mut p_data,
                    results[1].parse::<i32>().unwrap(),
                    update_rate,
                    clear_rate,
                );
                println!(
                    "Cycle completed {} steps",
                    results[1].parse::<i32>().unwrap()
                );
            }
            _ => {}
        }
    }
}

// fn generate_material_map(p_data: &mut rigid_body::RigidBodiesData) {
//     for material1 in &p_data.materials {
//         for material2 in &p_data.materials {
//             let s = format!("{} {}", material1.id, material2.id);

//             let effective = sphere::EffectiveMaterialPreCalc {
//                 eff_radius: 1.0 / (1.0 / material1.radius + 1.0 / material2.radius),

//                 eff_youngs_mod: 1.0
//                     / ((1.0 - material1.poisson_ratio * material1.poisson_ratio)
//                         / material1.youngs_mod
//                         + (1.0 - material2.poisson_ratio * material2.poisson_ratio)
//                             / material2.youngs_mod),
//             };

//             p_data.sphere_material_map.insert(s, effective);
//         }
//     }
// }

fn relax(d_data: &mut domain::DomainData, p_data: &mut rigid_body::RigidBodiesData) {
    let dt = 0.00000001;

    print_vtp(p_data, -1);
    let mut still_relaxing = true;
    let mut count = 0;
    while still_relaxing {
        grid::inital_integrate(p_data, dt);

        for i in 0..p_data.body_data.len() {
            p_data.velocity[i] = p_data.velocity[i] * 0.2;
            p_data.angular_moment[i] = p_data.angular_moment[i] * 0.2;
            p_data.is_collision[i] = false;

            p_data.force[i] = Vector3::new(0.0, 0.0, 0.0);
            p_data.torque[i] = Vector3::new(0.0, 0.0, 0.0);
            for j in 0..p_data.body_data[i].wrap.len() {
                p_data.body_data[i].wrap[j] = Vector3::zeros();
            }
        }

        // p_data.update_all_body_spheres_positions();

        grid::relax_boundaries_box(d_data, p_data);
        grid::update(d_data, p_data);

        grid::relax(d_data, p_data, 0.00000002);

        grid::final_integrate(p_data, dt);

        let mut check_is_relaxed = true;
        for i in 0..p_data.body_data.len() {
            if p_data.is_collision[i] {
                check_is_relaxed = false;
            }
        }
        still_relaxing = !check_is_relaxed;

        if count % 2000 == 0 {
            println!("Still Relaxing");
            print_vtp(p_data, count);
        }
        count += 1;
    }

    //Set lees Edwards boundary condition velocity here, Also calculate volume fraction
    let mut volume = 0.0;
    let mut rng = rand::thread_rng();
    for i in 0..p_data.body_data.len() {
        volume += p_data.body_data[i].volume;

        let x_velocity: f64 =
            (p_data.position[i][1] - d_data.domain[1] * 0.5) * d_data.lees_edwards_boundary;

        // let y: f64 = rng.gen();
        // let z: f64 = rng.gen();
        // println!("Init Velocity: {}", x_velocity);
        p_data.velocity[i] = Vector3::new(x_velocity, 0.0, 0.0);
        p_data.angular_moment[i] = Vector3::zeros();
        // let x: f64 = rng.gen::<f64>();
        // let y: f64 = rng.gen::<f64>();
        // let z: f64 = rng.gen::<f64>();
        // p_data.velocity[i] = Vector3::new(
        //     (x - 0.5) * d_data.lees_edwards_boundary * d_data.domain[0] * 0.6,
        //     (y - 0.5) * d_data.lees_edwards_boundary * d_data.domain[0] * 0.6,
        //     (z - 0.5) * d_data.lees_edwards_boundary * d_data.domain[0] * 0.6,
        // );
    }

    let solid_fraction = volume / (d_data.domain[0] * d_data.domain[1] * d_data.domain[2]);

    p_data.volume_fraction = solid_fraction;

    println!("Volume Fraction: {}", solid_fraction);
    println!("Finished Relaxing");
}

fn cycle(
    d_data: &mut domain::DomainData,
    p_data: &mut rigid_body::RigidBodiesData,
    total_cycles: i32,
    update_rate: i32,
    clear_rate: i32,
) {
    let dt = calculate_delta_time(p_data);

    let mut f_data = rigid_body::ForceData {
        particle_indexes: Vec::new(),
        force: Vec::new(),
        del: Vec::new(),
        forcedata: Vec::new(),
    };

    let mut kinetic_tensor = Matrix3::zeros();
    let mut collision_tensor = Matrix3::zeros();
    let mut average_reset_count = 0;

    let mut ledisplace = 0.0;

    let mut change_direction = -1.0;

    for cycle_count in 0..total_cycles {
        //Update velocity and position based on forces
        // grid::euler_integration(p_data, dt);

        grid::inital_integrate(p_data, dt);

        //Resets if a particle is in collision, and resets forces to zero

        // if cycle_count % 5000 == 0 {
        //     change_direction = -change_direction;
        // }
        for i in 0..p_data.body_data.len() {
            p_data.is_collision[i] = false;
            p_data.force[i] = Vector3::new(0.0, 0.0, 0.0);
            // p_data.torque[i] = Vector3::new(
            //     0.000000001 * change_direction,
            //     0.000000001 * change_direction,
            //     0.000000001 * change_direction,
            // );
            p_data.torque[i] = Vector3::zeros();

            // if p_data.omega[i].norm() < 20.0 {
            //     println!("{}", p_data.omega[i].norm());
            // }

            for j in 0..p_data.body_data[i].wrap.len() {
                p_data.body_data[i].wrap[j] = Vector3::zeros();
            }
        }
        f_data.forcedata.clear();

        //Boundary Conditions
        grid::lees_edwards_boundaries(d_data, p_data, dt, ledisplace);
        // grid::relax_boundaries_box(d_data, p_data);
        // grid::simple_boundary(d_data, p_data);

        // grid::update(d_data, p_data);
        grid::simple_collisions(d_data, p_data, &mut f_data, dt, ledisplace);

        grid::final_integrate(p_data, dt);

        //Brute Force Collision Detection, this Updates the forces on each particle
        // grid::_simp_collisions(d_data, p_data, &mut f_data, dt, ledisplace);

        //calculates the kinetic stress tensor
        kinetic_tensor =
            calculations::calc_kinetic_tensor(p_data, d_data, kinetic_tensor, average_reset_count);
        collision_tensor = calculations::calc_collision_tensor(
            &f_data,
            d_data,
            collision_tensor,
            average_reset_count,
        );

        average_reset_count += 1;

        //Print statments to terminal and prints the VTP, and Stress data
        if cycle_count % update_rate == 0 {
            println!("CYC {} completed", cycle_count);
            // println!("{}", p_data.velocity[0][0]);
            println!("{:?}", (kinetic_tensor + collision_tensor));

            print::print_vtp(p_data, cycle_count);
        }

        //Resets the averaging of the kinetic tensor
        if cycle_count % clear_rate == 0 {
            print::print_stress(kinetic_tensor, collision_tensor, cycle_count);

            kinetic_tensor = Matrix3::zeros();
            collision_tensor = Matrix3::zeros();
            average_reset_count = 0;
        }

        ledisplace += dt * d_data.lees_edwards_boundary * d_data.domain[1];
        ledisplace -= (ledisplace / d_data.domain[0]).floor() * d_data.domain[0];
    }
}

fn calculate_delta_time(p_data: &rigid_body::RigidBodiesData) -> f64 {
    //Checks each particles Size for the smallest delta time the simulation should use
    let mut dt: f64 = 0.001;
    for i in 0..p_data.body_data.len() {
        let g = p_data.youngs_mod[i] / (2.0 * (1.0 + p_data.poisson_ratio[i]));

        // println!("g {}", g);
        let alpha = 0.1631 * p_data.poisson_ratio[i] + 0.876605;

        // println!("alpha {}", alpha);
        for j in 0..p_data.body_data[i].radius.len() {
            let delta = PI * p_data.body_data[i].radius[j] / alpha
                * (p_data.body_data[i].density / g).sqrt();
            // println!("delta {}", delta);
            dt = delta.min(dt);
        }
    }

    println!("Useing {} for delta time", dt * 0.5);
    //Fractional Factor set to 0.5 here,
    return dt * 0.5;
}
