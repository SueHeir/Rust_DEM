mod calculations;
pub(crate) mod grid;
mod print;
use core::f64::consts::PI;
use std::collections::HashMap;

use nalgebra::{Matrix3, Vector3};

use crate::domain;
use crate::sphere;

pub fn handle_commands(
    command_stack: Vec<String>,
    mut d_data: domain::DomainData,
    mut p_data: sphere::ParticleData,
    mut bond_data: HashMap<(usize, usize), sphere::BondData>,
) {
    // generate_material_map(&mut p_data);

    for command in command_stack {
        let results: Vec<&str> = command.split_whitespace().collect();

        match &results[0][0..3] {
            // "REL" => relax(&mut d_data, &mut p_data),
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
                    &mut bond_data,
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

fn cycle(
    d_data: &mut domain::DomainData,
    p_data: &mut sphere::ParticleData,
    bond_data: &mut HashMap<(usize, usize), sphere::BondData>,
    total_cycles: i32,
    update_rate: i32,
    _clear_rate: i32,
) {
    let dt = calculate_delta_time(p_data);

    let mut change = 1.0;

    for cycle_count in 0..total_cycles {
        //Update velocity and position based on forces
        // grid::euler_integration(p_data, dt);

        grid::inital_integrate(p_data, d_data, dt);

        //Boundary Conditions
        // grid::lees_edwards_boundaries(d_data, p_data, dt, ledisplace);
        if cycle_count % 20 == 0 {
            change *= -1.0;
        }
        //Resets if a particle is in collision, and resets forces to zero
        for i in 0..p_data.radius.len() {
            p_data.is_collision[i] = false;
            p_data.force[i] = Vector3::new(0.0, 0.0, 0.0);
            p_data.torque[i] = Vector3::new(0.0, 0.0, 0.0);
        }

        //Brute Force Collision Detection, this Updates the forces on each particle
        grid::simp_collisions(p_data, bond_data, dt);
        grid::simp_bonds(p_data, bond_data, dt);
        grid::simp_wall(p_data, dt);
        bond_data.retain(|key, value| !value.broken);

        // grid::update(d_data, p_data);
        // grid::collisions(d_data, p_data, &mut f_data, dt, ledisplace);

        grid::final_integrate(p_data, d_data, dt);
        if cycle_count % update_rate == 0 {
            print::print_positions(p_data);
            for i in 0..p_data.radius.len() {
                // println!("Omega {:?}", p_data.omega[i]);
                println!(
                    "Position {:?} {:?} {:?}",
                    p_data.position[i][0], p_data.position[i][1], p_data.position[i][2]
                );
            }
        }
    }
}

fn calculate_delta_time(p_data: &sphere::ParticleData) -> f64 {
    //Checks each particles Size for the smallest delta time the simulation should use
    let mut dt: f64 = 0.001;
    for i in 0..p_data.radius.len() {
        let g = p_data.youngs_mod[i] / (2.0 * (1.0 + p_data.poisson_ratio[i]));
        let alpha = 0.1631 * p_data.poisson_ratio[i] + 0.876605;
        let delta = PI * p_data.radius[i] / alpha * (p_data.density[i] / g).sqrt();
        dt = delta.min(dt);
    }

    println!("Useing {} for delta time", dt * 0.01);
    //Fractional Factor set to 0.5 here,
    return dt * 0.01;
}
