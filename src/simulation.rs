mod calculations;
pub(crate) mod grid;
mod print;
use core::f64::consts::PI;

use nalgebra::{Matrix3, Vector3};

use crate::domain;
use crate::sphere;

pub fn handle_commands(
    command_stack: Vec<String>,
    mut d_data: domain::DomainData,
    mut p_data: sphere::ParticleData,
) {
    generate_material_map(&mut p_data);

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

fn generate_material_map(p_data: &mut sphere::ParticleData) {
    for material1 in &p_data.materials {
        for material2 in &p_data.materials {
            let s = format!("{} {}", material1.id, material2.id);

            let effective = sphere::EffectiveMaterialPreCalc {
                eff_radius: 1.0 / (1.0 / material1.radius + 1.0 / material2.radius),

                eff_youngs_mod: 1.0
                    / ((1.0 - material1.poisson_ratio * material1.poisson_ratio)
                        / material1.youngs_mod
                        + (1.0 - material2.poisson_ratio * material2.poisson_ratio)
                            / material2.youngs_mod),
            };

            p_data.sphere_material_map.insert(s, effective);
        }
    }
}

fn relax(d_data: &mut domain::DomainData, p_data: &mut sphere::ParticleData) {
    // print_vtp(p_data, -1);
    let mut still_relaxing = true;
    let mut count = 0;
    while still_relaxing {
        for i in 0..p_data.radius.len() {
            // println!("{:?}", p_data.position[i]);
            p_data.position[i] += p_data.velocity[i] * 0.00013;
            // println!("{:?}", p_data.velocity[i]);
            p_data.velocity[i] = p_data.velocity[i] * 0.5;
            p_data.is_collision[i] = false;
        }
        grid::update(d_data, p_data);
        grid::relax_boundaries_box(d_data, p_data);

        grid::relax(d_data, p_data, 0.02);
        still_relaxing = !grid::is_relaxed(d_data, p_data, 1.0);
        if count % 2000 == 0 {
            println!("Still Relaxing");
        }
        count += 1;
        // print_vtp(p_data, count);
    }

    //Set lees Edwards boundary condition velocity here, Also calculate volume fraction
    let mut volume = 0.0;
    for i in 0..p_data.radius.len() {
        volume += 4.0 / 3.0 * PI * p_data.radius[i].powi(3);

        let x_velocity: f64 =
            (p_data.position[i][1] - d_data.domain[1] * 0.5) * d_data.lees_edwards_boundary;

        // let y: f64 = rng.gen();
        // let z: f64 = rng.gen();
        p_data.velocity[i] = Vector3::new(x_velocity, 0.0, 0.0);

        //p_data.velocity[i] = Eigen::Vector3d((randf()-0.5)*d_data.lees_edwards_boundary*d_data.domain(0),(randf()-0.5)*d_data.lees_edwards_boundary*d_data.domain(0),(randf()-0.5)*d_data.lees_edwards_boundary*d_data.domain(0));
    }

    let solid_fraction = volume / (d_data.domain[0] * d_data.domain[1] * d_data.domain[2]);

    p_data.volume_fraction = solid_fraction;

    println!("Volume Fraction: {}", solid_fraction);
    println!("Finished Relaxing");
}

fn cycle(
    d_data: &mut domain::DomainData,
    p_data: &mut sphere::ParticleData,
    total_cycles: i32,
    update_rate: i32,
    clear_rate: i32,
) {
    let dt = calculate_delta_time(p_data);

    let mut f_data = sphere::ForceData {
        particle_indexes: Vec::new(),
        force: Vec::new(),
        del: Vec::new(),
        forcedata: Vec::new(),
    };

    let mut kinetic_tensor = Matrix3::zeros();
    let mut collision_tensor = Matrix3::zeros();
    let mut average_reset_count = 0;

    let mut ledisplace = 0.0;

    for cycle_count in 0..total_cycles {
        //Update velocity and position based on forces
        // grid::euler_integration(p_data, dt);

        grid::inital_integrate(p_data, dt);

        //Boundary Conditions
        grid::lees_edwards_boundaries(d_data, p_data, dt, ledisplace);

        //Resets if a particle is in collision, and resets forces to zero
        for i in 0..p_data.radius.len() {
            p_data.is_collision[i] = false;
            p_data.force[i] = Vector3::new(0.0, 0.0, 0.0);
        }

        f_data.forcedata.clear();

        grid::update(d_data, p_data);
        grid::collisions(d_data, p_data, &mut f_data, dt, ledisplace);

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

        //Print statments to terminal and prints the VTP, and Stress data to fikkk,k,mles
        if cycle_count % update_rate == 0 {
            // println!("CYC {} completed", cycle_count);
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

fn calculate_delta_time(p_data: &sphere::ParticleData) -> f64 {
    //Checks each particles Size for the smallest delta time the simulation should use
    let mut dt: f64 = 0.001;
    for i in 0..p_data.radius.len() {
        let g = p_data.youngs_mod[i] / (2.0 * (1.0 + p_data.poisson_ratio[i]));
        let alpha = 0.1631 * p_data.poisson_ratio[i] + 0.876605;
        let delta = PI * p_data.radius[i] / alpha * (p_data.density[i] / g).sqrt();
        dt = delta.min(dt);
    }

    println!("Useing {} for delta time", dt * 0.5);
    //Fractional Factor set to 0.5 here,
    return dt * 0.5;
}
