use nalgebra::Vector2;

use crate::domain;
use crate::sphere;

pub fn relax_boundaries_box(d_data: &mut domain::DomainData, p_data: &mut sphere::ParticleData) {
    for i in 0..p_data.radius.len() {
        if p_data.position[i][1] > d_data.domain[1] {
            // p_data.position[i][1] = d_data.domain[1] - p_data.radius[i]
            p_data.position[i][1] -= d_data.domain[1];
        }
        // if particle is less than domain move to end of domain
        // Also apply velocity change to particle for shearing
        else if p_data.position[i][1] <= 0.0 {
            // p_data.position[i][1] = p_data.radius[i]
            p_data.position[i][1] += d_data.domain[1];
        }
        // std::cout << distb << std::endl;//X boundary condition
        // if particles is greater than domain move to beginning of domain
        if p_data.position[i][0] > d_data.domain[0] {
            // p_data.position[i][0] = d_data.domain[0] - p_data.radius[i]
            p_data.position[i][0] -= d_data.domain[0];
        }
        // if particle is less than domain move to end of domain
        else if p_data.position[i][0] <= 0.0 {
            // p_data.position[i][0] = p_data.radius[i]
            p_data.position[i][0] += d_data.domain[0];
        }

        // Z boundary condition
        // if particles is greater than domain move to beginning of domain
        if p_data.position[i][2] > d_data.domain[2] {
            // p_data.position[i][2] = d_data.domain[2] - p_data.radius[i]
            p_data.position[i][2] -= d_data.domain[2];
        }
        // if particle is less than domain move to end of domain
        else if p_data.position[i][2] <= 0.0 {
            // p_data.position[i][2] = p_data.radius[i]
            p_data.position[i][2] += d_data.domain[2];
        }
    }
}

pub fn update(d_data: &mut domain::DomainData, p_data: &mut sphere::ParticleData) {
    //Clears all boxes of particles in the box
    for i in 0..d_data.g_data.len() {
        for j in 0..d_data.g_data[i].len() {
            for k in 0..d_data.g_data[i][j].len() {
                d_data.g_data[i][j][k].real.clear();
                d_data.g_data[i][j][k].ghost.clear();

                for index in 0..p_data.radius.len() {
                    if d_data.g_data[i][j][k].is_position_in_box(p_data.position[index]) {
                        d_data.g_data[i][j][k].real.push(index.try_into().unwrap())
                    } else if d_data.g_data[i][j][k].is_sphere_aabb_in_radius_enlarged_box(
                        p_data.position[index],
                        p_data.radius[index],
                        p_data.max_radius,
                    ) {
                        d_data.g_data[i][j][k].ghost.push(index.try_into().unwrap())
                    } else if ((i == 0 || i == d_data.g_data.len() - 1)
                        || (j == 0 || j == d_data.g_data[i].len() - 1)
                        || (k == 0 || k == d_data.g_data[i][j].len() - 1))
                        && d_data.g_data[i][j][k].is_periodic_sphere(
                            p_data.position[index],
                            p_data.radius[index],
                            p_data.max_radius,
                            d_data,
                        )
                    {
                        d_data.g_data[i][j][k].ghost.push(index.try_into().unwrap())
                    }
                }
            }
        }
    }
}
pub fn relax(d_data: &mut domain::DomainData, p_data: &mut sphere::ParticleData, relax_rate: f64) {
    for box_i in 0..d_data.collision_boxes[0] {
        for box_j in 0..d_data.collision_boxes[1] {
            for box_k in 0..d_data.collision_boxes[2] {
                for ii in 0..d_data.g_data[box_i as usize][box_j as usize][box_k as usize]
                    .real
                    .len()
                {
                    for jj in ii + 1
                        ..d_data.g_data[box_i as usize][box_j as usize][box_k as usize]
                            .real
                            .len()
                    {
                        let i = d_data.g_data[box_i as usize][box_j as usize][box_k as usize].real
                            [ii] as usize;
                        let j = d_data.g_data[box_i as usize][box_j as usize][box_k as usize].real
                            [jj] as usize;

                        // println!("i {}, j {}", i,j);
                        let delta_position = p_data.position[j] - p_data.position[i];

                        let distance = delta_position.norm();

                        //println!("Distance: {} Radius+radius: {}", distance, p_data.radius[i] + p_data.radius[j]);

                        if distance < p_data.radius[i] + p_data.radius[j] {
                            p_data.is_collision[i] = true;
                            p_data.is_collision[j] = true;

                            let normalized_delta = delta_position / distance;

                            //println!("{:?}",normalized_delta);

                            p_data.velocity[i] -= (relax_rate) * normalized_delta;
                            p_data.velocity[j] += (relax_rate) * normalized_delta;
                        }
                    }
                }

                for ii in 0..d_data.g_data[box_i as usize][box_j as usize][box_k as usize]
                    .real
                    .len()
                {
                    for jj in 0..d_data.g_data[box_i as usize][box_j as usize][box_k as usize]
                        .ghost
                        .len()
                    {
                        let i = d_data.g_data[box_i as usize][box_j as usize][box_k as usize].real
                            [ii] as usize;
                        let j = d_data.g_data[box_i as usize][box_j as usize][box_k as usize].ghost
                            [jj] as usize;

                        // println!("i {}, j {}", i,j);
                        let mut p1 = p_data.position[i];
                        let mut p2 = p_data.position[j];

                        let r1 = p_data.radius[i];
                        let r2 = p_data.radius[j];

                        if p1[1] - r1 + d_data.domain[1] <= p2[1] + r2 {
                            p1[1] += d_data.domain[1];
                        } else if p2[1] - r2 + d_data.domain[1] <= p1[1] + r1 {
                            p2[1] += d_data.domain[1];
                        }
                        if p1[0] - r1 + d_data.domain[0] <= p2[0] + r2 {
                            p1[0] += d_data.domain[0];
                        } else if p2[0] - r2 + d_data.domain[0] <= p1[0] + r1 {
                            p2[0] += d_data.domain[0];
                        }

                        if p1[2] - r1 + d_data.domain[2] <= p2[2] + r2 {
                            p1[2] += d_data.domain[2];
                        } else if p2[2] - r2 + d_data.domain[2] <= p1[2] + r1 {
                            p2[2] += d_data.domain[2];
                        }

                        let delta_position = p2 - p1;
                        // let delta_position = p_data.position[j] - p_data.position[i];

                        let distance = delta_position.norm();

                        //println!("Distance: {} Radius+radius: {}", distance, p_data.radius[i] + p_data.radius[j]);

                        if distance < p_data.radius[i] + p_data.radius[j] {
                            p_data.is_collision[i] = true;
                            p_data.is_collision[j] = true;

                            let normalized_delta = delta_position / distance;

                            //println!("{:?}",normalized_delta);

                            p_data.velocity[i] -= (relax_rate * 0.5) * normalized_delta;
                            p_data.velocity[j] += (relax_rate * 0.5) * normalized_delta;
                        }
                    }
                }
            }
        }
    }
}

pub fn is_relaxed(
    d_data: &domain::DomainData,
    p_data: &sphere::ParticleData,
    radius_percentage: f64,
) -> bool {
    println!("is_Relaxed");
    for box_i in 0..d_data.collision_boxes[0] {
        for box_j in 0..d_data.collision_boxes[1] {
            for box_k in 0..d_data.collision_boxes[2] {
                for ii in 0..d_data.g_data[box_i as usize][box_j as usize][box_k as usize]
                    .real
                    .len()
                {
                    for jj in ii + 1
                        ..d_data.g_data[box_i as usize][box_j as usize][box_k as usize]
                            .real
                            .len()
                    {
                        let i = d_data.g_data[box_i as usize][box_j as usize][box_k as usize].real
                            [ii] as usize;
                        let j = d_data.g_data[box_i as usize][box_j as usize][box_k as usize].real
                            [jj] as usize;
                        let delta_position = p_data.position[j] - p_data.position[i];

                        let distance = delta_position.norm();

                        if distance < radius_percentage * (p_data.radius[i] + p_data.radius[j]) {
                            return false;
                        }
                    }
                }

                for ii in 0..d_data.g_data[box_i as usize][box_j as usize][box_k as usize]
                    .real
                    .len()
                {
                    for jj in 0..d_data.g_data[box_i as usize][box_j as usize][box_k as usize]
                        .ghost
                        .len()
                    {
                        let i = d_data.g_data[box_i as usize][box_j as usize][box_k as usize].real
                            [ii] as usize;
                        let j = d_data.g_data[box_i as usize][box_j as usize][box_k as usize].ghost
                            [jj] as usize;

                        let mut p1 = p_data.position[i];
                        let mut p2 = p_data.position[j];

                        let r1 = p_data.radius[i];
                        let r2 = p_data.radius[j];

                        if p1[1] - r1 + d_data.domain[1] <= p2[1] + r2 {
                            p1[1] += d_data.domain[1];
                        } else if p2[1] - r2 + d_data.domain[1] <= p1[1] + r1 {
                            p2[1] += d_data.domain[1];
                        }
                        if p1[0] - r1 + d_data.domain[0] <= p2[0] + r2 {
                            p1[0] += d_data.domain[0];
                        } else if p2[0] - r2 + d_data.domain[0] <= p1[0] + r1 {
                            p2[0] += d_data.domain[0];
                        }

                        if p1[2] - r1 + d_data.domain[2] <= p2[2] + r2 {
                            p1[2] += d_data.domain[2];
                        } else if p2[2] - r2 + d_data.domain[2] <= p1[2] + r1 {
                            p2[2] += d_data.domain[2];
                        }

                        let delta_position = p1 - p2;

                        let distance = delta_position.norm();

                        if distance < radius_percentage * (p_data.radius[i] + p_data.radius[j]) {
                            return false;
                        }
                    }
                }
            }
        }
    }
    return true;
}

pub fn _simp_collisions(
    d_data: &domain::DomainData,
    p_data: &mut sphere::ParticleData,
    f_data: &mut sphere::ForceData,
    _dt: f64,
    ledisplace: f64,
) {
    for i in 0..p_data.radius.len() {
        for j in i + 1..p_data.radius.len() {
            let mut p1 = p_data.position[i];
            let mut p2 = p_data.position[j];
            let mut v1 = p_data.velocity[i];
            let mut v2 = p_data.velocity[j];

            let r1 = p_data.radius[i];
            let r2 = p_data.radius[j];
            if p1[1] - r1 + d_data.domain[1] <= p2[1] + r2 {
                p1[1] += d_data.domain[1];
                v1[0] += d_data.lees_edwards_boundary * d_data.domain[1];
                p1[0] += ledisplace;
                if p1[0] > d_data.domain[0] {
                    p1[0] -= d_data.domain[0];
                }
            } else if p2[1] - r2 + d_data.domain[1] <= p1[1] + r1 {
                p2[1] += d_data.domain[1];
                v2[0] += d_data.lees_edwards_boundary * d_data.domain[1];
                p2[0] += ledisplace;
                if p2[0] > d_data.domain[0] {
                    p2[0] -= d_data.domain[0];
                }
            }
            if p1[0] - r1 + d_data.domain[0] <= p2[0] + r2 {
                p1[0] += d_data.domain[0];
            } else if p2[0] - r2 + d_data.domain[0] <= p1[0] + r1 {
                p2[0] += d_data.domain[0];
            }

            if p1[2] - r1 + d_data.domain[2] <= p2[2] + r2 {
                p1[2] += d_data.domain[2];
            } else if p2[2] - r2 + d_data.domain[2] <= p1[2] + r1 {
                p2[2] += d_data.domain[2];
            }

            let delta_position = p2 - p1;

            let distance = delta_position.norm();

            if distance < p_data.radius[i] + p_data.radius[j] {
                p_data.is_collision[i] = true;
                p_data.is_collision[j] = true;

                let normalized_delta = delta_position / distance;

                let distance_delta = (p_data.radius[i] + p_data.radius[j]) - distance;

                // let eff = p_data.sphere_material_map.get(&s).unwrap();
                // let effective_radius = eff.eff_radius;
                // let effective_youngs = eff.eff_youngs_mod;

                let effective_radius = 1.0 / (1.0 / p_data.radius[i] + 1.0 / p_data.radius[j]);

                let effective_youngs = 1.0
                    / ((1.0 - p_data.poisson_ratio[i] * p_data.poisson_ratio[i])
                        / p_data.youngs_mod[i]
                        + (1.0 - p_data.poisson_ratio[j] * p_data.poisson_ratio[j])
                            / p_data.youngs_mod[j]);
                let normal_force = 4.0 / 3.0
                    * effective_youngs
                    * effective_radius.sqrt()
                    * distance_delta.powf(3.0 / 2.0);

                let contact_stiffness =
                    2.0 * effective_youngs * (effective_radius * distance_delta).sqrt();
                let reduced_mass =
                    p_data.mass[i] * p_data.mass[j] / (p_data.mass[i] + p_data.mass[j]);

                let delta_veloctiy = v2 - v1;
                let f_dot = normalized_delta.dot(&delta_veloctiy);
                let _v_r = f_dot * normalized_delta;

                // std::cout << beta <<F_dot << " " << reduced_mass << " "<< contact_stiffness << "\n";

                let dissipation_force = 2.0
                    * 0.9128709292
                    * p_data.beta
                    * (contact_stiffness * reduced_mass).sqrt()
                    * f_dot;

                // if(distance_delta > 0.00001)
                // std::cout << dissipation_force << " " << normal_force << " "<< distance_delta <<"\n";

                p_data.force[i] -= (normal_force - dissipation_force) * normalized_delta;
                p_data.force[j] += (normal_force - dissipation_force) * normalized_delta;

                f_data
                    .force
                    .push((normal_force - dissipation_force) * normalized_delta);
                f_data.particle_indexes.push(Vector2::new(i, j));
                f_data.del.push(delta_position);

                // if(distance_delta > 10e-7){
                //     std::cout << distance_delta << " " << switched << "\n";
                // }

                // if(switched){
                //     if(delta_veloctiy.norm() > 0.1){
                //         std::cout << delta_veloctiy[0] << " " << delta_veloctiy[1] << " " << delta_veloctiy[2] << "\n";
                //         std::cout << delta_veloctiy.norm() << "\n";
                //         std::cout << error << "\n";

                //     }

                // }
            }
        }
    }
}

pub fn collisions(
    d_data: &domain::DomainData,
    p_data: &mut sphere::ParticleData,
    f_data: &mut sphere::ForceData,
    _dt: f64,
    ledisplace: f64,
) {
    for box_i in 0..d_data.collision_boxes[0] {
        for box_j in 0..d_data.collision_boxes[1] {
            for box_k in 0..d_data.collision_boxes[2] {
                for ii in 0..d_data.g_data[box_i as usize][box_j as usize][box_k as usize]
                    .real
                    .len()
                {
                    for jj in ii + 1
                        ..d_data.g_data[box_i as usize][box_j as usize][box_k as usize]
                            .real
                            .len()
                    {
                        let i = d_data.g_data[box_i as usize][box_j as usize][box_k as usize].real
                            [ii] as usize;
                        let j = d_data.g_data[box_i as usize][box_j as usize][box_k as usize].real
                            [jj] as usize;

                        // println!("i {}, j {}", i,j);
                        let delta_position = p_data.position[j] - p_data.position[i];

                        let distance = delta_position.norm();

                        if distance < p_data.radius[i] + p_data.radius[j] {
                            p_data.is_collision[i] = true;
                            p_data.is_collision[j] = true;

                            let normalized_delta = delta_position / distance;

                            let distance_delta = (p_data.radius[i] + p_data.radius[j]) - distance;

                            // let eff = p_data.sphere_material_map.get(&s).unwrap();
                            // let effective_radius = eff.eff_radius;
                            // let effective_youngs = eff.eff_youngs_mod;

                            let effective_radius =
                                1.0 / (1.0 / p_data.radius[i] + 1.0 / p_data.radius[j]);

                            let effective_youngs = 1.0
                                / ((1.0 - p_data.poisson_ratio[i] * p_data.poisson_ratio[i])
                                    / p_data.youngs_mod[i]
                                    + (1.0 - p_data.poisson_ratio[j] * p_data.poisson_ratio[j])
                                        / p_data.youngs_mod[j]);

                            let contact_stiffness =
                                2.0 * effective_youngs * (effective_radius * distance_delta).sqrt();

                            let normal_force = 2.0 / 3.0 * distance_delta * contact_stiffness;
                            let reduced_mass =
                                p_data.mass[i] * p_data.mass[j] / (p_data.mass[i] + p_data.mass[j]);

                            let delta_veloctiy = p_data.velocity[j] - p_data.velocity[i];
                            let f_dot = normalized_delta.dot(&delta_veloctiy);
                            let v_r_n = f_dot * normalized_delta;

                            let dissipation_force = 2.0
                                * 0.91287092917
                                * p_data.beta
                                * (contact_stiffness * reduced_mass).sqrt()
                                * v_r_n.norm()
                                * v_r_n.dot(&normalized_delta).signum();

                            // println!("{} {} {}", distance_delta, normal_force, dissipation_force);
                            p_data.force[i] -=
                                (normal_force - dissipation_force) * normalized_delta;
                            p_data.force[j] +=
                                (normal_force - dissipation_force) * normalized_delta;

                            let force_length_matrix = ((normal_force - dissipation_force)
                                * normalized_delta)
                                * delta_position.transpose();
                            f_data.forcedata.push(force_length_matrix);
                        }
                    }
                }

                for ii in 0..d_data.g_data[box_i as usize][box_j as usize][box_k as usize]
                    .real
                    .len()
                {
                    for jj in 0..d_data.g_data[box_i as usize][box_j as usize][box_k as usize]
                        .ghost
                        .len()
                    {
                        let i = d_data.g_data[box_i as usize][box_j as usize][box_k as usize].real
                            [ii] as usize;
                        let j = d_data.g_data[box_i as usize][box_j as usize][box_k as usize].ghost
                            [jj] as usize;
                        let mut p1 = p_data.position[i];
                        let mut p2 = p_data.position[j];
                        let mut v1 = p_data.velocity[i];
                        let mut v2 = p_data.velocity[j];

                        let r1 = p_data.radius[i];
                        let r2 = p_data.radius[j];
                        if p1[1] - r1 + d_data.domain[1] <= p2[1] + r2 {
                            p1[1] += d_data.domain[1];
                            v1[0] += d_data.lees_edwards_boundary * d_data.domain[1];
                            p1[0] += ledisplace;
                            if p1[0] > d_data.domain[0] {
                                p1[0] -= d_data.domain[0];
                            }
                        } else if p2[1] - r2 + d_data.domain[1] <= p1[1] + r1 {
                            p2[1] += d_data.domain[1];
                            v2[0] += d_data.lees_edwards_boundary * d_data.domain[1];
                            p2[0] += ledisplace;
                            if p2[0] > d_data.domain[0] {
                                p2[0] -= d_data.domain[0];
                            }
                        }
                        if p1[0] - r1 + d_data.domain[0] <= p2[0] + r2 {
                            p1[0] += d_data.domain[0];
                        } else if p2[0] - r2 + d_data.domain[0] <= p1[0] + r1 {
                            p2[0] += d_data.domain[0];
                        }

                        if p1[2] - r1 + d_data.domain[2] <= p2[2] + r2 {
                            p1[2] += d_data.domain[2];
                        } else if p2[2] - r2 + d_data.domain[2] <= p1[2] + r1 {
                            p2[2] += d_data.domain[2];
                        }

                        let delta_position = p2 - p1;

                        let distance = delta_position.norm();

                        if distance < p_data.radius[i] + p_data.radius[j] {
                            p_data.is_collision[i] = true;
                            p_data.is_collision[j] = true;

                            let normalized_delta = delta_position / distance;

                            let distance_delta = (p_data.radius[i] + p_data.radius[j]) - distance;

                            let effective_radius =
                                1.0 / (1.0 / p_data.radius[i] + 1.0 / p_data.radius[j]);

                            let effective_youngs = 1.0
                                / ((1.0 - p_data.poisson_ratio[i] * p_data.poisson_ratio[i])
                                    / p_data.youngs_mod[i]
                                    + (1.0 - p_data.poisson_ratio[j] * p_data.poisson_ratio[j])
                                        / p_data.youngs_mod[j]);

                            let contact_stiffness =
                                2.0 * effective_youngs * (effective_radius * distance_delta).sqrt();

                            let normal_force = 2.0 / 3.0 * distance_delta * contact_stiffness;
                            let reduced_mass =
                                p_data.mass[i] * p_data.mass[j] / (p_data.mass[i] + p_data.mass[j]);

                            let delta_veloctiy = v2 - v1;
                            let f_dot = normalized_delta.dot(&delta_veloctiy);
                            let v_r_n = f_dot * normalized_delta;

                            let dissipation_force = 2.0
                                * 0.9128709292
                                * p_data.beta
                                * (contact_stiffness * reduced_mass).sqrt()
                                * v_r_n.norm()
                                * v_r_n.dot(&normalized_delta).signum();
                            // println!("{} {} {}", distance_delta, normal_force, dissipation_force);
                            p_data.force[i] -=
                                (normal_force - dissipation_force) * normalized_delta * 0.5;
                            p_data.force[j] +=
                                (normal_force - dissipation_force) * normalized_delta * 0.5;

                            let force_length_matrix =
                                ((normal_force - dissipation_force) * normalized_delta * 0.5)
                                    * delta_position.transpose();
                            f_data.forcedata.push(force_length_matrix);
                        }
                    }
                }
            }
        }
    }
}

pub fn _euler_integration(p_data: &mut sphere::ParticleData, dt: f64) {
    for i in 0..p_data.radius.len() {
        p_data.velocity[i] += dt * p_data.force[i] / p_data.mass[i];
        p_data.position[i] += p_data.velocity[i] * dt;
    }
}

pub fn inital_integrate(p_data: &mut sphere::ParticleData, dt: f64) {
    for i in 0..p_data.radius.len() {
        p_data.velocity[i] += 0.5 * dt * p_data.force[i] / p_data.mass[i];
        p_data.position[i] += p_data.velocity[i] * dt;
    }
}

pub fn final_integrate(p_data: &mut sphere::ParticleData, dt: f64) {
    for i in 0..p_data.radius.len() {
        p_data.velocity[i] += 0.5 * dt * p_data.force[i] / p_data.mass[i];
    }
}

pub fn lees_edwards_boundaries(
    d_data: &domain::DomainData,
    p_data: &mut sphere::ParticleData,
    _dt: f64,
    ledisplace: f64,
) {
    for i in 0..p_data.radius.len() {
        // Y boundary condition
        // if particles is greater than domain move to beginning of domain
        // Also apply velocity change to particle for shearing
        if p_data.position[i][1] > d_data.domain[1] {
            p_data.position[i][1] -= d_data.domain[1];
            p_data.velocity[i][0] -= d_data.lees_edwards_boundary * d_data.domain[1];
            p_data.position[i][0] -= ledisplace;
            if p_data.position[i][0] <= 0.0 {
                p_data.position[i][0] += d_data.domain[0];
            }
        }
        // if particle is less than domain move to end of domain
        // Also apply velocity change to particle for shearing
        else if p_data.position[i][1] <= 0.0 {
            p_data.position[i][1] += d_data.domain[1];
            p_data.velocity[i][0] += d_data.lees_edwards_boundary * d_data.domain[1];
            p_data.position[i][0] += ledisplace;
            if p_data.position[i][0] > d_data.domain[0] {
                p_data.position[i][0] -= d_data.domain[0];
            }
        }
        // std::cout << distb << std::endl;//X boundary condition
        // if particles is greater than domain move to beginning of domain
        if p_data.position[i][0] > d_data.domain[0] {
            p_data.position[i][0] -= d_data.domain[0];
        }
        // if particle is less than domain move to end of domain
        else if p_data.position[i][0] <= 0.0 {
            p_data.position[i][0] += d_data.domain[0];
        }

        // Z boundary condition
        // if particles is greater than domain move to beginning of domain
        if p_data.position[i][2] > d_data.domain[2] {
            p_data.position[i][2] -= d_data.domain[2];
        }
        // if particle is less than domain move to end of domain
        else if p_data.position[i][2] <= 0.0 {
            p_data.position[i][2] += d_data.domain[2];
        }
    }
}
