use std::process;

use nalgebra::Vector2;
use nalgebra::Vector3;

use crate::domain;
use crate::sphere;
use crate::sphere::ContactStruct;
use crate::sphere::EffectiveMaterialPreCalc;
use noisy_float::prelude::*;

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

pub fn simp_collisions(
    d_data: &domain::DomainData,
    p_data: &mut sphere::ParticleData,
    c_data: &mut sphere::ContactData,
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
                // v1[0] += d_data.lees_edwards_boundary * d_data.domain[1];
                // p1[0] += ledisplace;
                // if p1[0] > d_data.domain[0] {
                //     p1[0] -= d_data.domain[0];
                // }
            } else if p2[1] - r2 + d_data.domain[1] <= p1[1] + r1 {
                p2[1] += d_data.domain[1];
                // v2[0] += d_data.lees_edwards_boundary * d_data.domain[1];
                // p2[0] += ledisplace;
                // if p2[0] > d_data.domain[0] {
                //     p2[0] -= d_data.domain[0];
                // }
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
                let (c_index, swap_index) = get_collision_struct(c_data, i, j);
                let mut c_i = 0;
                let mut c_j = 1;

                if swap_index < 0 {
                    c_j = 0;
                    c_i = 1;
                }

                p_data.is_collision[i] = true;
                p_data.is_collision[j] = true;

                let normalized_delta_position = delta_position / distance;

                let distance_delta = (p_data.radius[i] + p_data.radius[j]) - distance;
                // let eff = p_data.sphere_material_map.get(&s).unwrap();
                // let effective_radius = eff.eff_radius;
                // let effective_youngs = eff.eff_youngs_mod;
                let effective_radius =
                    p_data.radius[i] * p_data.radius[j] / (p_data.radius[i] + p_data.radius[j]);

                c_data.contact_struct[c_index].contact_radius =
                    (distance_delta * effective_radius).sqrt();

                let effective_youngs = 1.0
                    / ((1.0 - p_data.poisson_ratio[i] * p_data.poisson_ratio[i])
                        / p_data.youngs_mod[i]
                        + (1.0 - p_data.poisson_ratio[j] * p_data.poisson_ratio[j])
                            / p_data.youngs_mod[j]);

                let effective_shear_modulus = 1.0
                    / ((1.0 - p_data.poisson_ratio[i] * p_data.poisson_ratio[i])
                        / p_data.shear_mod[i]
                        + (1.0 - p_data.poisson_ratio[j] * p_data.poisson_ratio[j])
                            / p_data.shear_mod[j]);

                let contact_stiffness =
                    2.0 * effective_youngs * (effective_radius * distance_delta).sqrt();

                c_data.contact_struct[c_index].contact_vector[c_i] = delta_position / 2.0;
                c_data.contact_struct[c_index].contact_vector[c_j] = -delta_position / 2.0;

                let reduced_mass =
                    p_data.mass[i] * p_data.mass[j] / (p_data.mass[i] + p_data.mass[j]);
                let velocity_rel = v2 - v1
                    + p_data.omega[i].cross(&c_data.contact_struct[c_index].contact_vector[c_i])
                    - p_data.omega[j].cross(&c_data.contact_struct[c_index].contact_vector[c_j]);

                let velocity_rel_normal =
                    velocity_rel.dot(&normalized_delta_position) * normalized_delta_position;

                c_data.contact_struct[c_index].previous_force_magnitude[c_i] =
                    c_data.contact_struct[c_index].force_magnitude[c_i];
                c_data.contact_struct[c_index].previous_force_magnitude[c_j] =
                    c_data.contact_struct[c_index].force_magnitude[c_j];

                c_data.contact_struct[c_index].previous_tangential_force[c_i] =
                    c_data.contact_struct[c_index].tangential_force[c_i];
                c_data.contact_struct[c_index].previous_tangential_force[c_j] =
                    c_data.contact_struct[c_index].tangential_force[c_j];

                c_data.contact_struct[c_index].previous_tangential_force_magnitude[c_i] =
                    c_data.contact_struct[c_index].tangential_force_magnitude[c_i];

                let normal_force = -4.0 / 3.0
                    * effective_youngs
                    * (effective_radius).sqrt()
                    * distance_delta.powf(3.0 / 2.0);

                let normal_dissipation_force = (5.0 / 6.0).sqrt()
                    * p_data.beta
                    * (2.0 * contact_stiffness * reduced_mass).sqrt()
                    * velocity_rel_normal.norm()
                    * velocity_rel_normal.dot(&normalized_delta_position).signum();

                let mut normal_force_total = normal_force + normal_dissipation_force;

                c_data.contact_struct[c_index].force_magnitude[c_i] = normal_force_total;
                c_data.contact_struct[c_index].force_magnitude[c_j] = -normal_force_total;

                normal_force_total = -normal_force_total;

                if velocity_rel_normal.dot(&normalized_delta_position).signum() > 0.0 {
                    // switch from loading to unloading
                    if c_data.contact_struct[c_index].k_loading == 0 {
                        c_data.contact_struct[c_index].k_loading = 1;

                        c_data.contact_struct[c_index].t_star =
                            c_data.contact_struct[c_index].previous_tangential_force_magnitude[c_i];
                    }
                    // switch from reloading to unloading
                    if c_data.contact_struct[c_index].k_loading == 2 {
                        c_data.contact_struct[c_index].k_loading = 1;

                        c_data.contact_struct[c_index].t_star =
                            c_data.contact_struct[c_index].previous_tangential_force_magnitude[c_i];
                    }
                } else if velocity_rel_normal.dot(&normalized_delta_position).signum() < 0.0 {
                    // switch from unloading to reloading
                    if c_data.contact_struct[c_index].k_loading == 1 {
                        c_data.contact_struct[c_index].k_loading = 2;

                        c_data.contact_struct[c_index].t_double_star =
                            c_data.contact_struct[c_index].previous_tangential_force_magnitude[c_i];
                    }
                }

                let velocity_rel_tangential = velocity_rel - velocity_rel_normal;
                let tangential_direction = velocity_rel_tangential / velocity_rel_tangential.norm();

                let tangential_displacement_vector = velocity_rel_tangential * _dt; //questionable

                // let undm = velocity_rel.dot(&normalized_delta_position) * _dt;
                // let und = undm * normalized_delta_position;

                // let td = velocity_rel * _dt - und;

                // let mut sd = (p_data.omega[i] * p_data.radius[i]
                //     + p_data.omega[j] * p_data.radius[j])
                //     .cross(&normalized_delta_position);

                // sd = td - sd * _dt;

                let d_old = c_data.contact_struct[c_index].tangential_displacement;
                c_data.contact_struct[c_index].tangential_displacement =
                    tangential_displacement_vector;
                let d_sum = d_old.dot(&c_data.contact_struct[c_index].tangential_displacement);

                if d_sum < 0.0 {
                    // println!("change directions");
                    c_data.contact_struct[c_index].change_tang_direction *= -1.0;
                }

                let tangential_displacement = c_data.contact_struct[c_index]
                    .tangential_displacement
                    .norm()
                    * c_data.contact_struct[c_index].change_tang_direction; //questionable

                c_data.contact_struct[c_index].sd += tangential_displacement;

                let delta_force_normal = -c_data.contact_struct[c_index].force_magnitude[c_i]
                    + c_data.contact_struct[c_index].previous_force_magnitude[c_i];

                c_data.contact_struct[c_index].sum_normal_force += delta_force_normal;
                // println!("{}", delta_force_normal);
                let mut theta = 1.0;

                if p_data.friction != 0.0 {
                    //No Slip solution

                    if c_data.contact_struct[c_index].k_loading == 0 {
                        // println!(
                        //     "k=0 {} {} {}",
                        //     c_data.contact_struct[c_index].previous_tangential_force[c_i].norm(),
                        //     delta_force_normal,
                        //     normal_force_total
                        // );
                        let theda_cubed = (1.0
                            - (c_data.contact_struct[c_index].previous_tangential_force_magnitude
                                [c_i]
                                + p_data.friction * delta_force_normal)
                                / (p_data.friction * (normal_force_total)));
                        theta = theda_cubed.abs().powf(1.0 / 3.0) * theda_cubed.signum();
                        // println!("0 {} {}", theta, theda_cubed);
                    }
                    if c_data.contact_struct[c_index].k_loading == 1 {
                        // println!(
                        //     "{} {} {} {}",
                        //     c_data.contact_struct[c_index].previous_tangential_force[c_i]
                        //         .norm(),
                        //     delta_force_normal,
                        //     normal_force_total,
                        //     p_data.friction
                        // );

                        let theda_cubed = 1.0
                            - ((c_data.contact_struct[c_index].t_star
                                - c_data.contact_struct[c_index]
                                    .previous_tangential_force_magnitude[c_i]
                                + 2.0 * p_data.friction * delta_force_normal)
                                / (2.0 * p_data.friction * (normal_force_total)));
                        theta = theda_cubed.abs().powf(1.0 / 3.0) * theda_cubed.signum();
                        // println!("1 {} {}", theda, theda_cubed);
                    }
                    if c_data.contact_struct[c_index].k_loading == 2 {
                        let theda_cubed = (1.0
                            - ((c_data.contact_struct[c_index]
                                .previous_tangential_force_magnitude[c_i]
                                - c_data.contact_struct[c_index].t_double_star)
                                + 2.0 * p_data.friction * delta_force_normal)
                                / (2.0 * p_data.friction * (normal_force_total)));
                        theta = theda_cubed.abs().powf(1.0 / 3.0) * theda_cubed.signum();
                        // println!("2 {}", theda)
                    }

                    // println!(
                    //     "{} {} {} {} {}",
                    //     effective_shear_modulus,
                    //     c_data.contact_struct[c_index].contact_radius,
                    //     theda,
                    //     tangential_displacement,
                    //     delta_force_normal
                    // );
                    // println!("theta = {}", theta);

                    let mut delta_tangential_force = 8.0
                        * effective_shear_modulus
                        * c_data.contact_struct[c_index].contact_radius
                        * c_data.contact_struct[c_index].sd;

                    // println!(
                    //     "theda = 1 conditio , force = {} < {}",
                    //     delta_tangential_force,
                    //     p_data.friction * delta_force_normal.abs()
                    // );
                    if delta_tangential_force.abs()
                        < p_data.friction * c_data.contact_struct[c_index].sum_normal_force
                    {
                        theta = 1.0;
                    }

                    delta_tangential_force = 8.0
                        * effective_shear_modulus
                        * c_data.contact_struct[c_index].contact_radius
                        * theta
                        * tangential_displacement
                        + (-1.0f64).powi(c_data.contact_struct[c_index].k_loading)
                            * p_data.friction
                            * delta_force_normal
                            * (1.0 - theta);

                    // println!(
                    //     "{} {} {} {} {} ",
                    //     effective_shear_modulus,
                    //     c_data.contact_struct[c_index].contact_radius,
                    //     theda,
                    //     tangential_displacement,
                    //     (-1.0f64).powi(c_data.contact_struct[c_index].k_loading)
                    //         * p_data.friction
                    //         * delta_force_normal
                    //         * (1.0 - theda)
                    // );
                    let mut tangential_stiffness = 8.0
                        * effective_shear_modulus
                        * c_data.contact_struct[c_index].contact_radius
                        * theta
                        + (-1.0f64).powi(c_data.contact_struct[c_index].k_loading)
                            * p_data.friction
                            * delta_force_normal
                            * (1.0 - theta)
                            / tangential_displacement;
                    tangential_stiffness = tangential_stiffness.abs();
                    // println!(
                    //     "{} {} {} {}",
                    //     p_data.beta, tangential_stiffness, reduced_mass, velocity_rel_tangential
                    // );

                    let tangential_dissipation_force = (5.0 / 6.0).sqrt()
                        * p_data.beta
                        * (2.0 * tangential_stiffness * reduced_mass).sqrt()
                        * velocity_rel_tangential.norm();
                    //* velocity_rel_normal.dot(&normalized_delta_position).signum();

                    // println!(
                    //     "before adding delta tf {} {} {}",
                    //     c_data.contact_struct[c_index].tangential_force[c_i],
                    //     velocity_rel_tangential,
                    //     tangential_displacement,
                    // );

                    // println!("{} {}", tangential_direction, delta_tangential_force);
                    // println!(
                    //     "tangential stuff {} {} {} {}",
                    //     delta_tangential_force,
                    //     delta_tangential_force * tangential_direction,
                    //     tangential_dissipation_force,
                    //     delta_tangential_force * tangential_direction
                    //         - tangential_dissipation_force,
                    // );

                    // c_data.contact_struct[c_index].tangential_force_magnitude[c_i] =
                    //     c_data.contact_struct[c_index].previous_tangential_force_magnitude[c_i]
                    //         + delta_tangential_force;

                    c_data.contact_struct[c_index].tangential_force_magnitude[c_i] +=
                        delta_tangential_force - tangential_dissipation_force;

                    c_data.contact_struct[c_index].tangential_force_magnitude[c_j] +=
                        -delta_tangential_force + tangential_dissipation_force;

                    if c_data.contact_struct[c_index].tangential_force_magnitude[c_i].abs()
                        > p_data.friction * normal_force_total.abs()
                    {
                        // println!("force = nu * delta P");
                        c_data.contact_struct[c_index].tangential_force_magnitude[c_i] = p_data
                            .friction
                            * normal_force_total.abs()
                            * c_data.contact_struct[c_index].tangential_force_magnitude[c_i]
                                .signum()
                    }
                    if c_data.contact_struct[c_index].tangential_force_magnitude[c_j].abs()
                        > p_data.friction * normal_force_total.abs()
                    {
                        c_data.contact_struct[c_index].tangential_force_magnitude[c_j] = p_data
                            .friction
                            * normal_force_total.abs()
                            * c_data.contact_struct[c_index].tangential_force_magnitude[c_j]
                                .signum()
                    }

                    c_data.contact_struct[c_index].tangential_force[c_i] =
                        c_data.contact_struct[c_index].tangential_force_magnitude[c_i]
                            * tangential_displacement_vector;
                    c_data.contact_struct[c_index].tangential_force[c_j] =
                        c_data.contact_struct[c_index].tangential_force_magnitude[c_j]
                            * tangential_displacement_vector;

                    println!(
                        "{} {} {} {} {} {}",
                        c_data.contact_struct[c_index].k_loading,
                        distance_delta,
                        normal_force_total,
                        // delta_force_normal,
                        tangential_displacement,
                        theta,
                        // c_data.contact_struct[c_index].tangential_force_magnitude[c_i],
                        c_data.contact_struct[c_index].tangential_force_magnitude[c_i],
                    );

                    // println!("{}", c_data.contact_struct[c_index].tangential_force[c_i]);

                    // println!(
                    //     "before friction {}",
                    //     c_data.contact_struct[c_index].tangential_force[c_i]
                    // );

                    p_data.torque[i] += c_data.contact_struct[c_index].tangential_force[c_i]
                        .cross(&c_data.contact_struct[c_index].contact_vector[c_i]);
                    p_data.torque[j] += c_data.contact_struct[c_index].tangential_force[c_j]
                        .cross(&c_data.contact_struct[c_index].contact_vector[c_j]);
                }
                // let delta_veloctiy = v2 - v1;
                // let f_dot = normalized_delta.dot(&delta_veloctiy);
                // let f_dot = normalized_delta.dot(&velocity_rel_normal);

                // let v_r_n = f_dot * normalized_delta;

                //  println!("{} {} {}", distance_delta, normal_force, dissipation_force);

                // println!(
                //     "{} {}",
                //     normal_force_total, c_data.contact_struct[c_index].tangential_force[c_i]
                // );

                p_data.force[i] += c_data.contact_struct[c_index].force_magnitude[c_i]
                    * normalized_delta_position
                    - c_data.contact_struct[c_index].tangential_force[c_i];

                p_data.force[j] += c_data.contact_struct[c_index].force_magnitude[c_j]
                    * normalized_delta_position
                    - c_data.contact_struct[c_index].tangential_force[c_j];

                // println!("{} {}", p_data.omega[i][2], p_data.omega[j][2])

                // let contact_area_cubed = 3.0 * normal_force * effective_radius / 4.0 * effective_youngs;

                // let tangential_force = p_data.friction * normal_force;

                // let force_length_matrix = ((normal_force - dissipation_force) * normalized_delta)
                //     * delta_position.transpose();
                // f_data.forcedata.push(force_length_matrix);

                // println!(
                //     "{} {}",
                //     force_tangential_i.norm(),
                //     (normal_force - dissipation_force)
                // );
            } else {
                remove_collision_struct(c_data, i, j);
            }
        }
    }
}

pub fn _euler_integration(p_data: &mut sphere::ParticleData, dt: f64) {
    for i in 0..p_data.radius.len() {
        p_data.velocity[i] += dt * p_data.force[i] / p_data.mass[i];
        p_data.position[i] += p_data.velocity[i] * dt;

        p_data.omega[i] += dt * p_data.torque[i] / p_data.moment_of_inertia[i];
        // p_data.quaterion[i] += p_data.omega[i] * dt;

        if p_data.is_collision[i] == false {
            p_data.torque[i] = Vector3::zeros();
        }
    }
}

pub fn inital_integrate(p_data: &mut sphere::ParticleData, dt: f64) {
    for i in 0..p_data.radius.len() {
        p_data.velocity[i] += 0.5 * dt * p_data.force[i] / p_data.mass[i];
        p_data.position[i] += p_data.velocity[i] * dt;

        p_data.omega[i] += 0.5 * dt * p_data.torque[i] / p_data.moment_of_inertia[i];
        // p_data.quaterion[i] += p_data.omega[i] * dt;
    }
}

pub fn final_integrate(p_data: &mut sphere::ParticleData, dt: f64) {
    for i in 0..p_data.radius.len() {
        p_data.velocity[i] += 0.5 * dt * p_data.force[i] / p_data.mass[i];

        if p_data.is_collision[i] == false {
            p_data.torque[i] = Vector3::zeros();
        }

        p_data.omega[i] += 0.5 * dt * p_data.torque[i] / p_data.moment_of_inertia[i];

        // p_data.pervious_tangential_force[i] = p_data.torque[i];
    }
}

pub fn remove_collision_struct(c_data: &mut sphere::ContactData, i: usize, j: usize) {
    let mut index = -1;
    for a in 0..c_data.particle_indexes.len() {
        if c_data.particle_indexes[a][0] == i && c_data.particle_indexes[a][1] == j {
            index = a as i32;
            // println!("{}", index);
        }

        if c_data.particle_indexes[a][0] == j && c_data.particle_indexes[a][1] == i {
            index = a as i32;
            // println!("{}", index);
        }
    }
    if index >= 0 {
        c_data.contact_struct.remove(index as usize);
        c_data.particle_indexes.remove(index as usize);
        println!("removed particle collision")
    }
}
pub fn get_collision_struct(c_data: &mut sphere::ContactData, i: usize, j: usize) -> (usize, i32) {
    for a in 0..c_data.particle_indexes.len() {
        // i and j index match
        if c_data.particle_indexes[a][0] == i && c_data.particle_indexes[a][1] == j {
            return (a, 1);
        }
        // i and j index don't match
        if c_data.particle_indexes[a][0] == j && c_data.particle_indexes[a][1] == i {
            return (a, -1);
        }
    }

    // not found, need to make a new collision index
    c_data.particle_indexes.push(Vector2::new(i, j));

    let mut contact_struct = sphere::ContactStruct {
        force_magnitude: Vector2::new(0.0, 0.0),
        previous_force_magnitude: Vector2::new(0.0, 0.0),
        tangential_force_magnitude: Vector2::new(0.0, 0.0),
        previous_tangential_force_magnitude: Vector2::new(0.0, 0.0),
        force: Vector2::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0)),
        previous_force: Vector2::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0)),
        tangential_force: Vector2::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0)),
        previous_tangential_force: Vector2::new(
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0),
        ),
        contact_vector: Vector2::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0)),
        contact_radius: 0.0,
        t_star: 0.0,
        t_double_star: 0.0,
        k_loading: 0,
        sd: 0.0,
        sum_normal_force: 0.0,
        change_tang_direction: 1.0,
        tangential_displacement: Vector3::new(0.0, 0.0, 0.0),
    };
    c_data.contact_struct.push(contact_struct);
    return (c_data.particle_indexes.len() - 1, 0);
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

pub fn _collisions(
    d_data: &domain::DomainData,
    p_data: &mut sphere::ParticleData,
    f_data: &mut sphere::ContactData,
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

                            // let force_length_matrix = ((normal_force - dissipation_force)
                            //     * normalized_delta)
                            //     * delta_position.transpose();
                            // f_data.forcedata.push(force_length_matrix);
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

                            // let force_length_matrix =
                            //     ((normal_force - dissipation_force) * normalized_delta * 0.5)
                            //         * delta_position.transpose();
                            // f_data.forcedata.push(force_length_matrix);
                        }
                    }
                }
            }
        }
    }
}
