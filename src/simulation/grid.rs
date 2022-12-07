use std::process::exit;

use nalgebra::AbstractRotation;
use nalgebra::Quaternion;
use nalgebra::UnitQuaternion;
use nalgebra::Vector2;
use nalgebra::Vector3;
use nalgebra::Vector4;

use crate::domain;
use crate::rigid_body;

pub fn relax_boundaries_box(
    d_data: &mut domain::DomainData,
    p_data: &mut rigid_body::RigidBodiesData,
) {
    for i in 0..p_data.body_data.len() {
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

        p_data.update_one_body_spheres_positions(i);

        //Check Individual Spheres and apply a wrap int to them
        for j in 0..p_data.body_data[i].position.len() {
            if p_data.body_data[i].position[j][1] > d_data.domain[1] {
                // p_data.position[i][1] = d_data.domain[1] - p_data.radius[i]
                p_data.body_data[i].position[j][1] -= d_data.domain[1];
                p_data.body_data[i].wrap[j][1] = -1;
            }
            // if particle is less than domain move to end of domain
            // Also apply velocity change to particle for shearing
            else if p_data.body_data[i].position[j][1] <= 0.0 {
                // p_data.position[i][1] = p_data.radius[i]
                p_data.body_data[i].position[j][1] += d_data.domain[1];
                p_data.body_data[i].wrap[j][1] = 1;
            }
            // std::cout << distb << std::endl;//X boundary condition
            // if particles is greater than domain move to beginning of domain
            if p_data.body_data[i].position[j][0] > d_data.domain[0] {
                // p_data.position[i][0] = d_data.domain[0] - p_data.radius[i]
                p_data.body_data[i].position[j][0] -= d_data.domain[0];
                p_data.body_data[i].wrap[j][0] = -1;
            }
            // if particle is less than domain move to end of domain
            else if p_data.body_data[i].position[j][0] <= 0.0 {
                // p_data.position[i][0] = p_data.radius[i]
                p_data.body_data[i].position[j][0] += d_data.domain[0];
                p_data.body_data[i].wrap[j][0] = 1;
            }

            // Z boundary condition
            // if particles is greater than domain move to beginning of domain
            if p_data.body_data[i].position[j][2] > d_data.domain[2] {
                // p_data.position[i][2] = d_data.domain[2] - p_data.radius[i]
                p_data.body_data[i].position[j][2] -= d_data.domain[2];
                p_data.body_data[i].wrap[j][2] = -1;
            }
            // if particle is less than domain move to end of domain
            else if p_data.body_data[i].position[j][2] <= 0.0 {
                // p_data.position[i][2] = p_data.radius[i]
                p_data.body_data[i].position[j][2] += d_data.domain[2];
                p_data.body_data[i].wrap[j][2] = 1;
            }
        }
    }
}

pub fn update(d_data: &mut domain::DomainData, p_data: &mut rigid_body::RigidBodiesData) {
    //Clears all boxes of particles in the box
    for i in 0..d_data.g_data.len() {
        for j in 0..d_data.g_data[i].len() {
            for k in 0..d_data.g_data[i][j].len() {
                d_data.g_data[i][j][k].real.clear();
                d_data.g_data[i][j][k].ghost.clear();

                for index in 0..p_data.body_data.len() {
                    for sphere_index in 0..p_data.body_data[index].position.len() {
                        if d_data.g_data[i][j][k]
                            .is_position_in_box(p_data.body_data[index].position[sphere_index])
                        {
                            d_data.g_data[i][j][k]
                                .real
                                .push((index.try_into().unwrap(), sphere_index.try_into().unwrap()))
                        } else if d_data.g_data[i][j][k].is_sphere_aabb_in_radius_enlarged_box(
                            p_data.body_data[index].position[sphere_index],
                            p_data.body_data[index].radius[sphere_index],
                            p_data.max_radius,
                        ) {
                            d_data.g_data[i][j][k]
                                .ghost
                                .push((index.try_into().unwrap(), sphere_index.try_into().unwrap()))
                        } else if ((i == 0 || i == d_data.g_data.len() - 1)
                            || (j == 0 || j == d_data.g_data[i].len() - 1)
                            || (k == 0 || k == d_data.g_data[i][j].len() - 1))
                            && d_data.g_data[i][j][k].is_periodic_sphere(
                                p_data.body_data[index].position[sphere_index],
                                p_data.body_data[index].radius[sphere_index],
                                p_data.max_radius,
                                d_data,
                            )
                        {
                            d_data.g_data[i][j][k]
                                .ghost
                                .push((index.try_into().unwrap(), sphere_index.try_into().unwrap()))
                        }
                    }
                }
            }
        }
    }
}

pub fn relax(
    d_data: &mut domain::DomainData,
    p_data: &mut rigid_body::RigidBodiesData,
    _relax_rate: f64,
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
                        let (i_body_int, i_int) =
                            d_data.g_data[box_i as usize][box_j as usize][box_k as usize].real[ii];
                        let (j_body_int, j_int) =
                            d_data.g_data[box_i as usize][box_j as usize][box_k as usize].real[jj];

                        let i_body = i_body_int as usize;
                        let j_body = j_body_int as usize;
                        let i = i_int as usize;
                        let j = j_int as usize;

                        if i_body == j_body {
                            continue;
                        }

                        let mut p1 = p_data.body_data[i_body].position[i];
                        let mut p2 = p_data.body_data[j_body].position[j];

                        let r1 = p_data.body_data[i_body].radius[i];
                        let r2 = p_data.body_data[j_body].radius[j];

                        let mut wrap1 = p_data.body_data[i_body].wrap[i];
                        let mut wrap2 = p_data.body_data[j_body].wrap[j];

                        if p1[1] - r1 + d_data.domain[1] <= p2[1] + r2 {
                            p1[1] += d_data.domain[1];
                            wrap1[1] += 1;
                        } else if p2[1] - r2 + d_data.domain[1] <= p1[1] + r1 {
                            p2[1] += d_data.domain[1];
                            wrap2[1] += 1;
                        }
                        if p1[0] - r1 + d_data.domain[0] <= p2[0] + r2 {
                            p1[0] += d_data.domain[0];
                            wrap1[0] += 1;
                        } else if p2[0] - r2 + d_data.domain[0] <= p1[0] + r1 {
                            p2[0] += d_data.domain[0];
                            wrap2[0] += 1;
                        }

                        if p1[2] - r1 + d_data.domain[2] <= p2[2] + r2 {
                            p1[2] += d_data.domain[2];
                            wrap1[2] += 1;
                        } else if p2[2] - r2 + d_data.domain[2] <= p1[2] + r1 {
                            p2[2] += d_data.domain[2];
                            wrap2[2] += 1;
                        }
                        let delta_position = p2 - p1;

                        let distance = delta_position.norm();

                        // println!(
                        //     "Distance: {} Radius+radius: {}",
                        //     distance,
                        //     p_data.body_data[i_body].radius[i] + p_data.body_data[j_body].radius[j]
                        // );

                        if distance
                            < p_data.body_data[i_body].radius[i]
                                + p_data.body_data[j_body].radius[j]
                        {
                            let contact_point =
                                p_data.body_data[i_body].position[i] + delta_position / 2.0;
                            let body_delta_position_i = contact_point - p_data.position[i_body]
                                + Vector3::new(
                                    -wrap1[0] as f64 * d_data.domain[0],
                                    -wrap1[1] as f64 * d_data.domain[1],
                                    -wrap1[2] as f64 * d_data.domain[2],
                                );

                            let body_delta_position_j = contact_point - p_data.position[j_body]
                                + Vector3::new(
                                    -wrap2[0] as f64 * d_data.domain[0],
                                    -wrap2[1] as f64 * d_data.domain[1],
                                    -wrap2[2] as f64 * d_data.domain[2],
                                );
                            p_data.is_collision[i_body] = true;
                            p_data.is_collision[j_body] = true;

                            let normalized_delta = delta_position / distance;

                            // let normalized_body_delta_i = body_delta_position_i.normalize();
                            // let normalized_body_delta_j = body_delta_position_j.normalize();

                            // let normal_force_i = normalized_body_delta_i
                            //     * (-normalized_delta).dot(&normalized_body_delta_i)
                            //     / 1.0;
                            // let normal_force_j = normalized_body_delta_j
                            //     * normalized_delta.dot(&normalized_body_delta_j)
                            //     / 1.0;

                            let torque_i = body_delta_position_i.cross(&(-normalized_delta));
                            let torque_j = body_delta_position_j.cross(&(normalized_delta));

                            p_data.torque[i_body] += torque_i;
                            p_data.torque[j_body] += torque_j;

                            //println!("{:?}",normalized_delta);

                            p_data.force[i_body] += -normalized_delta;
                            p_data.force[j_body] += normalized_delta;
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
                        let (i_body_int, i_int) =
                            d_data.g_data[box_i as usize][box_j as usize][box_k as usize].real[ii];
                        let (j_body_int, j_int) =
                            d_data.g_data[box_i as usize][box_j as usize][box_k as usize].ghost[jj];

                        let i_body = i_body_int as usize;
                        let j_body = j_body_int as usize;
                        let i = i_int as usize;
                        let j = j_int as usize;

                        if i_body == j_body {
                            continue;
                        }

                        // println!("i {}, j {}", i,j);
                        let mut p1 = p_data.body_data[i_body].position[i];
                        let mut p2 = p_data.body_data[j_body].position[j];

                        let r1 = p_data.body_data[i_body].radius[i];
                        let r2 = p_data.body_data[j_body].radius[j];

                        let mut wrap1 = p_data.body_data[i_body].wrap[i];
                        let mut wrap2 = p_data.body_data[j_body].wrap[j];

                        if p1[1] - r1 + d_data.domain[1] <= p2[1] + r2 {
                            p1[1] += d_data.domain[1];
                            wrap1[1] += 1;
                        } else if p2[1] - r2 + d_data.domain[1] <= p1[1] + r1 {
                            p2[1] += d_data.domain[1];
                            wrap2[1] += 1;
                        }
                        if p1[0] - r1 + d_data.domain[0] <= p2[0] + r2 {
                            p1[0] += d_data.domain[0];
                            wrap1[0] += 1;
                        } else if p2[0] - r2 + d_data.domain[0] <= p1[0] + r1 {
                            p2[0] += d_data.domain[0];
                            wrap2[0] += 1;
                        }

                        if p1[2] - r1 + d_data.domain[2] <= p2[2] + r2 {
                            p1[2] += d_data.domain[2];
                            wrap1[2] += 1;
                        } else if p2[2] - r2 + d_data.domain[2] <= p1[2] + r1 {
                            p2[2] += d_data.domain[2];
                            wrap2[2] += 1;
                        }

                        let delta_position = p2 - p1;
                        // let delta_position = p_data.position[j] - p_data.position[i];

                        let distance = delta_position.norm();

                        //println!("Distance: {} Radius+radius: {}", distance, p_data.radius[i] + p_data.radius[j]);

                        if distance
                            < p_data.body_data[i_body].radius[i]
                                + p_data.body_data[j_body].radius[j]
                        {
                            let contact_point = p1 + delta_position / 2.0;

                            if p_data.body_data[i_body].wrap[i][0].abs() > 1
                                || p_data.body_data[i_body].wrap[i][1].abs() > 1
                                || p_data.body_data[i_body].wrap[i][2].abs() > 1
                            {
                                println!("{}", p_data.body_data[i_body].wrap[i])
                            }
                            let body_delta_position_i = contact_point - p_data.position[i_body]
                                + Vector3::new(
                                    -wrap1[0] as f64 * d_data.domain[0],
                                    -wrap1[1] as f64 * d_data.domain[1],
                                    -wrap1[2] as f64 * d_data.domain[2],
                                );

                            if p_data.body_data[j_body].wrap[j][0].abs() > 1
                                || p_data.body_data[j_body].wrap[j][1].abs() > 1
                                || p_data.body_data[j_body].wrap[j][2].abs() > 1
                            {
                                println!("{}", p_data.body_data[j_body].wrap[j])
                            }

                            let body_delta_position_j = contact_point - p_data.position[j_body]
                                + Vector3::new(
                                    -wrap2[0] as f64 * d_data.domain[0],
                                    -wrap2[1] as f64 * d_data.domain[1],
                                    -wrap2[2] as f64 * d_data.domain[2],
                                );
                            p_data.is_collision[i_body] = true;
                            p_data.is_collision[j_body] = true;

                            let normalized_delta = delta_position / distance;

                            // let normalized_body_delta_i = body_delta_position_i.normalize();
                            // let normalized_body_delta_j = body_delta_position_j.normalize();

                            // let normal_force_i = normalized_body_delta_i
                            //     * (-normalized_delta).dot(&normalized_body_delta_i)
                            //     / 1.0;
                            // let normal_force_j = normalized_body_delta_j
                            //     * normalized_delta.dot(&normalized_body_delta_j)
                            //     / 1.0;

                            let torque_i = body_delta_position_i.cross(&(-normalized_delta));
                            let torque_j = body_delta_position_j.cross(&(normalized_delta));

                            p_data.torque[i_body] += torque_i * 0.5;
                            p_data.torque[j_body] += torque_j * 0.5;

                            //println!("{:?}",normalized_delta);

                            p_data.force[i_body] += -normalized_delta * 0.5;
                            p_data.force[j_body] += normalized_delta * 0.5;
                        }
                    }
                }
            }
        }
    }
}

pub fn collisions(
    d_data: &domain::DomainData,
    p_data: &mut rigid_body::RigidBodiesData,
    f_data: &mut rigid_body::ForceData,
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
                        let (i_body_int, i_int) =
                            d_data.g_data[box_i as usize][box_j as usize][box_k as usize].real[ii];
                        let (j_body_int, j_int) =
                            d_data.g_data[box_i as usize][box_j as usize][box_k as usize].real[jj];

                        let i_body = i_body_int as usize;
                        let j_body = j_body_int as usize;
                        let i = i_int as usize;
                        let j = j_int as usize;

                        if i_body == j_body {
                            continue;
                        }

                        let mut p1 = p_data.body_data[i_body].position[i];
                        let mut p2 = p_data.body_data[j_body].position[j];

                        let mut v1 = p_data.velocity[i_body];
                        let mut v2 = p_data.velocity[j_body];

                        let r1 = p_data.body_data[i_body].radius[i];
                        let r2 = p_data.body_data[j_body].radius[j];

                        let mut wrap1 = p_data.body_data[i_body].wrap[i];
                        let mut wrap2 = p_data.body_data[j_body].wrap[j];

                        // if p1[1] - r1 + d_data.domain[1] <= p2[1] + r2 {
                        //     p1[1] += d_data.domain[1];
                        //     wrap1[1] += 1;
                        //     // v1[0] += d_data.lees_edwards_boundary * d_data.domain[1];
                        //     // p1[0] += ledisplace;
                        //     // if p1[0] > d_data.domain[0] {
                        //     //     p1[0] -= d_data.domain[0];
                        //     //     wrap1[0] -= 1;
                        //     // }
                        // } else if p2[1] - r2 + d_data.domain[1] <= p1[1] + r1 {
                        //     p2[1] += d_data.domain[1];
                        //     wrap2[1] += 1;
                        //     // v2[0] += d_data.lees_edwards_boundary * d_data.domain[1];
                        //     // p2[0] += ledisplace;
                        //     // if p2[0] > d_data.domain[0] {
                        //     //     p2[0] -= d_data.domain[0];
                        //     //     wrap2[0] -= 1;
                        //     // }
                        // }
                        // if p1[0] - r1 + d_data.domain[0] <= p2[0] + r2 {
                        //     p1[0] += d_data.domain[0];
                        //     wrap1[0] += 1;
                        // } else if p2[0] - r2 + d_data.domain[0] <= p1[0] + r1 {
                        //     p2[0] += d_data.domain[0];
                        //     wrap2[0] += 1;
                        // }

                        // if p1[2] - r1 + d_data.domain[2] <= p2[2] + r2 {
                        //     p1[2] += d_data.domain[2];
                        //     wrap1[2] += 1;
                        // } else if p2[2] - r2 + d_data.domain[2] <= p1[2] + r1 {
                        //     p2[2] += d_data.domain[2];
                        //     wrap2[2] += 1;
                        // }
                        let delta_position = p2 - p1;

                        let distance = delta_position.norm();

                        // println!(
                        //     "Distance: {} Radius+radius: {}",
                        //     distance,
                        //     p_data.body_data[i_body].radius[i] + p_data.body_data[j_body].radius[j]
                        // );

                        if distance
                            < p_data.body_data[i_body].radius[i]
                                + p_data.body_data[j_body].radius[j]
                        {
                            let contact_point = p1 + delta_position / 2.0;
                            let body_delta_position_i = contact_point - p_data.position[i_body]
                                + Vector3::new(
                                    -wrap1[0] as f64 * d_data.domain[0],
                                    -wrap1[1] as f64 * d_data.domain[1],
                                    -wrap1[2] as f64 * d_data.domain[2],
                                );

                            let body_delta_position_j = contact_point - p_data.position[j_body]
                                + Vector3::new(
                                    -wrap2[0] as f64 * d_data.domain[0],
                                    -wrap2[1] as f64 * d_data.domain[1],
                                    -wrap2[2] as f64 * d_data.domain[2],
                                );

                            // if body_delta_position_i.norm() > d_data.domain[2]
                            //     || body_delta_position_j.norm() > d_data.domain[2]
                            // {
                            //     println!("wrap error");
                            //     exit(1);
                            // }

                            let normalized_delta = delta_position / distance;

                            let effective_radius = 1.0
                                / (1.0 / p_data.body_data[i_body].radius[i]
                                    + 1.0 / p_data.body_data[j_body].radius[j]);

                            let distance_delta = (p_data.body_data[i_body].radius[i]
                                + p_data.body_data[j_body].radius[j])
                                - distance;

                            let effective_youngs = 1.0
                                / ((1.0 - p_data.poisson_ratio[i] * p_data.poisson_ratio[i])
                                    / p_data.youngs_mod[i]
                                    + (1.0 - p_data.poisson_ratio[j] * p_data.poisson_ratio[j])
                                        / p_data.youngs_mod[j]);

                            let contact_stiffness =
                                2.0 * effective_youngs * (effective_radius * distance_delta).sqrt();

                            let reduced_mass = p_data.body_data[i_body].mass
                                * p_data.body_data[j_body].mass
                                / (p_data.body_data[i_body].mass + p_data.body_data[j_body].mass);

                            let delta_veloctiy = p_data.velocity[j] - p_data.velocity[i]
                                + p_data.omega[i].cross(&body_delta_position_i)
                                + p_data.omega[j].cross(&body_delta_position_j);

                            let normal_force = -4.0 / 3.0
                                * effective_youngs
                                * (effective_radius).sqrt()
                                * distance_delta.powf(3.0 / 2.0);

                            let normal_dissipation_force = (5.0 as f64 / 6.0 as f64).sqrt()
                                * p_data.beta
                                * (2.0 * contact_stiffness * reduced_mass).sqrt()
                                * delta_veloctiy.norm()
                                * delta_veloctiy.dot(&normalized_delta).signum();

                            let mut normal_force_total = normal_force + normal_dissipation_force;

                            // println!("beta: {}, contact_stiffness: {}, reduced_mass: {}, delta_velocity: {}", p_data.beta, contact_stiffness, reduced_mass,delta_veloctiy);

                            let force_i = normal_force_total * normalized_delta;
                            let force_j = -normal_force_total * normalized_delta;

                            p_data.is_collision[i_body] = true;
                            p_data.is_collision[j_body] = true;

                            let normalized_body_delta_i = body_delta_position_i.normalize();
                            let normalized_body_delta_j = body_delta_position_j.normalize();

                            let normal_force_i =
                                normalized_body_delta_i * force_i.dot(&normalized_body_delta_i);
                            let normal_force_j =
                                normalized_body_delta_j * force_j.dot(&normalized_body_delta_j);

                            let tangent_force_i = force_i - normal_force_i;
                            let tangent_force_j = force_j - normal_force_j;

                            let torque_i = body_delta_position_i.cross(&tangent_force_i);
                            let torque_j = body_delta_position_j.cross(&tangent_force_j);

                            p_data.torque[i_body] += torque_i;
                            p_data.torque[j_body] += torque_j;

                            //

                            p_data.force[i_body] += normal_force_i;
                            p_data.force[j_body] += normal_force_j;

                            // println!("force {} {}", p_data.force[i_body], p_data.force[j_body]);
                            // println!("torque {} {}", p_data.torque[i_body], p_data.torque[j_body]);
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
                        let (i_body_int, i_int) =
                            d_data.g_data[box_i as usize][box_j as usize][box_k as usize].real[ii];
                        let (j_body_int, j_int) =
                            d_data.g_data[box_i as usize][box_j as usize][box_k as usize].ghost[jj];

                        let i_body = i_body_int as usize;
                        let j_body = j_body_int as usize;
                        let i = i_int as usize;
                        let j = j_int as usize;

                        if i_body == j_body {
                            continue;
                        }

                        // println!("i {}, j {}", i,j);
                        let mut p1 = p_data.body_data[i_body].position[i];
                        let mut p2 = p_data.body_data[j_body].position[j];

                        let mut v1 = p_data.velocity[i_body];
                        let mut v2 = p_data.velocity[j_body];

                        let r1 = p_data.body_data[i_body].radius[i];
                        let r2 = p_data.body_data[j_body].radius[j];

                        let mut wrap1 = p_data.body_data[i_body].wrap[i];
                        let mut wrap2 = p_data.body_data[j_body].wrap[j];

                        if p1[1] - r1 + d_data.domain[1] <= p2[1] + r2 {
                            p1[1] += d_data.domain[1];
                            wrap1[1] += 1;
                            // v1[0] += d_data.lees_edwards_boundary * d_data.domain[1];
                            // p1[0] += ledisplace;
                            // if p1[0] > d_data.domain[0] {
                            //     p1[0] -= d_data.domain[0];
                            //     wrap1[0] -= 1;
                            // }
                        } else if p2[1] - r2 + d_data.domain[1] <= p1[1] + r1 {
                            p2[1] += d_data.domain[1];
                            wrap2[1] += 1;
                            // v2[0] += d_data.lees_edwards_boundary * d_data.domain[1];
                            // p2[0] += ledisplace;
                            // if p2[0] > d_data.domain[0] {
                            //     p2[0] -= d_data.domain[0];
                            //     wrap2[0] -= 1;
                            // }
                        }
                        if p1[0] - r1 + d_data.domain[0] <= p2[0] + r2 {
                            p1[0] += d_data.domain[0];
                            wrap1[0] += 1;
                        } else if p2[0] - r2 + d_data.domain[0] <= p1[0] + r1 {
                            p2[0] += d_data.domain[0];
                            wrap2[0] += 1;
                        }

                        if p1[2] - r1 + d_data.domain[2] <= p2[2] + r2 {
                            p1[2] += d_data.domain[2];
                            wrap1[2] += 1;
                        } else if p2[2] - r2 + d_data.domain[2] <= p1[2] + r1 {
                            p2[2] += d_data.domain[2];
                            wrap2[2] += 1;
                        }

                        let delta_position = p2 - p1;
                        // let delta_position = p_data.position[j] - p_data.position[i];

                        let distance = delta_position.norm();

                        //println!("Distance: {} Radius+radius: {}", distance, p_data.radius[i] + p_data.radius[j]);

                        if distance
                            < p_data.body_data[i_body].radius[i]
                                + p_data.body_data[j_body].radius[j]
                        {
                            let contact_point = p1 + delta_position / 2.0;
                            let body_delta_position_i = contact_point - p1
                                + Vector3::new(
                                    -wrap1[0] as f64 * d_data.domain[0],
                                    -wrap1[1] as f64 * d_data.domain[1],
                                    -wrap1[2] as f64 * d_data.domain[2],
                                );

                            let body_delta_position_j = contact_point - p2
                                + Vector3::new(
                                    -wrap2[0] as f64 * d_data.domain[0],
                                    -wrap2[1] as f64 * d_data.domain[1],
                                    -wrap2[2] as f64 * d_data.domain[2],
                                );

                            let normalized_delta = delta_position / distance;

                            let effective_radius = 1.0
                                / (1.0 / p_data.body_data[i_body].radius[i]
                                    + 1.0 / p_data.body_data[j_body].radius[j]);

                            let distance_delta = (p_data.body_data[i_body].radius[i]
                                + p_data.body_data[j_body].radius[j])
                                - distance;

                            let effective_youngs = 1.0
                                / ((1.0 - p_data.poisson_ratio[i] * p_data.poisson_ratio[i])
                                    / p_data.youngs_mod[i]
                                    + (1.0 - p_data.poisson_ratio[j] * p_data.poisson_ratio[j])
                                        / p_data.youngs_mod[j]);

                            let contact_stiffness =
                                2.0 * effective_youngs * (effective_radius * distance_delta).sqrt();

                            let normal_force = 2.0 / 3.0 * distance_delta * contact_stiffness;
                            let reduced_mass = p_data.body_data[i_body].mass
                                * p_data.body_data[j_body].mass
                                / (p_data.body_data[i_body].mass + p_data.body_data[j_body].mass);

                            let delta_veloctiy = v2 - v1
                                + p_data.omega[i].cross(&body_delta_position_i)
                                - p_data.omega[j].cross(&body_delta_position_j);

                            let f_dot = normalized_delta.dot(&delta_veloctiy);
                            let v_r_n = f_dot * normalized_delta;

                            let dissipation_force = 2.0
                                * 0.91287092917
                                * p_data.beta
                                * (contact_stiffness * reduced_mass).sqrt()
                                * v_r_n.norm()
                                * v_r_n.dot(&normalized_delta).signum();

                            let force_i = -(normal_force - dissipation_force) * normalized_delta;
                            let force_j = (normal_force - dissipation_force) * normalized_delta;

                            p_data.is_collision[i_body] = true;
                            p_data.is_collision[j_body] = true;

                            let normalized_body_delta_i = body_delta_position_i.normalize();
                            let normalized_body_delta_j = body_delta_position_j.normalize();

                            let normal_force_i =
                                normalized_body_delta_i * force_i.dot(&normalized_body_delta_i);
                            let normal_force_j =
                                normalized_body_delta_j * force_j.dot(&normalized_body_delta_j);

                            let tangent_force_i = force_i - normal_force_i;
                            let tangent_force_j = force_j - normal_force_j;

                            let torque_i = body_delta_position_i.cross(&tangent_force_i);
                            let torque_j = body_delta_position_j.cross(&tangent_force_j);

                            p_data.torque[i_body] += torque_i * 0.5;
                            p_data.torque[j_body] += torque_j * 0.5;

                            //println!("{:?}",normalized_delta);

                            p_data.force[i_body] += normal_force_i * 0.5;
                            p_data.force[j_body] += normal_force_j * 0.5;
                        }
                    }
                }
            }
        }
    }
}

pub fn simple_collisions(
    d_data: &domain::DomainData,
    p_data: &mut rigid_body::RigidBodiesData,
    f_data: &mut rigid_body::ForceData,
    _dt: f64,
    ledisplace: f64,
) {
    for i_body in 0..p_data.body_data.len() {
        for j_body in i_body + 1..p_data.body_data.len() {
            for i in 0..p_data.body_data[i_body].position.len() {
                for j in 0..p_data.body_data[j_body].position.len() {
                    let mut p1 = p_data.body_data[i_body].position[i];
                    let mut p2 = p_data.body_data[j_body].position[j];

                    let mut v1 = p_data.body_data[i_body].velocity[i];
                    let mut v2 = p_data.body_data[j_body].velocity[j];

                    let r1 = p_data.body_data[i_body].radius[i];
                    let r2 = p_data.body_data[j_body].radius[j];

                    let mut wrap1 = p_data.body_data[i_body].wrap[i];
                    let mut wrap2 = p_data.body_data[j_body].wrap[j];

                    if p1[1] - r1 + d_data.domain[1] <= p2[1] + r2 {
                        p1[1] += d_data.domain[1];
                        wrap1[1] += 1;
                        v1[0] += d_data.lees_edwards_boundary * d_data.domain[1];
                        p1[0] += ledisplace;
                        if p1[0] > d_data.domain[0] {
                            p1[0] -= d_data.domain[0];
                            wrap1[0] -= 1;
                        }
                    } else if p2[1] - r2 + d_data.domain[1] <= p1[1] + r1 {
                        p2[1] += d_data.domain[1];
                        wrap2[1] += 1;
                        v2[0] += d_data.lees_edwards_boundary * d_data.domain[1];
                        p2[0] += ledisplace;
                        if p2[0] > d_data.domain[0] {
                            p2[0] -= d_data.domain[0];
                            wrap2[0] -= 1;
                        }
                    }
                    if p1[0] - r1 + d_data.domain[0] <= p2[0] + r2 {
                        p1[0] += d_data.domain[0];
                        wrap1[0] += 1;
                    } else if p2[0] - r2 + d_data.domain[0] <= p1[0] + r1 {
                        p2[0] += d_data.domain[0];
                        wrap2[0] += 1;
                    }

                    if p1[2] - r1 + d_data.domain[2] <= p2[2] + r2 {
                        p1[2] += d_data.domain[2];
                        wrap1[2] += 1;
                    } else if p2[2] - r2 + d_data.domain[2] <= p1[2] + r1 {
                        p2[2] += d_data.domain[2];
                        wrap2[2] += 1;
                    }
                    let delta_position = p2 - p1;

                    let distance = delta_position.norm();

                    // println!(
                    //     "Distance: {} Radius+radius: {}",
                    //     distance,
                    //     p_data.body_data[i_body].radius[i] + p_data.body_data[j_body].radius[j]
                    // );

                    if distance
                        < p_data.body_data[i_body].radius[i] + p_data.body_data[j_body].radius[j]
                    {
                        let contact_point = p1 + delta_position * 0.5;
                        let body_delta_position_i = contact_point - p_data.position[i_body]
                            + Vector3::new(
                                -wrap1[0] as f64 * d_data.domain[0] - wrap1[1] as f64 * ledisplace,
                                -wrap1[1] as f64 * d_data.domain[1],
                                -wrap1[2] as f64 * d_data.domain[2],
                            );

                        let body_delta_position_j = contact_point - p_data.position[j_body]
                            + Vector3::new(
                                -wrap2[0] as f64 * d_data.domain[0] - wrap2[1] as f64 * ledisplace,
                                -wrap2[1] as f64 * d_data.domain[1],
                                -wrap2[2] as f64 * d_data.domain[2],
                            );

                        // if body_delta_position_i.norm() > p_data.body_data[i_body].radius[i] * 2.001
                        //     || body_delta_position_j.norm()
                        //         > p_data.body_data[i_body].radius[i] * 2.001
                        // {
                        //     println!("wtf lever arm");
                        // }

                        let normalized_delta = delta_position.normalize();

                        let effective_radius = 1.0
                            / (1.0 / p_data.body_data[i_body].radius[i]
                                + 1.0 / p_data.body_data[j_body].radius[j]);

                        let distance_delta = (p_data.body_data[i_body].radius[i]
                            + p_data.body_data[j_body].radius[j])
                            - distance;

                        let effective_youngs = 1.0
                            / ((1.0 - p_data.poisson_ratio[i_body] * p_data.poisson_ratio[i_body])
                                / p_data.youngs_mod[i_body]
                                + (1.0
                                    - p_data.poisson_ratio[j_body] * p_data.poisson_ratio[j_body])
                                    / p_data.youngs_mod[j_body]);

                        let contact_stiffness =
                            2.0 * effective_youngs * (effective_radius * distance_delta).sqrt();

                        let reduced_mass = p_data.body_data[i_body].mass
                            * p_data.body_data[j_body].mass
                            / (p_data.body_data[i_body].mass + p_data.body_data[j_body].mass);

                        let delta_velocity = v2 - v1;

                        // (p_data.velocity[j_body]
                        //     + p_data.omega[j_body].cross(&body_delta_position_j))
                        //     - (p_data.velocity[i_body]
                        //         + p_data.omega[i_body].cross(&body_delta_position_i));
                        // let delta_velocity = p_data.body_data[j_body].velocity[j]
                        //     - p_data.body_data[i_body].velocity[i];

                        // if delta_velocity.norm() > d_data.lees_edwards_boundary * d_data.domain[1] {
                        //     // println!(
                        //     //     "i {} v1{} wrap1 {}",
                        //     //     p_data.body_data[i_body].velocity[i], v1, wrap1
                        //     // );
                        //     // println!(
                        //     //     "j{} v2{}  wrap1 {}",
                        //     //     p_data.body_data[j_body].velocity[j], v2, wrap2
                        //     // );

                        //     // exit(1);
                        // }

                        let normal_force = -4.0 / 3.0
                            * effective_youngs
                            * (effective_radius).sqrt()
                            * distance_delta.powf(3.0 / 2.0);

                        let normal_dissipation_force = (5.0 as f64 / 6.0 as f64).sqrt()
                            * p_data.beta
                            * (2.0 * contact_stiffness * reduced_mass).sqrt()
                            * delta_velocity.norm()
                            * (delta_velocity.dot(&normalized_delta).signum());

                        let mut normal_force_total = normal_force + normal_dissipation_force;

                        if normal_force_total.signum() != normal_force.signum() {
                            normal_force_total = normal_force;
                            // println!("beta: {}, contact_stiffness: {}, reduced_mass: {}, delta_velocity: {}", p_data.beta, contact_stiffness, reduced_mass,delta_velocity.norm());
                        }

                        // println!("{} {}", normal_force, normal_dissipation_force);

                        //

                        let force_i = normal_force_total * normalized_delta;
                        let force_j = -normal_force_total * normalized_delta;

                        p_data.is_collision[i_body] = true;
                        p_data.is_collision[j_body] = true;

                        // let normalized_body_delta_i = body_delta_position_i.normalize();
                        // let normalized_body_delta_j = body_delta_position_j.normalize();

                        // let normal_force_i =
                        //     normalized_body_delta_i * force_i.dot(&normalized_body_delta_i);
                        // let normal_force_j =
                        //     normalized_body_delta_j * force_j.dot(&normalized_body_delta_j);

                        // let tangent_force_i = force_i - normal_force_i;
                        // let tangent_force_j = force_j - normal_force_j;

                        // if force_i.norm() < tangent_force_i.norm() {
                        //     println!("wtf tangent force i");

                        //     println!(
                        //         "f{} n{} t{} ",
                        //         force_i.norm(),
                        //         normal_force_i.norm(),
                        //         tangent_force_i.norm()
                        //     );
                        // }
                        // if force_j.norm() < tangent_force_j.norm() {
                        //     println!("wtf tangent force j");
                        //     println!(
                        //         "f{} n{} t{} ",
                        //         force_j.norm(),
                        //         normal_force_j.norm(),
                        //         tangent_force_j.norm()
                        //     );
                        //}

                        let torque_i = body_delta_position_i.cross(&force_i);
                        let torque_j = body_delta_position_j.cross(&force_j);

                        p_data.torque[i_body] += torque_i;
                        p_data.torque[j_body] += torque_j;

                        //

                        p_data.force[i_body] += force_i;
                        p_data.force[j_body] += force_j;

                        let force_length_matrix_i = (force_i) * body_delta_position_i.transpose();
                        f_data.forcedata.push(force_length_matrix_i);
                        let force_length_matrix_j = (force_j) * body_delta_position_j.transpose();
                        f_data.forcedata.push(force_length_matrix_j);
                    }
                }
            }
        }
    }
}

// pub fn euler_integration(p_data: &mut rigid_body::RigidBodiesData, dt: f64) {
//     for i in 0..p_data.body_data.len() {
//         p_data.velocity[i] += dt * p_data.force[i] / p_data.body_data[i].mass;
//         p_data.position[i] += p_data.velocity[i] * dt;

//         let global_interia_tensor = p_data.quaternion[i].to_rotation_matrix()
//             * p_data.moment_of_inertia_tensor[i]
//             * p_data.quaternion[i].to_rotation_matrix().transpose();

//         let temp_i =
//             p_data.torque[i] - p_data.omega[i].cross(&(global_interia_tensor * p_data.omega[i]));

//         p_data.omega[i] += dt
//             * global_interia_tensor
//                 .try_inverse()
//                 .expect("MMOI not inveratble")
//             * (temp_i);

//         p_data.update_rotation_with_omega(i, dt);
//     }
//     p_data.update_all_body_spheres_positions();
// }

pub fn inital_integrate(p_data: &mut rigid_body::RigidBodiesData, dt: f64) {
    for i in 0..p_data.body_data.len() {
        p_data.velocity[i] += 0.5 * dt * p_data.force[i] / p_data.body_data[i].mass;
        p_data.position[i] += p_data.velocity[i] * dt;

        p_data.angular_moment[i] += 0.5 * dt * p_data.torque[i];

        // println!("am{}", p_data.angular_moment[i][2]);
        p_data.omega[i] = angmom_to_omega(
            p_data.angular_moment[i],
            p_data.ex_space[i],
            p_data.ey_space[i],
            p_data.ez_space[i],
            p_data.diagonal_inertia[i],
            p_data.omega[i],
        );

        // let global_interia_tensor = p_data.quaternion[i].to_rotation_matrix()
        //     * p_data.moment_of_inertia_tensor[i]
        //     * p_data.quaternion[i].to_rotation_matrix().transpose();

        // p_data.omega[i] =
        //     global_interia_tensor.try_inverse().expect("cant inverse") * p_data.angular_moment[i];

        (p_data.quaternion[i], p_data.omega[i]) = richardson(
            p_data.quaternion[i],
            p_data.angular_moment[i],
            p_data.omega[i],
            p_data.diagonal_inertia[i],
            0.5 * dt,
        );

        // if i == 0 {
        //     println!("after {}", p_data.quaternion[i]);
        // }

        p_data.update_space_exyz_with_q(i);

        // let temp_i =
        //     p_data.torque[i] - p_data.omega[i].cross(&(global_interia_tensor * p_data.omega[i]));

        // p_data.omega[i] += 0.5
        //     * dt
        //     * global_interia_tensor
        //         .try_inverse()
        //         .expect("MMOI not inveratble")
        //     * (temp_i);

        // p_data.update_rotation_with_omega(i, dt);
    }
    p_data.update_all_body_spheres_positions();
}

pub fn angmom_to_omega(
    m: Vector3<f64>,
    ex: Vector3<f64>,
    ey: Vector3<f64>,
    ez: Vector3<f64>,
    idiag: Vector3<f64>,
    mut w: Vector3<f64>,
) -> Vector3<f64> {
    let mut wbody: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);

    if idiag[0] != 0.0 {
        wbody[0] = (m[0] * ex[0] + m[1] * ex[1] + m[2] * ex[2]) / idiag[0];
        // println!("idiag{}", idiag[0]);
    }
    if idiag[1] != 0.0 {
        wbody[1] = (m[0] * ey[0] + m[1] * ey[1] + m[2] * ey[2]) / idiag[1];
        // println!("idiag{}", idiag[1]);
    }
    if idiag[2] != 0.0 {
        wbody[2] = (m[0] * ez[0] + m[1] * ez[1] + m[2] * ez[2]) / idiag[2];
        // println!("idiag{}", idiag[2]);
    }

    w[0] = wbody[0] * ex[0] + wbody[1] * ey[0] + wbody[2] * ez[0];
    w[1] = wbody[0] * ex[1] + wbody[1] * ey[1] + wbody[2] * ez[1];
    w[2] = wbody[0] * ex[2] + wbody[1] * ey[2] + wbody[2] * ez[2];

    return w;
}

fn mq_to_omega(
    m: Vector3<f64>,
    q: UnitQuaternion<f64>,
    moments: Vector3<f64>,
    mut w: Vector3<f64>,
) -> Vector3<f64> {
    let rot = q.to_rotation_matrix();

    let mut wbody = rot.transpose() * m;

    if moments[0] == 0.0 {
        wbody[0] = 0.0;
    } else {
        wbody[0] /= moments[0];
    }
    if moments[1] == 0.0 {
        wbody[1] = 0.0;
    } else {
        wbody[1] /= moments[1];
    }
    if moments[2] == 0.0 {
        wbody[2] = 0.0;
    } else {
        wbody[2] /= moments[2];
    }
    w = rot * wbody;

    return w;
}

fn richardson(
    mut q: UnitQuaternion<f64>,
    m: Vector3<f64>,
    mut w: Vector3<f64>,
    moments: Vector3<f64>,
    dtq: f64,
) -> (UnitQuaternion<f64>, Vector3<f64>) {
    // full update from dq/dt = 1/2 w q

    //     // full update from dq/dt = 1/2 w q

    //   double wq[4];
    //   MathExtra::vecquat(w,q,wq);

    let mut wq: Quaternion<f64> = Quaternion::new(
        -w[0] * q[0] - w[1] * q[1] - w[2] * q[2],
        q[3] * w[0] + w[1] * q[2] - w[2] * q[1],
        q[3] * w[1] + w[2] * q[0] - w[0] * q[2],
        q[3] * w[2] + w[0] * q[1] - w[1] * q[0],
    );
    //   double qfull[4];
    //   qfull[0] = q[0] + dtq * wq[0];
    //   qfull[1] = q[1] + dtq * wq[1];
    //   qfull[2] = q[2] + dtq * wq[2];
    //   qfull[3] = q[3] + dtq * wq[3];
    //   MathExtra::qnormalize(qfull);

    let qfull_vec = q.as_vector() + dtq * wq.as_vector();
    let qfull_non_unit = Quaternion::new(qfull_vec[3], qfull_vec[0], qfull_vec[1], qfull_vec[2]);
    let mut qfull = UnitQuaternion::from_quaternion(qfull_non_unit);

    //   // 1st half update from dq/dt = 1/2 w q

    //   double qhalf[4];
    //   qhalf[0] = q[0] + 0.5*dtq * wq[0];
    //   qhalf[1] = q[1] + 0.5*dtq * wq[1];
    //   qhalf[2] = q[2] + 0.5*dtq * wq[2];
    //   qhalf[3] = q[3] + 0.5*dtq * wq[3];
    //   MathExtra::qnormalize(qhalf);

    let qhalf_vec = q.as_vector() + 0.5 * dtq * wq.as_vector();
    let qhalf_non_unit = Quaternion::new(qhalf_vec[3], qhalf_vec[0], qhalf_vec[1], qhalf_vec[2]);
    let mut qhalf = UnitQuaternion::from_quaternion(qhalf_non_unit);

    //   // re-compute omega at 1/2 step from m at 1/2 step and q at 1/2 step
    //   // recompute wq

    //   MathExtra::mq_to_omega(m,qhalf,moments,w);
    //   MathExtra::vecquat(w,qhalf,wq);

    w = mq_to_omega(m, qhalf, moments, w);

    wq = Quaternion::new(
        -w[0] * qhalf[0] - w[1] * qhalf[1] - w[2] * qhalf[2],
        qhalf[3] * w[0] + w[1] * qhalf[2] - w[2] * qhalf[1],
        qhalf[3] * w[1] + w[2] * qhalf[0] - w[0] * qhalf[2],
        qhalf[3] * w[2] + w[0] * qhalf[1] - w[1] * qhalf[0],
    );

    //   // 2nd half update from dq/dt = 1/2 w q

    // qhalf[0] += 0.5 * dtq * wq[0];
    // qhalf[1] += 0.5 * dtq * wq[1];
    // qhalf[2] += 0.5 * dtq * wq[2];
    // qhalf[3] += 0.5 * dtq * wq[3];
    // MathExtra::qnormalize(qhalf);

    let temp = qhalf.as_vector() + 0.5 * dtq * wq.as_vector();
    let temp_non_unit = Quaternion::new(temp[3], temp[0], temp[1], temp[2]);
    qhalf = UnitQuaternion::from_quaternion(temp_non_unit);

    // // corrected Richardson update

    // q[0] = 2.0 * qhalf[0] - qfull[0];
    // q[1] = 2.0 * qhalf[1] - qfull[1];
    // q[2] = 2.0 * qhalf[2] - qfull[2];
    // q[3] = 2.0 * qhalf[3] - qfull[3];
    // MathExtra::qnormalize(q);
    let q_temp = 2.0 * qhalf.as_vector() - qfull.as_vector();
    let q_temp_non_unit = Quaternion::new(q_temp[3], q_temp[0], q_temp[1], q_temp[2]);
    let finstuff = UnitQuaternion::from_quaternion(q_temp_non_unit);

    // println!("{:?}", finstuff);
    return (finstuff, w);
}

pub fn final_integrate(p_data: &mut rigid_body::RigidBodiesData, dt: f64) {
    for i in 0..p_data.body_data.len() {
        p_data.velocity[i] += 0.5 * dt * p_data.force[i] / p_data.body_data[i].mass;

        p_data.angular_moment[i] += 0.5 * dt * p_data.torque[i];

        p_data.omega[i] = angmom_to_omega(
            p_data.angular_moment[i],
            p_data.ex_space[i],
            p_data.ey_space[i],
            p_data.ez_space[i],
            p_data.diagonal_inertia[i],
            p_data.omega[i],
        );

        // let global_interia_tensor = p_data.quaternion[i].to_rotation_matrix()
        //     * p_data.moment_of_inertia_tensor[i]
        //     * p_data.quaternion[i].to_rotation_matrix().transpose();

        // p_data.omega[i] =
        //     global_interia_tensor.try_inverse().expect("cant inverse") * p_data.angular_moment[i];

        // let global_interia_tensor = p_data.quaternion[i].to_rotation_matrix()
        //     * p_data.moment_of_inertia_tensor[i]
        //     * p_data.quaternion[i].to_rotation_matrix().transpose();

        // let temp_i =
        //     p_data.torque[i] - p_data.omega[i].cross(&(global_interia_tensor * p_data.omega[i]));

        // p_data.omega[i] += 0.5
        //     * dt
        //     * global_interia_tensor
        //         .try_inverse()
        //         .expect("MMOI not inveratble")
        //     * (temp_i);
    }
}

pub fn lees_edwards_boundaries(
    d_data: &domain::DomainData,
    p_data: &mut rigid_body::RigidBodiesData,
    _dt: f64,
    ledisplace: f64,
) {
    for i in 0..p_data.body_data.len() {
        if p_data.position[i][1] > d_data.domain[1] {
            // p_data.position[i][1] = d_data.domain[1] - p_data.radius[i]
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
            // p_data.position[i][1] = p_data.radius[i]
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

        p_data.update_one_body_spheres_positions(i);

        //Check Individual Spheres and apply a wrap int to them
        for j in 0..p_data.body_data[i].position.len() {
            if p_data.body_data[i].position[j][1] > d_data.domain[1] {
                p_data.body_data[i].position[j][1] -= d_data.domain[1];
                p_data.body_data[i].wrap[j][1] += -1;

                p_data.body_data[i].velocity[j][0] -=
                    d_data.lees_edwards_boundary * d_data.domain[1];
                p_data.body_data[i].position[j][0] -= ledisplace;
                if p_data.body_data[i].position[j][0] <= 0.0 {
                    p_data.body_data[i].position[j][0] += d_data.domain[0];
                    p_data.body_data[i].wrap[j][0] += 1
                }
            }
            // if particle is less than domain move to end of domain
            // Also apply velocity change to particle for shearing
            else if p_data.body_data[i].position[j][1] <= 0.0 {
                p_data.body_data[i].position[j][1] += d_data.domain[1];
                p_data.body_data[i].wrap[j][1] += 1;
                p_data.body_data[i].velocity[j][0] +=
                    d_data.lees_edwards_boundary * d_data.domain[1];
                p_data.body_data[i].position[j][0] += ledisplace;
                if p_data.body_data[i].position[j][0] > d_data.domain[0] {
                    p_data.body_data[i].position[j][0] -= d_data.domain[0];
                    p_data.body_data[i].wrap[j][0] += -1;
                }
            }
            // std::cout << distb << std::endl;//X boundary condition
            // if particles is greater than domain move to beginning of domain
            if p_data.body_data[i].position[j][0] > d_data.domain[0] {
                p_data.body_data[i].position[j][0] -= d_data.domain[0];
                p_data.body_data[i].wrap[j][0] += -1;
            }
            // if particle is less than domain move to end of domain
            else if p_data.body_data[i].position[j][0] <= 0.0 {
                p_data.body_data[i].position[j][0] += d_data.domain[0];
                p_data.body_data[i].wrap[j][0] += 1;
            }

            // Z boundary condition
            // if particles is greater than domain move to beginning of domain
            if p_data.body_data[i].position[j][2] > d_data.domain[2] {
                p_data.body_data[i].position[j][2] -= d_data.domain[2];
                p_data.body_data[i].wrap[j][2] += -1;
            }
            // if particle is less than domain move to end of domain
            else if p_data.body_data[i].position[j][2] <= 0.0 {
                p_data.body_data[i].position[j][2] += d_data.domain[2];
                p_data.body_data[i].wrap[j][2] += 1;
            }
        }
    }
}

pub fn simple_boundary(d_data: &mut domain::DomainData, p_data: &mut rigid_body::RigidBodiesData) {
    for i in 0..p_data.body_data.len() {
        if p_data.position[i][1] > d_data.domain[1] {
            p_data.position[i][1] = d_data.domain[1];
            p_data.velocity[i][1] = 0.0;
        }
        // if particle is less than domain move to end of domain
        // Also apply velocity change to particle for shearing
        else if p_data.position[i][1] <= 0.0 {
            p_data.position[i][1] = 0.0;
            p_data.velocity[i][1] = 0.0;
        }
        // std::cout << distb << std::endl;//X boundary condition
        // if particles is greater than domain move to beginning of domain
        if p_data.position[i][0] > d_data.domain[0] {
            p_data.position[i][0] = d_data.domain[0];
            p_data.velocity[i][0] = 0.0;
        }
        // if particle is less than domain move to end of domain
        else if p_data.position[i][0] <= 0.0 {
            p_data.position[i][0] = 0.0;
            p_data.velocity[i][0] = 0.0;
        }

        // Z boundary condition
        // if particles is greater than domain move to beginning of domain
        if p_data.position[i][2] > d_data.domain[2] {
            p_data.position[i][2] = d_data.domain[2];
            p_data.velocity[i][2] = 0.0;
        }
        // if particle is less than domain move to end of domain
        else if p_data.position[i][2] <= 0.0 {
            p_data.position[i][2] = 0.0;
            p_data.velocity[i][2] = 0.0;
        }

        p_data.update_one_body_spheres_positions(i);
    }
}
