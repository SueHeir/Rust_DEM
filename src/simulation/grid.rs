use core::prelude::v1;
use std::collections::HashMap;
use std::f64::consts::PI;

use nalgebra::Quaternion;
use nalgebra::UnitQuaternion;
use nalgebra::Vector3;

use crate::domain;
use crate::sphere;

pub fn simp_collisions(
    p_data: &mut sphere::ParticleData,
    bond_data: &mut HashMap<(usize, usize), sphere::BondData>,
    _dt: f64,
) {
    for i in 0..p_data.radius.len() {
        for j in i + 1..p_data.radius.len() {
            if bond_data.contains_key(&(i, j)) {
                return;
            }
            let p1 = p_data.position[i];
            let p2 = p_data.position[j];
            let v1 = p_data.velocity[i];
            let v2 = p_data.velocity[j];
            let r1 = p_data.radius[i];
            let r2 = p_data.radius[j];

            let delta_position = p2 - p1;

            let distance = delta_position.norm();

            if distance < r1 + r2 {
                p_data.is_collision[i] = true;
                p_data.is_collision[j] = true;

                let normalized_delta = delta_position / distance;

                let distance_delta = (p_data.radius[i] + p_data.radius[j]) - distance;

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
                let v_r_n = f_dot * normalized_delta;

                let dissipation_force = 2.0
                    * 0.91287092917
                    * p_data.beta
                    * (contact_stiffness * reduced_mass).sqrt()
                    * v_r_n.norm()
                    * v_r_n.dot(&normalized_delta).signum();

                p_data.force[i] -= (normal_force - dissipation_force) * normalized_delta;
                p_data.force[j] += (normal_force - dissipation_force) * normalized_delta;

                // println!("Collision")
            }
        }
    }
}

pub fn simp_bonds(
    p_data: &mut sphere::ParticleData,
    bond_data: &mut HashMap<(usize, usize), sphere::BondData>,
    _dt: f64,
) {
    for ((i, j), b_data) in bond_data.iter_mut() {
        let x_p1 = p_data.position[*i];
        let x_p2 = p_data.position[*j];

        let x_rel = x_p2 - x_p1;

        let vel1 = p_data.velocity[*i];
        let vel2 = p_data.velocity[*j];

        let wrot1 = p_data.omega[*i];
        let wrot2 = p_data.omega[*j];

        let w_rel = wrot2 - wrot1;

        let nom = x_rel.magnitude();

        let n12 = x_rel.normalize();

        let r1_rel = 0.5 * nom * n12;
        let r2_rel = -0.5 * nom * n12;

        let v1_w = wrot1.cross(&r1_rel);
        let v2_w = wrot2.cross(&r2_rel);

        let vc1 = vel1 + v1_w;
        let vc2 = vel2 + v2_w;

        let vn12 = (vc2 - vc1).dot(&n12);

        let vn_rel = vn12 * n12;
        let vt_rel = (vc2 - vc1) - vn_rel;

        let wn12 = w_rel.dot(&n12);

        let wn_rel = wn12 * n12;

        let wt_rel = w_rel - wn_rel;

        let a_bnd = PI * p_data.radius[*i].powi(2);
        let j_bnd = 0.5 * PI * p_data.radius[*i].powi(4);

        // let estart = 1.0
        //     / ((1.0 - p_data.poisson_ratio[*i] * p_data.poisson_ratio[*i]) / p_data.youngs_mod[*i]
        //         + (1.0 - p_data.poisson_ratio[*j] * p_data.poisson_ratio[*j])
        //             / p_data.youngs_mod[*j]);
        // let rstar = 1.0 / (1.0 / p_data.radius[*i] + 1.0 / p_data.radius[*j]);

        let lbond = p_data.radius[*i] + p_data.radius[*j];

        let ovlp = nom - lbond;

        let tem_ft = b_data.tangential_force_sum + b_data.incremental_tangential_force;
        let tem_mn = b_data.normal_moment_sum + b_data.incremental_normal_moment;
        let tem_mt = b_data.tangential_moment_sum + b_data.incremental_tangential_moment;

        let len_ft = tem_ft.magnitude();
        let len_mn = tem_mn.magnitude();
        let len_mt = tem_mt.magnitude();

        let norm_ft = tem_ft.dot(&n12);
        let norm_mn = tem_mn.dot(&n12);
        let norm_mt = tem_mt.dot(&n12);

        let tem_ft1 = tem_ft - norm_ft * n12;
        let tem_mn1 = norm_mn * n12;
        let tem_mt1 = tem_mt - norm_mt * n12;

        let len_ft1 = tem_ft1.magnitude();
        let len_mn1 = tem_mn1.magnitude();
        let len_mt1 = tem_mt1.magnitude();

        if len_ft1 != 0.0 {
            b_data.tangential_force_sum = len_ft / len_ft1 * tem_ft1;
            b_data.incremental_tangential_force = Vector3::zeros();
        }

        if len_mn1 != 0.0 {
            b_data.normal_moment_sum = len_mn / len_mn1 * tem_mn1;
            b_data.incremental_normal_moment = Vector3::zeros();
        }

        if len_mt1 != 0.0 {
            b_data.tangential_moment_sum = len_mt / len_mt1 * tem_mt1;
            b_data.incremental_tangential_moment = Vector3::zeros();
        }

        let fbnd_norm = b_data.n_stiff * a_bnd / lbond * ovlp;

        let thres = p_data.radius[*i];
        let mut thres3 = thres * 0.0001;
        let thres1 = thres * (b_data.mt_stiff * j_bnd / lbond);
        let thres2 = thres * (b_data.mn_stiff * (0.5 * j_bnd) / lbond);
        thres3 = thres3 * b_data.t_stiff * a_bnd / lbond;

        b_data.normal_force_sum = fbnd_norm * n12;

        b_data.incremental_tangential_force += b_data.t_stiff * a_bnd / lbond * vt_rel * _dt;
        if b_data.incremental_tangential_force.magnitude() > thres3 {
            b_data.tangential_force_sum += b_data.incremental_tangential_force;
            b_data.incremental_tangential_force = Vector3::zeros();
        }

        b_data.incremental_normal_moment += b_data.mt_stiff * j_bnd / lbond * wn_rel * _dt;
        if b_data.incremental_normal_moment.magnitude() > thres1 {
            b_data.normal_moment_sum += b_data.incremental_normal_moment;
            b_data.incremental_normal_moment = Vector3::zeros();
        }

        b_data.incremental_tangential_moment +=
            b_data.mn_stiff * (0.5 * j_bnd) / lbond * wt_rel * _dt;
        if b_data.incremental_tangential_moment.magnitude() > thres2 {
            b_data.tangential_moment_sum += b_data.incremental_tangential_moment;
            b_data.incremental_tangential_moment = Vector3::zeros();
        }

        let ftan = b_data.tangential_force_sum + b_data.incremental_tangential_force;
        let mom_nor = b_data.normal_moment_sum + b_data.incremental_normal_moment;
        let mom_tan = b_data.tangential_moment_sum + b_data.incremental_tangential_moment;

        b_data.potential_energy_comp_ext =
            lbond / (2.0 * b_data.n_stiff * a_bnd) * fbnd_norm.powi(2);

        if b_data.t_stiff != 0.0 {
            b_data.potential_energy_tang_deform =
                lbond / (2.0 * b_data.t_stiff * a_bnd) * ftan.magnitude_squared();
        } else {
            b_data.potential_energy_tang_deform = 0.0;
        }

        if b_data.mt_stiff != 0.0 {
            b_data.potential_energy_rot_deform =
                lbond / (2.0 * b_data.mt_stiff * j_bnd) * mom_nor.magnitude_squared();
        } else {
            b_data.potential_energy_rot_deform = 0.0;
        }

        if b_data.mn_stiff != 0.0 {
            b_data.pontential_energy_bend_deform =
                lbond / (2.0 * b_data.mn_stiff * (0.5 * j_bnd)) * mom_tan.magnitude_squared();
        } else {
            b_data.pontential_energy_bend_deform = 0.0;
        }

        let mstar = 0.5 * p_data.mass[*i];
        let moi_star = 0.5 * p_data.diagonal_inertia[*i][0];

        let damp_n = 2.0 * b_data.e_n * (mstar * b_data.n_stiff * a_bnd / lbond).sqrt();
        let damp_t = 2.0 * b_data.e_t * (mstar * b_data.t_stiff * a_bnd / lbond).sqrt();
        let damp_tor = 2.0 * b_data.e_tor * (moi_star * b_data.mt_stiff * j_bnd / lbond).sqrt();
        let damp_bend =
            2.0 * b_data.e_bend * (moi_star * b_data.mn_stiff * (0.5 * j_bnd) / lbond).sqrt();

        let fdamp_n = damp_n * vn_rel;
        let fdamp_t = damp_t * vt_rel;
        let mdamp_n = damp_tor * wn_rel;
        let mdamp_t = damp_bend * wt_rel;

        let fn_bnd = fbnd_norm;

        let ft_bnd = ftan.magnitude();
        let mn_bnd = mom_nor.magnitude();
        let mt_bnd = mom_tan.magnitude();

        b_data.current_sigma = fn_bnd / a_bnd + (2.0 * mt_bnd.abs() * p_data.radius[*i]) / j_bnd;
        b_data.current_tau = ft_bnd.abs() / a_bnd + mn_bnd.abs() * p_data.radius[*i] / j_bnd;

        if b_data.current_sigma > b_data.sigma_max || b_data.current_tau > b_data.tau_max {
            b_data.broken = true;
            println!("Broken Bond");
            if b_data.current_sigma > b_data.sigma_max {
                println!("because sigma")
            }
            if b_data.current_tau > b_data.tau_max {
                println!("because tau")
            }
            if ovlp < 0.0 {
                p_data.force[*i] += fbnd_norm * n12;
                p_data.force[*j] -= fbnd_norm * n12;
            }
        }

        let ftot_tan = ftan + fdamp_t;

        let torq1 = r1_rel.cross(&ftot_tan);
        let torq2 = r2_rel.cross(&-ftot_tan);

        p_data.force[*i] += b_data.normal_force_sum + ftan + fdamp_n + fdamp_t;
        p_data.force[*j] -= b_data.normal_force_sum + ftan + fdamp_n + fdamp_t;

        // println!("{:?}", b_data.normal_force_sum + ftan + fdamp_n + fdamp_t);

        p_data.torque[*i] += mom_nor + mom_tan + torq1 + mdamp_n + mdamp_t;
        p_data.torque[*j] += -mom_nor - mom_tan + torq2 - mdamp_n - mdamp_t;

        // println!("{:?}", mom_nor + mom_tan + torq1 + mdamp_n + mdamp_t);
    }
}

pub fn inital_integrate(p_data: &mut sphere::ParticleData, dt: f64) {
    for i in 0..p_data.radius.len() {
        p_data.velocity[i] += 0.5 * dt * p_data.force[i] / p_data.mass[i];
        p_data.position[i] += p_data.velocity[i] * dt;

        p_data.angular_moment[i] += 0.5 * dt * p_data.torque[i];

        p_data.omega[i] = angmom_to_omega(
            p_data.angular_moment[i],
            p_data.ex_space[i],
            p_data.ey_space[i],
            p_data.ez_space[i],
            p_data.diagonal_inertia[i],
            p_data.omega[i],
        );

        (p_data.quaternion[i], p_data.omega[i]) = richardson(
            p_data.quaternion[i],
            p_data.angular_moment[i],
            p_data.omega[i],
            p_data.diagonal_inertia[i],
            0.5 * dt,
        );

        p_data.update_space_exyz_with_q(i);
    }
}

pub fn final_integrate(p_data: &mut sphere::ParticleData, dt: f64) {
    for i in 0..p_data.radius.len() {
        p_data.velocity[i] += 0.5 * dt * p_data.force[i] / p_data.mass[i];

        p_data.angular_moment[i] += 0.5 * dt * p_data.torque[i];

        p_data.omega[i] = angmom_to_omega(
            p_data.angular_moment[i],
            p_data.ex_space[i],
            p_data.ey_space[i],
            p_data.ez_space[i],
            p_data.diagonal_inertia[i],
            p_data.omega[i],
        );
    }
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

fn mq_to_omega(m: Vector3<f64>, q: UnitQuaternion<f64>, moments: Vector3<f64>) -> Vector3<f64> {
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
    let w = rot * wbody;

    return w;
}

fn richardson(
    q: UnitQuaternion<f64>,
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
    let qfull = UnitQuaternion::from_quaternion(qfull_non_unit);

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

    w = mq_to_omega(m, qhalf, moments);

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
