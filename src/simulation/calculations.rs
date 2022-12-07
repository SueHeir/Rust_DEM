use crate::{domain, rigid_body};
use nalgebra::{Matrix3, Vector3};

// pub fn _get_kinetic_energy(p_data: &mut rigid_body::RigidBodiesData) -> f64 {
//     let mut sum = 0.0;
//     // Energy is defined as 1/2 m * velocity dot velocity
//     for i in 0..p_data.radius.len() {
//         sum += 1.0 / 2.0 * p_data.mass[i] * p_data.velocity[i].dot(&p_data.velocity[i]);
//     }
//     return sum;
// }

pub fn calc_kinetic_tensor(
    p_data: &rigid_body::RigidBodiesData,
    d_data: &domain::DomainData,
    kinetic_tensor: Matrix3<f64>,
    average_reset_count: i32,
) -> Matrix3<f64> {
    let n_particles = p_data.body_data.len();
    // Get average velocity
    let mut average_velocity = Vector3::new(0.0, 0.0, 0.0);
    for i in 0..p_data.body_data.len() {
        average_velocity += p_data.velocity[i];
    }
    average_velocity = average_velocity.scale(1.0 / n_particles as f64);

    // println!("{:?}", average_velocity);

    // Get average kinetic tensor
    let mut temp_kinetic_tensor = Matrix3::zeros();

    for i in 0..p_data.body_data.len() {
        temp_kinetic_tensor += p_data.body_data[i].mass
            * ((p_data.velocity[i] - average_velocity)
                * (p_data.velocity[i] - average_velocity).transpose());
    }
    temp_kinetic_tensor = temp_kinetic_tensor.scale(1.0 / (d_data.domain_volume));

    // println!("{:?}", temp_kinetic_tensor);

    // average this frames kinetic tensor with prevous tensors
    let kt: Matrix3<f64> = (kinetic_tensor.scale(average_reset_count as f64) + temp_kinetic_tensor)
        .scale(1.0 / (average_reset_count + 1) as f64) as Matrix3<f64>;
    // println!("{:?}", kt);

    return kt;
}

pub fn calc_collision_tensor(
    f_data: &rigid_body::ForceData,
    d_data: &domain::DomainData,
    collision_tensor: Matrix3<f64>,
    average_reset_count: i32,
) -> Matrix3<f64> {
    let mut temp_collision_tensor = Matrix3::zeros();

    for i in 0..f_data.forcedata.len() {
        temp_collision_tensor += f_data.forcedata[i];
    }
    temp_collision_tensor = temp_collision_tensor.scale(1.0 / d_data.domain_volume);

    // average this frames kinetic tensor with prevous tensors
    let ct: Matrix3<f64> = (collision_tensor.scale(average_reset_count as f64)
        + temp_collision_tensor)
        .scale(1.0 / (average_reset_count + 1) as f64) as Matrix3<f64>;
    // println!("{:?}", kt);

    return ct;
}
