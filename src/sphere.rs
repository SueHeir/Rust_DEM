use nalgebra::{Matrix3, Quaternion, UnitQuaternion, UnitVector3, Vector2, Vector3, Vector4};
use std::collections::HashMap;

macro_rules! pub_struct {
    ($name:ident {$($field:ident: $t:ty,)*}) => {
        #[derive(Debug, Clone, PartialEq)] // ewww
        pub struct $name {
            $(pub $field: $t),*
        }
    }
}

pub_struct!(Material {
    radius: f64,
    mass: f64,
    youngs_mod: f64,
    poisson_ratio: f64,
    density: f64,
    id: i32,
});

pub_struct!(EffectiveMaterialPreCalc {
    eff_radius: f64,
    eff_youngs_mod: f64,
});

pub_struct!( ParticleData {

    radius: Vec<f64>,
    max_radius: f64,
    mass: Vec<f64>,
    youngs_mod: Vec<f64>,
    poisson_ratio: Vec<f64>,
    density: Vec<f64>,
    position: Vec<Vector3<f64>>,
    velocity: Vec<Vector3<f64>>,
    force: Vec<Vector3<f64>>,




    is_collision: Vec<bool>,

    materials: Vec<Material>,

    sphere_material: Vec<usize>,

    sphere_material_map: HashMap<String,EffectiveMaterialPreCalc>,

    restitution_coefficient: f64,
    beta: f64,
    friction: f64,
    volume_fraction: f64,

    quaternion: Vec<UnitQuaternion<f64>>,
    omega: Vec<Vector3<f64>>,
    torque: Vec<Vector3<f64>>,
    angular_moment: Vec<Vector3<f64>>,
    ex_space: Vec<Vector3<f64>>,
    ey_space: Vec<Vector3<f64>>,
    ez_space: Vec<Vector3<f64>>,
    diagonal_inertia: Vec<Vector3<f64>>,


});

impl ParticleData {
    pub fn exyz_to_q(
        &self,
        ex: Vector3<f64>,
        ey: Vector3<f64>,
        ez: Vector3<f64>,
    ) -> UnitQuaternion<f64> {
        // squares of quaternion components

        let q0sq = 0.25 * (ex[0] + ey[1] + ez[2] + 1.0);
        let q1sq = q0sq - 0.5 * (ey[1] + ez[2]);
        let q2sq = q0sq - 0.5 * (ex[0] + ez[2]);
        let q3sq = q0sq - 0.5 * (ex[0] + ey[1]);

        // some component must be greater than 1/4 since they sum to 1
        // compute other components from it

        let mut q = [0.0, 0.0, 0.0, 0.0];

        if q0sq >= 0.25 {
            q[0] = q0sq.sqrt();
            q[1] = (ey[2] - ez[1]) / (4.0 * q[0]);
            q[2] = (ez[0] - ex[2]) / (4.0 * q[0]);
            q[3] = (ex[1] - ey[0]) / (4.0 * q[0]);
        } else if q1sq >= 0.25 {
            q[1] = q1sq.sqrt();
            q[0] = (ey[2] - ez[1]) / (4.0 * q[1]);
            q[2] = (ey[0] + ex[1]) / (4.0 * q[1]);
            q[3] = (ex[2] + ez[0]) / (4.0 * q[1]);
        } else if q2sq >= 0.25 {
            q[2] = q2sq.sqrt();
            q[0] = (ez[0] - ex[2]) / (4.0 * q[2]);
            q[1] = (ey[0] + ex[1]) / (4.0 * q[2]);
            q[3] = (ez[1] + ey[2]) / (4.0 * q[2]);
        } else if q3sq >= 0.25 {
            q[3] = q3sq.sqrt();
            q[0] = (ex[1] - ey[0]) / (4.0 * q[3]);
            q[1] = (ez[0] + ex[2]) / (4.0 * q[3]);
            q[2] = (ez[1] + ey[2]) / (4.0 * q[3]);
        }
        let quat = UnitQuaternion::new_normalize(Quaternion::from_parts(
            q[0],
            Vector3::new(q[1], q[2], q[3]),
        ));
        return quat;
    }

    pub fn update_all_rotation_with_omega(&mut self, step_size: f64) {
        for i in 0..self.omega.len() {
            self.update_rotation_with_omega(i, step_size);
        }
    }

    pub fn update_rotation_with_omega(&mut self, index: usize, step_size: f64) {
        let step_angle = step_size * self.omega[index].norm();

        if step_angle == 0.0 {
            return;
        }

        let step_rotation_axis: UnitVector3<f64> = UnitVector3::new_normalize(self.omega[index]);

        let rotation_quaternion = UnitQuaternion::from_axis_angle(&step_rotation_axis, step_angle);

        self.quaternion[index] = rotation_quaternion * self.quaternion[index];
    }

    pub fn update_all_space_exyz_with_q(&mut self) {
        for i in 0..self.ex_space.len() {
            self.update_space_exyz_with_q(i);
        }
    }

    pub fn update_space_exyz_with_q(&mut self, index: usize) {
        let qqq = self.quaternion[index].as_vector();
        let q = Vector4::new(qqq[3], qqq[0], qqq[1], qqq[2]);

        self.ex_space[index][0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
        self.ex_space[index][1] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
        self.ex_space[index][2] = 2.0 * (q[1] * q[3] - q[0] * q[2]);

        self.ey_space[index][0] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
        self.ey_space[index][1] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
        self.ey_space[index][2] = 2.0 * (q[2] * q[3] + q[0] * q[1]);

        self.ez_space[index][0] = 2.0 * (q[1] * q[3] + q[0] * q[2]);
        self.ez_space[index][1] = 2.0 * (q[2] * q[3] - q[0] * q[1]);
        self.ez_space[index][2] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
        // let binding = self.quaternion[index].to_rotation_matrix();
        // let r = binding.matrix();

        // self.ex_space[index] = Vector3::new(r[0], r[3], r[6]);
        // self.ey_space[index] = Vector3::new(r[1], r[4], r[7]);
        // self.ez_space[index] = Vector3::new(r[2], r[5], r[8]);
    }
}

pub_struct!( ForceData {
    particle_indexes: Vec<Vector2<usize>>,
    force: Vec<Vector3<f64>>,
    del: Vec<Vector3<f64>>,
    forcedata:  Vec<Matrix3<f64>>,
});

pub_struct!(BondData {
    //* BONd nstiff tstiff Mnstiff Mtstiff  ibond sigma_max tau_max on/off
    //bond 4.3e10 1.6e10 4.3e10 1.6e10 1 2.2e8 2.2e8 on //bdp 0.8763 0.8763 0.8763 0.8763
    n_stiff: f64,
    t_stiff: f64,
    mn_stiff: f64,
    mt_stiff: f64,
    sigma_max: f64,
    tau_max: f64,
    //* Bond damping 'BDP e_n e_t e_tor e_bend'
    //bdp 0.8763 0.8763 0.8763 0.8763
    e_n: f64,
    e_t: f64,
    e_tor: f64,
    e_bend: f64,

    normal_force_sum: Vector3<f64>,
    tangential_force_sum: Vector3<f64>,
    normal_moment_sum: Vector3<f64>,
    tangential_moment_sum: Vector3<f64>,

    incremental_normal_moment: Vector3<f64>,
    incremental_tangential_moment: Vector3<f64>,

    potential_energy_comp_ext: f64,
    potential_energy_tang_deform: f64,
    potential_energy_rot_deform: f64,
    pontential_energy_bend_deform: f64,

    incremental_tangential_force: Vector3<f64>,

    resistance_normal_force: f64,

    current_sigma: f64,
    current_tau: f64,

    broken: bool,


});
