use nalgebra::{
    Matrix3, Matrix4, Quaternion, UnitQuaternion, UnitVector3, Vector2, Vector3, Vector4,
};
use rand::Rng;
use std::{collections::HashMap, process::exit};

macro_rules! pub_struct {
    ($name:ident {$($field:ident: $t:ty,)*}) => {
        #[derive(Debug, Clone, PartialEq)] // ewww
        pub struct $name {
            $(pub $field: $t),*
        }
    }
}

pub_struct!(Material {
    youngs_mod: f64,
    poisson_ratio: f64,
    density: f64,
    id: i32,
});

#[derive(Debug, Clone, PartialEq)]
pub struct RigidBodySphereData {
    pub radius: Vec<f64>,
    pub position: Vec<Vector3<f64>>,
    pub local_position: Vec<Vector3<f64>>,
    pub velocity: Vec<Vector3<f64>>,

    pub wrap: Vec<Vector3<i8>>,

    pub density: f64,
    pub mass: f64,
    pub volume: f64,

    pub local_inertia_tensor: Matrix3<f64>,
    pub diagonal_inertia: Vector3<f64>,
    pub ex_space: Vector3<f64>,
    pub ey_space: Vector3<f64>,
    pub ez_space: Vector3<f64>,
}

impl RigidBodySphereData {
    //You need position and radius already loaded to call this
    pub fn initialize(&mut self, density: f64) {
        let repeat_count = 100_000_000;
        let mut rng = rand::thread_rng();

        let (aabb_min, aabb_max) = self.get_aabb();

        let length = aabb_max - aabb_min;

        let mut points_inside_particle: Vec<Vector3<f64>> = Vec::new();
        for i in 0..repeat_count {
            let rand_value = Vector3::new(rng.gen::<f64>(), rng.gen::<f64>(), rng.gen::<f64>());
            let rand_position = aabb_min
                + Vector3::new(
                    length[0] * rand_value[0],
                    length[1] * rand_value[1],
                    length[2] * rand_value[2],
                );

            for (j, pos) in self.position.iter().enumerate() {
                let distance = rand_position - pos;
                if distance.norm() < self.radius[j] {
                    points_inside_particle.push(rand_position);
                    break;
                }
            }
        }

        let mut center_of_mass: Vector3<f64> = Vector3::zeros();
        for pos in points_inside_particle.iter() {
            center_of_mass += pos;
        }
        center_of_mass = center_of_mass / points_inside_particle.len() as f64;

        println!("Center of mass {}", center_of_mass);

        for pos in self.position.iter_mut() {
            *pos -= center_of_mass;
        }
        for pos in self.local_position.iter_mut() {
            *pos -= center_of_mass;
        }

        println!("Local Position {:?}", self.local_position);
        println!("Position {:?}", self.position);

        self.volume = (points_inside_particle.len() as f64) / (repeat_count as f64)
            * length[0]
            * length[1]
            * length[2];

        println!(
            "Volume {:?} {} {} {}",
            self.volume, length[0], length[1], length[2]
        );

        self.density = density;

        self.mass = self.volume * density;

        println!("Mass {:?}", self.mass);

        let point_mass = self.mass / points_inside_particle.len() as f64;

        for posi in points_inside_particle.iter() {
            let pos = posi - center_of_mass;

            self.local_inertia_tensor[0] += point_mass * (pos[1].powi(2) + pos[2].powi(2));
            self.local_inertia_tensor[4] += point_mass * (pos[0].powi(2) + pos[2].powi(2));
            self.local_inertia_tensor[8] += point_mass * (pos[0].powi(2) + pos[1].powi(2));
            self.local_inertia_tensor[1] += -point_mass * pos[0] * pos[1];
            self.local_inertia_tensor[3] += -point_mass * pos[0] * pos[1];
            self.local_inertia_tensor[5] += -point_mass * pos[1] * pos[2];
            self.local_inertia_tensor[7] += -point_mass * pos[1] * pos[2];
            self.local_inertia_tensor[6] += -point_mass * pos[0] * pos[2];
            self.local_inertia_tensor[2] += -point_mass * pos[0] * pos[2];
        }

        println!("Inertia Tensor {:?}", self.local_inertia_tensor);

        self.diagonal_inertia = self.local_inertia_tensor.diagonal();
    }

    pub fn get_aabb(&self) -> (Vector3<f64>, Vector3<f64>) {
        let mut aabb_min = Vector3::new(
            self.position[0][0],
            self.position[0][1],
            self.position[0][2],
        );
        let mut aabb_max = Vector3::new(
            self.position[0][0],
            self.position[0][1],
            self.position[0][2],
        );
        for (i, pos) in self.position.iter().enumerate() {
            if aabb_min[0] > pos[0] - self.radius[i] {
                aabb_min[0] = pos[0] - self.radius[i]
            }
            if aabb_min[1] > pos[1] - self.radius[i] {
                aabb_min[1] = pos[1] - self.radius[i]
            }
            if aabb_min[2] > pos[2] - self.radius[i] {
                aabb_min[2] = pos[2] - self.radius[i]
            }

            if aabb_max[0] < pos[0] + self.radius[i] {
                aabb_max[0] = pos[0] + self.radius[i]
            }
            if aabb_max[1] < pos[1] + self.radius[i] {
                aabb_max[1] = pos[1] + self.radius[i]
            }
            if aabb_max[2] < pos[2] + self.radius[i] {
                aabb_max[2] = pos[2] + self.radius[i]
            }
        }

        return (aabb_min, aabb_max);
    }

    pub fn calculate_eigen_system(&mut self) {
        //-----------------------------------------
        // calculate eigenvalues and eigenvectors
        //-----------------------------------------

        let (ierror, evectors, evalues) = self.get_jacobi();
        self.diagonal_inertia = evalues;
        if ierror {
            println!("Insufficient Jacobi rotations for rigid body");
            exit(1);
        }

        println!("{:?}", self.diagonal_inertia);
        self.diagonal_inertia = evalues;

        self.ex_space[0] = *evectors.index((0, 0));
        self.ex_space[1] = *evectors.index((1, 0));
        self.ex_space[2] = *evectors.index((2, 0));
        self.ey_space[0] = *evectors.index((0, 1));
        self.ey_space[1] = *evectors.index((1, 1));
        self.ey_space[2] = *evectors.index((2, 1));
        self.ez_space[0] = *evectors.index((0, 2));
        self.ez_space[1] = *evectors.index((1, 2));
        self.ez_space[2] = *evectors.index((2, 2));

        // if any principal moment < scaled EPSILON, set to 0.0

        let mut max = self.diagonal_inertia.max();
        const EPSILON: f64 = 1.0e-7;

        if self.diagonal_inertia[0] < EPSILON * max {
            self.diagonal_inertia[0] = 0.0;
        }
        if self.diagonal_inertia[1] < EPSILON * max {
            self.diagonal_inertia[1] = 0.0;
        }
        if self.diagonal_inertia[2] < EPSILON * max {
            self.diagonal_inertia[2] = 0.0;
        }

        // enforce 3 evectors as a right-handed coordinate system
        // flip 3rd evector if needed

        // double ez[3];
        // vectorCross3D(ex_space_, ey_space_, ez);
        let ez = self.ex_space.cross(&self.ey_space);
        let dot = ez.dot(&self.ez_space);
        if dot < 0.0 {
            self.ez_space = self.ez_space.scale(-1.0);
        }

        println!(
            "ex:{} ey:{} ez:{}",
            self.ex_space, self.ey_space, self.ez_space
        );

        println!("{:?}", self.diagonal_inertia);
    }

    pub fn get_jacobi(&mut self) -> (bool, Matrix3<f64>, Vector3<f64>) {
        const MAXJACOBI: i32 = 100;

        //    i,j,k;
        let mut tresh;
        let mut tau;
        let mut theta;
        let mut t;
        let mut sm;
        let mut s;
        let mut h;
        let mut g;
        let mut c;

        let mut evectors = Matrix3::identity();
        let mut matrix = self.local_inertia_tensor;

        let mut evalues = self.diagonal_inertia;
        let mut b: Vector3<f64> = self.diagonal_inertia;

        let mut z: Vector3<f64> = Vector3::zeros();

        for iter in 1..MAXJACOBI {
            sm = 0.0;
            for i in 0..2 {
                for j in i + 1..3 {
                    sm += matrix.index((i, j)).abs();
                }
            }
            if (sm == 0.0) {
                return (false, evectors, evalues);
            }

            if iter < 4 {
                tresh = 0.2 * sm / (3.0 * 3.0);
            }
            //change 3 to 3.0
            else {
                tresh = 0.0;
            }

            for i in 0..2 {
                for j in (i + 1)..3 {
                    g = 100.0 * (matrix.index((i, j))).abs();
                    if (iter > 4
                        && (evalues[i]).abs() + g == (evalues[i]).abs()
                        && (evalues[j]).abs() + g == (evalues[j]).abs())
                    {
                        matrix[i + j * 3] = 0.0;
                    } else if ((matrix.index((i, j))).abs() > tresh) {
                        h = evalues[j] - evalues[i];
                        if ((h).abs() + g == (h).abs()) {
                            t = (matrix.index((i, j))) / h;
                        } else {
                            theta = 0.5 * h / (matrix.index((i, j)));
                            t = 1.0 / ((theta).abs() + (1.0 + theta * theta).sqrt());
                            if (theta < 0.0) {
                                t = -t;
                            }
                        }
                        c = 1.0 / (1.0 + t * t).sqrt();
                        s = t * c;
                        tau = s / (1.0 + c);
                        h = t * matrix.index((i, j));
                        z[i] -= h;
                        z[j] += h;
                        evalues[i] -= h;
                        evalues[j] += h;
                        matrix[i + j * 3] = 0.0;
                        for k in 0..i {
                            matrix = self.roation_jacobi(matrix, k, i, k, j, s, tau);
                        }
                        for k in i + 1..j {
                            matrix = self.roation_jacobi(matrix, i, k, k, j, s, tau);
                        }
                        for k in j + 1..3 {
                            matrix = self.roation_jacobi(matrix, i, k, j, k, s, tau);
                        }
                        for k in 0..3 {
                            evectors = self.roation_jacobi(evectors, k, i, k, j, s, tau);
                        }
                    }
                }
            }

            b += z;
            evalues += z;
            z = Vector3::zeros();
        }

        return (true, evectors, evalues);
    }

    pub fn roation_jacobi(
        &self,
        mut matrix: Matrix3<f64>,
        i: usize,
        j: usize,
        k: usize,
        l: usize,
        s: f64,
        tau: f64,
    ) -> Matrix3<f64> {
        let g = *matrix.index((i, j));
        let h = *matrix.index((k, l));
        matrix[i + j * 3] = g - s * (h + g * tau);
        matrix[k + l * 3] = h + s * (g - h * tau);
        return matrix;
    }
}

pub_struct!( RigidBodiesData {

    max_radius: f64,
    body_data: Vec<RigidBodySphereData>,

    youngs_mod: Vec<f64>,
    poisson_ratio: Vec<f64>,

    position: Vec<Vector3<f64>>,
    velocity: Vec<Vector3<f64>>,
    force: Vec<Vector3<f64>>,

    moment_of_inertia_tensor: Vec<Matrix3<f64>>,
    diagonal_inertia: Vec<Vector3<f64>>,
    quaternion: Vec<UnitQuaternion<f64>>,
    omega: Vec<Vector3<f64>>,
    torque: Vec<Vector3<f64>>,
    angular_moment: Vec<Vector3<f64>>,
    ex_space: Vec<Vector3<f64>>,
    ey_space: Vec<Vector3<f64>>,
    ez_space: Vec<Vector3<f64>>,

    is_collision: Vec<bool>,

    materials: Vec<Material>,

    sphere_material: Vec<usize>,



    restitution_coefficient: f64,
    beta: f64,
    friction: f64,
    volume_fraction: f64,


});

impl RigidBodiesData {
    pub fn update_all_body_spheres_positions(&mut self) {
        for i in 0..self.body_data.len() {
            self.update_one_body_spheres_positions(i);
        }
    }

    pub fn update_one_body_spheres_positions(&mut self, index: usize) {
        for i in 0..self.body_data[index].local_position.len() {
            let temp = self.position[index]
                + self.quaternion[index].transform_vector(&self.body_data[index].local_position[i]);
            self.body_data[index].position[i] = temp;

            let local_level_arm = self.body_data[index].position[i] - self.position[index];

            self.body_data[index].velocity[i] =
                self.velocity[index] + self.omega[index].cross(&local_level_arm);
        }
    }

    pub fn update_all_rotation_with_omega(&mut self, step_size: f64) {
        for i in 0..self.body_data.len() {
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
        for i in 0..self.body_data.len() {
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
