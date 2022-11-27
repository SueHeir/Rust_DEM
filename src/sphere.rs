use nalgebra::{Matrix3, Vector2, Vector3, Vector4};
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
    shear_mod: Vec<f64>,
    poisson_ratio: Vec<f64>,
    density: Vec<f64>,
    position: Vec<Vector3<f64>>,
    velocity: Vec<Vector3<f64>>,
    force: Vec<Vector3<f64>>,

    moment_of_inertia: Vec<f64>,
    quaternion: Vec<Vector4<f64>>,
    omega: Vec<Vector3<f64>>,
    torque: Vec<Vector3<f64>>,



    is_collision: Vec<bool>,

    materials: Vec<Material>,

    sphere_material: Vec<usize>,

    sphere_material_map: HashMap<String,EffectiveMaterialPreCalc>,

    restitution_coefficient: f64,
    beta: f64,
    friction: f64,
    volume_fraction: f64,


});

pub_struct!( ContactData {
    particle_indexes: Vec<Vector2<usize>>,
    contact_struct: Vec<ContactStruct>,


});

pub_struct!( ContactStruct {
    force_magnitude: Vector2<f64>,
    previous_force_magnitude: Vector2<f64>,
    previous_force: Vector2<Vector3<f64>>,
    force: Vector2<Vector3<f64>>,
    tangential_force: Vector2<Vector3<f64>>,
    previous_tangential_force: Vector2<Vector3<f64>>,
    tangential_force_magnitude: Vector2<f64>,
    previous_tangential_force_magnitude: Vector2<f64>,
    sd: f64,
    sum_normal_force: f64,
    change_tang_direction: f64,
    tangential_displacement: Vector3<f64>,
    contact_vector: Vector2<Vector3<f64>>,
    contact_radius: f64,
    t_star: f64,
    t_double_star: f64,
    k_loading: i32, }
);
