use nalgebra::{Matrix3, Vector2, Vector3};
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


});

pub_struct!( ForceData {
    particle_indexes: Vec<Vector2<usize>>,
    force: Vec<Vector3<f64>>,
    del: Vec<Vector3<f64>>,
    forcedata:  Vec<Matrix3<f64>>,
});
