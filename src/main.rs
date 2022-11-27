mod simulation;
use nalgebra::{Vector3, Vector4};

use ex::fs::File;
use rand::prelude::*;
use std::{
    collections::HashMap,
    env,
    f64::consts::PI,
    io::{prelude::*, BufReader},
    process,
};

mod domain;
mod sphere;

fn lines_from_file(filename: String) -> Vec<String> {
    let filenamestring: String = filename.clone();
    let file = match File::open(filename) {
        Ok(file) => file,
        Err(err) => {
            println!("Error: {}", err);
            std::process::exit(1);
        }
    };

    println!("Opened File: {}", filenamestring);

    let buf = BufReader::new(file);
    buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect()
}

fn main() {
    let mut args: Vec<String> = env::args().collect();

    let mut command_stack: Vec<String> = Vec::<String>::new();

    println!("DEM code written by Elizabeth Suehr");

    if args.len() < 2 {
        println!("Please enter file to read as argument");
        process::exit(1);
    }

    let filename = args.remove(1);
    let lines = lines_from_file(filename);

    let mut p_data = sphere::ParticleData {
        max_radius: 0.0,
        radius: Vec::<f64>::new(),
        mass: Vec::<f64>::new(),
        youngs_mod: Vec::<f64>::new(),
        shear_mod: Vec::<f64>::new(),
        poisson_ratio: Vec::<f64>::new(),
        density: Vec::<f64>::new(),
        position: Vec::<Vector3<f64>>::new(),
        velocity: Vec::<Vector3<f64>>::new(),
        force: Vec::<Vector3<f64>>::new(),
        moment_of_inertia: Vec::<f64>::new(),
        quaternion: Vec::<Vector4<f64>>::new(),
        omega: Vec::<Vector3<f64>>::new(),
        torque: Vec::<Vector3<f64>>::new(),
        is_collision: Vec::<bool>::new(),
        materials: Vec::<sphere::Material>::new(),
        sphere_material: Vec::<usize>::new(),
        sphere_material_map: HashMap::new(),
        restitution_coefficient: 0.95,
        beta: 0.0,
        friction: 0.35,
        volume_fraction: 0.0,
    };
    let mut d_data = domain::DomainData {
        domain: Vector3::new(1.0, 1.0, 1.0),
        domain_volume: 1.0,
        collision_boxes: Vector3::new(1, 1, 1),
        g_data: Vec::new(),
        lees_edwards_boundary: 1.0,
    };

    let mut argument = 0;
    for line in lines {
        if line.len() <= 3 {
            continue;
        }

        let results: Vec<&str> = line.split_whitespace().collect();

        if argument == 0 {
            if &results[0][0..3] != "STA" {
                println!("First Command must be START");
                process::exit(1);
            }
        }

        match &results[0][0..3] {
            "STA" => {
                println!("{}", line);

                d_data.domain = Vector3::new(
                    results[1].parse::<f64>().unwrap(),
                    results[2].parse::<f64>().unwrap(),
                    results[3].parse::<f64>().unwrap(),
                );
                d_data.collision_boxes = Vector3::new(
                    results[4].parse::<i32>().unwrap(),
                    results[5].parse::<i32>().unwrap(),
                    results[6].parse::<i32>().unwrap(),
                );

                d_data.domain_volume = d_data.domain[0] * d_data.domain[1] * d_data.domain[2];

                let mut g_data: Vec<Vec<Vec<domain::Box>>> = Vec::new();

                for i in 0..d_data.collision_boxes[0] {
                    g_data.push(Vec::new());
                    for j in 0..d_data.collision_boxes[1] {
                        g_data[i as usize].push(Vec::new());
                        for k in 0..d_data.collision_boxes[2] {
                            let len = Vector3::new(
                                d_data.domain[0] / d_data.collision_boxes[0] as f64,
                                d_data.domain[1] / d_data.collision_boxes[1] as f64,
                                d_data.domain[2] / d_data.collision_boxes[2] as f64,
                            );

                            let the_box = domain::Box {
                                position: Vector3::new(i, j, k),
                                real: Vec::<i32>::new(),
                                ghost: Vec::<i32>::new(),
                                lo: Vector3::new(
                                    i as f64 * len[0],
                                    j as f64 * len[1],
                                    k as f64 * len[2],
                                ),
                                hi: Vector3::new(
                                    i as f64 * len[0] + len[0],
                                    j as f64 * len[1] + len[1],
                                    k as f64 * len[2] + len[2],
                                ),
                            };
                            g_data[i as usize][j as usize].push(the_box);
                        }
                    }
                }
                d_data.g_data = g_data;
            }
            "DAM" => {
                println!("{}", line);
                p_data.restitution_coefficient = results[1].parse::<f64>().unwrap();

                let log_e = p_data.restitution_coefficient.ln();
                let beta = -log_e / (PI * PI + log_e * log_e).sqrt();

                p_data.beta = beta;
            }
            "LEB" => {
                println!("{}", line);
                d_data.lees_edwards_boundary = results[1].parse::<f64>().unwrap();
            }
            "MAT" => {
                println!("{}", line);

                let material = sphere::Material {
                    radius: results[2].parse::<f64>().unwrap(),
                    mass: results[3].parse::<f64>().unwrap() * PI * 4.0 / 3.0
                        * results[2].parse::<f64>().unwrap().powi(3),
                    youngs_mod: results[4].parse::<f64>().unwrap(),
                    poisson_ratio: results[5].parse::<f64>().unwrap(),
                    density: results[3].parse::<f64>().unwrap(),
                    id: results[1].parse::<i32>().unwrap(),
                };

                p_data.materials.push(material);

                println!("Material Loaded: {:?} ", p_data.materials.last());

                let mut max_radius = p_data.max_radius;
                for material in &p_data.materials {
                    max_radius = material.radius.max(max_radius);
                }

                p_data.max_radius = max_radius;
            }
            "RGP" => {
                println!("{}", line);

                let mut rng = rand::thread_rng();

                let num_particles: i32 = results[1].parse::<i32>().unwrap();

                for (_i, material) in p_data.materials.iter().enumerate() {
                    if material.id == results[2].parse::<i32>().unwrap() {
                        println!("Generating particles with Material: {:?}", material);
                        for _j in 0..num_particles {
                            p_data.sphere_material.push(material.id as usize);

                            p_data.radius.push(material.radius);

                            p_data.mass.push(material.mass);
                            p_data
                                .moment_of_inertia
                                .push(2.0 / 5.0 * material.mass * material.radius.powi(2));

                            p_data.density.push(material.density);
                            p_data.youngs_mod.push(material.youngs_mod);
                            p_data
                                .shear_mod
                                .push(material.youngs_mod / (2.0 * (1.0 + material.poisson_ratio)));
                            p_data.poisson_ratio.push(material.poisson_ratio);

                            let x: f64 = rng.gen::<f64>();
                            let y: f64 = rng.gen::<f64>();
                            let z: f64 = rng.gen::<f64>();

                            let vx: f64 = rng.gen::<f64>();
                            let vy: f64 = rng.gen::<f64>();
                            let vz: f64 = rng.gen::<f64>();

                            p_data.position.push(Vector3::new(
                                d_data.domain.x * x,
                                d_data.domain.y * y,
                                d_data.domain.z * z,
                            ));
                            p_data.velocity.push(Vector3::new(
                                vx * 0.5 - 0.25,
                                vy * 0.5 - 0.25,
                                vz * 0.5 - 0.25,
                            ));
                            p_data.force.push(Vector3::new(0.0, 0.0, 0.0));

                            p_data.quaternion.push(Vector4::new(0.0, 0.0, 0.0, 0.0));
                            p_data.omega.push(Vector3::new(0.0, 0.0, 0.0));
                            p_data.torque.push(Vector3::new(0.0, 0.0, 0.0));

                            p_data.is_collision.push(false);
                        }
                    }
                }

                println!(
                    "Generating {} randomonly located spheres",
                    p_data.radius.len()
                );
            }
            "FOR" => {
                println!("{}", line);

                for (_i, material) in p_data.materials.iter().enumerate() {
                    if material.id == results[1].parse::<i32>().unwrap() {
                        println!("Generating two particle for force check");

                        p_data.sphere_material.push(material.id as usize);

                        p_data.radius.push(material.radius);

                        p_data.mass.push(material.mass);
                        p_data
                            .moment_of_inertia
                            .push(2.0 / 5.0 * material.mass * material.radius.powi(2));
                        p_data.density.push(material.density);
                        p_data.youngs_mod.push(material.youngs_mod);
                        p_data
                            .shear_mod
                            .push(material.youngs_mod / (2.0 * (1.0 + material.poisson_ratio)));
                        p_data.poisson_ratio.push(material.poisson_ratio);

                        p_data.position.push(Vector3::new(
                            d_data.domain.x * 0.4,
                            d_data.domain.y * 0.5,
                            d_data.domain.z * 0.5,
                        ));
                        p_data.velocity.push(Vector3::new(0.05, 0.0, 0.0));
                        p_data.force.push(Vector3::new(0.0, 0.0, 0.0));

                        p_data.quaternion.push(Vector4::new(0.0, 0.0, 0.0, 0.0));
                        p_data.omega.push(Vector3::new(0.0, 0.0, 0.0));
                        p_data.torque.push(Vector3::new(0.0, 0.0, 0.0));

                        p_data.is_collision.push(false);

                        p_data.sphere_material.push(material.id as usize);

                        p_data.radius.push(material.radius);

                        p_data.mass.push(material.mass);
                        p_data
                            .moment_of_inertia
                            .push(2.0 / 5.0 * material.mass * material.radius.powi(2));
                        p_data.density.push(material.density);
                        p_data.youngs_mod.push(material.youngs_mod);
                        p_data
                            .shear_mod
                            .push(material.youngs_mod / (2.0 * (1.0 + material.poisson_ratio)));
                        p_data.poisson_ratio.push(material.poisson_ratio);

                        p_data.position.push(Vector3::new(
                            d_data.domain.x * 0.6,
                            d_data.domain.y * 0.5 - 0.5e-4,
                            d_data.domain.z * 0.5,
                        ));
                        p_data.velocity.push(Vector3::new(-0.05, 0.0, 0.0));
                        p_data.force.push(Vector3::new(0.0, 0.0, 0.0));

                        p_data.quaternion.push(Vector4::new(0.0, 0.0, 0.0, 0.0));
                        p_data.omega.push(Vector3::new(0.0, 0.0, 0.0));
                        p_data.torque.push(Vector3::new(0.0, 0.0, 0.0));

                        p_data.is_collision.push(false);
                    }
                }

                println!("Generating {} located spheres", p_data.radius.len());
            }
            "REL" => {
                println!("{}", line);
                command_stack.push(line);
            }
            "CYC" => {
                println!("{}", line);
                command_stack.push(line);
            }
            _ => {}
        }

        argument += 1;
    }

    simulation::handle_commands(command_stack, d_data, p_data);
}
