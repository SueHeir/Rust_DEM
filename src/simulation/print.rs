use nalgebra::Matrix3;
use std::{fs::File, io::Write};

use crate::rigid_body;
// Prints positions of particles to a vtp file in the vtp folder (no check is done for opening file, must include folder or no printing)
pub fn print_vtp(p_data: &mut rigid_body::RigidBodiesData, count: i32) {
    let filename = format!("./vtp/{}CYCLE.vtp", count);
    let mut file = File::create(filename).unwrap();

    // Write a &str in the file (ignoring the result).
    write!(&mut file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n<PolyData>\n").unwrap();

    let mut sum: i32 = 0;
    for i in 0..p_data.body_data.len() {
        sum += p_data.body_data[i].radius.len() as i32;
    }

    write!(&mut file, "<Piece NumberOfPoints=\"{}\">\n", sum).unwrap();

    write!(&mut file, "<Points>").unwrap();

    write!(
        &mut file,
        "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"
    )
    .unwrap();
    for i in 0..p_data.body_data.len() {
        for j in 0..p_data.body_data[i].radius.len() {
            writeln!(
                &mut file,
                "{} {} {}",
                p_data.body_data[i].position[j][0],
                p_data.body_data[i].position[j][1],
                p_data.body_data[i].position[j][2]
            )
            .unwrap();
        }
    }
    write!(&mut file, "</DataArray>\n").unwrap();
    write!(&mut file, "</Points>\n").unwrap();

    write!(&mut file, "<PointData Scalars=\"\" Vectors=\"\">\n").unwrap();
    write!(
        &mut file,
        "<DataArray type=\"Float32\" Name=\"Radius\" format=\"ascii\">\n"
    )
    .unwrap();
    for i in 0..p_data.body_data.len() {
        for j in 0..p_data.body_data[i].radius.len() {
            writeln!(&mut file, "{}", p_data.body_data[i].radius[j]).unwrap();
        }
    }

    write!(&mut file, "</DataArray>").unwrap();

    write!(
        &mut file,
        "<DataArray type=\"Float32\" Name=\"Collision\" format=\"ascii\">\n"
    )
    .unwrap();
    for i in 0..p_data.body_data.len() {
        for j in 0..p_data.body_data[i].radius.len() {
            if p_data.is_collision[i] {
                writeln!(&mut file, "1",).unwrap();
            } else if !p_data.is_collision[i] {
                writeln!(&mut file, "0",).unwrap();
            }
        }
    }

    write!(&mut file, "</DataArray>").unwrap();

    write!(
        &mut file,
        "</PointData>\n</Piece>\n</PolyData>\n</VTKFile>\n"
    )
    .unwrap();
}

pub fn print_stress(kinetic_tensor: Matrix3<f64>, collision_tensor: Matrix3<f64>, count: i32) {
    let mut file = match File::options()
        .read(true)
        .write(true)
        .append(true)
        .open("stress.txt")
    {
        Ok(file) => file,
        Err(err) => {
            println!("Error: {}", err);
            std::process::exit(1);
        }
    };

    writeln!(
        &mut file,
        "{} {} {} {} {} {} {} {} {} {}",
        count,
        kinetic_tensor.index((0, 0)),
        kinetic_tensor.index((1, 0)),
        kinetic_tensor.index((1, 1)),
        collision_tensor.index((0, 0)),
        collision_tensor.index((1, 0)),
        collision_tensor.index((1, 1)),
        kinetic_tensor.index((0, 0)).abs() + collision_tensor.index((0, 0)).abs(),
        kinetic_tensor.index((1, 0)).abs() + collision_tensor.index((1, 0)).abs(),
        kinetic_tensor.index((1, 1)).abs() + collision_tensor.index((1, 1)).abs()
    )
    .ok();
}
