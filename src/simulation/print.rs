use nalgebra::Matrix3;
use std::{collections::HashMap, fs::File, io::Write};

use crate::sphere::{self, BondData};
// Prints positions of particles to a vtp file in the vtp folder (no check is done for opening file, must include folder or no printing)
pub fn print_vtp(p_data: &mut sphere::ParticleData, count: i32) {
    let filename = format!("./vtp/{}CYCLE.vtp", count);
    let mut file = File::create(filename).unwrap();

    // Write a &str in the file (ignoring the result).
    write!(&mut file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n<PolyData>\n").unwrap();

    write!(
        &mut file,
        "<Piece NumberOfPoints=\"{}\">\n",
        p_data.radius.len()
    )
    .unwrap();

    write!(&mut file, "<Points>").unwrap();

    write!(
        &mut file,
        "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"
    )
    .unwrap();
    for i in 0..p_data.radius.len() {
        writeln!(
            &mut file,
            "{} {} {}",
            p_data.position[i][0], p_data.position[i][1], p_data.position[i][2]
        )
        .unwrap();
    }
    write!(&mut file, "</DataArray>\n").unwrap();
    write!(&mut file, "</Points>\n").unwrap();

    write!(&mut file, "<PointData Scalars=\"\" Vectors=\"\">\n").unwrap();
    write!(
        &mut file,
        "<DataArray type=\"Float32\" Name=\"Radius\" format=\"ascii\">\n"
    )
    .unwrap();
    for i in 0..p_data.radius.len() {
        writeln!(&mut file, "{}", p_data.radius[i]).unwrap();
    }

    write!(&mut file, "</DataArray>").unwrap();

    write!(
        &mut file,
        "</PointData>\n</Piece>\n</PolyData>\n</VTKFile>\n"
    )
    .unwrap();
}

// pub fn print_stress(kinetic_tensor: Matrix3<f64>, collision_tensor: Matrix3<f64>, count: i32)
// {

//      let mut file = match File::options()
//                             .read(true)
//                             .write(true)
//                             .append(true)
//                             .open("stress.txt") {
//         Ok(file) => {
//             file
//         }
//         Err(err) => {

//             println!("Error: {}", err);
//             std::process::exit(1);
//         }
//     };

//     writeln!(&mut file,"{} {} {} {} {} {} {} {} {} {}",count,kinetic_tensor.index((0,0)),kinetic_tensor.index((1,0)), kinetic_tensor.index((1,1)),collision_tensor.index((0,0)),collision_tensor.index((1,0)), collision_tensor.index((1,1)),kinetic_tensor.index((0,0)).abs() + collision_tensor.index((0,0)).abs(),kinetic_tensor.index((1,0)).abs() + collision_tensor.index((1,0)).abs(), kinetic_tensor.index((1,1)).abs() + collision_tensor.index((1,1)).abs()).ok();
//    }

pub fn print_positions(
    p_data: &mut sphere::ParticleData,
    bond_data: &mut HashMap<(usize, usize), sphere::BondData>,
) {
    let mut file = match File::options()
        .read(true)
        .write(true)
        .append(true)
        .open("data.txt")
    {
        Ok(file) => file,
        Err(err) => {
            println!("Error: {}", err);
            std::process::exit(1);
        }
    };
    for i in 0..p_data.radius.len() {
        writeln!(
            &mut file,
            "{} {} {}",
            p_data.position[i][0], p_data.position[i][1], p_data.position[i][2]
        );
    }

    let mut file = match File::options()
        .read(true)
        .write(true)
        .append(true)
        .open("bond_data.txt")
    {
        Ok(file) => file,
        Err(err) => {
            println!("Error: {}", err);
            std::process::exit(1);
        }
    };
    for ((i, j), b_data) in bond_data.iter() {
        writeln!(
            &mut file,
            "{} {} {} {} {} {}",
            i, j, b_data.current_sigma, b_data.current_tau, b_data.sigma_max, b_data.tau_max
        );
    }
}
