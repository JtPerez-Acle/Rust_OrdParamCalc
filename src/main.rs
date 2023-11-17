use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;

mod math_procss;
use crate::math_procss::*;

// Function to extract molecule data from a line
fn extract_molecule_data(line: &str, index: usize) -> Result<(String, f64, f64, f64), io::Error> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 7 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Insufficient data in line."));
    }
    let molecule_name = format!("{}_{}", parts[1], index);  // Append the index to the molecule name
    let x: f64 = parts[3].parse().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to parse x coordinate"))?;
    let y: f64 = parts[4].parse().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to parse y coordinate"))?;
    let z: f64 = parts[5].parse().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to parse z coordinate"))?;
    Ok((molecule_name, x, y, z))
}

fn handle_new_frame(line: &str, molecule_order: &mut Vec<String>, molecule_values: &mut HashMap<String, Vec<(f64, f64, f64)>>, frame_count: &mut usize, central_atom: &str) -> f64 {
    if molecule_order.len() > 0 {
        let mut total_order_param_frame = 0.0;
        let mut order_params_frame: Vec<f64> = Vec::new();

        for (index, key) in molecule_order.iter().enumerate() {
            let parts: Vec<&str> = key.split("_").collect();
            let molecule_name = parts[0];
            if molecule_name.contains(central_atom) {
                if index > 0 && index < molecule_order.len() - 1 {
                    let key_prev = &molecule_order[index - 1];
                    let key_next = &molecule_order[index + 1];

        
                    let coords_prev = molecule_values.get(key_prev).unwrap_or(&vec![(0.0, 0.0, 0.0)])[0];
                    let coords = molecule_values.get(key).unwrap_or(&vec![(0.0, 0.0, 0.0)])[0];
                    let coords_next = molecule_values.get(key_next).unwrap_or(&vec![(0.0, 0.0, 0.0)])[0];
        
                    // Check the distance between i-1 and i+1 atoms
                    if distance(coords_prev, coords_next) <= 0.3 {
                        // Calculate the vector for the bond of interest using neighboring atoms
                        let vector_bond = (coords_next.0 - coords_prev.0, coords_next.1 - coords_prev.1, coords_next.2 - coords_prev.2);

                        // Reference vector (z-axis)
                        let vector_ref = (0.0, 0.0, 1.0);
        
                        // Calculate the angle between the bond of interest and the reference vector
                        let angle = angle_between(vector_bond, vector_ref);
                        
                        // println!("Index: {}, Key: {}, Prev: {}, Next: {} ", index, key, key_prev, key_next);
                        // println!("Coords Prev: {:?}, Coords: {:?}, Coords Next: {:?}", coords_prev, coords, coords_next);

                        // Remove comments to check for bug and issues while processing
                        // Calculate the order parameter
                        let order_param: f64 = order_parameter(angle);
                        // Checks the value of order_param and prints it to console
                        // print!("{} ", order_param);
                        total_order_param_frame += order_param;
                        order_params_frame.push(order_param);
                        println!("Frame {}: Calculated order parameters count: {} ", frame_count, order_params_frame.len());
                        // println!("Frame {}: Total order parameter: {}", frame_count, total_order_param_frame);

                    } else {
                        break;
                    }
                }
            }
        }

        let average_order_param_frame = if order_params_frame.len() > 0 {
            total_order_param_frame / order_params_frame.len() as f64
        } else {
            println!("Warning: No order parameters were calculated for frame {}. This might indicate an issue with the input data or selected atom.", frame_count);
            0.0
        };
        println!("Frame {}: Average order parameter: {}", frame_count, average_order_param_frame);

        molecule_values.clear();
        molecule_order.clear();
    }

    let current_frame_time: f64 = line.split("t=").nth(1).unwrap_or("0.0").trim().split_whitespace().next().unwrap_or("0.0").parse().unwrap_or(0.0);
    *frame_count += 1;
    current_frame_time
}

fn process_molecule_data(line: &str, index: usize, molecule_order: &mut Vec<String>, molecule_values: &mut HashMap<String, Vec<(f64, f64, f64)>>) {
    if let Ok((molecule_name, x, y, z)) = extract_molecule_data(&line, index) {
        molecule_values.entry(molecule_name.clone()).or_insert_with(Vec::new).push((x, y, z));
        molecule_order.push(molecule_name.clone());
    }
}

fn main() -> io::Result<()> {
    // Path to the .gro file
    let path: &Path = Path::new("/home/tomas/O_Data/TTA_tta_2adn_400ns_R1.gro");

    let mut molecule_values: HashMap<String, Vec<(f64, f64, f64)>> = HashMap::new();
    let mut molecule_order: Vec<String> = Vec::new();

    // User-specified range
    let mut range_input = String::new();
    print!("Enter the range for t (MIN 0 MAX 400000; 1 FRAME = 10): ");
    io::stdout().flush()?;
    io::stdin().read_line(&mut range_input)?;
    let range_parts: Vec<f64> = range_input.trim().split_whitespace().map(|s| s.parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, format!("Failed to parse float: {}", e)))).collect::<Result<Vec<f64>, _>>()?;
    if range_parts.len() != 2 {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid time range format."));
    }
    let (range_start, range_end) = (range_parts[0], range_parts[1]);
    let mut frame_count = 0;

    // Ask the user for the central molecule
    let mut input = String::new();
    print!("Enter the central atom (e.g., CAL): ");
    io::stdout().flush()?;
    io::stdin().read_line(&mut input)?;
    let central_atom = input.trim();
    let mut current_frame_time: f64 = 0.0; // Track the current frame's time
    let mut is_processing_frame = false;

    for (index, line) in io::BufReader::new(File::open(&path)?).lines().enumerate() {
        let line: String = line?;
    
        if line.contains("t=") {
            // Start of a new frame
            current_frame_time = handle_new_frame(&line, &mut molecule_order, &mut molecule_values, &mut frame_count, &central_atom);
            if current_frame_time >= range_start && current_frame_time < range_end + 10.0 { 
                is_processing_frame = true;
            } else {
                is_processing_frame = false;
            }
        }
        
        if is_processing_frame {
            process_molecule_data(&line, index, &mut molecule_order, &mut molecule_values);
        }
       
    }
    

    Ok(())
}
