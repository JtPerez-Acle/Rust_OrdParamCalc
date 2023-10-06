use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;

// Function to compute dot product of two vectors
fn dot_product((x1, y1, z1): (f64, f64, f64), (x2, y2, z2): (f64, f64, f64)) -> f64 {
    x1 * x2 + y1 * y2 + z1 * z2
}

// Function to compute magnitude of a vector - inc Debug message
fn magnitude((x, y, z): (f64, f64, f64)) -> f64 {
    if x.is_infinite() || y.is_infinite() || z.is_infinite() || x.is_nan() || y.is_nan() || z.is_nan() {
        println!("Invalid values: x={}, y={}, z={}", x, y, z);
        return 0.0;
    }
    let result = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
    if result.is_nan() {
        println!("NaN result for values: x={}, y={}, z={}", x, y, z);
    }
    result
}

fn angle_between(v1: (f64, f64, f64), v2: (f64, f64, f64)) -> f64 {
    let dot = dot_product(v1, v2);
    let mag_v1 = magnitude(v1);
    let mag_v2 = magnitude(v2);
    if mag_v1 == 0.0 || mag_v2 == 0.0 {
        println!("Warning: Zero magnitude vector encountered. Check your data.");
        return 0.0;
    }
    let value = dot / (mag_v1 * mag_v2);
    let clamped_value = value.max(-1.0).min(1.0);  // Clamp the value to the range [-1, 1]
    clamped_value.acos().to_degrees()
}


fn order_parameter(angle_deg: f64) -> f64 {
    let angle_rad = angle_deg.to_radians();  // Convert the angle from degrees to radians
    1.5 * angle_rad.cos().powi(2) - 0.5
}


fn median(values: &mut Vec<f64>) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        (values[mid - 1] + values[mid]) / 2.0
    } else {
        values[mid]
    }
}


// Function to extract molecule data from a line
fn extract_molecule_data(line: &str) -> Result<(String, i32, f64, f64, f64), io::Error> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 7 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Insufficient data in line."));
    }
    let molecule_name = parts[1].to_string();
    let molecule_number: i32 = parts[2].parse().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to parse molecule number"))?;
    let x: f64 = parts[3].parse().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to parse x coordinate"))?;
    let y: f64 = parts[4].parse().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to parse y coordinate"))?;
    let z: f64 = parts[5].parse().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to parse z coordinate"))?;
    Ok((molecule_name, molecule_number, x, y, z))
}

fn main() -> io::Result<()> {
    // Path to the .gro file
    let path = Path::new("/home/jt/Desktop/TTA_tta_2adn_400ns_R1.gro");
    let file = File::open(&path)?;

    let mut molecule_values: HashMap<String, Vec<(f64, f64, f64)>> = HashMap::new();
    let mut molecule_order: Vec<String> = Vec::new();

    // User-specified range
    let mut range_input = String::new();
    print!("Enter the range for t (e.g., 0.00000 4000.00000): ");
    io::stdout().flush()?;
    io::stdin().read_line(&mut range_input)?;
    let range_parts: Vec<f64> = range_input.trim().split_whitespace().map(|s| s.parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, format!("Failed to parse float: {}", e)))).collect::<Result<Vec<f64>, _>>()?;
    if range_parts.len() != 2 {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid time range format."));
    }
    let (range_start, range_end) = (range_parts[0], range_parts[1]);
    let mut is_processing_frame = false;
    let mut frame_count = 0;

    // Ask the user for the central molecule
    let mut input = String::new();
    print!("Enter the central atom (e.g., CAL): ");
    io::stdout().flush()?;
    io::stdin().read_line(&mut input)?;
    let central_atom = input.trim();
    let mut current_frame_time: f64 = 0.0; // Track the current frame's time

    for line in io::BufReader::new(File::open(&path)?).lines() {
        let line = line?;
    
        if line.contains("t=") {
            // End of previous frame
            if molecule_order.len() > 0 {
                let mut total_order_param_frame = 0.0;
                let mut order_params_frame: Vec<f64> = Vec::new();
    
                for (index, key) in molecule_order.iter().enumerate() {
                    if key.contains(central_atom) {
                        if index > 0 && index < molecule_order.len() - 1 {
                            let key_prev = &molecule_order[index - 1];
                            let key_next = &molecule_order[index + 1];
    
                            let coords_prev = molecule_values.get(key_prev).unwrap_or(&vec![(0.0, 0.0, 0.0)])[0];
                            let coords = molecule_values.get(key).unwrap_or(&vec![(0.0, 0.0, 0.0)])[0];
                            let coords_next = molecule_values.get(key_next).unwrap_or(&vec![(0.0, 0.0, 0.0)])[0];
    
                            // Calculate the vector for the bond of interest using neighboring atoms
                            let vector_bond = (coords_next.0 - coords_prev.0, coords_next.1 - coords_prev.1, coords_next.2 - coords_prev.2);
    
                            // Reference vector (z-axis)
                            let vector_ref = (0.0, 0.0, 1.0);
    
                            // Calculate the angle between the bond of interest and the reference vector
                            let angle = angle_between(vector_bond, vector_ref);
    
                            // Calculate the order parameter
                            let order_param = order_parameter(angle);
                            total_order_param_frame += order_param;
                            order_params_frame.push(order_param);
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
                println!("Frame {}: Median order parameter: {}", frame_count, median(&mut order_params_frame));
    
                molecule_values.clear();
                molecule_order.clear();
            }
    
            // Start of a new frame
            current_frame_time = line.split("t=").nth(1).unwrap_or("0.0").trim().split_whitespace().next().unwrap_or("0.0").parse().unwrap_or(0.0);
            if current_frame_time >= range_start && current_frame_time < range_end + 10.0 { // Correctly handle the frame's time range
                is_processing_frame = true;
                frame_count += 1;
            } else {
                is_processing_frame = false;
            }
        }
    
        if is_processing_frame {
            if let Ok((molecule_name, molecule_number, x, y, z)) = extract_molecule_data(&line) {
                let key = format!("{} {}", molecule_name, molecule_number);
                molecule_values.entry(key.clone()).or_insert_with(Vec::new).push((x, y, z));
                molecule_order.push(key.clone());
            }
        }
    }
    

    Ok(())
}


/* Calculates Order Parameter, dot point. Issue is it averages the X, Y and Z of an i atom, and their i+1 and i-1.use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;
use std::f64::consts::PI;

// Function to compute dot product of two vectors
fn dot_product((x1, y1, z1): (f64, f64, f64), (x2, y2, z2): (f64, f64, f64)) -> f64 {
    x1 * x2 + y1 * y2 + z1 * z2
}

// Function to compute magnitude of a vector
fn magnitude((x, y, z): (f64, f64, f64)) -> f64 {
    (x.powi(2) + y.powi(2) + z.powi(2)).sqrt()
}

// Function to compute the angle between two vectors
fn angle_between(v1: (f64, f64, f64), v2: (f64, f64, f64)) -> f64 {
    let dot = dot_product(v1, v2);
    let mag_v1 = magnitude(v1);
    let mag_v2 = magnitude(v2);
    (dot / (mag_v1 * mag_v2)).acos()
}

// Function to compute the order parameter
fn order_parameter(angle: f64) -> f64 {
    1.5 * angle.cos().powi(2) - 0.5
}

fn main() -> io::Result<()> {
    // Path to the .gro file
    let path = Path::new("/home/jt/Desktop/TTA_tta_2adn_400ns_R1.gro");
    let file = File::open(&path)?;
    let reader = io::BufReader::new(file);

    let mut molecule_values: HashMap<String, Vec<(f64, f64, f64)>> = HashMap::new();
    let mut molecule_order: Vec<String> = Vec::new();

    // User-specified range
    let mut range_input = String::new();
    print!("Enter the range for t (e.g., 0.00000 4000.00000): ");
    io::stdout().flush()?;
    io::stdin().read_line(&mut range_input)?;
    let range_parts: Vec<f64> = range_input.trim().split_whitespace().map(|s| s.parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, format!("Failed to parse float: {}", e)))).collect::<Result<Vec<f64>, _>>()?;
    if range_parts.len() != 2 {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid time range format."));
    }
    let (range_start, range_end) = (range_parts[0], range_parts[1]);
    let mut entry_count = 0;

    for line in reader.lines() {
        let line = line?;
        
        if line.contains("t=") {
            let time_value: f64 = line.split("t=").nth(1).unwrap_or("0.0").trim().split_whitespace().next().unwrap_or("0.0").parse().unwrap_or(0.0);
            if time_value == range_start {
                entry_count = 1;
            } else if time_value == range_end {
                break;
            }
        }
        
        if entry_count > 0 {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 7 {
                entry_count += 1;
                let molecule_name = parts[1];
                let molecule_number: i32 = parts[2].parse().unwrap_or(0);
                let x: f64 = parts[3].parse().unwrap_or(0.0);
                let y: f64 = parts[4].parse().unwrap_or(0.0);
                let z: f64 = parts[5].parse().unwrap_or(0.0);
                
                let key = format!("{} {}", molecule_name, molecule_number);
                molecule_values.entry(key.clone()).or_insert_with(Vec::new).push((x, y, z));
                molecule_order.push(key.clone());

                }
            }
        }

        // Calculate average coordinates for each molecule
    for (_key, coords_vec) in molecule_values.iter_mut() {
        let sum_coords = coords_vec.iter().fold((0.0, 0.0, 0.0), |(sum_x, sum_y, sum_z), &(x, y, z)| {
            (sum_x + x, sum_y + y, sum_z + z)
        });
        let count = coords_vec.len() as f64;
        *coords_vec = vec![(sum_coords.0 / count, sum_coords.1 / count, sum_coords.2 / count)];
    }


    // Ask the user for the central molecule and its number
    let mut input = String::new();
    print!("Enter the central molecule and its number (e.g., CAM 2): ");
    io::stdout().flush()?;
    io::stdin().read_line(&mut input)?;
    let parts: Vec<&str> = input.trim().split_whitespace().collect();
    let molecule_name = parts[0];
    let molecule_number: i32 = parts[1].parse().unwrap_or(0);

    let key = format!("{} {}", molecule_name, molecule_number);

    println!("Searching for molecule: {}", key);

    if let Some(position) = molecule_order.iter().position(|k| k == &key) {
        if position > 0 && position < molecule_order.len() - 1 {
            let key_prev = &molecule_order[position - 1];
            let key_next = &molecule_order[position + 1];

            let coords_prev = molecule_values.get(key_prev).unwrap_or(&vec![(0.0, 0.0, 0.0)])[0];
            let coords = molecule_values.get(&key).unwrap_or(&vec![(0.0, 0.0, 0.0)])[0];
            let coords_next = molecule_values.get(key_next).unwrap_or(&vec![(0.0, 0.0, 0.0)])[0];
            
            // Calculate the vector for the bond of interest using neighboring atoms (i-1 and i+1)
            let vector_bond = (coords_next.0 - coords_prev.0, coords_next.1 - coords_prev.1, coords_next.2 - coords_prev.2);
            
            // Reference vector (z-axis)
            let vector_ref = (0.0, 0.0, 1.0);
            
            // Calculate the angle between the bond of interest and the reference vector
            let angle = angle_between(vector_bond, vector_ref);
            
            // Calculate the order parameter
            let order_param = order_parameter(angle);
            println!("Coordinates for {}: {:?}", key_prev, coords_prev);
            println!("Coordinates for {}: {:?}", key, coords);
            println!("Coordinates for {}: {:?}", key_next, coords_next);
            println!("Vector for bond of interest ({} to {}): {:?}", key_prev, key_next, vector_bond);
            println!("Angle with reference vector (z-axis): {} radians or {} degrees", angle, angle * (180.0 / PI));
            println!("Order parameter: {}", order_param);


        } else {
            println!("Cannot find neighboring atoms for {}.", key);
        }
    } else {
        println!("Cannot find the central atom {}.", key);
    }

    Ok(())
}  
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;
use std::f64::consts::PI;

// Function to compute dot product of two vectors
fn dot_product((x1, y1, z1): (f64, f64, f64), (x2, y2, z2): (f64, f64, f64)) -> f64 {
    x1 * x2 + y1 * y2 + z1 * z2
}

// Function to compute magnitude of a vector
fn magnitude((x, y, z): (f64, f64, f64)) -> f64 {
    (x.powi(2) + y.powi(2) + z.powi(2)).sqrt()
}

// Function to compute the angle between two vectors
fn angle_between(v1: (f64, f64, f64), v2: (f64, f64, f64)) -> f64 {
    let dot = dot_product(v1, v2);
    let mag_v1 = magnitude(v1);
    let mag_v2 = magnitude(v2);
    (dot / (mag_v1 * mag_v2)).acos()
}

// Function to compute the order parameter
fn order_parameter(angle: f64) -> f64 {
    1.5 * angle.cos().powi(2) - 0.5
}

fn main() -> io::Result<()> {
    // Path to the .gro file
    let path = Path::new("/home/jt/Desktop/TTA_tta_2adn_400ns_R1.gro");
    let file = File::open(&path)?;
    let reader = io::BufReader::new(file);

    let mut molecule_values: HashMap<String, Vec<(f64, f64, f64)>> = HashMap::new();
    let mut molecule_order: Vec<String> = Vec::new();

    // User-specified range
    let mut range_input = String::new();
    print!("Enter the range for t (e.g., 0.00000 4000.00000): ");
    io::stdout().flush()?;
    io::stdin().read_line(&mut range_input)?;
    let range_parts: Vec<f64> = range_input.trim().split_whitespace().map(|s| s.parse::<f64>().map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, format!("Failed to parse float: {}", e)))).collect::<Result<Vec<f64>, _>>()?;
    if range_parts.len() != 2 {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid time range format."));
    }
    let (range_start, range_end) = (range_parts[0], range_parts[1]);
    let mut entry_count = 0;

    for line in reader.lines() {
        let line = line?;
        
        if line.contains("t=") {
            let time_value: f64 = line.split("t=").nth(1).unwrap_or("0.0").trim().split_whitespace().next().unwrap_or("0.0").parse().unwrap_or(0.0);
            if time_value == range_start {
                entry_count = 1;
            } else if time_value == range_end {
                break;
            }
        }
        
        if entry_count > 0 {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 7 {
                entry_count += 1;
                let molecule_name = parts[1];
                let molecule_number: i32 = parts[2].parse().unwrap_or(0);
                let x: f64 = parts[3].parse().unwrap_or(0.0);
                let y: f64 = parts[4].parse().unwrap_or(0.0);
                let z: f64 = parts[5].parse().unwrap_or(0.0);
                
                let key = format!("{} {}", molecule_name, molecule_number);
                molecule_values.entry(key.clone()).or_insert_with(Vec::new).push((x, y, z));
                molecule_order.push(key.clone());

                }
            }
        }

        // Calculate average coordinates for each molecule
    for (_key, coords_vec) in molecule_values.iter_mut() {
        let sum_coords = coords_vec.iter().fold((0.0, 0.0, 0.0), |(sum_x, sum_y, sum_z), &(x, y, z)| {
            (sum_x + x, sum_y + y, sum_z + z)
        });
        let count = coords_vec.len() as f64;
        *coords_vec = vec![(sum_coords.0 / count, sum_coords.1 / count, sum_coords.2 / count)];
    }


    // Ask the user for the central molecule and its number
    let mut input = String::new();
    print!("Enter the central molecule and its number (e.g., CAM 2): ");
    io::stdout().flush()?;
    io::stdin().read_line(&mut input)?;
    let parts: Vec<&str> = input.trim().split_whitespace().collect();
    let molecule_name = parts[0];
    let molecule_number: i32 = parts[1].parse().unwrap_or(0);

    let key = format!("{} {}", molecule_name, molecule_number);

    println!("Searching for molecule: {}", key);

    if let Some(position) = molecule_order.iter().position(|k| k == &key) {
        if position > 0 && position < molecule_order.len() - 1 {
            let key_prev = &molecule_order[position - 1];
            let key_next = &molecule_order[position + 1];

            let coords_prev = molecule_values.get(key_prev).unwrap_or(&vec![(0.0, 0.0, 0.0)])[0];
            let coords = molecule_values.get(&key).unwrap_or(&vec![(0.0, 0.0, 0.0)])[0];
            let coords_next = molecule_values.get(key_next).unwrap_or(&vec![(0.0, 0.0, 0.0)])[0];
            
            // Calculate the vector for the bond of interest using neighboring atoms (i-1 and i+1)
            let vector_bond = (coords_next.0 - coords_prev.0, coords_next.1 - coords_prev.1, coords_next.2 - coords_prev.2);
            
            // Reference vector (z-axis)
            let vector_ref = (0.0, 0.0, 1.0);
            
            // Calculate the angle between the bond of interest and the reference vector
            let angle = angle_between(vector_bond, vector_ref);
            
            // Calculate the order parameter
            let order_param = order_parameter(angle);
            println!("Coordinates for {}: {:?}", key_prev, coords_prev);
            println!("Coordinates for {}: {:?}", key, coords);
            println!("Coordinates for {}: {:?}", key_next, coords_next);
            println!("Vector for bond of interest ({} to {}): {:?}", key_prev, key_next, vector_bond);
            println!("Angle with reference vector (z-axis): {} radians or {} degrees", angle, angle * (180.0 / PI));
            println!("Order parameter: {}", order_param);


        } else {
            println!("Cannot find neighboring atoms for {}.", key);
        }
    } else {
        println!("Cannot find the central atom {}.", key);
    }

    Ok(())
}    
*/

/* Calculates the dot_product, angle, and order parameter of an I, I+1 and I+2 atoms.
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;
use std::f64::consts::PI;

// Function to compute dot product of two vectors
fn dot_product((x1, y1, z1): (f64, f64, f64), (x2, y2, z2): (f64, f64, f64)) -> f64 {
    x1 * x2 + y1 * y2 + z1 * z2
}

// Function to compute magnitude of a vector
fn magnitude((x, y, z): (f64, f64, f64)) -> f64 {
    (x.powi(2) + y.powi(2) + z.powi(2)).sqrt()
}

// Function to compute the angle between two vectors
fn angle_between(v1: (f64, f64, f64), v2: (f64, f64, f64)) -> f64 {
    let dot = dot_product(v1, v2);
    let mag_v1 = magnitude(v1);
    let mag_v2 = magnitude(v2);
    (dot / (mag_v1 * mag_v2)).acos()
}

// Function to compute the order parameter
fn order_parameter(angle: f64) -> f64 {
    1.5 * angle.cos().powi(2) - 0.5
}

fn main() -> io::Result<()> {
    // Path to the .gro file
    let path = Path::new("/home/jt/Desktop/TTA_tta_2adn_400ns_R1.gro");
    let file = File::open(&path)?;
    let reader = io::BufReader::new(file);

    let mut molecules = Vec::new();
    let mut molecule_values: HashMap<String, Vec<(f64, f64, f64)>> = HashMap::new();

    // User-specified range
    let mut range_input = String::new();
    print!("Enter the range (e.g., 0 15000): ");
    io::stdout().flush()?;
    io::stdin().read_line(&mut range_input)?;
    let range_parts: Vec<i32> = range_input.trim().split_whitespace().map(|s| s.parse().unwrap()).collect();
    let (range_start, range_end) = (range_parts[0], range_parts[1]);

    let mut entry_count = 0;

    for line in reader.lines() {
        let line = line?;
        
        if line.contains(&format!("step= {}", range_end)) {
            break;
        }

        if line.contains(&format!("step= {}", range_start)) || entry_count > 0 {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 7 {
                entry_count += 1;
                let molecule_name = parts[1];
                let molecule_number: i32 = parts[2].parse().unwrap_or(0);
                let x: f64 = parts[3].parse().unwrap_or(0.0);
                let y: f64 = parts[4].parse().unwrap_or(0.0);
                let z: f64 = parts[5].parse().unwrap_or(0.0);
                
                let key = format!("{} {}", molecule_name, molecule_number);
                molecule_values.entry(key.clone()).or_insert_with(Vec::new).push((x, y, z));
            }
        }
    }

    // Ask the user for 3 molecules and their numbers
    let mut input = String::new();
    print!("Enter 3 molecules and their numbers (e.g., CAM 1; CAA 3; CAQ 16): ");
    io::stdout().flush()?;
    io::stdin().read_line(&mut input)?;
    let selections: Vec<&str> = input.trim().split(';').collect();

    for selection in selections {
        let parts: Vec<&str> = selection.trim().split_whitespace().collect();
        if parts.len() == 2 {
            let molecule_name = parts[0];
            let molecule_number: i32 = parts[1].parse().unwrap_or(0);
            let key = format!("{} {}", molecule_name, molecule_number);
            
            if let Some(values) = molecule_values.get(&key) {
                println!("Molecule {}: Found {} times", key, values.len());
                for (x, y, z) in values {
                    println!("Values: x={}, y={}, z={}", x, y, z);
                }
                
                // Ask the user if they want to average the values
                let mut avg_input = String::new();
                print!("Do you want to average these values for {}? (yes/no): ", key);
                io::stdout().flush()?;
                io::stdin().read_line(&mut avg_input)?;
                let avg_choice = avg_input.trim().to_lowercase();
    
                if avg_choice == "yes" || avg_choice == "y" {
                    let total_values = values.len() as f64;
                    let avg_x: f64 = values.iter().map(|(x, _, _)| x).sum::<f64>() / total_values;
                    let avg_y: f64 = values.iter().map(|(_, y, _)| y).sum::<f64>() / total_values;
                    let avg_z: f64 = values.iter().map(|(_, _, z)| z).sum::<f64>() / total_values;
                    println!("Averaged Values for {}: x={}, y={}, z={}", key, avg_x, avg_y, avg_z);
                    molecules.push((molecule_name.to_string(), molecule_number, avg_x, avg_y, avg_z));
                } else {
                    let &(x, y, z) = &values[0];    
                    molecules.push((molecule_name.to_string(), molecule_number, x, y, z));
                }
            }
        }
    }

    // Ask the user which vectors represent the bond of interest and the bond of reference
    let mut bond_input = String::new();
    print!("Do you want to use the first and second molecules as the bond of interest and the third molecule as the reference point? (yes/no): ");
    io::stdout().flush()?;
    io::stdin().read_line(&mut bond_input)?;
    let bond_choice = bond_input.trim().to_lowercase();

    let (bond_of_interest, reference_point) = if bond_choice == "yes" || bond_choice == "y" {
        // First and second molecules are bond of interest, third molecule is reference point
        ((&molecules[0], &molecules[1]), &molecules[2])
    } else {
        // User will specify
        // TODO: Ask the user to specify which molecules form the bond of interest and which one is the reference point
        // For now, we'll default to the first choice
        ((&molecules[0], &molecules[1]), &molecules[2])
    };

    // Calculate the vector for the bond of interest
    let (name1, number1, x1, y1, z1) = bond_of_interest.0;
    let (name2, number2, x2, y2, z2) = bond_of_interest.1;
    let vector = (x2 - x1, y2 - y1, z2 - z1);
    println!("Vector for bond of interest between {} {} and {} {}: ({}, {}, {})", name1, number1, name2, number2, vector.0, vector.1, vector.2);

    // Calculate the angle and order parameter using the reference point
    let angle = angle_between((0.0, 0.0, -1.0), vector); //let angle = angle_between((ref_x - x1, ref_y - y1, ref_z - z1), vector);
    let order_param = order_parameter(angle);
    println!("Angle with reference point {}: {} radians", reference_point.0, angle);
    println!("Angle in degrees: {}", angle * (180.0 / PI));
    println!("Order parameter: {}", order_param);

    Ok(())
}
*/


/* -- THE FOLLOWING CODE WORKS, DONT MODIFY IT, IT WILL ASK THE RANGE, THE MOLECULES, IF AVERAGE (Y/N), SAVING AVERAGE VALUES AND CALCULATING ITS VECTORS

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;

fn main() -> io::Result<()> {
    // Path to the .gro file
    let path = Path::new("/home/jt/Desktop/TTA_tta_2adn_400ns_R1.gro");
    let file = File::open(&path)?;
    let reader = io::BufReader::new(file);

    let mut molecules = Vec::new();
    let mut molecule_values: HashMap<String, Vec<(f64, f64, f64)>> = HashMap::new();

    // User-specified range
    let mut range_input = String::new();
    print!("Enter the range (e.g., 0 15000): ");
    io::stdout().flush()?;
    io::stdin().read_line(&mut range_input)?;
    let range_parts: Vec<i32> = range_input.trim().split_whitespace().map(|s| s.parse().unwrap()).collect();
    let (range_start, range_end) = (range_parts[0], range_parts[1]);

    let mut entry_count = 0;

    for line in reader.lines() {
        let line = line?;
        
        if line.contains(&format!("step= {}", range_end)) {
            break;
        }

        if line.contains(&format!("step= {}", range_start)) || entry_count > 0 {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 7 {
                entry_count += 1;
                let molecule_name = parts[1];
                let molecule_number: i32 = parts[2].parse().unwrap_or(0);
                let x: f64 = parts[3].parse().unwrap_or(0.0);
                let y: f64 = parts[4].parse().unwrap_or(0.0);
                let z: f64 = parts[5].parse().unwrap_or(0.0);
                
                let key = format!("{} {}", molecule_name, molecule_number);
                molecule_values.entry(key.clone()).or_insert_with(Vec::new).push((x, y, z));
            }
        }
    }

    // Ask the user for 3 molecules and their numbers
    let mut input = String::new();
    print!("Enter 3 molecules and their numbers (e.g., CAM 1; CAA 3; CAQ 16): ");
    io::stdout().flush()?;
    io::stdin().read_line(&mut input)?;
    let selections: Vec<&str> = input.trim().split(';').collect();

    for selection in selections {
        let parts: Vec<&str> = selection.trim().split_whitespace().collect();
        if parts.len() == 2 {
            let molecule_name = parts[0];
            let molecule_number: i32 = parts[1].parse().unwrap_or(0);
            let key = format!("{} {}", molecule_name, molecule_number);
            
            if let Some(values) = molecule_values.get(&key) {
                println!("Molecule {}: Found {} times", key, values.len());
                for (x, y, z) in values {
                    println!("Values: x={}, y={}, z={}", x, y, z);
                }
                
                // Ask the user if they want to average the values
                let mut avg_input = String::new();
                print!("Do you want to average these values for {}? (yes/no): ", key);
                io::stdout().flush()?;
                io::stdin().read_line(&mut avg_input)?;
                let avg_choice = avg_input.trim().to_lowercase();
    
                if avg_choice == "yes" || avg_choice == "y" {
                    let total_values = values.len() as f64;
                    let avg_x: f64 = values.iter().map(|(x, _, _)| x).sum::<f64>() / total_values;
                    let avg_y: f64 = values.iter().map(|(_, y, _)| y).sum::<f64>() / total_values;
                    let avg_z: f64 = values.iter().map(|(_, _, z)| z).sum::<f64>() / total_values;
                    println!("Averaged Values for {}: x={}, y={}, z={}", key, avg_x, avg_y, avg_z);
                    molecules.push((molecule_name.to_string(), molecule_number, avg_x, avg_y, avg_z));
                } else {
                    let &(x, y, z) = &values[0];    
                    molecules.push((molecule_name.to_string(), molecule_number, x, y, z));
                }
            }
        }
    }

    // Calculate and display the vectors between the selected molecules
    for i in 0..molecules.len() {
        for j in i + 1..molecules.len() {
            let (name1, number1, x1, y1, z1) = &molecules[i];
            let (name2, number2, x2, y2, z2) = &molecules[j];
            let dx = x2 - x1;
            let dy = y2 - y1;
            let dz = z2 - z1;
            println!("Vector between {} {} and {} {}: ({}, {}, {})", name1, number1, name2, number2, dx, dy, dz);
        }
    }

    Ok(())
}
*/