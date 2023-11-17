// Function to compute dot product of two vectors
pub fn dot_product((x1, y1, z1): (f64, f64, f64), (x2, y2, z2): (f64, f64, f64)) -> f64 {
    x1 * x2 + y1 * y2 + z1 * z2
}

pub fn normalize(v: (f64, f64, f64)) -> (f64, f64, f64) {
    let magnitude = (v.0.powf(2.0) + v.1.powf(2.0) + v.2.powf(2.0)).sqrt();
    if magnitude == 0.0 {
        println!("Warning: Zero magnitude vector encountered. Check your data.");
        return (0.0, 0.0, 0.0);
    }
    (v.0 / magnitude, v.1 / magnitude, v.2 / magnitude)
}

// Calculate the angle between two normalized vectors
pub fn angle_between(v1: (f64, f64, f64), v2: (f64, f64, f64)) -> f64 {
    let v1_normalized = normalize(v1);
    let v2_normalized = normalize(v2);

    let dot = dot_product(v1_normalized, v2_normalized);
    dot.acos()
}

// Get the order parameter using the documented function
pub fn order_parameter(angle_deg: f64) -> f64 {
    1.5 * (angle_deg.cos().powf(2.0)) - 0.5
}

// Calculate the distance between i-1 and i+1 to check if it's greater than 0.3nm
pub fn distance((x1, y1, z1): (f64, f64, f64), (x2, y2, z2): (f64, f64, f64)) -> f64 {
    ((x2 - x1).powi(2) + (y2 - y1).powi(2) + (z2 - z1).powi(2)).sqrt()
}