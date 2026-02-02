
/// In place cross-product  
/// $$ \vec{a} \times \vec{b} = \vec{c} $$
#[inline(always)]
pub fn cross(
    a: &[f64; 3], 
    b: &[f64; 3], 
    c: &mut [f64; 3]
) {
    c[0] = a[1]*b[2]- a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2]; 
    c[2] = a[0]*b[1] - a[1]*b[0];
}

/// Custom min/max finding function  
/// 
/// # Returns: 
/// `(min, max)`: the maximum and minimum values in the slice
pub fn min_and_max(arr: &[f64]) -> Option<(f64, f64)> {
    let mut min: f64 = 0.0;
    let mut max: f64 = 0.0;

    if arr.len() < 1 {
        return None;
    }

    for a in arr.iter() {
        if *a < min {
            min = *a;
        }
        if *a > max {
            max = *a;
        }
    }
    Some((min, max))
}

/// Sort a slice using pre-computed indices
/// 
/// # Arguments:
/// - `data`:       slice to be sorted
/// - `scratch`:    pre-allocated buffer
/// - `indices`:    pre-computed index for each value in `data`
pub fn sort_by_indices<T>(data: &mut [T], scratch: &mut [T], indices: &[usize]) 
    where T: Copy { 
    scratch.copy_from_slice(data);
    for (i, &idx) in indices.iter().enumerate() {
        data[i] = scratch[idx];
    }
}


/// Take the vector magnitude (norm) of a slice 
#[inline(always)]
pub fn mag(slice: &[f64]) -> f64 {
    let mut sumsq: f64 = 0.0;
    for &value in slice.iter() {
        sumsq += value*value;
    }
    sumsq.sqrt()
}


/// Compute the distance between two points in 3D space
/// 
/// # Arguments 
/// - `a`: start vector
/// - `b`: end vector
/// 
/// # Returns
/// magnitude of the distance between the two points
#[inline(always)]
pub fn distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    mag(&vec_distance(a, b))
}

/// Distance vector from a to b
#[inline(always)]
pub fn vec_distance(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    let mut vec_dist: [f64; 3] = [0.0; 3];
    for i in 0..3 {
        vec_dist[i] = b[i] - a[i];
    }
    vec_dist
}