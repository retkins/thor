
/// In place cross-product  
/// $$ \vec{a} \times \vec{b} = \vec{c} $$
#[inline(always)]
pub fn cross(
    a: (f64, f64, f64), 
    b: (f64, f64, f64), 
    c: &mut (f64, f64, f64)
) {
    c.0 = a.1*b.2 - a.2*b.1;
    c.1 = a.2*b.0 - a.0*b.2; 
    c.2 = a.0*b.1 - a.1*b.0;
}

/// Custom min/max finding function
/// Returns the bounds of the array (min, max)
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


/// Sort a slice `data` using a temporary pre-allocated buffer `scratch`, 
/// using pre-calculated `indices` for the sort operation
pub fn sort_by_indices<T>(data: &mut [T], scratch: &mut [T], indices: &[usize]) 
    where T: Copy { 
    scratch.copy_from_slice(data);
    for (i, &idx) in indices.iter().enumerate() {
        data[i] = scratch[idx];
    }
}