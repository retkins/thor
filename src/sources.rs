/// Octree sources
/// Represent monopole/dipole/quadrapole sources for a collection
/// of source points acting on a target point 
/// 

/// Monopole approximation of the effect of a group of source points
/// at a target point. Approximates the group of sources as a 
/// single 'super source'
pub fn monopole(
    source_centroids: (&[f64], &[f64], &[f64]), 
    volumes: &[f64], 
    jdensity: (&[f64], &[f64], &[f64]), 
    target: (f64, f64, f64)
) -> (f64, f64, f64) {
    // TODO: length-check guards
    let m: usize = volumes.len();

    let mut mx = vec![0.0; m];
    let mut my = vec![0.0; m];
    let mut mz = vec![0.0; m];
    let (x, y, z) = source_centroids;
    let (jx, jy, jz) = jdensity;

    // Compute current density moments 
    for i in 0..m { 
        mx[i] = volumes[i]*jx[i];
        my[i] = volumes[i]*jy[i];
        mz[i] = volumes[i]*jz[i];
    }

    // Compute moment-weighted centroid location
    let mx_total: f64 = mx.iter().sum();
    let my_total: f64 = my.iter().sum();
    let mz_total: f64 = mz.iter().sum();
    let mnet_total = (mx_total*mx_total + my_total*my_total + mz_total*mz_total).sqrt();

    (0.0, 0.0, 0.0)
}