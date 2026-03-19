#![allow(non_snake_case)]

use crate::MU0_4PI;
use crate::vec3::Vec3;

/// Compute the magnetic field at target points (x, y, z) using a direct (O(N^2)) Biot-Savart summation
///
/// This is the 'scalar' version of this function; i.e. it uses only one thread.
///
/// # Arguments
/// - `centx`, `centy`, `centz`: (m) locations of source element centroids in 3D space
/// - `vol`:                     (m^3) volume of each source element
/// - `jx`, `jy`, `jz`:          (A/m^2) current density vector of each source element
/// - `x`, `y`, `z`:             (m) location of each target point
/// - `bx`, `by`, `bz`:          (T) magnetic flux density at each target point
pub fn bfield_direct(
    centx: &[f64],
    centy: &[f64],
    centz: &[f64],
    vol: &[f64],
    jx: &[f64],
    jy: &[f64],
    jz: &[f64],
    x: &[f64],
    y: &[f64],
    z: &[f64],
    bx: &mut [f64],
    by: &mut [f64],
    bz: &mut [f64],
) -> Result<(), ()> {
    // TODO: length checks on input arrays

    // Outer loop over source elements
    for i in 0..centx.len() {
        // Hoist invariants out of the loop explicitly
        let centxi: f64 = centx[i];
        let centyi: f64 = centy[i];
        let centzi: f64 = centz[i];
        let vol_mu0_4pi: f64 = vol[i] * MU0_4PI;
        let jxi = jx[i];
        let jyi = jy[i];
        let jzi = jz[i];

        // Inner loop over target points
        for (((xj, yj), zj), ((bxj, byj), bzj)) in x
            .iter()
            .zip(y.iter())
            .zip(z.iter())
            .zip(bx.iter_mut().zip(by.iter_mut()).zip(bz.iter_mut()))
        {
            // Vector from the element centroid to the target point: r'
            let rx: f64 = xj - centxi;
            let ry: f64 = yj - centyi;
            let rz: f64 = zj - centzi;
            let r: f64 = (rx * rx + ry * ry + rz * rz).sqrt();

            // J x r'
            let jxrpx: f64 = jyi * rz - jzi * ry;
            let jxrpy: f64 = jzi * rx - jxi * rz;
            let jxrpz: f64 = jxi * ry - jyi * rx;

            // Null out the singularity around the element centroid
            // This avoids `jmp` instructions and enables auto-vectorization of the inner loop
            // Hard-coded a tolerance of 0.1mm (TODO: update)
            let mask = if r > 1e-4 { 1.0 } else { 0.0 };
            let r3 = r * r * r + (1.0 - mask);
            let constant = vol_mu0_4pi * mask / r3;

            // Accumulation
            *bxj += constant * jxrpx;
            *byj += constant * jxrpy;
            *bzj += constant * jxrpz;
        }
    }
    Ok(())
}

// Older version of the above that uses explicit indices for the inner loop
#[allow(unused)]
fn bfield_direct_old(
    centx: &[f64],
    centy: &[f64],
    centz: &[f64],
    vol: &[f64],
    jx: &[f64],
    jy: &[f64],
    jz: &[f64],
    x: &[f64],
    y: &[f64],
    z: &[f64],
    bx: &mut [f64],
    by: &mut [f64],
    bz: &mut [f64],
) {
    let m: usize = centx.len();
    let n: usize = x.len();

    for i in 0..m {
        let centxi: f64 = centx[i];
        let centyi: f64 = centy[i];
        let centzi: f64 = centz[i];
        let vol_mu0_4pi: f64 = vol[i] * MU0_4PI;
        let jxi = jx[i];
        let jyi = jy[i];
        let jzi = jz[i];

        for j in 0..n {
            let rx: f64 = x[j] - centxi;
            let ry: f64 = y[j] - centyi;
            let rz: f64 = z[j] - centzi;
            let r: f64 = (rx * rx + ry * ry + rz * rz).sqrt();

            let jxrpx: f64 = jyi * rz - jzi * ry;
            let jxrpy: f64 = jzi * rx - jxi * rz;
            let jxrpz: f64 = jxi * ry - jyi * rx;

            let mask = if r > 1e-4 { 1.0 } else { 0.0 };
            let r3 = r * r * r + (1.0 - mask); // avoid div by zero
            let constant = vol_mu0_4pi * mask / r3;
            bx[j] += constant * jxrpx;
            by[j] += constant * jxrpy;
            bz[j] += constant * jxrpz;
        }
    }
}

use crate::sources::tet4::h_field_tet4;

/// Compute the magnetic field using the direct tetrahedral integration method
pub fn hfield_direct_tet(
    nodes_flat: &[f64],
    vol: &[f64],
    jdensity_flat: &[f64],
    x: &[f64],
    y: &[f64],
    z: &[f64],
    hx: &mut [f64],
    hy: &mut [f64],
    hz: &mut [f64],
) -> Result<(), ()> {
    let n_sources = vol.len();
    let n_targets = x.len();
    let mut f = vec![Vec3([0.0; 3]); n_targets];

    // TODO: length checks
    for i in 0..n_sources {
        let node_slice = &nodes_flat[12 * i..12 * i + 12];
        let nodes = [
            Vec3([node_slice[0], node_slice[1], node_slice[2]]),
            Vec3([node_slice[3], node_slice[4], node_slice[5]]),
            Vec3([node_slice[6], node_slice[7], node_slice[8]]),
            Vec3([node_slice[9], node_slice[10], node_slice[11]]),
        ];

        let jdensity = Vec3([
            jdensity_flat[3 * i],
            jdensity_flat[3 * i + 1],
            jdensity_flat[3 * i + 2],
        ]);
        h_field_tet4(&nodes, &jdensity, (x, y, z), &mut f, (hx, hy, hz));
        f.fill(Vec3([0.0; 3]));
    }

    Ok(())
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_direct() {}
}
