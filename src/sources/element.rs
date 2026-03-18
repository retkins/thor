use crate::math::{cross, dot3, mag3, unit_vector};

// Compute the edge integral for a finite element source
//
// Includes a small finite value to prevent singularities when the target is close to the edge
//
// Reference: Eq 52 of CRYO-06-034
pub fn edge_integral(x: f64, y: f64, z: f64) -> f64 {
    let r: f64 = mag3(x, y, z);

    if r <= 1e-8 {
        return 0.0;
    }

    // r could be close to -x, resulting in a singularity
    let mut result = if (x + r).abs() > 1e-8 {
        (x + r).ln()
    } else {
        0.0
    };

    // prevent singularities when target is on the edge
    let za = z.abs();
    if (y.abs() > 1e-12 * r) && (za > 1e-12 * r) {
        result += za / y * ((x * za / (y * r)).atan() - (x / y).atan());
    }

    result
}

// Compute the unit vectors of an edge csys given 3 nodes
// Nodes should be in right hand rule order
// middle node is the 'anchor', xhat points from node 1 to node 0
// yhat points away from the center of the face
pub fn edge_csys(x: &[f64; 3], y: &[f64; 3], z: &[f64; 3]) -> ([f64; 3], [f64; 3], [f64; 3]) {
    // first edge defines x'-axis
    let xp_hat = unit_vector(&[x[1], y[1], z[1]], &[x[0], y[0], z[0]]);

    // second edge defines the next in plane direction to find the out of plane z
    // it is not returned because it is not guaranteed to be orthogonal to x and z
    let yp_hat_temp = unit_vector(&[x[2], y[2], z[2]], &[x[1], y[1], z[1]]);

    // unit vector normal to surface
    let mut zp_hat = [0.0; 3];
    cross(&xp_hat, &yp_hat_temp, &mut zp_hat);
    let norm = mag3(zp_hat[0], zp_hat[1], zp_hat[2]);
    zp_hat[0] /= norm;
    zp_hat[1] /= norm;
    zp_hat[2] /= norm;

    // other in-plane direction, not necessarily aligned with an edge
    let mut yp_hat = [0.0; 3];
    cross(&xp_hat, &zp_hat, &mut yp_hat);
    return (xp_hat, yp_hat, zp_hat);
}

use crate::vec3::Vec3;
/// Fast version of the above
#[inline(always)]
pub fn edge_csys_fast(node1: &Vec3, node2: &Vec3, node3: &Vec3) -> (Vec3, Vec3, Vec3) {
    let xp_hat = Vec3(unit_vector(node2.to_slice(), node1.to_slice()));
    let yp_hat_temp = unit_vector(node3.to_slice(), node1.to_slice());
    let mut zp_hat = xp_hat.cross(&Vec3(yp_hat_temp));
    let norm_inv = 1.0 / zp_hat.mag();
    zp_hat *= norm_inv;
    let yp_hat = xp_hat.cross(&zp_hat);
    (xp_hat, yp_hat, zp_hat)
}

// transform a global xyz position into a local coordinate system
pub fn transform(g: &[f64; 3], xhat: &[f64; 3], yhat: &[f64; 3], zhat: &[f64; 3]) -> [f64; 3] {
    [dot3(xhat, g), dot3(yhat, g), dot3(zhat, g)]
}

// transform a global xyz position into a local coordinate system
#[inline(always)]
pub fn transform_fast(g: &Vec3, xhat: &Vec3, yhat: &Vec3, zhat: &Vec3) -> Vec3 {
    Vec3([xhat.dot(g), yhat.dot(g), zhat.dot(g)])
}
