#![allow(unused)]
use crate::MU0_4PI;
use crate::{
    archive::octree::{SourceNode, SourceOctree},
    math::{cross, distance, vec_distance},
};

// Compute the bfield contributions from a leaf node at a target point
// by directly integrating each source in the leaf
//
// A leaf may contain multipole sources
// TODO: from testing, leaf node should likely only have one source point
pub fn bfield_leaf(
    centroids: (&[f64], &[f64], &[f64]),
    vj: (&[f64], &[f64], &[f64]),
    target: &[f64; 3],
) -> [f64; 3] {
    // Length checks for auto-vectorization
    let n = centroids.0.len();
    assert_eq!(n, centroids.1.len());
    assert_eq!(n, centroids.2.len());
    assert_eq!(n, vj.0.len());
    assert_eq!(n, vj.1.len());
    assert_eq!(n, vj.2.len());

    let mut b = [0.0; 3];
    let x = target[0];
    let y = target[1];
    let z = target[2];
    let (centx, centy, centz) = centroids;
    let (vjx, vjy, vjz) = vj;

    // Outer loop over sources
    for i in 0..n {
        let rx: f64 = x - centx[i];
        let ry: f64 = y - centy[i];
        let rz: f64 = z - centz[i];
        let r: f64 = (rx * rx + ry * ry + rz * rz).sqrt();

        let jxrpx: f64 = vjy[i] * rz - vjz[i] * ry;
        let jxrpy: f64 = vjz[i] * rx - vjx[i] * rz;
        let jxrpz: f64 = vjx[i] * ry - vjy[i] * rx;

        let mask = if r > 1e-4 { 1.0 } else { 0.0 };
        let r3 = r * r * r + (1.0 - mask); // avoid div by zero
        let constant = MU0_4PI * mask / r3;
        b[0] += constant * jxrpx;
        b[1] += constant * jxrpy;
        b[2] += constant * jxrpz;
    }

    b
}

// Direct calculation of leaf contributions for only a single source in the leaf
#[inline(always)]
fn bfield_leaf_single(centroid: &[f64; 3], vj: &[f64; 3], target: &[f64; 3]) -> [f64; 3] {
    let r = vec_distance(centroid, target);
    let rmag = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
    let mut b = [0.0; 3];
    if rmag > 1e-4 {
        let invrmag3 = (1.0 / rmag).powi(3);
        cross(vj, &r, &mut b);
        for i in 0..b.len() {
            b[i] *= MU0_4PI * invrmag3;
        }
        b
    } else {
        b
    }
}

/// Recursively compue the contributions of a node and all child nodes at a target point
///
/// # Arguments
/// - `tree`:               the octree to traverse, containing biot-savart sources
/// - `current_index`:      the current node
/// - `target`:             (m) the location in 3D space at which to calculate fields
/// - `theta`:              Barnes-Hut angle opening parameter
///
/// # Returns
/// (T) magnetic flux density at the target point
pub fn bfield_node(
    tree: &SourceOctree,
    current_index: u32,
    target: &[f64; 3],
    theta: f64,
) -> [f64; 3] {
    // Accumulator
    let mut b = [0.0; 3];

    match tree.nodes[current_index as usize] {
        // For branch nodes, we perform the BH acceptance test
        // if we pass, recursion stops here and we compute the 'super-source' contribution of the node
        SourceNode::Branch {
            level: _,
            size,
            children,
            centroid,
            vj,
        } => {
            let d = distance(&centroid, target);
            if d * theta > size {
                // Accept node, compute contribution
                let r = vec_distance(&centroid, target);
                let mut _b = [0.0; 3];
                cross(&vj, &r, &mut _b);
                let invr3 = 1.0 / d.powi(3);
                b[0] += _b[0] * MU0_4PI * invr3;
                b[1] += _b[1] * MU0_4PI * invr3;
                b[2] += _b[2] * MU0_4PI * invr3;
                b
            } else {
                for child in children {
                    if child > 0 {
                        let _b = bfield_node(tree, child, target, theta);
                        b[0] += _b[0];
                        b[1] += _b[1];
                        b[2] += _b[2];
                    }
                }
                b
            }
        }
        // Leaves require direct calculation of contribution from discrete sources
        SourceNode::Leaf {
            level: _,
            source_range,
            ..
        } => {
            let (i, j) = (source_range.0 as usize, source_range.1 as usize);
            let _b = bfield_leaf(
                (
                    &tree.sources.xg[i..j],
                    &tree.sources.yg[i..j],
                    &tree.sources.zg[i..j],
                ),
                (
                    &tree.sources.vjx[i..j],
                    &tree.sources.vjy[i..j],
                    &tree.sources.vjz[i..j],
                ),
                target,
            );
            // let centroid = [tree.sources.xg[i], tree.sources.yg[i], tree.sources.zg[i]];
            // let vj = [tree.sources.vjx[i], tree.sources.vjy[i], tree.sources.vjz[i]];
            // let _b = bfield_leaf_single(&centroid, &vj, target);
            b[0] += _b[0];
            b[1] += _b[1];
            b[2] += _b[2];
            b
        }
    }
}

/// Compute the magnetic field at (x,y,z) using the barnes-hut (octree) method
///
/// # Arguments
/// - `centx`, `centy`, `centz`: (m) locations of source element centroids in 3D space
/// - `vol`:                     (m^3) volume of each source element
/// - `jx`, `jy`, `jz`:          (A/m^2) current density vector of each source element
/// - `x`, `y`, `z`:             (m) location of each target point
/// - `bx`, `by`, `bz`:          (T) magnetic flux density at each target point
/// - `theta`:                   Barnes-Hut angle-opening parameter (recommended < 0.5)
/// - `leaf_threshold`:          number of source points in each leaf (recommended = 1)
pub fn bfield_octree(
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
    theta: f64,
    leaf_threshold: u32,
) -> Result<(), ()> {
    // Build the source octree
    let max_depth: u8 = 21;
    let tree = SourceOctree::from_source_points(
        (centx, centy, centz),
        vol,
        (jx, jy, jz),
        max_depth,
        leaf_threshold,
    )
    .unwrap();
    let n = x.len();
    // TODO: length checks
    for i in 0..n {
        let centroid = [x[i], y[i], z[i]];
        let b = bfield_node(&tree, 0, &centroid, theta);
        bx[i] += b[0];
        by[i] += b[1];
        bz[i] += b[2];
    }
    Ok(())
}
