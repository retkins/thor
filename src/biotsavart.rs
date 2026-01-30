#![allow(non_snake_case, unused)]


use crate::{MU0_4PI}; 
use crate::math::{cross, distance, vec_distance};
use crate::octree;
use crate::octree::{SourceOctree, Node};

/// Compute the magnetic field at target points (x, y, z) using a direct (O(N^2)) Biot-Savart summation
pub fn bfield_direct(
    centx: &[f64], centy: &[f64], centz: &[f64],
    vol: &[f64], 
    jx: &[f64], jy: &[f64], jz: &[f64], 
    x: &[f64], y: &[f64], z: &[f64], 
    Bx: &mut [f64], By: &mut [f64], Bz: &mut [f64], 
) {

    // TODO: length checks on input arrays

    // Outer loop over source elements
    for i in 0..centx.len(){
        
        // Hoist invariants out of the loop explicitly
        let centxi: f64 = centx[i]; 
        let centyi: f64 = centy[i];
        let centzi: f64 = centz[i];
        let vol_mu0_4pi: f64 = vol[i]*MU0_4PI;
        let jxi = jx[i]; 
        let jyi = jy[i]; 
        let jzi = jz[i];

        // Inner loop over target points
        for (((xj, yj), zj), ((bxj, byj), bzj)) in x.iter().zip(y.iter()).zip(z.iter()).zip(Bx.iter_mut().zip(By.iter_mut()).zip(Bz.iter_mut())) {

            // Vector from the element centroid to the target point: r'
            let rx: f64 = xj - centxi; 
            let ry: f64 = yj - centyi; 
            let rz: f64 = zj - centzi; 
            let r: f64 = (rx*rx + ry*ry + rz*rz).sqrt(); 

            // J x r'
            let jxrpx: f64 = jyi*rz - jzi*ry; 
            let jxrpy: f64 = jzi*rx - jxi*rz; 
            let jxrpz: f64 = jxi*ry - jyi*rx; 

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
}


// Older version of the above that uses explicit indices for the inner loop
fn bfield_direct_old(
    centx: &[f64], centy: &[f64], centz: &[f64],
    vol: &[f64], 
    jx: &[f64], jy: &[f64], jz: &[f64], 
    x: &[f64], y: &[f64], z: &[f64], 
    Bx: &mut [f64], By: &mut [f64], Bz: &mut [f64], 
) {

    let m: usize = centx.len();
    let n: usize = x.len();

    for i in 0..m {
        
        let centxi: f64 = centx[i]; 
        let centyi: f64 = centy[i];
        let centzi: f64 = centz[i];
        let vol_mu0_4pi: f64 = vol[i]*MU0_4PI;
        let jxi = jx[i]; 
        let jyi = jy[i]; 
        let jzi = jz[i];

        for j in 0..n {

            let rx: f64 = x[j] - centxi; 
            let ry: f64 = y[j] - centyi; 
            let rz: f64 = z[j] - centzi; 
            let r: f64 = (rx*rx + ry*ry + rz*rz).sqrt(); 

            let jxrpx: f64 = jyi*rz - jzi*ry; 
            let jxrpy: f64 = jzi*rx - jxi*rz; 
            let jxrpz: f64 = jxi*ry - jyi*rx; 

            let mask = if r > 1e-4 { 1.0 } else { 0.0 };
            let r3 = r * r * r + (1.0 - mask); // avoid div by zero
            let constant = vol_mu0_4pi * mask / r3;
            Bx[j] += constant * jxrpx; 
            By[j] += constant * jxrpy; 
            Bz[j] += constant * jxrpz;
        }
    }
}


/// Compute the bfield contributions from a leaf node at a target point 
/// by directly integrating each source in the leaf
/// 
/// A leaf may contain multipole sources
/// TODO: from testing, leaf node should likely only have one source point
fn bfield_leaf(
    centroids: (&[f64], &[f64], &[f64]), 
    vj: (&[f64], &[f64], &[f64]), 
    target: &[f64; 3]
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
            let r: f64 = (rx*rx + ry*ry + rz*rz).sqrt(); 

            let jxrpx: f64 = vjy[i]*rz - vjz[i]*ry; 
            let jxrpy: f64 = vjz[i]*rx - vjx[i]*rz; 
            let jxrpz: f64 = vjx[i]*ry - vjy[i]*rx; 

            let mask = if r > 1e-4 { 1.0 } else { 0.0 };
            let r3 = r * r * r + (1.0 - mask); // avoid div by zero
            let constant = MU0_4PI * mask / r3;
            b[0] += constant * jxrpx; 
            b[1] += constant * jxrpy; 
            b[2] += constant * jxrpz;
        }

    return b;
}

// Direct calculation of leaf contributions for only a single source in the leaf
#[inline(always)]
fn bfield_leaf_single(
    centroid: &[f64; 3], 
    vj: &[f64; 3], 
    target: &[f64; 3]
) -> [f64; 3] {
    let r = vec_distance(centroid, target); 
    let rmag = (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]).sqrt();
    let mut b = [0.0; 3]; 
    if rmag > 1e-4 {
        let invrmag3 = (1.0/rmag).powi(3);
        cross(&vj, &r, &mut b); 
        for i in 0..b.len() {
            b[i] *= MU0_4PI * invrmag3;
        }
        return b;
    }
    else {
        return b;
    }


}

// Recursively call to compute the contributions of a node and all child nodes at a target point
fn bfield_node(
    tree: &SourceOctree, 
    current_index: u32, 
    target: &[f64; 3], 
    theta: f64
) -> [f64; 3] {

    // Accumulator
    let mut b = [0.0; 3];

    match tree.nodes[current_index as usize] {

        // For interior nodes, we perform the BH acceptance test
        // if we pass, recursion stops here and we compute the 'super-source' contribution of the node
        Node::Interior { level: _, size, children, centroid, vj } => {
            let d = distance(&centroid, target);
            if d*theta > size {
                // Accept node, compute contribution
                let r = vec_distance(&centroid, target);
                let mut _b = [0.0; 3];
                cross(&vj, &r, &mut _b);
                let invr3 = 1.0 / d.powi(3);
                b[0] += _b[0]*MU0_4PI * invr3;
                b[1] += _b[1]*MU0_4PI * invr3;
                b[2] += _b[2]*MU0_4PI * invr3;
                return b;
            } 
            else {
                for child in children {
                    if child > 0 {
                        let _b = bfield_node(tree, child, target, theta);
                        b[0] += _b[0];
                        b[1] += _b[1]; 
                        b[2] += _b[2];
                    }
                }
                return b;
            }
        },
        // Leaves require direct calculation of contribution from discrete sources
        Node::Leaf { level: _, source_range } => {
            let (i, _) = (source_range.0 as usize, source_range.1 as usize);
            // let _b = bfield_leaf(
            //     (&tree.sources.xg[i..j], &tree.sources.yg[i..j], &tree.sources.zg[i..j]), 
            //     (&tree.sources.vjx[i..j], &tree.sources.vjy[i..j], &tree.sources.vjz[i..j]), 
            //     target
            // );
            let centroid = [tree.sources.xg[i], tree.sources.yg[i], tree.sources.zg[i]];
            let vj = [tree.sources.vjx[i], tree.sources.vjy[i], tree.sources.vjz[i]];
            let _b = bfield_leaf_single(&centroid, &vj, target);
            b[0] += _b[0];
            b[1] += _b[1]; 
            b[2] += _b[2];
            b
        }
    }
}

/// Compute the magnetic field at (x,y,z) using the barnes-hut (octree) method
pub fn bfield_octree(
    centx: &[f64], centy: &[f64], centz: &[f64],
    vol: &[f64], 
    jx: &[f64], jy: &[f64], jz: &[f64], 
    x: &[f64], y: &[f64], z: &[f64], 
    Bx: &mut [f64], By: &mut [f64], Bz: &mut [f64], 
    theta: f64, leaf_threshold: u32
) {
    // Build the source octree 
    let max_depth: u8 = 21; 
    let tree = octree::SourceOctree::from_source_points(
        (&centx, &centy, &centz), vol, (&jx, &jy, &jz), max_depth, leaf_threshold
    ).unwrap();
    let n = x.len(); 
    // TODO: length checks
    for i in 0..n {
        let centroid = [x[i], y[i], z[i]];
        let b = bfield_node(&tree, 0, &centroid, theta);
        Bx[i] += b[0]; 
        By[i] += b[1]; 
        Bz[i] += b[2];
    }
}

#[cfg(test)]
mod tests{

    #[test]
    fn test_direct() {
        

    }
}
