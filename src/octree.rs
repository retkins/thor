
use crate::morton;
use crate::math::{min_and_max, sort_by_indices};

/// Determines the location and extent of a collection of source points
#[derive(Debug)]
struct BoundingBox {
    xc: f64,
    yc: f64,
    zc: f64,
    side_length: f64,
    xbounds: (f64, f64),
    ybounds: (f64, f64),
    zbounds: (f64, f64),
}


// Get the maximum side length of a bounding box cube that encloses all x,y,z bounds
fn side_length_from_bounds(xbounds: (f64, f64), ybounds: (f64, f64), zbounds: (f64, f64)) -> f64 {
    let xrange: f64 = xbounds.1 - xbounds.0;
    let yrange: f64 = ybounds.1 - ybounds.0;
    let zrange: f64 = zbounds.1 - zbounds.0;

    let xymax = match xrange > yrange {
        true => xrange,
        false => yrange,
    };

    match xymax > zrange {
        true => xymax,
        false => zrange,
    }
}


impl BoundingBox {
    fn from_centroids(centroids: (&[f64], &[f64], &[f64])) -> Option<Self> {
        // TODO: length check

        let xbounds: Option<(f64, f64)> = min_and_max(centroids.0);
        let ybounds: Option<(f64, f64)> = min_and_max(centroids.1);
        let zbounds: Option<(f64, f64)> = min_and_max(centroids.2);

        if xbounds.is_none() | ybounds.is_none() | zbounds.is_none() {
            return None;
        }

        let (mut xb, mut yb, mut zb) = (xbounds.unwrap(), ybounds.unwrap(), zbounds.unwrap());
        let side_length = side_length_from_bounds(xb, yb, zb);
        let xc: f64 = xb.0 + 0.5 * side_length;
        let yc: f64 = yb.0 + 0.5 * side_length;
        let zc: f64 = zb.0 + 0.5 * side_length;

        xb.0 = xc - 0.5 * side_length;
        xb.1 = xc + 0.5 * side_length;
        yb.0 = yc - 0.5 * side_length;
        yb.1 = yc + 0.5 * side_length;
        zb.0 = zc - 0.5 * side_length;
        zb.1 = zc + 0.5 * side_length;

        Some(Self {
            xc: xc,
            yc: yc,
            zc: zc,
            side_length: side_length,
            xbounds: xb,
            ybounds: yb,
            zbounds: zb,
        })
    }
}


/// A collection of sources that can be sorted by morton codes
pub struct Sources {
    codes: Vec<u64>, // Morton code for each source point
    xg: Vec<f64>,    // Geometric location of the source point
    yg: Vec<f64>,
    zg: Vec<f64>,
    vjx: Vec<f64>, // Current density moment product
    vjy: Vec<f64>,
    vjz: Vec<f64>,
    bbox: Option<BoundingBox>,
}


impl Sources {
    /// Basic constructor; zeroed with initial size
    pub fn new(n: usize) -> Self {
        Self {
            codes: Vec::with_capacity(n),
            xg: Vec::with_capacity(n),
            yg: Vec::with_capacity(n),
            zg: Vec::with_capacity(n),
            vjx: Vec::with_capacity(n),
            vjy: Vec::with_capacity(n),
            vjz: Vec::with_capacity(n),
            bbox: None,
        }
    }

    /// Add source points (elements) to the octree
    ///
    /// Note that this involves a data copy. We're going to be sorting the arrays
    /// and need to own the memory. (?)
    pub fn from_source_points(
        centroids: (&[f64], &[f64], &[f64]),
        volumes: &[f64],
        jdensity: (&[f64], &[f64], &[f64]),
        max_depth: u8
    ) -> Option<Self> {
        // TODO: length check on all vectors
        let n: usize = centroids.0.len();

        let _bbox: Option<BoundingBox> = BoundingBox::from_centroids(centroids);
        if _bbox.is_none() {
            return None;
        }
        let bbox = _bbox.unwrap();

        let (x, y, z) = centroids;
        let (jx, jy, jz) = jdensity;
        let mut codes: Vec<u64> = Vec::with_capacity(n);
        let mut xg: Vec<f64> = Vec::with_capacity(n);
        let mut yg: Vec<f64> = Vec::with_capacity(n);
        let mut zg: Vec<f64> = Vec::with_capacity(n);
        let mut vjx: Vec<f64> = Vec::with_capacity(n);
        let mut vjy: Vec<f64> = Vec::with_capacity(n);
        let mut vjz: Vec<f64> = Vec::with_capacity(n);

        let scale: f64 = morton::calculate_scale_factor(max_depth as u32);
        let min_corner: (f64, f64, f64) = (bbox.xbounds.0, bbox.ybounds.0, bbox.zbounds.0);

        for i in 0..n {
            let pt: (f64, f64, f64) = (x[i], y[i], z[i]);
            let vol: f64 = volumes[i];

            codes.push(morton::encode(pt, scale, bbox.side_length, min_corner));
            xg.push(pt.0);
            yg.push(pt.1);
            zg.push(pt.2);
            vjx.push(vol * jx[i]);
            vjy.push(vol * jy[i]);
            vjz.push(vol * jz[i]);
        }

        return Some(Sources {
            codes: codes,
            xg: xg,
            yg: yg,
            zg: zg,
            vjx: vjx,
            vjy: vjy,
            vjz: vjz,
            bbox: Some(bbox),
        });
    }

    /// Sort the octree by morton order
    pub fn morton_sort(&mut self) {

        // Get the indices to sort by from the codes
        let mut indices: Vec<usize> = (0..self.codes.len()).collect();
        indices.sort_by(|&i, &j| self.codes[i].cmp(&self.codes[j]));

        // Now sort all of the other data arrays
        let mut scratch = vec![0.0; self.codes.len()]; 

        sort_by_indices(&mut self.xg, &mut scratch, &indices);
        sort_by_indices(&mut self.yg, &mut scratch, &indices);
        sort_by_indices(&mut self.zg, &mut scratch, &indices);
        sort_by_indices(&mut self.vjx, &mut scratch, &indices);
        sort_by_indices(&mut self.vjy, &mut scratch, &indices);
        sort_by_indices(&mut self.vjz, &mut scratch, &indices);

        // Now sort the codes themselves
        let mut scratch_codes: Vec<u64> = vec![0;self.codes.len()];
        sort_by_indices(&mut self.codes, &mut scratch_codes, &indices);
    }
}


/// Represents a single node in the octree
/// 
/// Contains both geometric and physics data. Future optimizations might break it up.
pub struct Node {
    // Tree structure
    n_children: u8,         
    first_child: u32,       // we check n_children first; index into the Nodes array
    level: u8, 

    // Indexes into the sources array 
    source_start: u32, 
    source_end: u32,

    // Geometry - this determines the acceptance criteria 
    center: [f64; 3], 
    size: f64,

    // Physics
    centroid: [f64; 3], 
    vj: [f64; 3]
}


pub struct SourceOctree {
    max_depth: u8, 
    // we use sources to build the octree and then do not need it again
    sources: Sources, 
    nodes: Vec<Node>, 
    leaf_threshold: u32
}

impl SourceOctree {

    /// Construct a new octree using source point definitions
    pub fn from_source_points(
        centroids: (&[f64], &[f64], &[f64]),
        volumes: &[f64],
        jdensity: (&[f64], &[f64], &[f64]),
        max_depth: u8
    ) -> Option<Self> {
        // TODO: length check

        let n_sources: usize = centroids.0.len(); 
        let sources = match Sources::from_source_points(centroids, volumes, jdensity, max_depth) {
            Some(s) => s, 
            None => {return None}
        };

        // Memory usage will be significantly more than n_sources, as 
        // every source point has a leaf node, but this is a good place to start
        let mut nodes: Vec<Node> = Vec::with_capacity(n_sources);
        let leaf_threshold: u32 = 32;

        make_nodes(&sources, &mut nodes, leaf_threshold);

        Some(Self {
            max_depth: max_depth, 
            sources: sources, 
            nodes: nodes, 
            leaf_threshold: leaf_threshold
        })
    }
}


// Build the tree
fn make_nodes(sources: &Sources, nodes: &mut Vec<Node>, leaf_threshold: u32) -> Result<(), ()> {

    let n_sources: usize = sources.codes.len(); 


    Ok(())
}   



#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_octree_build_1() {
        // test along the diagonal
        let n = 5;
        let x = vec![-0.2, -0.1, 0.0, 0.1, 0.2];
        let y = vec![-0.2, -0.1, 0.0, 0.1, 0.2];
        let z = vec![-0.2, -0.1, 0.0, 0.1, 0.2];
        let volumes = vec![1e-4; n];
        let jx = vec![1e8; n];
        let jy = vec![-1e8; n];
        let jz = vec![2e8; n];
        let max_depth: u8 = 21;

        let mut tree = Sources::from_source_points((&x, &y, &z), &volumes, (&jx, &jy, &jz), max_depth).unwrap();

        println!("Bounding box: {:?}", tree.bbox);
        println!("Max depth: {}", 21);
        println!("Morton codes:");
        for (i, code) in tree.codes.iter().enumerate() {
            println!("Point: ({:.2}, {:.2}, {:.2})", x[i], y[i], z[i]);
            println!("\tCode (decimal): {}", code);
            println!("\tCode (binary):  {:064b}", code);
            println!("\tCode (hex):     0x{:016x}", code);
        }
        tree.morton_sort();
        println!("After sorting: ");
        for (i, code) in tree.codes.iter().enumerate() {
            println!("Point: ({:.2}, {:.2}, {:.2})", x[i], y[i], z[i]);
            println!("\tCode (decimal): {}", code);
            println!("\tCode (binary):  {:064b}", code);
            println!("\tCode (hex):     0x{:016x}", code);
        }
    }

    #[test]
    fn test_octree_build_2() {
        // keep xy=0, increase z
        let n = 5;
        let x = vec![0.0; n];
        let y = vec![0.0; n];
        let z = vec![-0.2, -0.1, 0.0, 0.1, 0.2];
        let volumes = vec![1e-4; n];
        let jx = vec![1e8; n];
        let jy = vec![-1e8; n];
        let jz = vec![2e8; n];
        let max_depth: u8 = 21;

        let mut tree = Sources::from_source_points((&x, &y, &z), &volumes, (&jx, &jy, &jz), max_depth).unwrap();

        println!("Bounding box: {:?}", tree.bbox);
        println!("Max depth: {}", 21);
        println!("Morton codes:");
        for (i, code) in tree.codes.iter().enumerate() {
            println!("Point: ({:.2}, {:.2}, {:.2})", x[i], y[i], z[i]);
            println!("\tCode (decimal): {}", code);
            println!("\tCode (binary):  {:064b}", code);
            println!("\tCode (hex):     0x{:016x}", code);
        }

        tree.morton_sort();
        println!("After sorting: ");
        for (i, code) in tree.codes.iter().enumerate() {
            println!("Point: ({:.2}, {:.2}, {:.2})", x[i], y[i], z[i]);
            println!("\tCode (decimal): {}", code);
            println!("\tCode (binary):  {:064b}", code);
            println!("\tCode (hex):     0x{:016x}", code);
        }

    }

    #[test]
    fn test_tree_build() {
        let n: usize = 7;
        let x = vec![0.0, 0.4, 0.9, 0.6, 0.6, 0.1, 1.0];
        let y = vec![0.0, 0.2, 0.4, 0.6, 0.9, 0.7, 1.0];
        let z = vec![0.0, 0.6, 0.6, 0.6, 0.6, 0.6, 1.0];
        let volumes = vec![1e-4; n];
        let jx = vec![1e8; n];
        let jy = vec![2e7; n];
        let jz = vec![3e6; n];

        let max_depth: u8 = 21; 

        let mut tree = SourceOctree::from_source_points((&x, &y, &z), &volumes, (&jx, &jy, &jz), max_depth).unwrap();

        let i: usize = 1;
        println!("code:\n{:b}", tree.sources.codes[i]);
        for k in 0..max_depth {
            let prefix = tree.sources.codes[i] >> (3*(max_depth-k));
            println!("k = {}, (x,y,z) = ({:.1}, {:.1}, {:.1}), prefix={:03b}", k, x[i], y[i], z[i], prefix);
        }

    }
}
