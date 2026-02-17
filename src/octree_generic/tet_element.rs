
use crate::{
    math::sort_by_indices, morton, octree::BoundingBox, octree_generic::{
        CurrentSources, HFieldSolver, Sources
    }, sources::{hfield_tetrahedron, h_point}, vec3::Vec3
};


pub struct TetSources {
    pub nodes: Vec<[Vec3; 4]>,      // of the element, not the tree
    pub centroids: Vec<Vec3>,
    pub volumes: Vec<f64>,
    pub jdensity: Vec<Vec3>,
    pub bbox: BoundingBox
}

impl TetSources {

    /// Constructor 
    pub fn new(
        nodes_flat: &[f64], 
        centroids_flat: &[f64],
        vol: &[f64], 
        jdensity_flat: &[f64]
    ) -> Self {

        // TODO: length check on inputs 
        let n: usize = vol.len(); 

        let mut nodes: Vec<[Vec3; 4]> = vec![[Vec3([0.0;3]); 4]; n];
        let mut centroids: Vec<Vec3> = vec![Vec3([0.0;3]); n];
        let volumes: Vec<f64> = vol.to_vec();
        let mut jdensity: Vec<Vec3> = vec![Vec3([0.0;3]); n];

        for i in 0..n {
            let node_slice = &nodes_flat[12*i..12*i+12];
            for j in 0..4 {
                for k in 0..3 {
                    nodes[i][j][k] = node_slice[3*j+k];
                }
            }
            let j = 3 * i;
            centroids[i] = Vec3([centroids_flat[j], centroids_flat[j+1], centroids_flat[j+2]]);
            jdensity[i] = Vec3([jdensity_flat[j], jdensity_flat[j+1], jdensity_flat[j+2]]);
        }

        let bbox = BoundingBox::from_centroids_vec(&centroids);
        Self { nodes: nodes, centroids: centroids, volumes: volumes, jdensity: jdensity, bbox: bbox}
    }
}

impl Sources for TetSources {
    fn len(&self) -> usize {
        return self.volumes.len(); 
    }

    fn centroid(&self, i: usize) -> Vec3 {
        self.centroids[i]
    }
    
    fn moment(&self, i: usize) -> Vec3 {
        self.jdensity[i] * self.volumes[i]
    }

    fn sort(&mut self, indices: &[usize]) {
        let n = self.len();
        let mut scratch_nodes = vec![[Vec3([0.0;3]); 4]; n];
        sort_by_indices(&mut self.nodes, &mut scratch_nodes, indices);

        let mut scratch_vecs = vec![Vec3([0.0;3]); n];
        sort_by_indices(&mut self.centroids, &mut scratch_vecs, indices);
        sort_by_indices(&mut self.jdensity, &mut scratch_vecs, indices);

        let mut scratch_vol = vec![0.0; n];
        sort_by_indices(&mut self.volumes, &mut scratch_vol, indices);
    }

    fn bbox(&self) -> &crate::octree::BoundingBox {
        &self.bbox
    }

    fn encode(&mut self, max_depth: u8) -> (&BoundingBox, Vec<u64>) {
        let n = self.len();
        let mut codes: Vec<u64> = Vec::with_capacity(n);
        let bbox = self.bbox();
        let scale: f64 = morton::calculate_scale_factor(max_depth as u32);
        let min_corner: (f64, f64, f64) = bbox.min_corner();

        for i in 0..n {
            let pt: (f64, f64, f64) = (self.centroids[i][0], self.centroids[i][1], self.centroids[i][2]);
            codes.push(morton::encode(pt, scale, bbox.side_length, min_corner));
        }

        (bbox, codes)
    }

}

impl HFieldSolver for CurrentSources<TetSources> {
    fn h_field_branch(&self, centroid: &Vec3, vj: &Vec3, target: &Vec3) -> Vec3 {
        let radius = 0.0;
        h_point(centroid, vj, radius, target)
    }

    fn h_field_leaf(&self, start: usize, end: usize, target: &Vec3) -> Vec3 {
        let mut h = Vec3([0.0; 3]);
        for i in start..end {

            let nx = [self.0.nodes[i][0][0], self.0.nodes[i][1][0], self.0.nodes[i][2][0], self.0.nodes[i][3][0]];
            let ny = [self.0.nodes[i][0][1], self.0.nodes[i][1][1], self.0.nodes[i][2][1], self.0.nodes[i][3][1]];
            let nz = [self.0.nodes[i][0][2], self.0.nodes[i][1][2], self.0.nodes[i][2][2], self.0.nodes[i][3][2]];

            h += hfield_tetrahedron(&nx, &ny, &nz, &self.0.jdensity[i], target)
        }
        h
    }
}


// impl HFieldSolver for DipoleSources<PointSources> {
//     fn h_field_branch(&self, centroid: &Vec3, moment: &Vec3, target: &Vec3) -> Vec3 {
//         let radius = 0.0;   // far-field
//         h_point_dipole(centroid, moment, radius, target)
//     }

//     fn h_field_leaf(&self, start: usize, end: usize, target: &Vec3) -> Vec3 {
//         let mut h = Vec3([0.0;3]);
//         for i in start..end {
//             let centroid = Vec3::from_slice_tuple((&self.0.xg, &self.0.yg, &self.0.zg), i);
//             let moment = Vec3::from_slice_tuple((&self.0.vjx, &self.0.vjy, &self.0.vjz), i);
//             h += h_point_dipole(&centroid, &moment, self.0.r[i], target);
//         }
//         h
//     }
// }