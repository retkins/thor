
use crate::octree::BoundingBox;
use crate::octree_generic::{
    Sources, HFieldSolver, CurrentSources
};
use crate::sources::h_point;
use crate::math::sort_by_indices;
use crate::morton;
use std::f64::consts::PI;


pub struct PointSources {
    pub xg: Vec<f64>,
    pub yg: Vec<f64>,
    pub zg: Vec<f64>,
    pub r: Vec<f64>,
    pub vjx: Vec<f64>,
    pub vjy: Vec<f64>,
    pub vjz: Vec<f64>,
    bbox: BoundingBox
}

impl PointSources {
    pub fn new(x: &[f64], y: &[f64], z: &[f64], vol: &[f64], jx: &[f64], jy: &[f64], jz: &[f64]) -> Self {
        let n = x.len();
        let xg: Vec<f64> = x.to_vec();
        let yg: Vec<f64> = y.to_vec();
        let zg: Vec<f64> = z.to_vec();
        let mut r  = vec![0.0; n];
        let mut vjx  = vec![0.0; n];
        let mut vjy  = vec![0.0; n];
        let mut vjz  = vec![0.0; n];

        // compute radius of each point and current density moment
        for i in 0..n {
            r[i] = ((3.0*PI/4.0)*vol[i]).powi(3);
            vjx[i] = vol[i]*jx[i];
            vjy[i] = vol[i]*jy[i];
            vjz[i] = vol[i]*jz[i];
        }
        let bbox = BoundingBox::from_centroids((&xg, &yg, &zg)).unwrap();

        Self { xg: xg, yg: yg, zg: zg, r: r, vjx: vjx, vjy: vjy, vjz: vjz, bbox }
    }
}


impl Sources for PointSources {
    fn len(&self) -> usize {
        self.xg.len()
    }

    fn centroid(&self, i: usize) -> [f64; 3] {
        [self.xg[i], self.yg[i], self.zg[i]]
    }

    fn moment(&self, i: usize) -> [f64; 3] {
        [self.vjx[i], self.vjy[i], self.vjz[i]]
    }

    fn sort(&mut self, indices: &[usize]) {
        let mut scratch = vec![0.0; self.len()];
        sort_by_indices(&mut self.xg, &mut scratch, &indices);
        sort_by_indices(&mut self.yg, &mut scratch, &indices);
        sort_by_indices(&mut self.zg, &mut scratch, &indices);
        sort_by_indices(&mut self.r, &mut scratch, &indices);
        sort_by_indices(&mut self.vjx, &mut scratch, &indices);
        sort_by_indices(&mut self.vjy, &mut scratch, &indices);
        sort_by_indices(&mut self.vjz, &mut scratch, &indices);
        
    }

    fn bbox(&self) -> &BoundingBox {
        &self.bbox
    }

    fn encode(&mut self, max_depth: u8) -> (&BoundingBox, Vec<u64>) {
        let n = self.len();
        let mut codes: Vec<u64> = Vec::with_capacity(n);
        let bbox = self.bbox();
        let scale: f64 = morton::calculate_scale_factor(max_depth as u32);
        let min_corner: (f64, f64, f64) = bbox.min_corner();

        for i in 0..n {
            let pt: (f64, f64, f64) = (self.xg[i], self.yg[i], self.zg[i]);
            codes.push(morton::encode(pt, scale, bbox.side_length, min_corner));
        }

        (bbox, codes)
    }
}

impl HFieldSolver for CurrentSources<PointSources> {
    fn h_field_branch(&self, centroid: &[f64;3], vj: &[f64;3], target: &[f64;3]) -> [f64;3] {
        let radius = 0.0;
        h_point(centroid, vj, radius, target)
    }

    fn h_field_leaf(&self, start: usize, end: usize, target: &[f64;3]) -> [f64;3] {
        let mut h = [0.0; 3]; 

        for i in start..end {
            let centroid = [self.0.xg[i], self.0.yg[i], self.0.zg[i]];
            let vj = [self.0.vjx[i], self.0.vjy[i], self.0.vjz[i]];
            let _h = h_point(&centroid, &vj, self.0.r[i], target);
            h[0] += _h[0]; 
            h[1] += _h[1]; 
            h[2] += _h[2];
        }
        h
    }
}