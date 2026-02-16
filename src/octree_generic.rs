#![allow(unused)]
// Rebuild the octree implementation so that it's generic over different types of sources in the tree

use crate::octree::BoundingBox;

/// Previous node definition
pub enum Node {
    Branch {
        level: u8, 
        size: f64,              // or maybe f32?
        children: [u32; 8], 
        centroid: [f64; 3],     // of the cluster of sources it represents
        vj: [f64; 3],           
    },
    Leaf {
        level: u8,              // maybe unecessary but harmless
        source_range: (u32, u32), 
        centroid: [f64; 3]
    }
}


pub struct Octree<S> {
    pub nodes: Vec<Node>,
    pub codes: Vec<u64>, 
    pub bbox: BoundingBox,
    pub sources: S
}

impl <S: Sources> Octree<S> {
    fn build_from_sources(s: S) -> Self {}
}

impl <S: HFieldSolver> Octree<S> {
    fn h_field(&self, targets: (&[f64], &[f64], &[f64]), h: (&mut [f64], &mut [f64], &mut [f64])) -> () {}
}

impl <S: GradHFieldSolver> Octree<S> {
    // need to figure out proper layout for gradh, should it be 9 mutable slices? or something like:
    // gradh: &mut [[[f64; 3]; 3]]
    // perhaps it needs its own type?
    fn gradh_field(&self, targets: (&[f64], &[f64], &[f64]), gradh: (&mut [f64], &mut [f64], &mut [f64])) -> () {}
}

pub struct PointSources {
    pub xg: Vec<f64>,
    pub yg: Vec<f64>,
    pub zg: Vec<f64>,
    pub vjx: Vec<f64>,
    pub vjy: Vec<f64>,
    pub vjz: Vec<f64>,
}

pub struct TetSources {
    pub nodes: Vec<[[f64; 3]; 4]>,
    pub jdensity: Vec<[f64; 3]>,
    pub centroid: Vec<[f64; 3]>,
    pub volume: Vec<f64>,
}

pub trait Sources {
    fn len(&self) -> usize;
    fn centroid(&self, i: usize) -> [f64; 3];
    fn moment(&self, i: usize) -> [f64; 3];         // vj for current sources, dipole for dipole sources
    fn sort(&mut self, indices: &[usize]);
}

pub struct DipoleSources<S>(pub S);     // for magnetic materials

pub struct CurrentSources<S>(pub S);    // for current-carrying conductors

pub trait HFieldSolver {
    fn h_field_branch(&self, centroid: &[f64;3], vj: &[f64;3], target: &[f64;3]) -> [f64;3];
    fn h_field_leaf(&self, start: usize, end: usize, target: &[f64;3]) -> [f64;3];
}

pub trait DipoleHFieldSolver {
    fn h_field(&self, targets: (&[f64], &[f64], &[f64]), h: (&[f64], &[f64], &[f64]));
}

pub trait GradHFieldSolver {
    fn gradh_field(&self, targets: (&[f64], &[f64], &[f64]), h: (&[f64], &[f64], &[f64]));
}

pub trait AFieldSolver {
    fn a_field(&self, targets: (&[f64], &[f64], &[f64]), h: (&[f64], &[f64], &[f64]));
}

impl HFieldSolver for CurrentSources<PointSources> {
    fn h_field_branch(&self, centroid: &[f64;3], vj: &[f64;3], target: &[f64;3]) -> [f64;3] {

    }

    fn h_field_leaf(&self, start: usize, end: usize, target: &[f64;3]) -> [f64;3] {

    };
}

