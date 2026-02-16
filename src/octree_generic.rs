#![allow(unused)]
// Rebuild the octree implementation so that it's generic over different types of sources in the tree

// these will be moved into this file eventually
use crate::octree::{BoundingBox, size_at_level, get_range_in_same_node, update_centroid};
use crate::math::{mag};


pub mod point;
pub mod tet_element;

pub trait Sources {
    fn len(&self) -> usize;
    fn bbox(&self) -> &BoundingBox;
    fn centroid(&self, i: usize) -> [f64; 3];
    fn moment(&self, i: usize) -> [f64; 3];         // vj for current sources, dipole for dipole sources
    fn sort(&mut self, indices: &[usize]);
    fn encode(&mut self, max_depth: u8) -> (&BoundingBox, Vec<u64>);
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

/// Previous node definition
pub enum Node {
    Branch {
        level: u8, 
        size: f64,              // or maybe f32?
        children: [u32; 8], 
        centroid: [f64; 3],     // of the cluster of sources it represents
        moment: [f64; 3],           
    },
    Leaf {
        level: u8,              // maybe unecessary but harmless
        source_range: (u32, u32), 
        centroid: [f64; 3]
    }
}


/// Recursively add nodes to the octree
/// 
/// Returns: the index of the node that it created
fn add_node<S: Sources+Copy>(
    sources: &S,
    codes: &Vec<u64>,
    nodes: &mut Vec<Node>,
    max_depth: u8,
    leaf_threshold: u32,
    start: usize,
    end: usize,
    level: u8,
) -> u32 {

    let current_index: usize;

    // Create a leaf if the number of sources in the range is small enough (dead-end for recursion)
    if end - start <= leaf_threshold as usize {
        let mut cx = 0.0; 
        let mut cy = 0.0; 
        let mut cz = 0.0;
        for i in start..end {
            let cent = sources.centroid(i);
            cx += cent[0]; 
            cy += cent[1];
            cz += cent[2];
        }
        cx /= (end - start) as f64;
        cy /= (end - start) as f64;
        cz /= (end - start) as f64;
        nodes.push(
            // Do not descend level
            Node::Leaf { level: level, source_range: (start as u32, end as u32), centroid: [cx, cy, cz] }
        );
        current_index = nodes.len() - 1;
    } 
    
    // Otherwise, create an branch node
    else {
        
        // Size of the node (cube side length) for BH acceptance criteria
        let size: f64 = size_at_level(sources.bbox().side_length, level);

        // Initialize the branch node first, then recursive calls later fill it
        nodes.push(
            Node::Branch { level: level, size: size, children: [0; 8], centroid: [0.0; 3], moment: [0.0;3] }
        );

        current_index = nodes.len() - 1;            // index of the node just created

        // Now recurse to create its children
        let mut child_indices: [u32; 8] = [0; 8];
        let mut n_children: usize = 0;
        let mut cursor = start;

        while cursor < end {
            // We descend to next level here, hence `level+1`
            // index to sources array
            let child_end = get_range_in_same_node(&codes, level+1, max_depth, cursor);

            // index in nodes array
            let child_idx = add_node(sources, codes, nodes, max_depth, leaf_threshold, cursor, child_end, level+1);
            child_indices[n_children] = child_idx;
            n_children += 1;
            cursor = child_end;     // index in sources array
        } 

        // Update the previously-initialized node with indices to its children
        // and the volume current density product 
        let mut parent_centroid: [f64; 3] = [0.0; 3];
        let mut parent_moment: [f64; 3] = [0.0; 3];
        
        for idx in child_indices {
            if idx < 1 {
                // Skip empty nodes
                break;
            }

            let mut parent_mag: f64 = mag(&parent_moment);

            match nodes[idx as usize] {
                Node::Branch { level: _, size: _, children: _, centroid, moment } => {
                    update_centroid(&mut parent_centroid, parent_mag, &centroid, mag(&moment));
                    for k in 0..3 as usize {
                        parent_moment[k] += moment[k];
                    }
                }, 
                Node::Leaf { level: _, source_range , centroid} => {
                    for _i in source_range.0..source_range.1  {
                        let i = _i as usize;
                        let moment = &sources.moment(i);
                        let moment_mag = mag(&[moment[0], moment[1], moment[2]]);
                        
                        update_centroid(&mut parent_centroid, parent_mag, &sources.centroid(i), moment_mag);
                        for k in 0..3 as usize {
                            parent_moment[k] += moment[k];
                        }
                        parent_mag = mag(&parent_moment);
                    }
                }
            }
        }

        if let Node::Branch { children, centroid, moment, ..} = &mut nodes[current_index] {
            *children = child_indices; 
            *centroid = parent_centroid;
            *moment = parent_moment;
        }
    }

    return current_index as u32;
}



pub struct Octree<S> {
    pub nodes: Vec<Node>,
    pub codes: Vec<u64>, 
    pub bbox: BoundingBox,
    pub sources: S
}

impl <S: Sources> Octree<S> {
    fn build_from_sources(s: S) -> Self {
        let max_depth: u8 = 21; 
        let leaf_threshold: u32 = 1;
        let (bbox, codes) = s.encode(max_depth);
        
        // sort the sources by morton code
        let mut indices: Vec<usize> = (0..codes.len()).collect();
        indices.sort_by(|&i, &j| codes[i].cmp(&codes[j]));
        s.sort(&&indices);      // why double borrow here? 

        // make nodes 
        let mut nodes: Vec<Node> = Vec::with_capacity(s.len());
        let start = 0; 
        let end = s.len();
        let level: u8 = 0;
        add_node(s, codes, nodes, max_depth, leaf_threshold, start, end, level);

        Self {nodes: nodes, codes: codes, bbox: bbox, }
    }
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