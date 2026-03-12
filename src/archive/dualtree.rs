#[allow(unused)]
use crate::archive::biotsavart::bfield_leaf;
use crate::archive::octree::{
    BoundingBox, SourceNode, SourceOctree, Sources, get_range_in_same_node, size_at_level,
};
use crate::math::{cross, distance, mag, running_average, sort_by_indices, vec_distance};
use crate::{MU0_4PI, morton};

/// A collection of target points that can be sorted by morton codes
pub struct Targets {
    pub codes: Vec<u64>, // Morton code for each source point
    pub x: Vec<f64>,     // Geometric location of each target point
    pub y: Vec<f64>,
    pub z: Vec<f64>,
    pub bx: Vec<f64>, // Magnetic field at each location
    pub by: Vec<f64>,
    pub bz: Vec<f64>,
    bbox: Option<BoundingBox>,
    original_indices: Vec<usize>, // needed to properly return the magnetic field vectors
}

impl Targets {
    /// Create a new set of targets by using only their positions in 3D space
    pub fn from_target_points(targets: (&[f64], &[f64], &[f64]), max_depth: u8) -> Self {
        let n: usize = targets.0.len();
        assert_eq!(n, targets.1.len());
        assert_eq!(n, targets.2.len());

        let bbox: BoundingBox = BoundingBox::from_centroids(targets).unwrap();
        let mut codes: Vec<u64> = Vec::with_capacity(n);
        let x: Vec<f64> = targets.0.to_vec();
        let y: Vec<f64> = targets.1.to_vec();
        let z: Vec<f64> = targets.2.to_vec();
        let bx: Vec<f64> = vec![0.0; n];
        let by: Vec<f64> = vec![0.0; n];
        let bz: Vec<f64> = vec![0.0; n];
        let original_indices: Vec<usize> = vec![0; n];

        let scale: f64 = morton::calculate_scale_factor(max_depth as u32);
        let min_corner: (f64, f64, f64) = bbox.min_corner();
        let side_length: f64 = bbox.side_length;

        // Encode all targets
        for i in 0..n {
            let pt: (f64, f64, f64) = (x[i], y[i], z[i]);
            codes.push(morton::encode(pt, scale, side_length, min_corner));
        }

        Self {
            codes,
            x,
            y,
            z,
            bx,
            by,
            bz,
            bbox: Some(bbox),
            original_indices,
        }
    }

    /// Sort the octree by morton order
    pub fn morton_sort(&mut self) {
        // Get the indices to sort by from the codes
        let mut indices: Vec<usize> = (0..self.codes.len()).collect();
        indices.sort_by(|&i, &j| self.codes[i].cmp(&self.codes[j]));
        self.original_indices = indices.clone();

        // Now sort all of the other data arrays
        let mut scratch = vec![0.0; self.codes.len()];

        sort_by_indices(&mut self.x, &mut scratch, &indices);
        sort_by_indices(&mut self.y, &mut scratch, &indices);
        sort_by_indices(&mut self.z, &mut scratch, &indices);

        // Now sort the codes themselves
        let mut scratch_codes: Vec<u64> = vec![0; self.codes.len()];
        sort_by_indices(&mut self.codes, &mut scratch_codes, &indices);
    }

    /// Get the centroid of a particular Target
    pub fn centroid(&self, idx: u32) -> [f64; 3] {
        let i = idx as usize;
        [self.x[i], self.y[i], self.z[i]]
    }
}

/// Represents a node in the tree that represents a cluster of target points
pub enum TargetNode {
    Branch {
        level: u8,
        size: f64,
        children: [u32; 8],
        centroid: [f64; 3],
        b: [f64; 3],
        n_targets: u32,
    },
    Leaf {
        level: u8,
        target_range: (u32, u32),
        centroid: [f64; 3],
    },
}

/// Represents the target tree for dual-tree Barnes Hut
pub struct TargetOctree {
    pub targets: Targets,
    pub nodes: Vec<TargetNode>,
}

impl TargetOctree {
    /// Recursively add nodes to the octree
    ///
    /// Arguments:
    /// * `targets`:    the list of targets to add to the tree
    /// * `nodes`:      container for all nodes in the tree, modified by this function
    /// * `max_depth`:  maximum number of levels in the tree, 21 for problems with u64 morton codes
    /// * `start`:      the initial index in the range of morton codes to add
    /// * `end`:        the last index in the range of morton codes to add
    /// * `level`:      current level in the octree
    ///
    /// Returns:
    ///     the index of the node that it created
    fn add_node(
        targets: &Targets,
        nodes: &mut Vec<TargetNode>,
        max_depth: u8,
        leaf_threshold: u32,
        start: usize,
        end: usize,
        level: u8,
    ) -> u32 {
        let current_index: usize;

        // Create a leaf if the number of targets in the range is small enough (dead-end for recursion)
        if end - start <= leaf_threshold as usize {
            let mut cx = 0.0;
            let mut cy = 0.0;
            let mut cz = 0.0;
            for i in start..end {
                cx += targets.x[i];
                cy += targets.y[i];
                cz += targets.z[i];
            }
            cx /= (end - start) as f64;
            cy /= (end - start) as f64;
            cz /= (end - start) as f64;
            nodes.push(
                // Do not descend level
                TargetNode::Leaf {
                    level,
                    target_range: (start as u32, end as u32),
                    centroid: [cx, cy, cz],
                },
            );
            current_index = nodes.len() - 1;
        }
        // Make a branch node
        else {
            // Size of the node (cube side length) for BH acceptance criteria
            let size: f64 = size_at_level(targets.bbox.as_ref().unwrap().side_length, level);

            // Initialize the branch node first, then recursive calls later fill it
            nodes.push(TargetNode::Branch {
                level,
                size,
                children: [0; 8],
                centroid: [0.0; 3],
                n_targets: 0,
                b: [0.0; 3],
            });

            current_index = nodes.len() - 1; // index of the node just created

            // Now recurse to create its children
            let mut child_indices: [u32; 8] = [0; 8];
            let mut n_children: usize = 0;
            let mut cursor: usize = start;

            while cursor < end {
                // We descend to next level here, hence `level+1`

                // index to targets array
                // this is the end index of all child nodes that share the same morton code at this level
                let child_end =
                    get_range_in_same_node(&targets.codes, level + 1, max_depth, cursor);

                // index in nodes array
                let child_idx = Self::add_node(
                    targets,
                    nodes,
                    max_depth,
                    leaf_threshold,
                    cursor,
                    child_end,
                    level + 1,
                );
                child_indices[n_children] = child_idx;
                n_children += 1;
                cursor = child_end; // index in sources array
            }

            // Update the previously-initialized node with indices to its children
            // Take the simple average of the locations of each of the target points in this cluster
            let mut parent_centroid: [f64; 3] = [0.0; 3];
            let mut n_targets_parent: u32 = 0;

            for idx in child_indices {
                if idx < 1 {
                    // Skip empty nodes
                    break;
                }

                match nodes[idx as usize] {
                    TargetNode::Branch {
                        level: _,
                        size: _,
                        children: _,
                        centroid,
                        b: _,
                        n_targets,
                    } => {
                        n_targets_parent = running_average(
                            n_targets_parent,
                            &mut parent_centroid,
                            n_targets,
                            &centroid,
                        );
                    }
                    TargetNode::Leaf {
                        level: _,
                        target_range,
                        ..
                    } => {
                        for i in (target_range.0 as usize)..(target_range.1 as usize) {
                            let leaf_centroid = [targets.x[i], targets.y[i], targets.z[i]];
                            n_targets_parent = running_average(
                                n_targets_parent,
                                &mut parent_centroid,
                                1,
                                &leaf_centroid,
                            );
                        }
                    }
                }
            }

            if let TargetNode::Branch {
                level: _,
                size: _,
                children,
                centroid,
                b: _,
                n_targets,
            } = &mut nodes[current_index]
            {
                *children = child_indices;
                *centroid = parent_centroid;
                *n_targets = n_targets_parent;
            }
        }

        current_index as u32
    }

    // Build the octree by making recursive calls to add_node()
    fn build_tree(
        targets: &Targets,
        nodes: &mut Vec<TargetNode>,
        max_depth: u8,
        leaf_threshold: u32,
    ) {
        let start: usize = 0;
        let end: usize = targets.codes.len();
        let level: u8 = 0;
        let _ = Self::add_node(targets, nodes, max_depth, leaf_threshold, start, end, level);
    }

    // Create a target octree from target points
    pub fn from_target_points(
        targets: (&[f64], &[f64], &[f64]),
        max_depth: u8,
        leaf_threshold: u32,
    ) -> Self {
        let n_targets: usize = targets.0.len();
        let mut targets: Targets = Targets::from_target_points(targets, max_depth);
        targets.morton_sort();

        let mut nodes: Vec<TargetNode> = Vec::with_capacity(n_targets);
        TargetOctree::build_tree(&targets, &mut nodes, max_depth, leaf_threshold);

        Self { targets, nodes }
    }

    /// Update all targets with the contributions from their parent nodes, starting at the root  
    /// This is the "downwards pass"
    pub fn update_targets(&mut self) {
        self.propagate(0, [0.0; 3]);
    }

    // Recursively propagate all contributions to the targets in the leaves
    fn propagate(&mut self, node_index: u32, b_contribution: [f64; 3]) {
        match self.nodes[node_index as usize] {
            TargetNode::Branch {
                level: _,
                size: _,
                children,
                centroid: _,
                b,
                n_targets: _,
            } => {
                for child in children {
                    if child < 1 {
                        break;
                    } else {
                        let b_total = [
                            b_contribution[0] + b[0],
                            b_contribution[1] + b[1],
                            b_contribution[2] + b[2],
                        ];
                        self.propagate(child, b_total);
                    }
                }
            }
            TargetNode::Leaf {
                level: _,
                target_range,
                centroid: _,
            } => {
                for i in (target_range.0 as usize)..(target_range.1 as usize) {
                    self.targets.bx[i] += b_contribution[0];
                    self.targets.by[i] += b_contribution[1];
                    self.targets.bz[i] += b_contribution[2];
                }
            }
        }
    }
}

/// Compute the magnetic field generated by a current-carrying
/// finite element mesh using a dual (source/target) octree Barnes Hut approach
pub fn bfield_dualtree(
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
    theta_source: f64,
    theta_target: f64,
    leaf_threshold: u32,
) {
    let max_depth: u8 = 21;

    let source_tree: SourceOctree = SourceOctree::from_source_points(
        (centx, centy, centz),
        vol,
        (jx, jy, jz),
        max_depth,
        leaf_threshold,
    )
    .unwrap();

    let mut target_tree: TargetOctree =
        TargetOctree::from_target_points((x, y, z), max_depth, leaf_threshold);

    bfield_trees(
        &source_tree,
        &mut target_tree,
        0,
        0,
        theta_source,
        theta_target,
    );
    target_tree.update_targets();

    // Return the data in the order it was originally supplied to the function
    // Perhaps there's a more efficient way to do this?
    for sorted_i in 0..bx.len() {
        let original_i = target_tree.targets.original_indices[sorted_i];
        bx[original_i] = target_tree.targets.bx[sorted_i];
        by[original_i] = target_tree.targets.by[sorted_i];
        bz[original_i] = target_tree.targets.bz[sorted_i];
    }
}

// Recursively compute the interaction between source and target trees
fn bfield_trees(
    source_tree: &SourceOctree,
    target_tree: &mut TargetOctree,
    source_node_idx: u32,
    target_node_idx: u32,
    theta_source: f64,
    theta_target: f64,
) {
    let source_node = &source_tree.nodes[source_node_idx as usize];
    let target_node = &mut target_tree.nodes[target_node_idx as usize];

    match (source_node, target_node) {
        // Leaf to leaf interaction: direct B-S solve
        (
            &SourceNode::Leaf {
                level: _,
                source_range,
                centroid: _,
                ..
            },
            &mut TargetNode::Leaf {
                level: _,
                target_range,
                centroid: _,
            },
        ) => {
            bfield_leaf_leaf(
                source_range,
                target_range,
                &source_tree.sources,
                &mut target_tree.targets,
            );
        }

        // Source branch to target leaf
        (
            &SourceNode::Branch {
                level: _,
                size: src_size,
                children: src_children,
                centroid: src_centroid,
                vj,
            },
            &mut TargetNode::Leaf {
                level: tgt_level,
                target_range,
                centroid: tgt_centroid,
            },
        ) => {
            let d = distance(&src_centroid, &tgt_centroid);
            let side_length = target_tree.targets.bbox.as_ref().unwrap().side_length;
            let tgt_size: f64 = size_at_level(side_length, tgt_level);
            if src_size / d < theta_source && tgt_size / d < theta_target {
                bfield_branch_leaf(&src_centroid, &vj, target_range, &mut target_tree.targets);
            } else {
                for i in src_children {
                    if i < 1 {
                        break;
                    }
                    bfield_trees(
                        source_tree,
                        target_tree,
                        i,
                        target_node_idx,
                        theta_source,
                        theta_target,
                    );
                }
            }
        }

        // Source leaf to target branch
        (
            &SourceNode::Leaf {
                level: src_level,
                source_range,
                centroid: src_centroid,
            },
            &mut TargetNode::Branch {
                level: tgt_level,
                size: _,
                children: tgt_children,
                centroid: tgt_centroid,
                ref mut b,
                n_targets: _,
            },
        ) => {
            let d = distance(&src_centroid, &tgt_centroid);
            let src_size = size_at_level(
                source_tree.sources.bbox.as_ref().unwrap().side_length,
                src_level,
            );
            let side_length = target_tree.targets.bbox.as_ref().unwrap().side_length;
            let tgt_size: f64 = size_at_level(side_length, tgt_level);
            if src_size / d < theta_source && tgt_size / d < theta_target {
                //src_size.max(tgt_size) / d < theta
                let b_contribution =
                    bfield_leaf_branch(source_range, &source_tree.sources, &tgt_centroid);
                b[0] += b_contribution[0];
                b[1] += b_contribution[1];
                b[2] += b_contribution[2];
            } else {
                for i in tgt_children {
                    if i < 1 {
                        break;
                    }
                    bfield_trees(
                        source_tree,
                        target_tree,
                        source_node_idx,
                        i,
                        theta_source,
                        theta_target,
                    );
                }
            }
        }

        // Source branch to target branch
        (
            &SourceNode::Branch {
                level: _,
                size: src_size,
                children: src_children,
                centroid: src_centroid,
                vj,
            },
            &mut TargetNode::Branch {
                level: _,
                size: tgt_size,
                children: tgt_children,
                centroid: tgt_centroid,
                ref mut b,
                n_targets: _,
            },
        ) => {
            let d = distance(&src_centroid, &tgt_centroid);
            if src_size / d < theta_source && tgt_size / d < theta_target {
                let b_contribution = bfield_point_point(&src_centroid, &vj, &tgt_centroid);
                b[0] += b_contribution[0];
                b[1] += b_contribution[1];
                b[2] += b_contribution[2];
            } else if src_size > tgt_size {
                for i in src_children {
                    if i < 1 {
                        break;
                    }
                    bfield_trees(
                        source_tree,
                        target_tree,
                        i,
                        target_node_idx,
                        theta_source,
                        theta_target,
                    );
                }
            } else {
                for i in tgt_children {
                    if i < 1 {
                        break;
                    }
                    bfield_trees(
                        source_tree,
                        target_tree,
                        source_node_idx,
                        i,
                        theta_source,
                        theta_target,
                    );
                }
            }
        }
    }
}

// Compute the interaction between a source leaf and a target leaf
fn bfield_leaf_leaf(
    source_range: (u32, u32),
    target_range: (u32, u32),
    sources: &Sources,
    targets: &mut Targets,
) {
    let (is0, is1) = (source_range.0 as usize, source_range.1 as usize);
    let centroids = (
        &sources.xg[is0..is1],
        &sources.yg[is0..is1],
        &sources.zg[is0..is1],
    );
    let vj = (
        &sources.vjx[is0..is1],
        &sources.vjy[is0..is1],
        &sources.vjz[is0..is1],
    );

    for it in (target_range.0 as usize)..(target_range.1 as usize) {
        let target = [targets.x[it], targets.y[it], targets.z[it]];
        let b_contribution = bfield_leaf(centroids, vj, &target);
        targets.bx[it] += b_contribution[0];
        targets.by[it] += b_contribution[1];
        targets.bz[it] += b_contribution[2];
    }
}

// Compute the interaction between a single point to a single point
fn bfield_point_point(centroid: &[f64; 3], vj: &[f64; 3], pt: &[f64; 3]) -> [f64; 3] {
    let r = vec_distance(centroid, pt);
    let rmag = mag(&r);
    let mut jxr = [0.0; 3];
    if rmag > 1e-4 {
        cross(vj, &r, &mut jxr);
        jxr[0] *= MU0_4PI / rmag.powi(3);
        jxr[1] *= MU0_4PI / rmag.powi(3);
        jxr[2] *= MU0_4PI / rmag.powi(3);
    }
    jxr
}

// Compute the interaction of a source branch to a target leaf
fn bfield_branch_leaf(
    centroid: &[f64; 3],
    vj: &[f64; 3],
    target_range: (u32, u32),
    targets: &mut Targets,
) {
    for i in (target_range.0 as usize)..(target_range.1 as usize) {
        let pt = targets.centroid(i as u32);
        let b = bfield_point_point(centroid, vj, &pt);
        targets.bx[i] += b[0];
        targets.by[i] += b[1];
        targets.bz[i] += b[2];
    }
}

fn bfield_leaf_branch(
    source_range: (u32, u32),
    sources: &Sources,
    tgt_centroid: &[f64; 3],
) -> [f64; 3] {
    let (is0, is1) = (source_range.0 as usize, source_range.1 as usize);
    let centroids = (
        &sources.xg[is0..is1],
        &sources.yg[is0..is1],
        &sources.zg[is0..is1],
    );
    let vj = (
        &sources.vjx[is0..is1],
        &sources.vjy[is0..is1],
        &sources.vjz[is0..is1],
    );

    bfield_leaf(centroids, vj, tgt_centroid)
}
