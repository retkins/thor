/// Dual tree Barnes-Hut Octree 
/// 
/// Used to calculate the effect of M source points on N target points
/// using the Biot-Savart law for magnetic fields
/// 
/// Each source point has associated volume and current density.
/// The direct summation algorithm is: 
/// delta_B = mu0/4pi * volume * J x r' / |r'|^3


// Vector or point in 3D space 
// test both f32 and f64 versions - might be able to get away with lower precision
pub struct Vector<T> {
    x: T, 
    y: T, 
    z: T
}

// Source Point: individual component of the finite sum
pub struct SourcePoint<T> {
    location: Vector<T>,
    moment: Vector<T>
}

pub struct SourceNode<T> {
    geometric_centroid: Vector<T>,      // spatial position
    size: T,                            // width of each side (cube)
    centroid: Vector<T>,                // scaled centroid based on moments of subnodes
    moment: Vector<T>,                  // volume * current density
    children: [u32;8],                  // indices of child nodes; 0 means "no child in this octant"
    leaf: Option<SourcePoint<T>>
}

pub struct SourceTree<T>  {
    nodes: Vec<SourceNode<T>>, 
}


pub struct TargetPoint<T> {
    location: Vector<T>, 
    bfield: Vector<T>
}

pub struct TargetNode<T> {
    geometric_centroid: Vector<T>,      // spatial position
    size: T,                            // width of each side (cube) (spanning volume of node)
    centroid: Vector<T>,                // scaled centroid based on simple average of child node positions
    bfield: Vector<T>,                  // magnetic field associated with the node (applied to all childen in downward pass)
    children: [u32;8],                  // indices of child nodes; 0 means "no child in this octant"
    leaf: Option<TargetPoint<T>>
}


// pub struct SourceOctree {
//     xg: Vec<f64>, 
//     yg: Vec<f64>, 
//     zg: Vec<f64>, 
//     size: Vec<f64>, 
//     xc: Vec<f64>, 
//     yc: Vec<f64>, 
//     zc: Vec<f64>, 
//     vjx: Vec<f64>, 
//     vjy: Vec<f64>, 
//     vjz: Vec<f64> 
// }