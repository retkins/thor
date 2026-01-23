use std::f64::consts::PI;

/// Biot-Savart integration constant:  
/// $$\frac{\mu_0}{4\pi} [H/m]$$
pub const MU0_4PI: f64 = 1e-7;

/// Magnetic permeability of free space:  
/// $$\mu_0 = 4\pi \cdot 10^{-7} H/m$$
pub const MU0: f64 = 4.0*PI*MU0_4PI;

pub mod vector;

/// Direct Biot-Savart Law integration functions
pub mod direct;

/// Analytical expressions for the magnetic field under highly specific conditions
pub mod analytical;

/// Basic read/write file operations
pub mod io; 

/// Error types for `thor`
pub mod errors;

/// Expressions for multipole Biot-Savart sources
pub mod sources;

/// Low-level math expressions
pub mod math;

/// Compute Morton codes and related functions
/// 
/// These are generally taken from previous work on the C implementation
/// 
/// Notes: 
/// 1. Grids, etc. are based on a *cube* in 3d space
pub mod morton;


/// Octree build and traversal
pub mod octree;

/// Python bindings
#[cfg(feature="python")]
pub mod python;