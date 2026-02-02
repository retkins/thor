//! Lightning-fast magnetic field calculations using octrees and the Barnes-Hut algorithm  
//! 
//! This is the Rust API documentation. 
//! Main documentation, including a theory manual and the Python API, 
//! is hosted [here](https://freestatelabs.com/thor).


use std::f64::consts::PI;

/// Biot-Savart integration constant:  
/// $$\frac{\mu_0}{4\pi} [H/m]$$
pub const MU0_4PI: f64 = 1e-7;

/// Magnetic permeability of free space:  
/// $$\mu_0 = 4\pi \cdot 10^{-7} H/m$$
pub const MU0: f64 = 4.0*PI*MU0_4PI;

/// Expressions for multipole Biot-Savart sources
pub mod sources;

/// Magnetic field calculations
pub mod biotsavart;

/// Analytical expressions for the magnetic field under highly specific conditions
pub mod analytical;

/// Basic read/write file operations
pub mod io; 

/// Error types for `thor`
pub mod errors;

/// Low-level math and array operations
pub mod math;

/// Compute Morton codes and related functions
/// 
/// These are generally taken from previous work on the C implementation
/// 
/// Notes: 
/// 1. Grids, etc. are based on a *cube* in 3d space
pub mod morton;


/// Dual tree Barnes-Hut Octree
///
/// Used to calculate the effect of M source points on N target points
/// using the Biot-Savart law for magnetic fields
///
/// Each source point has associated volume and current density.
/// The direct summation algorithm is:
/// delta_B = mu0/4pi * volume * J x r' / |r'|^3
pub mod octree;

/// Python bindings
#[cfg(feature="python")]
pub mod python;

/// Parallel processing
#[cfg(feature="parallel")]
pub mod biotsavart_parallel;