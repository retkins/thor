use std::f64::consts::PI;

/// Biot-Savart integration constant:  
/// $$\frac{\mu_0}{4\pi} [H/m]$$
pub const MU0_4PI: f64 = 1e-7;

/// Magnetic permeability of free space:  
/// $$\mu_0 = 4\pi \cdot 10^{-7} H/m$$
pub const MU0: f64 = 4.0*PI*MU0_4PI;

pub mod vector;
pub mod direct;
pub mod analytical;
pub mod io; 
pub mod errors;
pub mod sources;
pub mod math;

/// Compute Morton codes and related functions
/// 
/// These are generally taken from previous work on the C implementation
/// 
/// Notes: 
/// 1. Grids, etc. are based on a *cube* in 3d space
pub mod morton;

/// Python bindings
#[cfg(feature="python")]
pub mod python;

