use std::f64::consts::PI;

const MU0_4PI: f64 = 1e-7;
const MU0: f64 = 4.0*PI*MU0_4PI;

pub mod direct;
pub mod analytical;