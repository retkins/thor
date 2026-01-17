#![allow(non_snake_case)]

use crate::MU0;

pub fn bfield_loop_axis(z: f64, I: f64, R: f64) -> f64 {
    let R2: f64 = R*R;
    return MU0 * I * R2 / (2.0 * (z*z + R2).powf(1.5));
}