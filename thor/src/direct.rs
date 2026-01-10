#![allow(non_snake_case)]

use crate::MU0_4PI; 

pub fn bfield_direct(
    x: &[f64], y: &[f64], z: &[f64], 
    centx: &[f64], centy: &[f64], centz: &[f64],
    vol: &[f64], 
    jx: &[f64], jy: &[f64], jz: &[f64], 
    Bx: &mut [f64], By: &mut [f64], Bz: &mut [f64], 
) {

    let m: usize = centx.len();
    let n: usize = x.len();

    for i in 0..m {
        
        let centxi: f64 = centx[i]; 
        let centyi: f64 = centy[i];
        let centzi: f64 = centz[i];
        let vol_mu0_4pi: f64 = vol[i]*MU0_4PI;
        let jxi = jx[i]; 
        let jyi = jy[i]; 
        let jzi = jz[i];

        for j in 0..n {

            let rx: f64 = x[j] - centxi; 
            let ry: f64 = y[j] - centyi; 
            let rz: f64 = z[j] - centzi; 
            let r: f64 = (rx*rx + ry*ry + rz*rz).sqrt(); 
            let rinv: f64 = 1.0 / r;
            let rinv3: f64 = rinv*rinv*rinv;

            let jxrpx: f64 = jyi*rz - jzi*ry; 
            let jxrpy: f64 = jzi*rx - jxi*rz; 
            let jxrpz: f64 = jxi*ry - jyi*rx; 

            let constant: f64 = vol_mu0_4pi*rinv3;
            Bx[j] += constant * jxrpx; 
            By[j] += constant * jxrpy; 
            Bz[j] += constant * jxrpz;
        }
    }
}

#[cfg(test)]
mod tests{

    #[test]
    fn test_direct() {
        

    }
}
