
/// In place cross-product  
/// $$ \vec{a} \times \vec{b} = \vec{c} $$
#[inline(always)]
pub fn cross(
    a: (f64, f64, f64), 
    b: (f64, f64, f64), 
    c: &mut (f64, f64, f64)
) {
    c.0 = a.1*b.2 - a.2*b.1;
    c.1 = a.2*b.0 - a.0*b.2; 
    c.2 = a.0*b.1 - a.1*b.0;
}