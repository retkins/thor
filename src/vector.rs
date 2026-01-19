
use std::ops::{Mul, Sub};

/// Vector or point in 3D space
pub struct Vector<T> 
    where T: Copy + Mul<Output=T> + Sub<Output = T> {
    pub x: T, 
    pub y: T, 
    pub z: T
}

impl<T:> Vector<T> 
    where T: Copy + Mul<Output=T> + Sub<Output = T> {
    pub fn new(x: T, y: T, z: T) ->Self {
        Self {
            x: x, 
            y: y, 
            z: z
        }
    }

    /// Cross product of two vectors, returning a new vector
    ///     c = a x b
    pub fn cross(&self, other: &Vector<T>) -> Self {
        let cx = self.y*other.z - self.z*other.y; 
        let cy = self.z*other.x - self.x*other.z; 
        let cz = self.x*other.y - self.y*other.x; 

        Vector::new(cx, cy, cz)
    }
}

/// Cross product of two vectors, overwriting a third vector in place 
pub fn cross<T>(a: &Vector<T>, b: &Vector<T>, c: &mut Vector<T>)
    where T: Copy + Mul<Output=T> + Sub<Output = T>  {
        c.x = a.y*b.z - a.z*b.y; 
        c.y = a.z*b.x - a.x*b.z; 
        c.z = a.x*b.y - a.y*b.x;
}