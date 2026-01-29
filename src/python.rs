#![allow(unused)]

use pyo3::prelude::*;
use numpy::{PyReadonlyArray1, PyReadwriteArray1};

use crate::biotsavart;

#[pyfunction]
fn test(x: f64) -> PyResult<f64> {
    Ok(2.0*x)
}

#[pyfunction]
fn gemv(
    alpha: f64, 
    a: PyReadonlyArray1<f64>, 
    b: PyReadonlyArray1<f64>, 
    beta: f64, 
    mut y: PyReadonlyArray1<f64>) -> PyResult<()> {
        let a = a.as_slice()?;
        let b = b.as_slice()?;
        let y = y.as_slice()?;
        // todo...
        Ok(())
}

#[pyfunction]
fn count_rs(
    a: PyReadonlyArray1<f64>, 
    mut b: PyReadwriteArray1<f64>
) -> PyResult<i64>{
    let a = a.as_slice()?; 
    let b = b.as_slice_mut()?;

    assert_eq!(a.len(), b.len());
    for (i, ai) in a.iter().enumerate() { 
        b[i] = 2.0 * ai; 
    }
    
    Ok(a.len() as i64)
}

#[pyfunction]
fn _bfield_direct(
    centx: PyReadonlyArray1<f64>, 
    centy: PyReadonlyArray1<f64>,
    centz: PyReadonlyArray1<f64>,
    vol: PyReadonlyArray1<f64>,
    jx: PyReadonlyArray1<f64>,
    jy: PyReadonlyArray1<f64>,
    jz: PyReadonlyArray1<f64>,
    x: PyReadonlyArray1<f64>,
    y: PyReadonlyArray1<f64>,
    z: PyReadonlyArray1<f64>,
    mut bx: PyReadwriteArray1<f64>, 
    mut by: PyReadwriteArray1<f64>,
    mut bz: PyReadwriteArray1<f64>
) -> PyResult<()> {

    biotsavart::bfield_direct(
        centx.as_slice()?, 
        centy.as_slice()?, 
        centz.as_slice()?, 
        vol.as_slice()?, 
        jx.as_slice()?, 
        jy.as_slice()?, 
        jz.as_slice()?, 
        x.as_slice()?, 
        y.as_slice()?, 
        z.as_slice()?, 
        bx.as_slice_mut()?, 
        by.as_slice_mut()?, 
        bz.as_slice_mut()?
    );
    Ok(())
}


#[pyfunction]
fn _bfield_octree(
    centx: PyReadonlyArray1<f64>, 
    centy: PyReadonlyArray1<f64>,
    centz: PyReadonlyArray1<f64>,
    vol: PyReadonlyArray1<f64>,
    jx: PyReadonlyArray1<f64>,
    jy: PyReadonlyArray1<f64>,
    jz: PyReadonlyArray1<f64>,
    x: PyReadonlyArray1<f64>,
    y: PyReadonlyArray1<f64>,
    z: PyReadonlyArray1<f64>,
    mut bx: PyReadwriteArray1<f64>, 
    mut by: PyReadwriteArray1<f64>,
    mut bz: PyReadwriteArray1<f64>, 
    theta: f64, 
    leaf_threshold: u32
) -> PyResult<()> {

    biotsavart::bfield_octree(
        centx.as_slice()?, 
        centy.as_slice()?, 
        centz.as_slice()?, 
        vol.as_slice()?, 
        jx.as_slice()?, 
        jy.as_slice()?, 
        jz.as_slice()?, 
        x.as_slice()?, 
        y.as_slice()?, 
        z.as_slice()?, 
        bx.as_slice_mut()?, 
        by.as_slice_mut()?, 
        bz.as_slice_mut()?, 
        theta, leaf_threshold
    );
    Ok(())
}

#[pymodule]
fn _thor<'py>(_py: Python, m: Bound<'py, PyModule>) -> PyResult<()> {
    // Circular filaments
    m.add_function(wrap_pyfunction!(count_rs, m.clone())?)?;
    m.add_function(wrap_pyfunction!(gemv, m.clone())?)?;
    m.add_function(wrap_pyfunction!(test, m.clone())?)?;
    m.add_function(wrap_pyfunction!(_bfield_direct, m.clone())?)?;
    m.add_function(wrap_pyfunction!(_bfield_octree, m.clone())?)?;
    Ok(())
}