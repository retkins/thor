#![allow(unused)]

use pyo3::prelude::*;
use numpy::{PyReadonlyArray1, PyReadwriteArray1};

use crate::biotsavart;

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
    mut bz: PyReadwriteArray1<f64>, 
    nthreads_requested: u32
) -> PyResult<()> {

    #[cfg(feature = "parallel")]
    if nthreads_requested != 1 {
        use crate::biotsavart_parallel;
        biotsavart_parallel::bfield_direct_parallel(
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
            nthreads_requested
        );
        return Ok(());
    }

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
    leaf_threshold: u32, 
    nthreads_requested: u32
) -> PyResult<()> {

    #[cfg(feature = "parallel")]
    if nthreads_requested != 1 {
        use crate::biotsavart_parallel;
        biotsavart_parallel::bfield_octree_parallel(
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
            theta, leaf_threshold, 
            nthreads_requested
        );
        return Ok(());
    }

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
    m.add_function(wrap_pyfunction!(_bfield_direct, m.clone())?)?;
    m.add_function(wrap_pyfunction!(_bfield_octree, m.clone())?)?;
    Ok(())
}