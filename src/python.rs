use pyo3::prelude::*;
use numpy::{PyReadonlyArray1, PyReadwriteArray1};

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

#[pymodule]
fn _thor<'py>(_py: Python, m: Bound<'py, PyModule>) -> PyResult<()> {
    // Circular filaments
    m.add_function(wrap_pyfunction!(count_rs, m.clone())?)?;
    m.add_function(wrap_pyfunction!(gemv, m.clone())?)?;
    m.add_function(wrap_pyfunction!(test, m.clone())?)?;
    Ok(())
}