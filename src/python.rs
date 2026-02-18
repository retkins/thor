#![allow(unused)]

use pyo3::prelude::*;
use numpy::{PyReadonlyArray1, PyReadwriteArray1};

use crate::biotsavart;
use crate::dualtree;
use crate::sources::bfield_hexahedron;

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

    // #[cfg(feature = "parallel")]
    // if nthreads_requested != 1 {
    //     use crate::biotsavart_parallel;
    //     biotsavart_parallel::bfield_octree_parallel(
    //         centx.as_slice()?, 
    //         centy.as_slice()?, 
    //         centz.as_slice()?, 
    //         vol.as_slice()?, 
    //         jx.as_slice()?, 
    //         jy.as_slice()?, 
    //         jz.as_slice()?, 
    //         x.as_slice()?, 
    //         y.as_slice()?, 
    //         z.as_slice()?, 
    //         bx.as_slice_mut()?, 
    //         by.as_slice_mut()?, 
    //         bz.as_slice_mut()?, 
    //         theta, leaf_threshold, 
    //         nthreads_requested
    //     );
    //     return Ok(());
    // }

    use crate::octree_generic::{Octree, point, HFieldSolver, CurrentSources};
    let mut sources: CurrentSources<point::PointSources>= CurrentSources(point::PointSources::new(
        centx.as_slice()?, 
        centy.as_slice()?, 
        centz.as_slice()?, 
        vol.as_slice()?, 
        jx.as_slice()?, 
        jy.as_slice()?, 
        jz.as_slice()?, 
    ));
    let tree: Octree<CurrentSources<point::PointSources>>= Octree::build_from_sources(sources);

    #[cfg(feature = "parallel")]
    tree.h_field_parallel(
        (x.as_slice()?, 
        y.as_slice()?, 
        z.as_slice()?), 
        (bx.as_slice_mut()?, 
        by.as_slice_mut()?, 
        bz.as_slice_mut()?), 
        theta,
        nthreads_requested
    );

    // original version:
    //
    // biotsavart::bfield_octree(
    //     centx.as_slice()?, 
    //     centy.as_slice()?, 
    //     centz.as_slice()?, 
    //     vol.as_slice()?, 
    //     jx.as_slice()?, 
    //     jy.as_slice()?, 
    //     jz.as_slice()?, 
    //     x.as_slice()?, 
    //     y.as_slice()?, 
    //     z.as_slice()?, 
    //     bx.as_slice_mut()?, 
    //     by.as_slice_mut()?, 
    //     bz.as_slice_mut()?, 
    //     theta, leaf_threshold
    // );
    Ok(())
}

#[pyfunction]
fn _bfield_dualtree(
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
    theta_source: f64, 
    theta_target: f64,
    leaf_threshold: u32, 
    nthreads_requested: u32
) -> PyResult<()> {
    dualtree::bfield_dualtree(
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
        theta_source, theta_target, 
        leaf_threshold
    );
    Ok(())
}

#[pyfunction]
fn _bfield_hexahedron(
    nx: PyReadonlyArray1<f64>, 
    ny: PyReadonlyArray1<f64>, 
    nz: PyReadonlyArray1<f64>, 
    jdensity: PyReadonlyArray1<f64>,
    target: PyReadonlyArray1<f64>
)-> PyResult<(f64, f64, f64)> {

    let b = bfield_hexahedron(
    nx.as_slice()?, 
    ny.as_slice()?,
    nz.as_slice()?,
    jdensity.as_slice()?,
    target.as_slice()?);

    Ok((b[0], b[1], b[2]))
}

#[pyfunction]
fn _hfield_tetrahedrons(
    nodes_flat: PyReadonlyArray1<f64>, 
    centroids_flat: PyReadonlyArray1<f64>, 
    vol: PyReadonlyArray1<f64>, 
    jdensity_flat: PyReadonlyArray1<f64>, 
    x: PyReadonlyArray1<f64>, 
    y: PyReadonlyArray1<f64>, 
    z: PyReadonlyArray1<f64>, 
    mut bx: PyReadwriteArray1<f64>,
    mut by: PyReadwriteArray1<f64>,
    mut bz: PyReadwriteArray1<f64>, 
    theta: f64, 
    nthreads_requested: u32
) -> PyResult<()> {

    use crate::octree_generic::{Octree, tet_element, HFieldSolver, CurrentSources};
    let mut sources: CurrentSources<tet_element::TetSources> = CurrentSources(tet_element::TetSources::new(
        nodes_flat.as_slice()?, 
        centroids_flat.as_slice()?, 
        vol.as_slice()?, 
        jdensity_flat.as_slice()?, 
    ));
    let tree: Octree<CurrentSources<tet_element::TetSources>>= Octree::build_from_sources(sources);

    tree.h_field_parallel(
        (x.as_slice()?, 
        y.as_slice()?, 
        z.as_slice()?), 
        (bx.as_slice_mut()?, 
        by.as_slice_mut()?, 
        bz.as_slice_mut()?), 
        theta,
        nthreads_requested
    );
    Ok(())
}

#[pyfunction]
pub fn _hfield_tetrahedrons_direct(
    nodes_flat: PyReadonlyArray1<f64>, 
    centroids_flat: PyReadonlyArray1<f64>, 
    vol: PyReadonlyArray1<f64>, 
    jdensity_flat: PyReadonlyArray1<f64>, 
    x: PyReadonlyArray1<f64>, 
    y: PyReadonlyArray1<f64>, 
    z: PyReadonlyArray1<f64>, 
    mut hx: PyReadwriteArray1<f64>,
    mut hy: PyReadwriteArray1<f64>,
    mut hz: PyReadwriteArray1<f64>, 
    nthreads_requested: u32
) -> PyResult<()> {

    use crate::biotsavart_parallel;
    biotsavart_parallel::hfield_direct_tet_parallel(
        nodes_flat.as_slice()?, 
        vol.as_slice()?, 
        jdensity_flat.as_slice()?,
        x.as_slice()?,
        y.as_slice()?,
        z.as_slice()?,
        hx.as_slice_mut()?,
        hy.as_slice_mut()?,
        hz.as_slice_mut()?,
        nthreads_requested
    );

    Ok(())
}


#[pyfunction]
fn _hfield_dipole(
    centx: PyReadonlyArray1<f64>, 
    centy: PyReadonlyArray1<f64>,
    centz: PyReadonlyArray1<f64>,
    vol: PyReadonlyArray1<f64>,
    mx: PyReadonlyArray1<f64>,
    my: PyReadonlyArray1<f64>,
    mz: PyReadonlyArray1<f64>,
    x: PyReadonlyArray1<f64>,
    y: PyReadonlyArray1<f64>,
    z: PyReadonlyArray1<f64>,
    mut hx: PyReadwriteArray1<f64>, 
    mut hy: PyReadwriteArray1<f64>,
    mut hz: PyReadwriteArray1<f64>, 
    theta: f64,
    leaf_threshold: u32, 
    nthreads_requested: u32
) ->PyResult<()> {
    use crate::octree_generic::{Octree, DipoleSources, point::PointSources};
    let sources = DipoleSources(PointSources::new_dipole(
        centx.as_slice()?, 
        centy.as_slice()?, 
        centz.as_slice()?, 
        vol.as_slice()?, 
        mx.as_slice()?, 
        my.as_slice()?, 
        mz.as_slice()?
    ));
let tree = Octree::build_from_sources(sources);

// Evaluate
tree.h_field((x.as_slice()?, 
        y.as_slice()?, 
        z.as_slice()?), (hx.as_slice_mut()?, 
        hy.as_slice_mut()?, 
        hz.as_slice_mut()?), theta);

Ok(())
}


#[pymodule]
fn _thor<'py>(_py: Python, m: Bound<'py, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_bfield_direct, m.clone())?)?;
    m.add_function(wrap_pyfunction!(_bfield_octree, m.clone())?)?;
    m.add_function(wrap_pyfunction!(_bfield_dualtree, m.clone())?)?;
    m.add_function(wrap_pyfunction!(_bfield_hexahedron, m.clone())?)?;
    m.add_function(wrap_pyfunction!(_hfield_tetrahedrons, m.clone())?)?;
    m.add_function(wrap_pyfunction!(_hfield_dipole, m.clone())?)?;
    m.add_function(wrap_pyfunction!(_hfield_tetrahedrons_direct, m.clone())?)?;
    
    Ok(())
}