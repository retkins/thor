use rayon::prelude::*;
use std::thread::available_parallelism;
use std::num::NonZeroUsize;

use crate::{biotsavart::{bfield_direct, bfield_node}, octree::SourceOctree};

fn get_nthreads(nthreads_requested: u32) -> usize {
    let nthreads: usize;
    let nthreads_available: usize = available_parallelism().unwrap_or(NonZeroUsize::MIN).get();
    if nthreads_requested as usize > nthreads_available || nthreads_requested == 0 {
        nthreads = nthreads_available;
    }
    else {
        nthreads = nthreads_requested as usize;
    }
    return nthreads;
}


pub fn bfield_direct_parallel(
    centx: &[f64], centy: &[f64], centz: &[f64],
    vol: &[f64], 
    jx: &[f64], jy: &[f64], jz: &[f64], 
    x: &[f64], y: &[f64], z: &[f64], 
    bx: &mut [f64], by: &mut [f64], bz: &mut [f64], 
    nthreads_requested: u32
) -> Result<(),()>{
    // TODO: length checks
    let n: usize = centx.len();
    let nthreads: usize = get_nthreads(nthreads_requested);
    let chunk_size: usize = (n / nthreads).max(1);

    // chunk the inputs 
    let _x = x.par_chunks(chunk_size);
    let _y = y.par_chunks(chunk_size);
    let _z = z.par_chunks(chunk_size);
    let _bx = bx.par_chunks_mut(chunk_size);
    let _by = by.par_chunks_mut(chunk_size);
    let _bz = bz.par_chunks_mut(chunk_size);

    (_x, _y, _z, _bx, _by, _bz)
        .into_par_iter()
        .try_for_each(|(_x, _y, _z, _bx, _by, _bz)| {bfield_direct(centx, centy, centz, vol, jx, jy, jz, _x, _y, _z, _bx, _by, _bz)})?;

    Ok(())
}


fn bfield_node_chunk(
    tree: &SourceOctree,
    x: &[f64], 
    y: &[f64], 
    z:  &[f64], 
    bx: &mut [f64], 
    by: &mut [f64], 
    bz: &mut [f64], 
    theta: f64
) -> Result<(), ()> {
    let n: usize = x.len();
    // TODO: length checks
    for i in 0..n {
        let centroid = [x[i], y[i], z[i]];
        let b = bfield_node(&tree, 0, &centroid, theta);
        bx[i] += b[0]; 
        by[i] += b[1]; 
        bz[i] += b[2];
    }
    Ok(())
}


pub fn bfield_octree_parallel(
    centx: &[f64], centy: &[f64], centz: &[f64],
    vol: &[f64], 
    jx: &[f64], jy: &[f64], jz: &[f64], 
    x: &[f64], y: &[f64], z: &[f64], 
    bx: &mut [f64], by: &mut [f64], bz: &mut [f64], 
    theta: f64, leaf_threshold: u32, 
    nthreads_requested: u32
) -> Result<(),()>{
    // TODO: length checks
    let n: usize = centx.len();
    let nthreads: usize = get_nthreads(nthreads_requested);
    let chunk_size: usize = (n / nthreads).max(1);

    // Build the source octree 
    let max_depth: u8 = 21; 
    let tree = SourceOctree::from_source_points(
        (&centx, &centy, &centz), vol, (&jx, &jy, &jz), max_depth, leaf_threshold
    ).unwrap();

    // chunk the inputs 
    let _x = x.par_chunks(chunk_size);
    let _y = y.par_chunks(chunk_size);
    let _z = z.par_chunks(chunk_size);
    let _bx = bx.par_chunks_mut(chunk_size);
    let _by = by.par_chunks_mut(chunk_size);
    let _bz = bz.par_chunks_mut(chunk_size);

    (_x, _y, _z, _bx, _by, _bz)
        .into_par_iter()
        .try_for_each(|(_x, _y, _z, _bx, _by, _bz)| {bfield_node_chunk(&tree, _x, _y, _z, _bx, _by, _bz, theta)})?;

    Ok(())
}