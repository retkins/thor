# Morton Code Octree for Biot Savart Law Calculation

This program uses the Biot-Savart law for finite source points (finite elements)
to calculate the magnetic field at a given location in 3D space:
```
Given: element e[i], location r[j]
dB(e[i], r[j]) = mu0/4pi * vol[i] * (J[i] x r[j]) / r'[i,j]^3
```

## Purpose
Refactor the naive pointer-chasing+regression octree to do the following:
- Encode positions using Morton codes 
- Convert all data from array of structs into struct of arrays architecture
- Traverse tree using DFS explicitly (no recursion)
- Use dual-tree traversal to evaluate source and target points via DFS

## Design Trades

### Finite Precision
Morton codes require finite precision (quantization) of the octree; this will be 
user-selectable. For example, the Morton code should fit into a 64-bit integer, 
given 21 bits for each spatial coordinate. For a 100m bounding box, this gives
`100m / (2^21) = 48um` resolution, which should be plenty for most engineering 
problems. 

### Struct of Arrays
For large structs that contain many scalar values, an array of structs approach
means that all values in the struct are more likely to be loaded together 
and available simultaneously in cache. If these are all the values used for a given 
calculation, then this approach can be efficient. However, if other values are 
used, then a struct of arrays approach might work better. 

## Method

### Initial Conditions

Calculate starts with the following: 
- `const double *x, *y, *z`: 3d coordinates of where to calculate fields (target pts)
- `double *Bx, *By, *Bz`: allocated arrays for magnetic fields vectors (targets)
- `const double *cx, *cy, *cz`: source element centroid 
- `const double *vol`: source element volumes
- `const double *Jx, *Jy, *Jz`: source element current density vectors
- `uint32_t n_sources, n_targets`: problem size, defining extent of the input arrays
- 
### Overview
1. Build two Morton octrees: one for the source points, and another for the target pts
2. Perform bottom-up pass on the source octree to evaluate the current-moment at each interior (branch) node
3. Perform dual tree depth-first search from source to target tree, subdividing the larger of the two nodes being compared, until the Barnes Hut criteria (theta) has been reached. Contributions from the calculation are stored in branch nodes
4. Perform top-down pass on the target tree, adding the branch node contributions to all children

### Determine the bounding box size
For each of the source *and* target points, compute the bounding box for each. 
- Loop through all values of `*x, *y, *z`
- Find size values: `xmin, xmax, ymin, ymax, zmin, zmax` 
- Bounding box size is determined by `(bbox_min, bbox_max) = (min(xmin, ymin, zmin), max(xmax, ymax, zmax))`

### Allocate memory for each tree
```c
typedef struct SourceTree {
    uint32_t n;             // number of source points in tree, size of all subsequent arrays
    uint64_t *morton;       // Morton code for each point
    double *x, *y, *z;      // location of each source point (to be morton-sorted)
    double *vol;            // source volume (to be morton-sorted)
    double *Jx, *Jy, *Jz    // current density vector (to be morton-sorted)

}
```

### Morton Code 

Assumptions: 
- Bounding box is a cube (all nodes are a cube). Box size is the side length.
- Box centroid is `xcent, ycent, zcent` -> same for all nodes
- Box corner dimensions are given by combinations of `(xmin, xmax), (ymin, ymax), (zmin, zmax)`
- `xmax - xmin = side length`, etc 
- `xcent = xmax - 1/2 side length`, `xcent = xmin + 1/2 side length`, etc
- All points are guaranteed to fit within the bounding box or are excluded
  
#### Quantize
To use a Morton code, all position information has to be quantized: 
1. Start with x,y,z coordinates (each a double float)
2. Choose a level of quantization; in this case we choose `L=21 -> 2^21 = 2.1M grid locations per axis`; `scale = (double)(1u << 21 - 1u)`
3. Normalize the position in the grid: `normalize(x) = (x - xmin) / (side_length) -> double`
4. Quantize the position into an integer: `quantize(x) = (uint32_t)floor(scale*normalize(x))`

#### Encode 
Starting with an (x,y,z) coordinate represented on a quantized grid as `(ix, iy, iz)`:

TODO: finish this section