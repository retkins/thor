pub struct TetSources {
    pub nodes: Vec<[[f64; 3]; 4]>,      // of the element, not the tree
    pub jdensity: Vec<[f64; 3]>,
    pub centroid: Vec<[f64; 3]>,
    pub volume: Vec<f64>,
}