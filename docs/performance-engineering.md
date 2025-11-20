# Performance Engineering Notes for `thor`

## Double Floating-Point Precision
Generally, all scientific calculations performed in `thor` use double-precision
(64-bit) floating point variables. It's possible that single-precision calculations
would work successfully, but the overhead of proving this for all problems was 
not deemed worth the effort. 