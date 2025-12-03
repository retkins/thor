# Portable SIMD 
*This code will eventually be moved to its own repository called `solar`.* 

Portable SIMD instrinsics for X86_64 and Apple ARM targets: 
* ARM NEON 
* AVX2 
* AVX512

## Background

This single-header library provides basic portable SIMD intrinsics for use across projects compiled for multiple targets. Simply add this header to your include directory and enjoy easy-to-read multi-platform SIMD intrinsics. 

This project is under development and only a fractional subset of SIMD operations are included. I plan to only ones I use, as I use them, so only basic SIMD operations will ever be included. 

### Roadmap 
- [ ] Double precision float operations for the following: 
  - [x] Load 
  - [ ] Store 
  - [ ] Set 
  - [ ] Add 
  - [ ] Subtract 
  - [ ] Multiply 
  - [ ] Fused Multiply-Add
  - [ ] Divide 
  - [ ] Inverse
  - [ ] Square root 
  - [ ] Inverse square root 
  - [ ] Norm (magnitude) of a 3-length vector 
  - [ ] Compare a <= b to generate a mask
  - [ ] Compare a >= b to generate a mask
  - [ ] Blend using mask
- [ ] Single-precision float operations for the above



## Data Structures 

This library defines two constants, which provide the vector length on the relevant platform:

| Name | Description | Value on NEON | Value on AVX2 | Value on AVX512 |
| --- | --- | --- | --- | --- | 
| `VENS` | Length of a single-precision (32-bit) vector | 4 | 8 | 16 | 
| `VEND` | Length of a double-precision (64-bit) vector | 2 | 4 | 8 | 

Supported vector types include:

| Name | Description | NEON Equivalent | AVX Equivalent | AVX512 Equivalent | 
| --- | --- | --- | --- | --- | 
| `vps` | **V**ector of **P**acked **S**ingle-precision Floats | `float32x4_t` | `__m256` | `__m512` | 
| `vpd` | **V**ector of **P**acked **D**ouble-precision Floats | `float64x2_t` | `__m256d` | `__m512d` | 
| `vpsu` | **V**ector of **P**acked **S**ingle-precision **U**nsigned Integers | `uint32x4_t` | `__m256i` | `__m512i` | 
| `vpdu` | **V**ector of **P**acked **D**ouble-precision **U**nsigned Integers | `uint64x2_t` | `__m256i` | `__m512i` | 

Note that AVX2/512 instrinsics do not distinguish between unsigned/signed and single/double-precision integers. 

## Operations 

### Load/Store/Set 

| Name | Description | NEON Instruction  | AVX2 Instruction | AVX512 Instruction | 
| --- | --- | --- | --- | --- | 
| `vpd loadpd(double a[VLEND])` | Load `VLEND` number of elements from an array into a vector | `vld1q_f64()` | `_mm256_loadu_pd()` |`_mm512_loadu_pd()`|
| `void storepd(double a[VLEND], vpd b)` | Store `VLEND` number of elements from a vector into an array | `vst1q_f64()` | `_mm256_storeu_pd` |`_mm512_storeu_pd()`|
| `vpd setpd(double a)` | Create a vector of `VLEND` elements with each set to the value of `a` | `vst1q_f64()` | `_mm256_set1_pd` |`_mm512_set1_pd()`|

### Arithmetic 
| Name | Description | NEON Instruction  | AVX2 Instruction | AVX512 Instruction | 
| --- | --- | --- | --- | --- | 
| `vpd addpd(vpd a, vpd b)` | Add two vectors of double-precision floats | `vaddq_f64()` | `_mm256_add_pd()` | `_mm512_add_pd()` | 
| `vpd subpd(vpd a, vpd b)` | Subtract `b` from `a` | `vsubq_f64()` | `_mm256_sub_pd()` | `_mm512_sub_pd()` | 
| `vpd mulpd(vpd a, vpd b)` | Multiply `a` by `b` | `vmulq_f64()` | `_mm256_mul_pd()` | `_mm512_mul_pd()` | 
| `vpd fmapd(vpd a, vpd b, vpd c)` | Fused multiply-add; return `a*b + c`| `vfmaq_f64()` | `_mm256_fmadd_pd()` | `_mm512_fmadd_pd()` |
| `vpd divpd(vpd a, vpd b)` | Divide `a` by `b` | `vdivq_f64()` | `_mm256_div_pd()` | `_mm512_div_pd()` | 

### Elementary Math Functions 


### Special Functions 

## Known Issues and Bugs 
- AVX intrinsics differentiate between aligned and unaligned vector loads; only unaligned loads are supported.