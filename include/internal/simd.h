/*   Portable SIMD intrinsics for thor

    Fairly brittle, generally speaking. Assumes the following:
    - You're either on an Apple M-series (ARM) machine, or an X86-64 machine
    - If on X86-64, you have access to AVX2 instruction sets at a minimum


*/

#ifndef SIMD 
#define SIMD

#ifdef __APPLE__ 
// ARM NEON intrinsics (128-bit/16-byte vectors)

#include <arm_neon.h>

// Vector stride lengths (for loops)
static const int VLENS = 4;
static const int VLEND = 2; 

// Types: define packed double and packed single precision floating point vectors
// Also define masks for each 
typedef float32x4_t vps;        // Vector of Packed Single [precision float]
typedef float64x2_t vpd;        // Vector of Packed Double [precision float]
typedef uint32x4_t vpsu;        // Vector of Packed Single [precision] unsigned int
typedef uint64x2_t vpdu;        // Vector of Packed Double [precision] unsigned int

// Load two doubles from `a`
static inline vpd loadpd(double a[2]) { 
    return vld1q_f64(a);
}

// Store two doubles from `b` into `a`
static inline void storepd(double a[2], vpd b) {
    vst1q_f64(a, b);
}

// Set a scalar to all indices in a vector 
static inline vpd setpd(double a) {
    return vdupq_n_f64(a);
}

// Add two vectors of doubles
static inline vpd addpd(vpd a, vpd b) {
    return vaddq_f64(a, b);
}

// Subtract `b` from `a`
static inline vpd subpd(vpd a, vpd b) {
    return vsubq_f64(a, b);
}

// Multiply two vectors of doubles 
static inline vpd mulpd(vpd a, vpd b) {
    return vmulq_f64(a, b);
}

// Fused multiply-add, return a*b + c 
static inline vpd fmapd(vpd a, vpd b, vpd c) {
    return vfmaq_f64(c, a, b); 
}

// Divide `a` by `b`
static inline vpd divpd(vpd a, vpd b) {
    return vdivq_f64(a, b);
}

// Take the inverse of `a` 
static inline vpd invpd(vpd a) {
    return vdivq_f64(vdupq_n_f64(1.0), a);
}

// Take the element-wise square root of `a`
static inline vpd sqrtpd(vpd a) {
    return vsqrtq_f64(a);
}

// Take the element-wise inverse square root of `a` 
static inline vpd invsqrtpd(vpd a) {
    return invpd(sqrtpd(a));
}

// Take the vector magnitude (norm) of a 3-length vector
static inline vpd mag3pd(vpd x, vpd y, vpd z) {
    vpd x2_plus_y2 = fmapd(x, x, mulpd(y,y));
    return sqrtpd(fmapd(z, z, x2_plus_y2));
}

// Compare a < b
static inline vpdu cmpltpd(vpd a, vpd b) {
    return vcltq_f64(a, b);
}

// Compare a > b
static inline vpdu cmpgtpd(vpd a, vpd b) {
    return vcgtq_f64(a, b);
}

// Compare a == b
static inline vpdu cmpeqpd(vpd a, vpd b) {
    return vceqq_f64(a, b);
}

// Blend two floating point vectors 
// if mask is 1, take from a; else take from b
static inline vpd blendpd(vpd a, vpd b, vpdu mask) {
    return vbslq_f64(mask, a, b);  
}


#else  
// X86-64

#if AVX512
// X86-64 intrinsics using AVX512 (512-bit/64-byte vectors)

#else 
// X86-64 intrinsics using AVX2 (256-bit/32-byte vectors)

#endif
#endif 

#endif
