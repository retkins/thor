# `thor`

Approximate Biot-Savart Law integration on finite element meshes in near-linear time 

## Abstract 

The Biot-Savart Law for magnetic field calculations is well understood to be a $O(N^2)$ algorithm: given a set of $N$ sources and $M$ target points, compute $N \cdot M$ interactions. This produces highly exact calculations at the expense of high numerical cost, particularly when $N$ and $M$ are large. 

This software package provides a means of performing the *approximate* calculation in $O(N \space log(N))$ time, which provides orders of magnitude increase in the speed of the calculation while maintaining acceptable error. 

## Problem Statement 

Given a finite element mesh of $N$ elements, each with a specified current density vector, compute either: 
1. The effect of those current-carrying sources on the magnetic flux density at a set of points $M$, or 
2. The self-field induced by this collection of elements on itself 

### Acceptable Error
In particular, this software package is concerned with *approximate* (sometimes called "low-accuracy") problems, where some small error is acceptable. Acceptable error for `thor` is defined as a deviation < 0.1% ($10^{-3}$) from the equivalent fully-integrated Biot Savart Law solution. 

### Problem Size 
This software package is concerned primarily with problems that are in the 1-10M source/target element range; i.e. up to $1 \cdot 10^{14}$ interactions. This problem size is on the upper bound of what is typically produced by commercial finite element programs for solution on standard workstation computers. 

### Benchmarking Configuration
All benchmarks are produced using the following workstation:
- AMD Ryzen 9 9950X, a 16-core high-end consumer-grade CPU released in late 2024
  - Base clock: 4.3 GHz (boost clock speed disabled)
  - Hyperthreading disabled 
  - 80 kB L1 cache per core; 16 MB and 64 MB L2 and L3 cache
- 64 GB DDR5 RAM. 

All performance-critical code is written in C99 and compiled with GCC 15.2 (released Aug 2025) using the following flags: `-O3 -march=native -funroll-loops` 
  - Of critical note, `-ffast-math` is *disabled*

## Biot-Savart Law for Magnetic Field Calculations

The general form of the Biot-Savart Law is as follows: 

Given an observation point in 3D space $r$, compute the effect of a solid conductor made of differential volume elements $dV$, each with current density $J$: 
$$ \vec{B}(r) = \frac{\mu_0}{4\pi} \int \int \int_V \frac{\vec{J} \times \vec{r'}}{{|\vec{r'}|}^3} dV$$

where $\mu_0 = 4 \pi \cdot 10^{-7} H/m$ is the magnetic permeability constant for free space, $\vec{B}$ is the magnetic flux density in units of $T$, and $\vec{J}$ and $\vec{r}$ have units of $A/m^2$ and $m$. 

Due to the complexity of the volume integral, only a few closed-form solutions are known for specific source geometries. 

### Point Source 
The point source formulation is the most simple and the one used for finite element integration. It assumes the current density is constant over the finite element, such that the source can be approximated as acting at the element centroid:
$$ \vec{B}(r) = \frac{\mu_0}{4\pi} \cdot V_e \cdot \frac{\vec{J} \times \vec{r'}}{{|\vec{r'}|}^3} $$

The volume integral has been replaced by a simple multiplication of the volume of the element by the remainder of the equation. In fact, this is essentially the point of the finite element process: given a complex physical phenomenon occurring over a generic 3D body, compute the physics using a collection of much smaller points over which the solution is well-posed.  

### Finite Length Wire Source

See: https://freestatelabs.github.io/Wired.jl/dev/theory/. 
(derived from https://web.mit.edu/6.013_book/www/chapter8/8.2.html)

### Ring Source: 
See: https://freestatelabs.github.io/Wired.jl/dev/theory/. 
(derived from https://ntrs.nasa.gov/citations/20010038494)

## Barnes-Hut and Octree N-body Calculations


