# ForcedBurgersTurbulence_1D
![untitled01](https://github.com/user-attachments/assets/3604d7db-33d7-49e3-922e-61e710b92dae)
## Overview
This repository contains a C++ implementation of a spectral method for solving the Forced Burgers’ Equation using Fourier-Galerkin discretization and FFTW for efficient spectral transforms. The code is optimized for high-performance computing with efficient memory management.
### Equation Formulation
The Forced Burgers’ Equation is a fundamental model for turbulence and nonlinear wave dynamics:

$$\frac{\partial u}{\partial t}+ u \frac{\partial u}{\partial x}=\nu \frac{\partial ^2 u}{\partial x ^2}+f(x,t)$$

where:
* $u=u(x,t)$: Velocity field,
* $\nu$: Kinematic viscosity coefficient,
* $f(x,t)$: External forcing term.
This equation is solved in a periodic domain $x \in [0,L]$ , using spectral methods.
### Forcing term
We injected a large-scale forcing term $f(x,t)$ into the system that is Gaussian white-in-time and acts on the wavenumber range $k\in [1,3]$.
## Numerical Methods
This code uses spectral Fourier-Galerkin methods, which leverage Fast Fourier Transforms (FFTs) for spatial differentiation:
* Pseudo-spectral method for computing nonlinear terms.
* Explicit time integration using a third-order Runge-Kutta scheme.
* Anti-aliasing via the 3/2-rule to reduce energy pile-up at high wavenumbers.
* FFTW for efficient Fourier transforms.
### 1. Fourier Discretization
The velocity field is represented in Fourier space as:

$$ u\left(x\right) = \sum_{k}{\hat{u}_k(t)e^{ikx}} $$

where:
* $\hat{u}_k$ are Fourier coefficients,
* $k$ are wavenumbers in a periodic domain.

Spatial derivatives are computed as:

$$\frac{\partial u}{\partial x} \rightarrow ik \hat{u}_k$$
$$\frac{\partial^2u}{\partial x^2} \rightarrow -k^2\hat{u}_k$$

### 2. Time Integration (Explicit compact RK3)
Time integration is performed using a third-order compact Runge-Kutta (RK3) scheme:

$$x(t+h)=x(t) + \frac{1}{9}(2K_1~+~3K_2~+~4K_3)$$

where:

$$K_1=h f(t,x)$$
$$K_2=h f(t+\frac{1}{2}h,~x+\frac{1}{2}K_1)$$
$$K_3=h f(t+\frac{3}{4}h,~x+\frac{3}{4}K_2)$$

And $K$ is computed right-hand side containing advection, diffusion, and forcing terms.
### 3. Aliasing Removal (3/2-Rule)
To avoid aliasing errors in nonlinear terms, the 3/2-rule is used:
* The spectral domain is padded with zeros before computing nonlinear terms.
* The result is then projected back onto the original wavenumber range.
## Technical Aspects
### 1. Memory Management
* FFTW uses `fftw_malloc` for aligned memory allocation, ensuring compatibility with SIMD vectorization.
* Raw pointers are used for FFTW arrays to optimize performance.
* All allocated memory is properly freed using `fftw_free`.
### 2. Usage of Pointers
* Fourier coefficients `u_hat` are stored as pointers using FFTW’s fftw_complex type.
* Functions modify `u_hat` in-place, reducing memory overhead.
### 3. FFTW Library Usage
FFTW `fftw3.h` is used for spectral transforms, including:
* Forward FFT: Converts real-space velocity field to Fourier space.
* Backward FFT: Converts spectral derivatives back to real space.
* In-place computation: FFTW plans are created once and reused.
### 4. Computational Efficiency
* FFTW pre-planning `FFTW_ESTIMATE` reduces computation time.
* Vectorized operations in spectral space accelerate performance.
* Explicit time integration avoids implicit matrix inversions.
## Installation & Compilation
### Requirements
* FFTW3: Install FFTW library (`fftw3.h` must be accessible).
* C++ Compiler: GCC, Clang, or MSVC with C++17 support.


### 📜 License
This project is licensed under the MIT License - see the [main.cpp](main.cpp) file for details.

