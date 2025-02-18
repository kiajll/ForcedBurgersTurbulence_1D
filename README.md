# ForcedBurgersTurbulence_1D
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

### 📜 License
This project is licensed under the MIT License - see the [main.cpp](main.cpp) file for details.

