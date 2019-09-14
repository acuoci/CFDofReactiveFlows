# Matlab4CFDofRF
Collection of codes in Matlab(R) for solving basic problems presented and discussed in the "Computational Fluid Dynamics of Reactive Flows" course (Politecnico di Milano)

## 1. Advection-diffusion equation in 1D
The advection-diffusion equation is solved on a 1D domain using the finite-difference method. Constant, uniform velocity and diffusion coefficients are assumed. The forward (or explicit) Euler method is adopted for the time discretization, while spatial derivatives are discretized using 2nd-order, centered schemes.
* Matlab script: [advection_diffusion_1d.m](codes/advection_diffusion_1d.m)
* Matlab live script: [advection_diffusion_1d_live.mlx](codes/advection_diffusion_1d_live.mlx)

## 2. Advection-diffusion equation in 2D
The advection-diffusion equation is solved on a 2D rectangular domain using the finite-difference method. Constant, uniform velocity components and diffusion coefficients are assumed. The forward (or explicit) Euler method is adopted for the time discretization, while spatial derivatives are discretized using 2nd-order, centered schemes.
* Matlab script: [advection_diffusion_2d.m](codes/advection_diffusion_2d.m)
* Matlab live script: [advection_diffusion_2d_live.mlx](codes/advection_diffusion_2d_live.mlx)

## 3. Poisson equation in 2D
The Poisson equation is solved on a 2D rectangular domain using the finite-difference method. A constant source term is initially adopted. Spatial derivatives are discretized using 2nd-order, centered schemes. Different methods are adopted for solving the equation: the Jacobi method, the Gauss-Siedler method, and the Successive Over-Relaxation (SOR) method
* Matlab script: [poisson_2d.m](codes/poisson_2d.m)
* Matlab live script: [poisson_2d_live.mlx](codes/poisson_2d_live.mlx)
