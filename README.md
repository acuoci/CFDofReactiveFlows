# Matlab4CFDofRF
Collection of codes in Matlab(R) for solving basic problems presented and discussed in the "Computational Fluid Dynamics of Reactive Flows" course (Politecnico di Milano)

## 1. Advection-diffusion equation in 1D
The advection-diffusion equation is solved on a 1D domain using the finite-difference method. Constant, uniform velocity and diffusion coefficients are assumed. The forward (or explicit) Euler method is adopted for the time discretization, while spatial derivatives are discretized using 2nd-order, centered schemes.
* Matlab script: [advection_diffusion_1d.m](codes/advection_diffusion_1d.m)
* Matlab live script: [advection_diffusion_1d_live.mlx](codes/advection_diffusion_1d_live.mlx)
