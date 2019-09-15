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

## 4. Navier-Stokes equations in 2D: vorticity-streamline formulation
The Navier-Stokes equations for an incompressible fluid are solved on a 2D rectangular domain according to the vorticity-streamline formulation. The vorticity advection-diffusion equation is solved using the forward Euler method and 2nd order, centered spatial discretizations. The streamline function Poisson equation is solved using the Successive Over-Relaxation method and 2nd order, centered discretization for the spatial derivatives. 
* Matlab script (square domain, uniform grid): [driven_cavity_2d_vorticity.m](codes/driven_cavity/driven_cavity_2d_vorticity.m)
* Matlab live script (square domain, uniform grid): [driven_cavity_2d_vorticity_live.mlx](codes/driven_cavity/driven_cavity_2d_vorticity_live.mlx)
* C++ code (square domain, uniform grid): [driven_cavity_2d_vorticity.cpp](codes/driven_cavity/driven_cavity_2d_vorticity.cpp)
* Matlab script (rectangular domain, non-uniform grid): [driven_cavity_2d_vorticity_nonuniform.m](codes/driven_cavity/driven_cavity_2d_vorticity_nonuniform.m)

## 5. Navier-Stokes equations in 2D: staggered grid and projection algorithm
The Navier-Stokes equations for an incompressible fluid are solved on a 2D rectangular domain meshed with a staggered grid. The momentum equations are solved using the forward Euler method and 2nd order, centered spatial discretizations. The projection algorithm is adopted for managing the coupling between pressure and velocity. In particular, the corresponding Poisson equation for pressure is solved using the Successive Over-Relaxation method and 2nd order, centered discretization for the spatial derivatives. 
* Matlab script (square domain, uniform grid): [driven_cavity_2d_staggered.m](codes/driven_cavity/driven_cavity_2d_staggered.m)
* Matlab script (rectangular domain, non-uniform grid): [driven_cavity_2d_staggered_nonuniform.m](codes/driven_cavity/driven_cavity_2d_staggered_nonuniform.m)
