# CFDofReactiveFlows
Collection of codes in Matlab(R) and C++ for solving basic problems presented and discussed in the "Computational Fluid Dynamics of Reactive Flows" course (Politecnico di Milano)

## 1. Advection-diffusion equation in 1D with the Finite Difference (FD) method
The advection-diffusion equation is solved on a 1D domain using the finite-difference method. Constant, uniform velocity and diffusion coefficients are assumed. The forward (or explicit) Euler method is adopted for the time discretization, while spatial derivatives are discretized using 2nd-order, centered schemes.
* Matlab script: [advection_diffusion_1d.m](codes/finite_difference/advection_diffusion_1d.m)
* Matlab live script: [advection_diffusion_1d_live.mlx](codes/finite_difference/advection_diffusion_1d_live.mlx)

## 2. Advection-diffusion equation in 2D with the Finite Difference (FD) method
The advection-diffusion equation is solved on a 2D rectangular domain using the finite-difference method. Constant, uniform velocity components and diffusion coefficients are assumed. The forward (or explicit) Euler method is adopted for the time discretization, while spatial derivatives are discretized using 2nd-order, centered schemes.
* Matlab script: [advection_diffusion_2d.m](codes/finite_difference/advection_diffusion_2d.m)
* Matlab live script: [advection_diffusion_2d_live.mlx](codes/finite_difference/advection_diffusion_2d_live.mlx)

## 3. Poisson equation in 2D
The Poisson equation is solved on a 2D rectangular domain using the finite-difference method. A constant source term is initially adopted. Spatial derivatives are discretized using 2nd-order, centered schemes. Different methods are adopted for solving the equation: the Jacobi method, the Gauss-Siedler method, and the Successive Over-Relaxation (SOR) method
* Matlab script: [poisson_2d.m](codes/finite_difference/poisson_2d.m)
* Matlab live script: [poisson_2d_live.mlx](codes/finite_difference/poisson_2d_live.mlx)

## 4. Navier-Stokes equations in 2D: vorticity-streamline formulation
The Navier-Stokes equations for an incompressible fluid are solved on a 2D rectangular domain according to the vorticity-streamline formulation. The vorticity advection-diffusion equation is solved using the forward Euler method and 2nd order, centered spatial discretizations. The streamline function Poisson equation is solved using the Successive Over-Relaxation method and 2nd order, centered discretization for the spatial derivatives. 
* Matlab script (square domain, uniform grid): [driven_cavity_2d_vorticity.m](codes/driven_cavity/driven_cavity_2d_vorticity.m)
* Matlab live script (square domain, uniform grid): [driven_cavity_2d_vorticity_live.mlx](codes/driven_cavity/driven_cavity_2d_vorticity_live.mlx)
* C++ code (square domain, uniform grid): [driven_cavity_2d_vorticity.cpp](codes/driven_cavity/driven_cavity_2d_vorticity.cpp)
* Matlab script (rectangular domain, non-uniform grid): [driven_cavity_2d_vorticity_nonuniform.m](codes/driven_cavity/driven_cavity_2d_vorticity_nonuniform.m)

## 5. Advection-diffusion equations in 1D with the Finite Volume (FV) method
The advection-diffusion equations are solved on a 1D domain using the finite volume method. Both explicit (forward) and implicit (backward) Euler methods are considered. Different discretization schemes for the advective term are implemented: centered, upwind, hybrid, power-law and QUICK.
* Diffusion equation with explicit (forward) Euler method and centered differencing scheme. Matlab script: [diffusion_1d.m](codes/finite_volume/diffusion_1d.m)
* Diffusion equation with implicit (backward) Euler method and centered differencing scheme. Matlab script: [diffusion_1d_implicit.m](codes/finite_volume/diffusion_1d_implicit.m)
* Steady-state advection-diffusion equation with implicit (backward) Euler method and several discretization schemes for the advective contribution (centered, upwind, hybrid, power-law). Matlab script: [steady_advection_diffusion_1d_implicit.m](codes/finite_volume/steady_advection_diffusion_1d_implicit.m)
* Steady-state advection-diffusion equation with implicit (backward) Euler method and QUICK scheme for the discretization of the advective contribution. Matlab script: [steady_advection_diffusion_1d_implicit_quick.m](codes/finite_volume/steady_advection_diffusion_1d_implicit_quick.m)

## 6. Navier-Stokes equations in 2D: staggered grid and projection algorithm
The Navier-Stokes equations for an incompressible fluid are solved on a 2D rectangular domain meshed with a staggered grid. The momentum equations are solved using the forward Euler method and 2nd order, centered spatial discretizations. The projection algorithm is adopted for managing the coupling between pressure and velocity. In particular, the corresponding Poisson equation for pressure is solved using the Successive Over-Relaxation method and 2nd order, centered discretization for the spatial derivatives. 
* Matlab script (square domain, uniform grid): [driven_cavity_2d_staggered.m](codes/driven_cavity/driven_cavity_2d_staggered.m)
* C++ code (square domain, uniform grid): [driven_cavity_2d_staggered.cpp](codes/driven_cavity/driven_cavity_2d_staggered.cpp)
* Matlab script (rectangular domain, non-uniform grid): [driven_cavity_2d_staggered_nonuniform.m](codes/driven_cavity/driven_cavity_2d_staggered_nonuniform.m)

#### Extensions/Modifications
* Addition of a **passive scalar** equation governed by the usual advection-diffusion equation without source terms. Passive scalar equation solved with the Finite Volume (FV) or Finite Difference (FD) methods. Matlab script (square domain, uniform grid, FV): [driven_cavity_2d_staggered_passivescalar_fv.m](codes/driven_cavity/driven_cavity_2d_staggered_passivescalar_fv.m). Matlab script (square domain, uniform grid, FD): [driven_cavity_2d_staggered_passivescalar_fd.m](codes/driven_cavity/driven_cavity_2d_staggered_passivescalar_fd.m). 
* Addition of **inflow and outflow boundaries** along the west and east sides, respectively. Matlab script (square domain, uniform grid): [driven_cavity_2d_staggered_inout.m](codes/driven_cavity/driven_cavity_2d_staggered_inout.m)
* Addition of inflow and outflow boundaries along the west and east sides, respectively (see above) and addition of a **passive scalar** equation. Matlab script (square domain, uniform grid): [driven_cavity_2d_staggered_inout_passivescalar.m](codes/driven_cavity/driven_cavity_2d_staggered_inout_passivescalar.m)
* Addition of a **rectangular obstacle** inside the computational domain. Matlab script (square domain, uniform grid): [driven_cavity_2d_staggered_obstacle.m](codes/driven_cavity/driven_cavity_2d_obstacle.m)
* Addition of **buoyancy** term in the momentum equation via Boussinesq approximation, together with a temperature equation. Matlab script (square domain, uniform grid): [driven_cavity_2d_staggered_buoyancy.m](codes/driven_cavity/driven_cavity_2d_staggered_buoyancy.m)
* Addition of **local residence time** equation. Matlab script (square domain, uniform grid): [driven_cavity_2d_staggered_tau.m](codes/driven_cavity/driven_cavity_2d_staggered_tau.m)
* Addition of automatic calculation of **Residence Time Distribution (RTD)** and Cumulative Distribution Function (CDF). Matlab script (square domain, uniform grid): [driven_cavity_2d_staggered_rtd.m](codes/driven_cavity/driven_cavity_2d_staggered_rtd.m)

## 7. Navier-Stokes equations in 2D: examples
### Taylor-Green vortex in 2D
The Taylor-Green vortex is an exact closed form solution of 2D, incompressible Navier-Stokes equations. This 2D decaying vortex defined in the square domain, 0-2pi, serves as a benchmark problem for testing and validation of incompressible Navier-Stokes codes. The implementation here proposed is based on the Finite Volume Method (FV) applied on a staggered mesh and coupled with the Porjection Algorithm. Matlab script (square domain, uniform grid, FV): [taylor_green_vortex_2d.m](codes/driven_cavity/taylor_green_vortex_2d.m).
