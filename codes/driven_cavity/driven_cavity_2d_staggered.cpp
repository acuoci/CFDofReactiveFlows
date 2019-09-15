// ----------------------------------------------------------------------- //
//                                    __  __  __       _  __   __          //
//      |\ /|  _  |_ |  _  |_  |__|  /   |_  |  \  _  (_ |__) | _          //
//      |   | (_| |_ | (_| |_)    |  \__ |   |__/ (_) |  |  \ |            //
//                                                                         //
// ----------------------------------------------------------------------- //
//                                                                         //
//   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       //
//   CRECK Modeling Group <http ://creckmodeling.chem.polimi.it>           //
//   Department of Chemistry, Materials and Chemical Engineering           //
//   Politecnico di Milano                                                 //
//   P.zza Leonardo da Vinci 32, 20133 Milano                              //
//                                                                         //
// ----------------------------------------------------------------------- //
//                                                                         //
//   This file is part of Matlab4CFDofRF framework.                        //
//                                                                         //
//   License                                                               //
//                                                                         //
//   Copyright(C) 2019 Alberto Cuoci                                       //
//   Matlab4CFDofRF is free software : you can redistribute it and/or      //
//   modify it under the terms of the GNU General Public License as        //
//   published by the Free Software Foundation, either version 3 of the    //
//   License, or (at your option) any later version.                       //
//                                                                         //
//   Matlab4CFDofRF is distributed in the hope that it will be useful,     //
//   but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the           //
//   GNU General Public License for more details.                          //
//                                                                         //
//   You should have received a copy of the GNU General Public License     //
//   along with Matlab4CRE.If not, see <http://www.gnu.org/licenses/>.     //
//                                                                         //
//------------------------------------------------------------------------ //
//                                                                         //
//  Code : 2D driven - cavity problem in pressure/velocity formulation     //
//        The code is adapted and extended from Tryggvason, Computational  //
//        Fluid Dynamics http ://www.nd.edu/~gtryggva/CFD-Course/          //
//                                                                         //
// ----------------------------------------------------------------------- //

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseLU>
#include <stdio.h>
#include <iostream>
#include <fstream>

// Function to write output files in Tecplot format
void WriteTecplotFile(	const std::string filename,
						const Eigen::VectorXd& x, const Eigen::VectorXd& y,
						const Eigen::MatrixXd& u, const Eigen::MatrixXd& v,
						const Eigen::MatrixXd& p);

int main()
{
	// Basic setup
	// Use only even numbers for number of cells
	const unsigned int nx = 200;						// number of cells along x
	const unsigned int ny = nx;							// number of cells along y
	const double L = 1.;								// length [m]
	const double h = 1. / static_cast<double>(nx);		// grid step along [m]
	const double nu = 1e-3;								// kinematic viscosity [m2/s]
	const double tau = 20;								// total time of simulation [s]


	// Parameters for SOR
	const unsigned int max_iterations = 10000;  // maximum number of iterations
	const double beta = 1.5;					// SOR coefficient
	const double max_error = 1e-5;				// error for convergence

	// Boundary conditions
	const double un = 1;       // north wall velocity [m/s]
	const double us = 0;       // south wall velocity [m/s]
	const double ve = 0;       // east wall velocity [m/s]
	const double vw = 0;       // west wall velocity [m/s]

	// Time step
	const double sigma = 0.5;										// safety factor for time step (stability)
	const double u2 = (un*un + us*us + ve*ve + vw*vw);				// squared velocity order of magnitude [m2/s2]
	const double dt_diff = h*h/4./nu;								// time step (diffusion stability) [s]
	const double dt_conv = 4*nu/u2;									// time step (convection stability) [s]
	const double dt = sigma*std::min(dt_diff, dt_conv);				// time step (stability) [s]
	const unsigned int nsteps = static_cast<unsigned int>(tau/dt);  // number of steps
	const double Re = un*L / nu;									// Reynolds' number

	// Output options
	const bool write_tecplot = false;

	// Summary
	std::cout << "Reynolds' number: " << Re << std::endl;
	std::cout << "Time step: " << dt << std::endl;
	std::cout << " - Diffusion: " << dt_diff << std::endl;
	std::cout << " - Convection: " << dt_conv << std::endl;

	// Memory allocation
	Eigen::MatrixXd u(nx + 1, ny + 2); u.setZero();				// x-velocity [m/s]
	Eigen::MatrixXd v(nx + 2, ny + 1); v.setZero();				// y-velocity [m/s]
	Eigen::MatrixXd p(nx + 2, ny + 2); p.setZero();				// pressure [Pa]
	Eigen::MatrixXd ut(nx + 1, ny + 2); ut.setZero();			// temporary x-velocity [m/s]
	Eigen::MatrixXd vt(nx + 2, ny + 1); vt.setZero();			// temporary y-velocity [m/s]
	Eigen::MatrixXd uu(nx + 1, ny + 1); uu.setZero();			// reconstructed x-velocity [m/s]
	Eigen::MatrixXd vv(nx + 1, ny + 1); vv.setZero();			// reconstructed y-velocity [m/s]
	Eigen::MatrixXd pp(nx + 1, ny + 1); pp.setZero();			// reconstructd pressure [Pa]

	// Grid construction
	Eigen::VectorXd x(nx+1);         // grid coordinates(x axis)
	Eigen::VectorXd y(ny+1);         // grid coordinates(y axis)
	for (unsigned int i = 0; i < nx+1; i++)
		x(i) = h*i;
	for (unsigned int j = 0; j < ny+1; j++)
		y(j) = h*j;

	// Coefficient for pressure equation
	Eigen::MatrixXd gamma(nx + 2, ny + 2); gamma.setConstant(1./4.);
	gamma(1,1) = 1./2.; 
	gamma(1,ny) = 1./2.; 
	gamma(nx,1) = 1./2.; 
	gamma(nx,ny) = 1./2.;
	for (unsigned int i = 2; i<ny; i++)
	{
		gamma(1, i) = 1./3.;
		gamma(nx, i) = 1./3.;
	}
	for (unsigned int i=2; i<nx; i++)
	{
		gamma(i, 1) = 1./3.;
		gamma(i, ny) = 1./3.;
	}

	// Time loop
	double t = 0;
	for (unsigned int istep = 1; istep <= nsteps; istep++)
	{
		// ------------------------------------------------------------------ //
		// Boundary conditions
		// ------------------------------------------------------------------ //
		for (unsigned int i = 0; i<nx+1; i++)
		{
			u(i, 0)    = 2.*us - u(i, 1);      // south wall
			u(i, ny+1) = 2.*un - u(i, ny);     // north wall
		}

		for (unsigned int i = 0; i < ny + 1; i++)
		{
			v(0, i) = 2.*vw - v(1, i);         // west wall
			v(nx+1, i) = 2.*ve - v(nx, i);     // east wall
		}

		// ------------------------------------------------------------------ //
		// Temporary velocity along x
		// ------------------------------------------------------------------ //
		for (unsigned int i=1; i<nx; i++)
			for (unsigned int j=1; j<ny+1; j++)
			{
				const double ue2 = 0.25*(u(i + 1, j) + u(i, j))*(u(i + 1, j) + u(i, j));
				const double uw2 = 0.25*(u(i, j) + u(i - 1, j))*(u(i, j) + u(i - 1, j));
				const double unv = 0.25*(u(i, j + 1) + u(i, j))*(v(i + 1, j) + v(i, j));
				const double usv = 0.25*(u(i, j) + u(i, j - 1))*(v(i + 1, j - 1) + v(i, j - 1));

				const double A = (ue2 - uw2 + unv - usv) / h;
				const double D = (nu/h/h)*(u(i + 1, j) + u(i - 1, j) + u(i, j + 1) + u(i, j - 1) - 4 * u(i, j));

				ut(i, j) = u(i, j) + dt*(-A + D);
			}

		// ------------------------------------------------------------------ //
		// Temporary velocity along x
		// ------------------------------------------------------------------ //
		for (unsigned int i = 1; i<nx+1; i++)
			for (unsigned int j = 1; j < ny; j++)
			{
				const double vn2 = 0.25*(v(i, j + 1) + v(i, j))*(v(i, j + 1) + v(i, j));
				const double vs2 = 0.25*(v(i, j) + v(i, j - 1))*(v(i, j) + v(i, j - 1));
				const double veu = 0.25*(u(i, j + 1) + u(i, j))*(v(i + 1, j) + v(i, j));
				const double vwu = 0.25*(u(i - 1, j + 1) + u(i - 1, j))*(v(i, j) + v(i - 1, j));
				const double A = (vn2 - vs2 + veu - vwu) / h;
				const double D = (nu/h/h)*(v(i + 1, j) + v(i - 1, j) + v(i, j + 1) + v(i, j - 1) - 4 * v(i, j));

				vt(i, j) = v(i, j) + dt*(-A + D);
			}

		// ------------------------------------------------------------------ //
		// Poisson equation(SOR)
		// ------------------------------------------------------------------ //
		unsigned int iter = 0;
		{
			for (iter = 1; iter <= max_iterations; iter++)
			{
				for (unsigned int i = 1; i < nx+1; i++)
					for (unsigned int j = 1; j < ny+1; j++)
					{
						const double delta = p(i + 1, j) + p(i - 1, j) + p(i, j + 1) + p(i, j - 1);
						const double S = (h / dt)*(ut(i, j) - ut(i - 1, j) + vt(i, j) - vt(i, j - 1));
						p(i, j) = beta*gamma(i, j)*(delta - S) + (1 - beta)*p(i, j);
					}

				// Estimate the error
				double epsilon = 0.0;
				for (unsigned int i = 1; i < nx + 1; i++)
					for (unsigned int j = 1; j < ny + 1; j++)
					{
						const double delta = p(i + 1, j) + p(i - 1, j) + p(i, j + 1) + p(i, j - 1);
						const double S = (h / dt)*(ut(i, j) - ut(i - 1, j) + vt(i, j) - vt(i, j - 1));
						epsilon += std::abs( p(i, j) - gamma(i, j)*(delta - S));
					}
				epsilon /= (nx*ny);

				// Check the error
				if (epsilon <= max_error) // stop if converged
					break;
			}
		}
		
		// ------------------------------------------------------------------ //
		// Correct the velocity
		// ------------------------------------------------------------------ //
		for (unsigned int i = 1; i < nx; i++)
			for (unsigned int j = 1; j < ny+1; j++)
				u(i,j) = ut(i,j) - (dt/h)*(p(i+1,j)-p(i,j));
		
		for (unsigned int i = 1; i < nx+1; i++)
			for (unsigned int j = 1; j < ny; j++)
				v(i, j) = vt(i, j) - (dt/h)*(p(i,j+1)-p(i,j));
		
		// ------------------------------------------------------------------ //
		// Interpolate velocities and pressure on the same FD grid
		// ------------------------------------------------------------------ //
		for (unsigned int i = 0; i < nx+1; i++)
			for (unsigned int j = 0; j < ny+1; j++)
			{
				uu(i,j)  = 0.5*(u(i,j+1) + u(i,j));
				vv(i,j)  = 0.5*(v(i+1,j) + v(i,j));
				pp(i, j) = 0.25*(p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1));
			}
		

		// Post-processing operations
		if (istep % 10 == 1)
		{
			std::cout << "Step: " << istep << " - Time: " << t << " - Poisson iterations: " << iter << std::endl;

			if (write_tecplot == true)
			{
				std::stringstream label; label << static_cast<int>(istep/10);
				std::string filename = "Solution.tec." + label.str();
				WriteTecplotFile(filename, x, y, uu, vv, pp);
			}

		}

		// Update time step
		t += dt;
	}

	// ------------------------------------------------------------------ //
	// Write velocity profiles along the centerlines for exp comparison
	// ------------------------------------------------------------------ //
	unsigned int iaxis = static_cast<unsigned int>(nx/2);
	unsigned int jaxis = static_cast<unsigned int>(ny/2);

	std::ofstream fileVertical("vertical.txt", std::ios::out);
	fileVertical.setf(std::ios::scientific);
	for (unsigned int i = 0; i < ny+1; i++)
		fileVertical << y(i) << " " << uu(iaxis, i) << std::endl;
	fileVertical.close();

	std::ofstream fileHorizontal("horizontal.txt", std::ios::out);
	fileHorizontal.setf(std::ios::scientific);
	for (unsigned int i = 0; i < nx+1; i++)
		fileHorizontal << x(i) << " " << vv(i, jaxis) << std::endl;
	fileHorizontal.close();

	// ------------------------------------------------------------------ //
	// Write Tecplot file
	// ------------------------------------------------------------------ //
	WriteTecplotFile("SolutionFinal.tec", x, y, uu, vv, pp);

	std::cout << "Calculations completed" << std::endl;
	std::cout << "Press enter to exit..." << std::endl;
	getchar();
	
	return 0;
}

void WriteTecplotFile(	const std::string filename, 
						const Eigen::VectorXd& x, const Eigen::VectorXd& y,
						const Eigen::MatrixXd& u, const Eigen::MatrixXd& v,
						const Eigen::MatrixXd& p)
{
	std::ofstream fTecplot(filename.c_str(), std::ios::out);
	fTecplot.setf(std::ios::scientific);
	fTecplot << "Title = Solution" << std::endl;
	fTecplot << "Variables = \"x\", \"y\", \"u\", \"v\", \"p\" " << std::endl;
	fTecplot << "Zone I = " << x.size() << ", J = " << y.size() << ", F = POINT" << std::endl;

	for (unsigned int i = 0; i < x.size(); i++)
		for (unsigned int j = 0; j < y.size(); j++)
			fTecplot << x(i) << " " << y(j) << " " << u(i, j) << " " << v(i, j) << " " << p(i, j) << std::endl;
	fTecplot.close();
}
