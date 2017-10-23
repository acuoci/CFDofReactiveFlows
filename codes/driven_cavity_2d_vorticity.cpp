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
//	 License                                                               //
//                                                                         //
//   Copyright(C) 2017 Alberto Cuoci                                       //
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
//  Code : 2D driven - cavity problem in vorticity/streamline formulation  //
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
						const Eigen::MatrixXd& omega, const Eigen::MatrixXd& psi);

// Function to set the sparse matrix solver as alternative to the SOR method
void SetSparseMatrixSolver(	const int nx, const int ny,
							Eigen::SparseLU< Eigen::SparseMatrix<double> >& solver);

int main()
{
	// Basic setup
	const unsigned int nx = 25;						    // number of grid points along x
	const unsigned int ny = nx;							// number of grid points along y
	const double h = 1. / static_cast<double>(nx - 1);	// grid step along[-]
	const double Re = 100;								// Reynolds number[-]
	const double tau = 10;								// total time of simulation[-]

	// Parameters for SOR
	const unsigned int max_iterations = 10000;  // maximum number of iterations
	const double beta = 1.9;					// SOR coefficient
	const double max_error = 0.0001;			// error for convergence

	// Data for reconstructing the velocity field
	const double L = 1;                    // length[m]
	const double nu = 1e-6;                // kinematic viscosity[m2/s]
	const double Uwall = nu*Re / L;        // wall velocity[m/s]

	// Time step
	const double sigma = 0.5;										// safety factor for time step(stability)
	const double dt_diff = h*h * Re / 4;							// time step(diffusion stability)
	const double dt_conv = 4 / Re;									// time step(convection stability)
	const double dt = sigma*std::min(dt_diff, dt_conv);				// time step(stability)
	const unsigned int nsteps = static_cast<unsigned int>(tau/dt);  // number of steps

	// Choosing the Poisson solver
	bool poisson_sor_solver = true;

	// Output options
	const bool write_tecplot = false;

	// Summary
	std::cout << "Time step: " << dt << std::endl;
	std::cout << " - Diffusion: " << dt_diff << std::endl;
	std::cout << " - Convection: " << dt_conv << std::endl;

	// Memory allocation
	Eigen::MatrixXd psi(nx, ny); psi.setZero();			// streamline function
	Eigen::MatrixXd omega(nx, ny); omega.setZero();		// vorticity
	Eigen::MatrixXd psio(nx, ny);  psio.setZero();		// streamline function at previous time
	Eigen::MatrixXd omegao(nx, ny);  omegao.setZero();  // vorticity at previous time
	Eigen::MatrixXd u(nx, ny); u.setZero();				// reconstructed dimensionless x - velocity
	Eigen::MatrixXd v(nx, ny); v.setZero();				// reconstructed dimensionless y - velocity
	Eigen::MatrixXd U(nx, ny); U.setZero();				// reconstructed x - velocity
	Eigen::MatrixXd V(nx, ny); V.setZero();				// reconstructed y - velocity

	// Grid construction
	Eigen::VectorXd x(nx);         // grid coordinates(x axis)
	Eigen::VectorXd y(ny);         // grid coordinates(y axis)
	for (unsigned int i = 0; i < nx; i++)
		x(i) = h*i;
	for (unsigned int j = 0; j < ny; j++)
			y(j) = h*j;

	// Sparse solver for Poisson equation
	Eigen::SparseLU< Eigen::SparseMatrix<double> > solver;
	Eigen::VectorXd b;
	Eigen::VectorXd solution;
	if (poisson_sor_solver == false)
	{
		b.resize(nx*ny);
		solution.resize(nx*ny);
		SetSparseMatrixSolver(nx, ny, solver);
	}

	// Time loop
	double t = 0;
	for (unsigned int istep = 1; istep <= nsteps; istep++)
	{
		// ------------------------------------------------------------------ //
		// Sparse LU factorization
		// ------------------------------------------------------------------ //
		if (poisson_sor_solver == false)
		{
			b.setZero();
			for (unsigned int i = 1; i < nx - 1; i++)
				for (unsigned int j = 1; j < ny - 1; j++)
				{
					const int index = j*nx + i;
					b(index) = h*h*omega(i, j);
				}

			//std::cout << b << std::endl;
			solution = solver.solve(b);
			for (unsigned int i = 0; i < nx; i++)
				for (unsigned int j = 0; j < ny; j++)
				{
					const int index = j*nx + i;
					psi(i, j) = solution(index);
				}
		}

		// ------------------------------------------------------------------ //
		// Poisson equation(SOR)
		// ------------------------------------------------------------------ //
		unsigned int iter = 0;
		if (poisson_sor_solver == true)
		{
			for (iter = 1; iter <= max_iterations; iter++)
			{
				psio = psi;
				for (unsigned int i = 1; i < nx - 1; i++)
					for (unsigned int j = 1; j < ny - 1; j++)
					{
						// solve for the stream function by SOR iteration
						psi(i, j) = 0.25*beta*(psi(i + 1, j) + psi(i - 1, j) + psi(i, j + 1) +
							psi(i, j - 1) + h*h*omega(i, j)) + (1.0 - beta)*psi(i, j);
					}

				// Estimate the error
				double epsilon = 0.0;
				for (unsigned int i = 0; i < nx; i++)
					for (unsigned int j = 0; j < ny; j++)
						epsilon += std::abs(psio(i, j) - psi(i, j));

				// Check the error
				if (epsilon <= max_error) // stop if converged
					break;
			}
		}
		
		// ------------------------------------------------------------------ //
		// Find vorticity on boundaries
		// ------------------------------------------------------------------ //

		for (unsigned int i = 1; i < nx - 1; i++)		// south
			omega(i, 0) = -2.0*psi(i, 1) / (h*h);

		for (unsigned int i = 1; i < nx - 1; i++)	// north
			omega(i, ny - 1) = -2.0*psi(i, ny - 2) / (h*h) - 2.0 / h*1.;

		for (unsigned int i = 1; i < ny - 1; i++)	// east
			omega(0, i) = -2.0*psi(1, i) / (h*h);

		for (unsigned int i = 1; i < ny - 1; i++)	// west
			omega(nx - 1, i) = -2.0*psi(nx - 2, i) / (h*h);

		// ------------------------------------------------------------------ //
		// Find new vorticity in interior points
		// ------------------------------------------------------------------ //
		omegao = omega;
		for (unsigned int i = 1; i < nx - 1; i++)
			for (unsigned int j = 1; j < ny - 1; j++)
			{
				omega(i, j) = omegao(i, j) + dt*(-0.25*((psi(i, j + 1) - psi(i, j - 1))*
							(omegao(i + 1, j) - omegao(i - 1, j)) - (psi(i + 1, j) - psi(i - 1, j))*
							(omegao(i, j + 1) - omegao(i, j - 1))) / (h*h) +
							1. / Re*(omegao(i + 1, j) + omegao(i - 1, j) + omegao(i, j + 1) +
								omegao(i, j - 1) - 4.*omegao(i, j)) / (h*h));
			}

		// Post-processing operations
		if (istep % 10 == 1)
		{
			std::cout << "Step: " << istep << " - Time: " << t << " - Poisson iterations: " << iter << std::endl;

			if (write_tecplot == true)
			{
				std::stringstream label; label << static_cast<int>(istep/10);
				std::string filename = "Solution.tec." + label.str();
				WriteTecplotFile(filename, x, y, u, v, omega, psi);
			}

		}

		// Update time step
		t += dt;

		// ------------------------------------------------------------------ //
		// Reconstruction of dimensionless velocity field
		// ------------------------------------------------------------------ //

		for (unsigned int i = 0; i < nx; i++)
			u(i, ny - 1) = 1.;
		for (unsigned int i = 1; i < nx - 1; i++)
			for (unsigned int j = 1; j < ny - 1; j++)
			{
				u(i, j) = (psi(i, j + 1) - psi(i, j - 1)) / 2. / h;
				v(i, j) = -(psi(i + 1, j) - psi(i - 1, j)) / 2. / h;
			}

		// ------------------------------------------------------------------ //
		// Reconstruction of velocity field
		// ------------------------------------------------------------------ //
		U = u*Uwall;
		V = v*Uwall;
	}

	// ------------------------------------------------------------------ //
	// Write velocity profiles along the centerlines for exp comparison
	// ------------------------------------------------------------------ //
	unsigned int iaxis = static_cast<unsigned int>(std::round((nx - 1) / 2));
	unsigned int jaxis = static_cast<unsigned int>(std::round((ny - 1) / 2));

	std::ofstream fileVertical("vertical.txt", std::ios::out);
	fileVertical.setf(std::ios::scientific);
	for (unsigned int i = 0; i < ny; i++)
		fileVertical << y(i) << " " << u(iaxis, i) << std::endl;
	fileVertical.close();

	std::ofstream fileHorizontal("horizontal.txt", std::ios::out);
	fileHorizontal.setf(std::ios::scientific);
	for (unsigned int i = 0; i < nx; i++)
		fileHorizontal << x(i) << " " << v(i, jaxis) << std::endl;
	fileHorizontal.close();

	// ------------------------------------------------------------------ //
	// Write Tecplot file
	// ------------------------------------------------------------------ //
	WriteTecplotFile("SolutionFinal.tec", x, y, u, v, omega, psi);

	std::cout << "Calculations completed" << std::endl;
	std::cout << "Press enter to exit..." << std::endl;
	getchar();
	
	return 0;
}

void WriteTecplotFile(	const std::string filename, 
						const Eigen::VectorXd& x, const Eigen::VectorXd& y,
						const Eigen::MatrixXd& u, const Eigen::MatrixXd& v,
						const Eigen::MatrixXd& omega, const Eigen::MatrixXd& psi)
{
	std::ofstream fTecplot(filename.c_str(), std::ios::out);
	fTecplot.setf(std::ios::scientific);
	fTecplot << "Title = Solution" << std::endl;
	fTecplot << "Variables = \"x\", \"y\", \"u\", \"v\", \"omega\", \"psi\" " << std::endl;
	fTecplot << "Zone I = " << x.size() << ", J = " << y.size() << ", F = POINT" << std::endl;

	for (unsigned int i = 0; i < x.size(); i++)
		for (unsigned int j = 0; j < y.size(); j++)
			fTecplot << x(i) << " " << y(j) << " " << u(i, j) << " " << v(i, j) << " " << omega(i, j) << " " << psi(i, j) << std::endl;
	fTecplot.close();
}

void SetSparseMatrixSolver(const int nx, const int ny, Eigen::SparseLU< Eigen::SparseMatrix<double> >& solver)
{
	typedef Eigen::Triplet<double> T;

	Eigen::SparseMatrix<double> A(nx*ny, nx*ny);
	std::vector<T> tripletList;
	tripletList.reserve(5 * nx*ny);

	// Internal points
	for (int i = 1; i < nx - 1; i++)
		for (int j = 1; j < ny - 1; j++)
		{
			const int index = j*nx + i;
			tripletList.push_back(T(index, index, 4.));
			tripletList.push_back(T(index, index + 1, -1.));
			tripletList.push_back(T(index, index - 1, -1.));
			tripletList.push_back(T(index, index + nx, -1.));
			tripletList.push_back(T(index, index - nx, -1.));
		}

	// South/North points
	for (int i = 0; i < nx; i++)
	{
		tripletList.push_back(T(i, i, 1.));
		tripletList.push_back(T((ny - 1)*nx + i, (ny - 1)*nx + i, 1.));
	}

	// East/West points
	for (int j = 1; j < ny - 1; j++)
	{
		tripletList.push_back(T(j*nx, j*nx, 1.));
		tripletList.push_back(T(j*nx + (nx - 1), j*nx + (nx - 1), 1.));
	}

	A.setFromTriplets(tripletList.begin(), tripletList.end());

	solver.analyzePattern(A);
	solver.factorize(A);
}


