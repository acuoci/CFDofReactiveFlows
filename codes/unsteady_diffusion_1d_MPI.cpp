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
//   Copyright(C) 2018 Alberto Cuoci                                       //
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
//  Code : unsteady 1D diffusion equation parallelized using MPI           //
//         dphi/dt=Gammad2phi/dt2                                          //
//         BC: phi(x=0) = 0; phi(x=1) = 0;                                 //
//         IC: phi(t=0)=phi0;                                              //
//                                                                         //
//         Adapted from:                                                   //
//         https://people.sc.fsu.edu/~jburkardt/cpp_src/                   //
// ----------------------------------------------------------------------- //

#include <cmath>
#include <stdio.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <mpi.h>


double boundary_condition(const double x, const double t);
double initial_condition(const double x, const double t);
double analytical_solution(const double x, const double t);
void update(const int id, const int n_procs, double* phi);

struct
{
	int		n_points;		// total number of grid points
	double	time_end;		// end time [s]
	double	length;			// length computational domain
	double	Gamma;			// diffusion coefficient [m2/s]
	double	cfl_max;		// max CFL
	double	x_delta;		// grid spacing [m]
	double	time_delta;		// time step [s]
	int		n_time_steps;	// total number of time steps
	double	cfl;			// current CFL

} data;

int main(int argc, char *argv[])
{
	int id;				// process identifier: 0=master
	int n_procs;		// total number of processes

	// Total number of grid points
	data.n_points = 1000;
	for (int i = 1; i < argc; ++i)
		if (argv[i] == "-npoints")
			data.n_points = std::strtol(argv[i + 1], NULL, 10);
	
	// Initialize MPI Communication
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	// Check number of points
	if (data.n_points % n_procs != 0)
	{
		std::cout << "ERROR!" << std::endl;
		std::cout << "  The number of grid points must be a multiple of number of processes!" << std::endl;
		return 0;
	}

	// Numerical data
	data.time_end = 1.0;		// end time [s]
	data.length = 1.0;			// length computational domain [m]
	data.Gamma = 0.002;			// diffusion coefficient [m2/s]
	data.cfl_max = 0.025;		// max CFL

	// Calculate grid step
	data.x_delta = data.length / static_cast<double>(data.n_points - 1);

	// Set time step
	data.time_delta = data.cfl_max * data.x_delta*data.x_delta / data.Gamma;

	// Set number of time steps and adjust time step
	data.n_time_steps = std::ceil(data.time_end / data.time_delta);
	data.time_delta = data.time_end / static_cast<double>(data.n_time_steps);

	//  Check the CFL condition and quit if it is too large.
	data.cfl = data.Gamma * data.time_delta / data.x_delta / data.x_delta;

	// Master: print CFL
	if (id == 0)
		std::cout << "  CFL stability criterion value = " << data.cfl << std::endl;

	// Master: print message on screen
	if (id == 0)
	{
		std::cout << std::endl;
		std::cout << "UnsteadyDiffusion1D_MPI:" << std::endl;
		std::cout << "  C++/MPI version" << std::endl;
		std::cout << "  Solve the 1D time-dependent diffusion equation." << std::endl;
		std::cout << std::endl;
		std::cout << "  Compute an approximate solution to the time dependent" << std::endl;
		std::cout << "  1D diffusion equation:" << std::endl;
		std::cout << std::endl;
		std::cout << "    dphi/dt = Gamma*d2phi/dx2" << std::endl;
		std::cout << std::endl;
		std::cout << "  for " << 0 << " = x_min < x < x_max = " << data.length << std::endl;
		std::cout << std::endl;
		std::cout << "  and " << 0. << " = time_start < t <= t_end = " << data.time_end << std::endl;
		std::cout << std::endl;
		std::cout << "  Boundary conditions are specified at x_min and x_max." << std::endl;
		std::cout << "  Initial conditions are specified at time_min." << std::endl;
		std::cout << std::endl;
		std::cout << "  The finite difference method is used to discretize the" << std::endl;
		std::cout << "  differential equation." << std::endl;
		std::cout << std::endl;
		std::cout << "  This uses " << data.n_points << " equally spaced points in X (delta_x = " << data.x_delta << " m)" << std::endl;
		std::cout << "  and " << data.n_time_steps << " equally spaced points in time (delta_t = " << data.time_delta << " s)" << std::endl;
		std::cout << std::endl;
		std::cout << "  Parallel execution is done using " << n_procs << " processors." << std::endl;
		std::cout << "  Domain decomposition is used." << std::endl;
		std::cout << "  Each processor works on " << data.n_points/n_procs << " nodes, " << std::endl;
		std::cout << "  and shares some information with its immediate neighbors." << std::endl;
	}

	// Master: start time registration
	double wtime;		// total CPU time
	if (id == 0)
		wtime = MPI_Wtime();

	// Solve the problem: computes the solution of the diffusion equation
	double *phi = new double[data.n_points / n_procs+2]; // unknown (local)
	update(id, n_procs, phi);
	
	// Master: close time registration and print exit message
	if (id == 0)
	{
		wtime = MPI_Wtime() - wtime;

		std::cout << std::endl;
		std::cout << "  Wall clock elapsed seconds = " << wtime << std::endl;
	}

	// Print solution (local)
	{
		std::string name = "solution." + std::to_string(id) + ".txt";

		std::ofstream f_solution(name.c_str(), std::ios::out);
		f_solution.setf(std::ios::scientific);
		for (int i = 1; i < (data.n_points / n_procs+2)-1; i++)
			f_solution << phi[i] << std::endl;

		f_solution.close();		
	}

	// Clean
	delete[] phi;
	
	// Close MPI environment
	MPI_Finalize();

	// Master: Reconstruct solution
	if (id == 0)
	{
		double* phi_sol = new double[data.n_points];

		int count = 0;
		for (int k = 0; k < n_procs; k++)
		{
			std::string name = "solution." + std::to_string(k) + ".txt";
			std::ifstream f_local(name.c_str(), std::ios::in);

			for (int i = 1; i < (data.n_points/n_procs+2) - 1; i++)
				f_local >> phi_sol[count++];

			f_local.close();
		}

		std::ofstream f_solution("solution.txt", std::ios::out);
		f_solution.setf(std::ios::scientific);

		double sum_error_squared = 0.;
		for (int i = 0; i < data.n_points; i++)
		{
			const double solution = analytical_solution(i * data.x_delta, data.time_end);
			const double error = phi_sol[i] - solution;
			sum_error_squared += error * error;
			f_solution	<< i * data.x_delta << " "
						<< phi_sol[i] << " "
						<< solution << " "
						<< error << std::endl;
		}

		f_solution.close();

		sum_error_squared /= static_cast<double>(data.n_points);
		std::cout << "  Absolute error: " << sum_error_squared << std::endl;
	}

	// Master: closure
	if (id == 0)
	{
		std::cout << std::endl;
		std::cout << "UnsteadyDiffusion1D_MPI:" << std::endl;
		std::cout << "  Normal end of execution." << std::endl;
		std::cout << std::endl;
	}

	return 0;
}

void update(const int id, const int n_procs, double* phi)
{
	const int n = data.n_points / n_procs;	// local number of points
	
	double *phi_new;						// updated unknown (local)
	double *x;								// grid coordinates (local)

	// Set the x coordinates of the n points
	// Ghost values of x are available as x[0] and x[n+1]
	x = new double[n + 2];
	for (int i = 0; i <= n + 1; i++)
		x[i] = ( static_cast<double>(id * n + i - 1) * data.length ) /
				 static_cast<double>(n_procs * n - 1);

	// Set the values of phi at the initial time
	double time = 0.;
	phi_new = new double[n + 2];
	phi[0] = 0.0;
	for (int i = 1; i <= n; i++)
		phi[i] = initial_condition(x[i], time);
	phi[n + 1] = 0.0;
	
	//  Compute the values of H at the next time, based on current data
	for (int j = 1; j <= data.n_time_steps; j++)
	{

		const double time_new = static_cast<double>(j)/
								static_cast<double>(data.n_time_steps) * data.time_end;
		
		// Message Passing Interface (MPI)
		{
			MPI_Status status;

			//  Send phi[1] to ID-1 (the first ID is excluded)
			if (id > 0)
			{
				const int tag = 1;
				MPI_Send(&phi[1], 1, MPI_DOUBLE, id - 1, tag, MPI_COMM_WORLD);
			}

			//  Receive phi[N+1] from ID+1 (the last ID is excluded)
			if (id < n_procs - 1)
			{
				const int tag = 1;
				MPI_Recv(&phi[n + 1], 1, MPI_DOUBLE, id + 1, tag, MPI_COMM_WORLD, &status);
			}

			//  Send phi[N] to ID+1 (the last ID is excluded)
			if (id < n_procs - 1)
			{
				const int tag = 2;
				MPI_Send(&phi[n], 1, MPI_DOUBLE, id + 1, tag, MPI_COMM_WORLD);
			}

			//  Receive phi[0] from ID-1 (the first ID is excluded)
			if (id > 0)
			{
				const int tag = 2;
				MPI_Recv(&phi[0], 1, MPI_DOUBLE, id - 1, tag, MPI_COMM_WORLD, &status);
			}
		}

		//  Update the unknown based on the four point stencil.
		for (int i = 1; i <= n; i++)
			phi_new[i] = phi[i] + data.cfl*(phi[i-1] - 2.*phi[i] + phi[i+1]);
		
		// phi at the boundaries was incorrectly computed using the differential equation.  
		// Replace that calculation by the boundary conditions.
		if (id == 0)
			phi_new[1] = boundary_condition(x[1], time_new);
		
		if (id == n_procs - 1)
			phi_new[n] = boundary_condition(x[n], time_new);
		
		//  Update time and unknown
		time = time_new;
		for (int i = 1; i <= n; i++)
			phi[i] = phi_new[i];
	}

	delete[] phi_new;
	delete[] x;

	return;
}

double boundary_condition(const double x, const double t)
{
	return 0.;
}

double initial_condition(const double x, const double t)
{
	return 100.;
}

double analytical_solution(const double x, const double t)
{
	const double phi0 = initial_condition(x, t);
	const double pi = std::acos(-1.);

	double sum = 0.;
	for (int n = 1; n <= 10000; n++)
	{
		const double _2_n_minus_1 = 2.*n - 1.;
		const double c = _2_n_minus_1 * pi / data.length;

		sum +=	std::sin(c*x)/(_2_n_minus_1)*
				std::exp(- std::pow(c,2.) * data.Gamma*t);
	}

	return 4.*phi0/pi*sum;
}
