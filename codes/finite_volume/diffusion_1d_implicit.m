% ----------------------------------------------------------------------- %
%                                     __  __  __       _  __   __         %
%        |\/|  _  |_ |  _  |_   |__| /   |_  |  \  _  (_ |__) |_          %
%        |  | (_| |_ | (_| |_)     | \__ |   |__/ (_) |  | \  |           %
%                                                                         %
% ----------------------------------------------------------------------- %
%                                                                         %
%   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       %
%   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            %
%   Department of Chemistry, Materials and Chemical Engineering           %
%   Politecnico di Milano                                                 %
%   P.zza Leonardo da Vinci 32, 20133 Milano                              %
%                                                                         %
% ----------------------------------------------------------------------- %
%                                                                         %
%   This file is part of Matlab4CFDofRF framework.                        %
%                                                                         %
%   License                                                               %
%                                                                         %
%   Copyright(C) 2019 Alberto Cuoci                                       %
%   Matlab4CFDofRF is free software: you can redistribute it and/or       %
%   modify it under the terms of the GNU General Public License as        %
%   published by the Free Software Foundation, either version 3 of the    %
%   License, or (at your option) any later version.                       %
%                                                                         %
%   Matlab4CFDofRF is distributed in the hope that it will be useful,     %
%   but WITHOUT ANY WARRANTY; without even the implied warranty of        %
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         %
%   GNU General Public License for more details.                          %
%                                                                         %
%   You should have received a copy of the GNU General Public License     %
%   along with Matlab4CRE. If not, see <http://www.gnu.org/licenses/>.    %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
%  The meaning of this code is to give a simple introduction to the       %
%  finite volume (FV) technique for discretizing transport equations      %
%  in space.                                                              %
%                                                                         %
%  Code: 1D diffusion equation with finite volume method and implicit     %
%        Euler method. A constant, uniform source term q is included.     %
%        The solution of the tridiagonal system of linear equations       %
%        is obtained by the standard solver available in MATLAB           %
%                                                                         %
%        rho*cp*dT/dt = k*d2T/dx2 + q                                     %
%        T(x=0)=T0, T(x=L)=TL, T(t=0)=Ti                                  %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% ----------------------------------------------------------------------- %
% User data
% ----------------------------------------------------------------------- %

rho = 1000.;    % density [kg/m3]
cp  = 1000.;    % specific heat [J/kg/K]
kappa = 0.5;    % thermal conductivity [W/m/K]
L = 0.02;       % length of computational domain [m]
tau = 1000;     % total time of simulation [s]
T0 = 100;       % temperature on the left side [C]
TL = 200;       % temperature on the right side [C]
Ti = 100;       % initial temperature [C]
q = 1e6;        % heat generation [W/m3]

% Numerical parameters
nx = 100;       % number of points [-]
dt = 1;         % time step [s]


% ----------------------------------------------------------------------- %
% Data processing
% ----------------------------------------------------------------------- %
h = L/(nx-1);           % step size [m]
alpha = kappa/rho/cp;   % thermal diffusivity [m2/s]
dt_diff=h^2/2/alpha;    % time step (diffusion stability) [s]
nsteps=tau/dt;          % number of time steps [-]
tc = L^2/alpha;         % characteristic diffusion time [s]
xi=(h/2):h:(L-h/2);     % grid internal points

fprintf('Max Time step [s]:       %f\n', dt_diff);
fprintf('Time step [s]:           %f\n', dt);
fprintf('Diffusivity [m2/s]:      %f\n', alpha);
fprintf('Characteristic time [s]: %f\n', tc);


% ----------------------------------------------------------------------- %
% Solution
% ----------------------------------------------------------------------- %
T = zeros(nx+1,1) + Ti;

% Matrix coefficients
Ap = 1/dt+2*alpha/h^2;
Ae = -alpha/h^2;
Aw = -alpha/h^2;

% Assembling tridiagonal matrix
n = nx+1;
b = zeros(n,1);
A = sparse(n, n);
A(1,1)=1; A(1,2)=1;
for i=2:n-1,   A(i,i-1)=Aw; end
for i=2:n-1,   A(i,i)=Ap;   end
for i=2:n-1,   A(i,i+1)=Ae; end
A(n,n)=1; A(n,n-1)=1;

% Advancing in time
for j=1:nsteps
    
    % Update boundary conditions (ghost points)
    b(1) = 2*T0; 
    b(n) = 2*TL;
    
    % Update RHS vector
    for i=2:n-1
        b(i) = 1/dt*T(i) + q/(rho*cp);
    end
    
    % Solve the system
    T = A\b;

    % On-the-fly post processing
    if (mod(j,25)==1)
        plot(xi,T(2:nx));
        drawnow;
    end
    
end

% ----------------------------------------------------------------------- %
% Data postprocessing
% ----------------------------------------------------------------------- %
x = 0:h:L;
TT = (T(1:nx)+T(2:nx+1))/2;
Ta = ((TL-T0)/L + q/2/kappa*(L-x)).*x + T0;     % analytical solution
plot(x,TT,'-', x,Ta,'o'); xlabel('x[m]'); ylabel('temperature [C]');

% Error estimation
error = norm(Ta-TT')/nx;
fprintf('Error: %e\n', error);
