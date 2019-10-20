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
%  Code: 1D steady-state advection-diffusion equation with finite volume  %
%        method and implicit Euler method                                 %
%        The solution of the tridiagonal system of linear equations       %
%        is obtained by the standard solver available in MATLAB           %
%                                                                         %
%        d/dx(rho*u*phi) = gamma*d2phi/dx2                                %
%        phi(x=0)=1, phi(x=L)=0                                           %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% ----------------------------------------------------------------------- %
% User data
% ----------------------------------------------------------------------- %

rho = 1.;           % density [kg/m3]
L = 1.0;            % length of computational domain [m]
gamma = 0.1;        % mass diffusion coefficient [kg/m/s]
u = 1;              % velocity [m/s]
phi0 = 1;           % left value of phi
phiL = 0;           % right value of phi

% Numerical parameters
nx = 20;            % number of points [-]
adv_scheme = 'CD';  % advection scheme: CD, UPWIND, HYBRID, POWER-LAW

% ----------------------------------------------------------------------- %
% Data processing
% ----------------------------------------------------------------------- %
h = L/(nx-1);           % step size [m]
Pe = rho*u*h/gamma;     % Peclet number [-]
xi=(h/2):h:(L-h/2);     % grid internal points

fprintf('Peclet number: %f\n', Pe);

% ----------------------------------------------------------------------- %
% Preparing system matrix
% ----------------------------------------------------------------------- %

Fe = rho*u;     Fw = rho*u;
De = gamma/h;   Dw = gamma/h;

if (strcmp(adv_scheme,'CD')) % Central differencing scheme
    Ae =  Fe/2 - De;
    Aw = -Fw/2 - Dw;
elseif (strcmp(adv_scheme,'UPWIND')) % Upwind scheme 
    Ae = -max(-Fe,0) -De;
    Aw = -max( Fw,0) -Dw;
elseif (strcmp(adv_scheme,'HYBRID')) % Hybrid
    Ae = -max( [-Fe, -Fe/2+De, 0] );
    Aw = -max( [ Fw,  Fw/2+Dw, 0] );
elseif (strcmp(adv_scheme,'POWER-LAW')) % Power-Law
    Ae = -( max(-Fe,0) + De*max(0,(1-0.1*abs(Pe))^5) );
    Aw = -( max( Fw,0) + Dw*max(0,(1-0.1*abs(Pe))^5) );
end

% Central coefficient
Ap =  -(Aw+Ae) + (Fe-Fw);

% Linear system
n = nx+1;
b = zeros(n,1);
A = sparse(n, n);
A(1,1)=1; A(1,2)=1;
for i=2:n-1,   A(i,i-1)=Aw; end
for i=2:n-1,   A(i,i)=Ap;   end
for i=2:n-1,   A(i,i+1)=Ae; end
A(n,n)=1; A(n,n-1)=1;

% ----------------------------------------------------------------------- %
% Solution
% ----------------------------------------------------------------------- %

% Boundary conditions (ghost points)
b(1) = 2*phi0; 
b(n) = 2*phiL;

% Update RHS vector
for i=2:n-1
    b(i) = 0;
end

% Solve the system (no iterations are needed)
phi = A\b;

% ----------------------------------------------------------------------- %
% Data postprocessing
% ----------------------------------------------------------------------- %
x = 0:h:L;
phiphi = (phi(1:nx)+phi(2:nx+1))/2;
phia = phi0+(phiL-phi0)*(exp(rho*u*x/gamma)-1)./(exp(rho*u*L/gamma)-1);     % analytical solution
plot(x,phiphi,'-', x,phia,'o'); xlabel('x[m]'); ylabel('phi');

% Error estimation
error = norm(phia-phiphi')/nx;
fprintf('Error: %e\n', error);