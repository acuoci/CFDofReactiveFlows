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
%        method and implicit Euler method. The QUICK (QUadratic           %
%        Upwind Interpolation for Convective Kinetics) discretization     %
%        is adopted for the advective term.                               %
%        The solution of the pentadiagonal system of linear equations     %
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

% ----------------------------------------------------------------------- %
% Data processing
% ----------------------------------------------------------------------- %
h = L/(nx-1);           % step size [m]
Pe = rho*u*h/gamma;     % Peclet number [-]
xi=(h/2):h:(L-h/2);     % grid internal points

fprintf('Peclet number: %f\n', Pe);

% ----------------------------------------------------------------------- %
% Preparin the linear system matrix (penta-diagonal)
% ----------------------------------------------------------------------- %
Fe = rho*u;     Fw = rho*u;
De = gamma/h;   Dw = gamma/h;
if (Fw>=0) alphaw=1; else alphaw=0; end
if (Fe>=0) alphae=1; else alphae=0; end

% Side coefficients
Ae  = - ( De - 3/8*alphae*Fe - 6/8*(1-alphae)*Fe - 1/8*(1-alphaw)*Fw ); 
Aee = -1/8*(1-alphae)*Fe;
Aw  = - ( Dw + 6/8*alphaw*Fw + 1/8*alphae*Fe + 3/8*(1-alphaw)*Fw ); 
Aww = 1/8*alphaw*Fw;

% Central coefficient
Ap =  -(Aww+Aw+Ae+Aee) + (Fe-Fw);

% Linear system
n = nx+3;
b = zeros(n,1);
A = sparse(n, n);

% Assembling for internal points only
for i=3:n-2,   A(i,i-2)=Aww; end
for i=3:n-2,   A(i,i-1)=Aw; end
for i=3:n-2,   A(i,i)=Ap; end
for i=3:n-2,   A(i,i+1)=Ae; end
for i=3:n-2,   A(i,i+2)=Aee; end

% Assembling for boundary points: west side
if (Fw>=0)
    A(1,1)=1; A(1,2)=-6; A(1,3)=-3;
    b(1)=-8*phi0; 
    A(2,2)=1; A(2,3)=1;
    b(2)=2*phi0; 
else
    A(1,1)=1; A(1,2)=-2; A(1,3)=1;
    b(1)=0;
    A(2,2)=3; A(2,3)=6; A(2,4)=-1;
    b(2)=8*phi0; 
end

% Assembling for boundary points: east side
if (Fe>=0)
    A(n-1,n-1)=3; A(n-1,n-2)=6; A(n-1,n-3)=-1;
    b(n-1)=8*phiL;
    A(n,n)=1; A(n,n-1)=-2; A(n,n-2)=1;
    b(n)=0;
else
    A(n-1,n-1)=1;  A(n-1, n-2)=1; 
    b(n-1)=2*phiL;
    A(n,n)=1; A(n,n-1)=-6; A(n,n-2)=-3;
    b(n)=-8*phiL;
end


% Update RHS vector
for i=3:n-2
    b(i) = 0;
end

% Solve the system
phi = A\b;

% ----------------------------------------------------------------------- %
% Data postprocessing
% ----------------------------------------------------------------------- %
x = [ 0, h/2:h:L-h/2, L];
phiphi = [ interp_w(Fw,phi,3), phi(3:nx+1)', interp_e(Fe,phi,n-2) ];
phia = phi0+(phiL-phi0)*(exp(rho*u*x/gamma)-1)./(exp(rho*u*L/gamma)-1);     % analytical solution
plot(x,phiphi,'-', x,phia,'o'); xlabel('x[m]'); ylabel('phi');

% Error estimation
error = norm(phia-phiphi')/nx;
fprintf('Error: %e\n', error);


% ----------------------------------------------------------------------- %
% Interpolation functions
% ----------------------------------------------------------------------- %
function phi_int = interp_w(Fw,phi,i)
    if (Fw>0)
        phi_int = 6/8*phi(i-1)+3/8*phi(i)-1/8*phi(i-2);
    else
        phi_int = 6/8*phi(i)+3/8*phi(i-1)-1/8*phi(i+1);
    end
end

function phi_int = interp_e(Fe,phi,i)
    if (Fe>0)
        phi_int = 6/8*phi(i)+3/8*phi(i+1)-1/8*phi(i-1);
    else
        phi_int = 6/8*phi(i+1)+3/8*phi(i)-1/8*phi(i+2);
    end
end
