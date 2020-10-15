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
%   Exact solutions from K. Masatsuka, "I do like CFD, Vol. 1", 2018      %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% User-defined data
%-------------------------------------------------------------------------%
nx=40;                  % number of grid points along x
ny=40;                  % number of grid points along y
nstep=150;              % number of time steps
lengthx=1.0;            % domain length along x [m]
lengthy=1.0;            % domain length along y [m]
D=0.050;                % diffusion coefficient [m2/s]
u=1;                    % velocity along x [m/s]
v=0.5;                  % velocity along y [m/s]
A = 1.25;               % arbitrary constant
B = 0.75;               % arbitrary constant
C = 1;                  % arbitrary constant


% Pre-processing of user-defined data
%-------------------------------------------------------------------------%
% Calculate grid steps
hx=lengthx/(nx-1);      % grid step along x [m]
hy=lengthy/(ny-1);      % grid step along y [m] 
x=0:hx:lengthx;         % x coordinates [m]
y=0:hy:lengthy;         % y coordinates [m]

% Numerical setup: time step (stability conditions)
sigma = 0.75;                       % safety coefficient
dt_diff  = 1/4*min(hx^2, hy^2)/D;   % diffusion [s]
dt_conv = 4*D/(u^2+v^2);            % convection [s]
dt = sigma*min(dt_diff, dt_conv);   % time step [s]
fprintf('Co=%f Di=%f Pe=%f \n', ...
    max(u,v)*dt/min(hx,hy), D*dt/min(hx^2,hy^2), max(hx,hy)*max(u,v)/D);

% Memory allocation
f=zeros(nx,ny);     % current numerical solution
f0=zeros(nx,ny);    % initial numerical solution
fan=zeros(nx,ny);   % analytical solution

% Initial condition
for i=1:nx
        for j=1:ny
            csi = u*x(i)+v*y(j);
            eta = v*x(i)-u*y(j);
            f0(i,j) = C*cos(B*pi*eta)*exp((1-sqrt(1+4*(A*pi*D)^2))/(2*D)*csi);
        end
end

% Advancing in time
%-------------------------------------------------------------------------%
t = 0.;
f = f0;
for m=1:nstep
    
    % Constant
    gamma = D*pi*pi*(u^2+v^2)*(B^2-A^2);
    
    % Analytical solution
    for i=1:nx
        for j=1:ny
            fan(i,j) = f0(i,j) * exp(-gamma*t);
        end
    end
    
    % Error (mean) between numerical and analytical solution
    error = sum(sum(abs(f-fan)))/(nx*ny);
    
    % Plot the current solution
    plot(x,f(:,ny/2), x,fan(:,ny/2));
    title('Solution along the horizontal axis at y=Ly/2');
    legend('numerical', 'analytical');
    xlabel('x-axis [m]'); ylabel('f value');
    pause(0.0001)
    
    % Boundary conditions (Dirichlet, constant in time)
    f(:,1)  = f0(:,1)*exp(-gamma*t);
    f(:,ny) = f0(:,ny)*exp(-gamma*t);
    f(1,:)  = f0(1,:)*exp(-gamma*t);
    f(nx,:) = f0(nx,:)*exp(-gamma*t);

    % Forward Euler method
    fo=f;
    for i=2:nx-1
        for j=2:ny-1
            f(i,j) = fo(i,j)...
                    -(0.5*dt*u/hx)*(fo(i+1,j)-fo(i-1,j))...
                    -(0.5*dt*v/hy)*(fo(i,j+1)-fo(i,j-1))...
                    +(D*dt/hx^2)*(fo(i+1,j)-2*fo(i,j)+fo(i-1,j))...
                    +(D*dt/hy^2)*(fo(i,j+1)-2*fo(i,j)+fo(i,j-1));
        end
    end   
    
    % New time step
    t=t+dt;
    
end

% Solution
figure;
pcolor(x, y, f'); 
colorbar; shading interp; colormap(jet); 
title('numerical solution');
xlabel('x-axis [m]'); ylabel('y-axis [m]');


% Difference
figure;
abserr = abs(f-fan);
pcolor(x, y, abserr'); 
colorbar; shading interp; colormap(jet); 
title('absolute error');
xlabel('x-axis [m]'); ylabel('y-axis [m]');
