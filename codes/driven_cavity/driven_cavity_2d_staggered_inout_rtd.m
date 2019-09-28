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
%  Code: 2D driven-cavity problem on a staggered grid with inclusion of   %
%        inlet and outlet sections and a passive tracer equation          %
%        The Residence Time Distribution (RTD) and the Cumulative         %
%        Distribution Function (CDF) are automatically calculated         %
%        The code is adapted and extended from Tryggvason, Computational  %
%        Fluid Dynamics http://www.nd.edu/~gtryggva/CFD-Course/           %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% ----------------------------------------------------------------------- %
% User data
% ----------------------------------------------------------------------- %

% Only even numbers of cells are acceptable
tau=5;       % total time of simulation [s]
alpha=0.005; % diffusion coefficient passive tracer [m2/s]
phi0=0;      % [TRACER] initial value of passive tracer

% Boundary conditions for passive tracer
% We assume that the passive tracer enters the domain through the inlet
% section with a fixed, user-defined value
phi_inlet=1;  % [TRACER] inlet value

% ----------------------------------------------------------------------- %
% Reading solution (from a previous steady-state calculation)
% ----------------------------------------------------------------------- %
[L, nx, ny, u, v, p, indices] = ReadSolution('solution.out');

% [INOUT] Inlet section (west side)
nin_start = indices(1);         % first cell index 
nin_end = indices(2);           % last cell index

% [INOUT] Outlet section (east side)
nout_start = indices(3);        % first cell index 
nout_end = indices(4);          % last cell index

% Process the grid
h=L/nx;                           % grid step (uniform grid) [m]

% [INOUT] Inlet/Outlet section areas
Ain = h*(nin_end-nin_start+1);      % inlet section area [m]
Aout = h*(nout_end-nout_start+1);   % outlet section area [m]

% Time step
umax = max([ max(max(abs(u))) max(max(abs(v)))]);
sigma = 0.50;                       % safety factor for time step (stability)
dt_diff=h^2/4/alpha;                % time step (diffusion stability) [s]
dt_conv=4*alpha/umax^2;             % time step (convection stability) [s]
dt=sigma*min(dt_diff, dt_conv);     % time step (stability) [s]
nsteps=tau/dt;                      % number of steps

fprintf('Time step: %f\n', dt);
fprintf(' - Diffusion:  %f\n', dt_diff);
fprintf(' - Convection: %f\n', dt_conv);

% Grid construction
x=0:h:L;                         % grid coordinates (x axis)
y=0:h:L;                         % grid coordinates (y axis)
[X,Y] = meshgrid(x,y);           % MATLAB grid

% ----------------------------------------------------------------------- %
% Memory allocation
% ----------------------------------------------------------------------- %

% Main fields (velocities and pressure)
phi=zeros(nx+2,ny+2)+phi0;  % [TRACER] passive tracer field

% Fields used only for graphical post-processing purposes
phiphi=zeros(nx+1,ny+1);    % [TRACER] passive tracer field

% ----------------------------------------------------------------------- %
% Solution over time
% ----------------------------------------------------------------------- %
t=0.0;
phiout_old = 0;     % [TRACER] passive tracer value at outlet section
Area = 0.;          % [TRACER] area below the RTD curve
for is=1:nsteps

    % ----------------------------------------------------------------------- %
    % [TRACER] passive tracer advection-diffusion equation
    % ----------------------------------------------------------------------- %

    % [TRACER] Boundary conditions (zero-gradient)
    phi(2:nx+1,1)=phi(2:nx+1,2);
    phi(2:nx+1,ny+2)=phi(2:nx+1,ny+1);
    phi(1,2:ny+1)=phi(2,2:ny+1);
    phi(nx+2,2:ny+1)=phi(nx+1,2:ny+1);
    
    % [TRACER] Boundary conditions at inlet section (Dirichlet)
    phi(1,nin_start:nin_end)=2*phi_inlet-phi(2,nin_start:nin_end);

    % [TRACER] Update passive tracer solution (tracer behaves like a passive scalar)
    phi = AdvectionDiffusion2DPassiveScalar( phi, u, v, nx, ny, h, dt, alpha );
    
    % Advance in time
    t=t+dt;
    
    % ----------------------------------------------------------------------- %
    % On the fly post-processing                                                   %
    % ----------------------------------------------------------------------- %
    if (mod(is,50)==1)
        
        % [TRACER] Adjust values on the corners (graphical purposes only)
        phi(1,1)=1/3*(phi(1,2)+phi(2,2)+phi(2,1));
        phi(nx+2,1)=1/3*(phi(nx+2,2)+phi(nx+1,1)+phi(nx+1,nx+1));
        phi(nx+2,1)=1/3*(phi(nx+2,2)+phi(nx+1,1)+phi(nx+1,2));
        phi(1,ny+2)=1/3*(phi(2,ny+1)+phi(1,ny+1)+phi(2,ny+1));
        phi(nx+2,ny+2)=1/3*(phi(nx+1,ny+2)+phi(nx+2,ny+1)+phi(nx+1,ny+1));

        % Field reconstruction
        phiphi(1:nx+1,1:ny+1)=0.25*(phi(1:nx+1,1:ny+1)+phi(2:nx+2,1:ny+1) + ...
                                    phi(1:nx+1,2:ny+2)+phi(2:nx+2,2:ny+2));

        % Outlet value
        phiout = mean(phi(nx+1,nout_start:nout_end));
        E = (phiout - phiout_old)/dt;
        Area = Area + E*dt;
        phiout_old = phiout;
        
        % Print on the screen
        fprintf( 'Step: %d - Time: %f - PhiOut: %f - Area: %f\n', ...
                 is, t, phiout, Area);
        
        % Surface: passive tracer
        subplot(221);
        surface(X,Y,phiphi');
        axis('square'); xlabel('x [m]');ylabel('y [m]'); title('passive tracer');
        
        % Plot CDF
        subplot(222);
        scatter(t,phiout, '.');
        axis('square'); xlabel('time [s]');ylabel('F'); title('CDF');
        ylim([0 1]); xlim([0 tau]);
        hold on;
        
        % Plot RTD
        subplot(223);
        scatter(t,E, '.');
        axis('square'); xlabel('time [s]');ylabel('E[1/s]'); title('RTD');
        ylim([0 Inf]); xlim([0 tau]);
        hold on;
    
        pause(0.01);
        
    end
 
end


% --------------------------------------------------------------------------------------
% Advection-diffusion equation for passive tracer
% --------------------------------------------------------------------------------------
function [phi] = AdvectionDiffusion2DPassiveScalar( phi, u, v, nx, ny, h, dt, alpha )
    
    phio = phi;
    for i=2:nx+1
            for j=2:ny+1
                
                % Velocity components on the cell faces
                ue = u(i,j);
                vn = v(i,j);
                uw = u(i-1,j);
                vs = v(i,j-1);
                
                % passive tracer on the faces (linear interpolation)
                phie = 0.50*(phio(i+1,j)+phio(i,j));
                phiw = 0.50*(phio(i-1,j)+phio(i,j));
                phin = 0.50*(phio(i,j+1)+phio(i,j));
                phis = 0.50*(phio(i,j-1)+phio(i,j));
                
                % Gradients of passive tracer on the faces (centered)
                dphi_dx_e = (phio(i+1,j)-phio(i,j))/h;
                dphi_dx_w = (phio(i,j)-phio(i-1,j))/h;
                dphi_dy_n = (phio(i,j+1)-phio(i,j))/h;
                dphi_dy_s = (phio(i,j)-phio(i,j-1))/h;
                
                % Convection and diffusion contributions
                convection = ue*phie*h - uw*phiw*h + ...
                             vn*phin*h - vs*phis*h;
                diffusion = alpha* ( dphi_dx_e*h -dphi_dx_w*h + ...
                                     dphi_dy_n*h -dphi_dy_s*h );
                 
                % Euler method
                phi(i,j)= phio(i,j) + dt/h^2 *( -convection + diffusion );
                
            end
    end

end

% --------------------------------------------------------------------------------------
% Read solution from a file
% --------------------------------------------------------------------------------------
function [L, nx, ny, u, v, p, indices] = ReadSolution(name)

    fsol = fopen(name, 'r');

    L = fscanf(fsol, '%f',1);
    nx = fscanf(fsol, '%d',1);
    
    indices(1) = fscanf(fsol, '%d',1);
    indices(2) = fscanf(fsol, '%d',1);
    indices(3) = fscanf(fsol, '%d',1);
    indices(4) = fscanf(fsol, '%d',1);
    
    % Memory allocation
    ny = nx;
    u = zeros(nx+1, ny+2);
    v = zeros(nx+2, ny+1);
    p = zeros(nx+2, ny+2);
   
    % Read fields
    for i=1:nx+1
        for j=1:ny+2
            u(i,j) = fscanf(fsol, '%e', 1);
        end
    end
    for i=1:nx+2
        for j=1:ny+1
            v(i,j) = fscanf(fsol, '%e', 1);
        end
    end
    for i=1:nx+2
        for j=1:ny+2
            p(i,j) = fscanf(fsol, '%e', 1);
        end
    end

    fclose(fsol);
    
end
