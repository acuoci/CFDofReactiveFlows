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
%   Copyright(C) 2020 Alberto Cuoci                                       %
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
%  Code: 2D driven-cavity problem based on a staggered grid               %
%        The code solves the RANS equations for turbulent conditions      %
%        Turbulence is modelled via the Prandtl k-l model*, requiring the %
%        solution of an additional transport equation for the turbulent   %
%        kinetic energy                                                   %
%      * https://www.cfd-online.com/Wiki/Prandtl%27s_one-equation_model   %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% ----------------------------------------------------------------------- %
% User data
% ----------------------------------------------------------------------- %

% Only even numbers of cells are acceptable
nx=16;              % number of (physical) cells along x
ny=nx;              % number of (physical) cells along y
L=1;                % length [m]
nu=0.01;            % kinematic viscosity [m2/s] 
tau=20;             % total time of simulation [s]
kappa0=0;           % initial value of turbulent kinetic energy [m2/s2]

% Boundary conditions (velocities)
un=10;      % north wall velocity [m/s]
us=0;       % south wall velocity [m/s]
ve=0;       % east wall velocity [m/s]
vw=0;       % west wall velocity [m/s]

% Boundary conditions for turbulent kinetic energy
kappan=0.;       % north wall turbulent kinetic energy [m2/s2]
kappas=0.;       % south wall turbulent kinetic energy [m2/s2]
kappae=0.;       % east wall turbulent kinetic energy [m2/s2]
kappaw=0.;       % west wall turbulent kinetic energy [m2/s2]

% Parameters for SOR
max_iterations=10000;   % maximum number of iterations
beta=1.5;               % SOR coefficient
max_error=1e-3;         % error for convergence

% ----------------------------------------------------------------------- %
% Data processing
% ----------------------------------------------------------------------- %
if (mod(nx,2)~=0 || mod(ny,2)~=0)
    error('Only even number of cells can be accepted (for graphical purposes only)');
end

% Grid step
h=L/nx;             % grid step [m]

% Time step
sigma = 0.5;                        % safety factor for time step (stability)
u2=(un^2+us^2+ve^2+vw^2);           % velocity measure [m2/s2]
dt_diff=h^2/4/nu;                   % time step (diffusion stability)
dt_conv=4*nu/u2;                    % time step (convection stability)
dt=sigma*min(dt_diff, dt_conv);     % time step (stability)
nsteps=tau/dt;                      % number of steps
Re = un*L/nu;                       % Reynolds' number

fprintf('Time step: %f\n', dt);
fprintf(' - Diffusion:  %f\n', dt_diff);
fprintf(' - Convection: %f\n', dt_conv);
fprintf('Reynolds number: %f\n', Re);

% Grid construction
x=0:h:L;                         % grid coordinates (x axis)
y=0:h:L;                         % grid coordinates (y axis)
[X,Y] = meshgrid(x,y);           % MATLAB grid

% ----------------------------------------------------------------------- %
% Memory allocation
% ----------------------------------------------------------------------- %

% Main fields (velocities and pressure)
u=zeros(nx+1,ny+2);
v=zeros(nx+2,ny+1);
p=zeros(nx+2,ny+2);
nut=zeros(nx+2,nx+2)+nu;
kappa=zeros(nx+2,ny+2)+kappa0;

% Temporary velocity fields
ut=zeros(nx+1,ny+2);
vt=zeros(nx+2,ny+1);

% Fields used only for graphical post-processing purposes
uu=zeros(nx+1,ny+1);
vv=zeros(nx+1,ny+1);
kk=zeros(nx+1,ny+1);

% Coefficient for pressure equation
c=zeros(nx+2,ny+2)+1/4;
c(2,3:ny)=1/3;c(nx+1,3:ny)=1/3;c(3:nx,2)=1/3;c(3:nx,ny+1)=1/3;
c(2,2)=1/2;c(2,ny+1)=1/2;c(nx+1,2)=1/2;c(nx+1,ny+1)=1/2;

% ----------------------------------------------------------------------- %
% Solution over time
% ----------------------------------------------------------------------- %
t=0.0;
for is=1:nsteps
    
    % Update the turbulent viscosity (Prandtl k-l model)
    Cmu = 0.09;
    nut = zeros(nx+2,ny+2);
    for i=2:nx+1
        for j=2:ny+1
            dx = (i-2)*h + h/2;
            dy = (j-2)*h + h/2;
            lw = min(dx,dy);
            nut(i,j) = Cmu*sqrt( max(kappa(i,j),0.) )*lw;
        end
    end
                
    % Boundary conditions (velocity)
    u(1:nx+1,1)=2*us-u(1:nx+1,2);
    u(1:nx+1,ny+2)=2*un-u(1:nx+1,ny+1);
    v(1,1:ny+1)=2*vw-v(2,1:ny+1);
    v(nx+2,1:ny+1)=2*ve-v(nx+1,1:ny+1);
    
    % Advection-diffusion equation (predictor)
    [ut, vt] = AdvectionDiffusion2D( ut, vt, u, v, nx, ny, h, dt, nu, nut);
    
    % Pressure equation (Poisson)
    [p, iter] = Poisson2D( p, ut, vt, c, nx, ny, h, dt, ...
                           beta, max_iterations, max_error);
    
    % Correct the velocity
    u(2:nx,2:ny+1)=ut(2:nx,2:ny+1)-(dt/h)*(p(3:nx+1,2:ny+1)-p(2:nx,2:ny+1));
    v(2:nx+1,2:ny)=vt(2:nx+1,2:ny)-(dt/h)*(p(2:nx+1,3:ny+1)-p(2:nx+1,2:ny));
    
    % Print time step on the screen
    if (mod(is,50)==1)
        nut_mean = mean(mean(nut(2:nx+1, 2:ny+1)));
        fprintf('Step: %d - Time: %f - Poisson iterations: %d - Mean nut: %f \n', ...
                is, t, iter, nut_mean);
    end
    
    % ----------------------------------------------------------------------- %
    % Turbulent kinetic energy: scalar advection-diffusion equation
    % ----------------------------------------------------------------------- %

    % Boundary conditions (Turbulent kinetic energy)
    kappa(2:nx+1,1)=2*kappas-kappa(2:nx+1,2);
    kappa(2:nx+1,ny+2)=2*kappan-kappa(2:nx+1,ny+1);
    kappa(1,2:ny+1)=2*kappaw-kappa(2,2:ny+1);
    kappa(nx+2,2:ny+1)=2*kappae-kappa(nx+1,2:ny+1);
    
    % Update passive scalar solution
    [kappa] = AdvectionDiffusion2DTurbulentKineticEnergy( kappa, u, v, nx, ny, h, dt, nu, nut);
  
    % Advance time
    t=t+dt;
 
    % ----------------------------------------------------------------------- %
    % Update graphical output
    % ----------------------------------------------------------------------- %
    if (mod(is,50)==1)
        
        % Post-processing (reconstructi on the corners of pressure cells)
        uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
        vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
        kk(1:nx+1,1:ny+1)=0.25*(kappa(1:nx+1,1:ny+1)+kappa(2:nx+2,1:ny+1) + ...
                                    kappa(1:nx+1,2:ny+2)+kappa(2:nx+2,2:ny+2)); 
        nutnut(1:nx+1,1:ny+1)=0.25*(nut(1:nx+1,1:ny+1)+nut(2:nx+2,1:ny+1) + ...
                                    nut(1:nx+1,2:ny+2)+nut(2:nx+2,2:ny+2));   
           
        % Contour: x-velocity
        subplot(231);
        contour(X,Y,uu');
        axis('square'); xlabel('x [m]');ylabel('y [m]'); title('x-velocity contour');

        % Contour: y-velocity
        subplot(232);
        contour(X,Y,vv');
        axis('square'); xlabel('x [m]');ylabel('y [m]'); title('y-velocity contour');

        % Plot: velocity components along the horizontal middle axis
        subplot(233);
        plot(x,uu(:, round(ny/2)));
        hold on;
        plot(x,vv(:, round(ny/2)));
        axis('square');
        xlabel('x [m]');ylabel('velocity [m/s]'); title('velocity along x-axis');
        legend('x-velocity', 'y-velocity');
        hold off;

        % Contour: turbulent kinetic energy
        subplot(234);
        contour(X,Y,kk');
        axis('square'); xlabel('x [m]');ylabel('y [m]'); title('turbulent kinetic energy');
        
        % Contour: passive scalar
        subplot(235);
        contour(X,Y,nutnut');
        axis('square'); xlabel('x [m]');ylabel('y [m]'); title('turbulent viscosity');

        % Plot: velocity components along the horizontal middle axis
        subplot(236);
        plot(x,kk(:, round(ny/2)));
        hold on;
        plot(y,kk(round(ny/2),:));
        axis('square');
        xlabel('x or y [m]');ylabel('turbulent kinetic energy'); title('turbulent kinetic energy (centerlines)');
        legend('along x', 'along y'); 
        hold off;

        pause(0.001)
    
    end
    
end


% --------------------------------------------------------------------------------------
% Poisson equation solver
% --------------------------------------------------------------------------------------
function [p, iter] = Poisson2D( p, ut, vt, gamma, nx, ny, h, dt, ...
                                beta, max_iterations, max_error)

    for iter=1:max_iterations
        
        for i=2:nx+1
            for j=2:ny+1
                
                delta = p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1);
                S = (h/dt)*(ut(i,j)-ut(i-1,j)+vt(i,j)-vt(i,j-1));
                p(i,j)=beta*gamma(i,j)*( delta-S )+(1-beta)*p(i,j);
                
            end
        end
        
        % Estimate the error
        epsilon=0.0; 
        for i=2:nx+1
            for j=2:ny+1
                delta = p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1);
                S = (h/dt)*(ut(i,j)-ut(i-1,j)+vt(i,j)-vt(i,j-1));              
                epsilon=epsilon+abs( p(i,j) - gamma(i,j)*( delta-S ) );
            end
        end
        epsilon = epsilon / (nx*ny);
        
        % Check the error
        if (epsilon <= max_error) % stop if converged
            break;
        end 
        
    end

end

% --------------------------------------------------------------------------------------
% Advection-diffusion equation
% --------------------------------------------------------------------------------------
function [ut, vt] = AdvectionDiffusion2D( ut, vt, u, v, nx, ny, h, dt, nu, nut)
                            
    % Temporary u-velocity
    for i=2:nx
        for j=2:ny+1 
            
            nutot = nu + nut(i,j);
            
            ue2 = 0.25*( u(i+1,j)+u(i,j) )^2;
            uw2 = 0.25*( u(i,j)+u(i-1,j) )^2;
            unv = 0.25*( u(i,j+1)+u(i,j) )*( v(i+1,j)+v(i,j) );
            usv = 0.25*( u(i,j)+u(i,j-1) )*( v(i+1,j-1)+v(i,j-1) );
            
            A = (ue2-uw2+unv-usv)/h;
            D = (nutot/h^2)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j));
            
            ut(i,j)=u(i,j)+dt*(-A+D);
            
        end
    end
    
    % Temporary v-velocity
    for i=2:nx+1
        for j=2:ny 
            
            nutot = nu + nut(i,j);
            
            vn2 = 0.25*( v(i,j+1)+v(i,j) )^2;
            vs2 = 0.25*( v(i,j)+v(i,j-1) )^2;
            veu = 0.25*( u(i,j+1)+u(i,j) )*( v(i+1,j)+v(i,j) );
            vwu = 0.25*( u(i-1,j+1)+u(i-1,j) )*( v(i,j)+v(i-1,j) );
            A = (vn2 - vs2 + veu - vwu)/h;
            D = (nutot/h^2)*(v(i+1,j)+v(i-1,j)+v(i,j+1)+v(i,j-1)-4*v(i,j));
            
            vt(i,j)=v(i,j)+dt*(-A+D);
            
        end
    end
    
end

% --------------------------------------------------------------------------------------
% Advection-diffusion equation for turbulent kinetic energy
% --------------------------------------------------------------------------------------
function [kappa] = AdvectionDiffusion2DTurbulentKineticEnergy( kappa, u, v, nx, ny, h, dt, nu, nut)

    % Prandtl number for k
    sigmak = 1.0;
    
    kappao = kappa;
    for i=2:nx+1
            for j=2:ny+1
                
                % Distance from the wall
                dx = (i-2)*h + h/2;
                dy = (j-2)*h + h/2;
                lw = min(dx,dy);
                
                % Total viscosity
                nutot = nu + nut(i,j)/sigmak;
                nutote = 0.5*( nutot + ( nu + nut(i+1,j)/sigmak) );
                nutotw = 0.5*( nutot + ( nu + nut(i-1,j)/sigmak) );
                nutotn = 0.5*( nutot + ( nu + nut(i,j+1)/sigmak) );
                nutots = 0.5*( nutot + ( nu + nut(i,j-1)/sigmak) );
                
                % Velocity components on the cell faces
                ue = u(i,j);
                vn = v(i,j);
                uw = u(i-1,j);
                vs = v(i,j-1);
                
                uen = (u(i,j)+u(i,j+1))/2;
                ues = (u(i,j)+u(i,j-1))/2;
                uwn = (u(i-1,j)+u(i-1,j+1))/2;
                uws = (u(i-1,j)+u(i-1,j-1))/2;
               
                vne = (v(i,j)+v(i+1,j))/2;
                vse = (v(i,j-1)+v(i+1,j-1))/2;  
                vnw = (v(i,j)+v(i-1,j))/2;
                vsw = (v(i,j-1)+v(i-1,j-1))/2;
                
                un = ( uen + uwn ) /2;
                us = ( ues + uws ) /2;
                ve = ( vne + vse ) /2;
                vw = ( vnw + vsw ) /2;
                
                % Scalar on the faces (linear interpolation)
                kappae = 0.50*(kappao(i+1,j)+kappao(i,j));
                kappaw = 0.50*(kappao(i-1,j)+kappao(i,j));
                kappan = 0.50*(kappao(i,j+1)+kappao(i,j));
                kappas = 0.50*(kappao(i,j-1)+kappao(i,j));
                
                % Gradients of scalar on the faces (centered)
                dkappa_dx_e = (kappao(i+1,j)-kappao(i,j))/h;
                dkappa_dx_w = (kappao(i,j)-kappao(i-1,j))/h;
                dkappa_dy_n = (kappao(i,j+1)-kappao(i,j))/h;
                dkappa_dy_s = (kappao(i,j)-kappao(i,j-1))/h;
                
                % Convection and diffusion contributions
                convection = ue*kappae*h - uw*kappaw*h + ...
                             vn*kappan*h - vs*kappas*h;
                         
                diffusion = nutote*dkappa_dx_e*h -nutotw*dkappa_dx_w*h + ...
                            nutotn*dkappa_dy_n*h -nutots*dkappa_dy_s*h ;
                                 
                % Production term: P=nut*S^2
                %                  S = sqrt(2*Sij*Sij)
                %                  Sij = 1/2*(dvi/dxj+dvj/dxi)
                du_over_dx = (ue-uw)/h;
                dv_over_dy = (vn-vs)/h;
                du_over_dy = (un-us)/h;
                dv_over_dx = (ve-vw)/h;
                Sxx = 0.5*(du_over_dx+du_over_dx);
                Sxy = 0.5*(du_over_dy+dv_over_dx);
                Syx = Sxy;
                Syy = 0.5*(dv_over_dy+dv_over_dy);
                S2 = 2*(Sxx^2+2*Sxy*Syx+Syy^2);
                P = nutot*S2;
                
                % Dissipation term: D=k^(3/2)/lw
                D = ( max(kappao(i,j),0)^1.5 )/lw;
                 
                % Euler method
                kappa(i,j)= kappao(i,j) + dt/h^2 *( -convection + diffusion ) + dt*(P-D);
                
            end
    end

end
