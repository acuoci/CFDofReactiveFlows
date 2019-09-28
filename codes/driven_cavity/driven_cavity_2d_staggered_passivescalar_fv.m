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
%  Code: 2D driven-cavity problem based on a staggered grid               %
%        with inclusion of a passive scalar equation (FV)                 %
%        The passive scalar equation is solved using the Finite           %
%        Volume Method                                                    %
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
nx=16;              % number of (physical) cells along x
ny=nx;              % number of (physical) cells along y
L=1;                % length [m]
nu=0.01;            % kinematic viscosity [m2/s] 
tau=20;             % total time of simulation [s]
alpha=0.01;         % diffusion coefficient passive scalar [m2/s]
phi0=0;             % initial value of passive scalar

% Boundary conditions (velocities)
un=1;       % north wall velocity [m/s]
us=0;       % south wall velocity [m/s]
ve=0;       % east wall velocity [m/s]
vw=0;       % west wall velocity [m/s]

% Boundary conditions for passive scalar
phin=0;       % north wall passive scalar
phis=1;       % south wall passive scalar
phie=0;       % east wall passive scalar
phiw=0;       % west wall passive scalar

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
phi=zeros(nx+2,ny+2)+phi0;

% Temporary velocity fields
ut=zeros(nx+1,ny+2);
vt=zeros(nx+2,ny+1);

% Fields used only for graphical post-processing purposes
uu=zeros(nx+1,ny+1);
vv=zeros(nx+1,ny+1);
phiphi=zeros(nx+1,ny+1);
omega=zeros(nx+1,ny+1);

% Coefficient for pressure equation
c=zeros(nx+2,ny+2)+1/4;
c(2,3:ny)=1/3;c(nx+1,3:ny)=1/3;c(3:nx,2)=1/3;c(3:nx,ny+1)=1/3;
c(2,2)=1/2;c(2,ny+1)=1/2;c(nx+1,2)=1/2;c(nx+1,ny+1)=1/2;

% ----------------------------------------------------------------------- %
% Solution over time
% ----------------------------------------------------------------------- %
t=0.0;
for is=1:nsteps
    
    % Boundary conditions (velocity)
    u(1:nx+1,1)=2*us-u(1:nx+1,2);
    u(1:nx+1,ny+2)=2*un-u(1:nx+1,ny+1);
    v(1,1:ny+1)=2*vw-v(2,1:ny+1);
    v(nx+2,1:ny+1)=2*ve-v(nx+1,1:ny+1);
    
    % Advection-diffusion equation (predictor)
    [ut, vt] = AdvectionDiffusion2D( ut, vt, u, v, nx, ny, h, dt, nu);
    
    % Pressure equation (Poisson)
    [p, iter] = Poisson2D( p, ut, vt, c, nx, ny, h, dt, ...
                           beta, max_iterations, max_error);
    
    % Correct the velocity
    u(2:nx,2:ny+1)=ut(2:nx,2:ny+1)-(dt/h)*(p(3:nx+1,2:ny+1)-p(2:nx,2:ny+1));
    v(2:nx+1,2:ny)=vt(2:nx+1,2:ny)-(dt/h)*(p(2:nx+1,3:ny+1)-p(2:nx+1,2:ny));
    
    % Print time step on the screen
    if (mod(is,50)==1)
        fprintf('Step: %d - Time: %f - Poisson iterations: %d\n', is, t, iter);
    end
    
    % ----------------------------------------------------------------------- %
    % Passive scalar advection-diffusion equation
    % ----------------------------------------------------------------------- %

    % Boundary conditions (passive scalar)
    phi(2:nx+1,1)=2*phis-phi(2:nx+1,2);
    phi(2:nx+1,ny+2)=2*phin-phi(2:nx+1,ny+1);
    phi(1,2:ny+1)=2*phiw-phi(2,2:ny+1);
    phi(nx+2,2:ny+1)=2*phie-phi(nx+1,2:ny+1);

    % Update passive scalar solution
    phi = AdvectionDiffusion2DPassiveScalar( phi, u, v, nx, ny, h, dt, alpha );
    
  
    % Advance time
    t=t+dt;
 
    % ----------------------------------------------------------------------- %
    % Update graphical output
    % ----------------------------------------------------------------------- %
    if (mod(is,50)==1)
        
        % Post-processing (reconstructi on the corners of pressure cells)
        uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
        vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
        phiphi(1:nx+1,1:ny+1)=0.25*(phi(1:nx+1,1:ny+1)+phi(2:nx+2,1:ny+1) + ...
                                    phi(1:nx+1,2:ny+2)+phi(2:nx+2,2:ny+2));
        omega(1:nx+1,1:ny+1)=(  u(1:nx+1,2:ny+2)-u(1:nx+1,1:ny+1)-...
                                v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1) )/(2*h);
           
        % Contour: vorticity
        subplot(231);
        contour(X,Y,omega',20);
        axis('square'); xlabel('x [m]');ylabel('y [m]'); title('vorticity contour');

        % Contour: x-velocity
        subplot(232);
        contour(X,Y,uu');
        axis('square'); xlabel('x [m]');ylabel('y [m]'); title('x-velocity contour');

        % Contour: y-velocity
        subplot(233);
        contour(X,Y,vv');
        axis('square'); xlabel('x [m]');ylabel('y [m]'); title('y-velocity contour');

        % Plot: velocity components along the horizontal middle axis
        subplot(234);
        plot(x,uu(:, round(ny/2)));
        hold on;
        plot(x,vv(:, round(ny/2)));
        axis('square');
        xlabel('x [m]');ylabel('velocity [m/s]'); title('velocity along x-axis');
        legend('x-velocity', 'y-velocity');
        hold off;

        % Contour: passive scalar
        subplot(235);
        contour(X,Y,phiphi');
        axis('square'); xlabel('x [m]');ylabel('y [m]'); title('passive scalar contour');

        % Plot: velocity components along the horizontal middle axis
        subplot(236);
        plot(x,phiphi(:, round(ny/2)));
        hold on;
        plot(y,phiphi(round(ny/2),:));
        axis('square');
        xlabel('x or y [m]');ylabel('passive scalar'); title('passive scalar (centerlines)');
        legend('along x', 'along y'); ylim([min([phin,phis,phie,phiw]), max([phin,phis,phie,phiw])]);
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
function [ut, vt] = AdvectionDiffusion2D( ut, vt, u, v, nx, ny, h, dt, nu)
                            
    % Temporary u-velocity
    for i=2:nx
        for j=2:ny+1 
            
            ue2 = 0.25*( u(i+1,j)+u(i,j) )^2;
            uw2 = 0.25*( u(i,j)+u(i-1,j) )^2;
            unv = 0.25*( u(i,j+1)+u(i,j) )*( v(i+1,j)+v(i,j) );
            usv = 0.25*( u(i,j)+u(i,j-1) )*( v(i+1,j-1)+v(i,j-1) );
            
            A = (ue2-uw2+unv-usv)/h;
            D = (nu/h^2)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j));
            
            ut(i,j)=u(i,j)+dt*(-A+D);
            
        end
    end
    
    % Temporary v-velocity
    for i=2:nx+1
        for j=2:ny 
            
            vn2 = 0.25*( v(i,j+1)+v(i,j) )^2;
            vs2 = 0.25*( v(i,j)+v(i,j-1) )^2;
            veu = 0.25*( u(i,j+1)+u(i,j) )*( v(i+1,j)+v(i,j) );
            vwu = 0.25*( u(i-1,j+1)+u(i-1,j) )*( v(i,j)+v(i-1,j) );
            A = (vn2 - vs2 + veu - vwu)/h;
            D = (nu/h^2)*(v(i+1,j)+v(i-1,j)+v(i,j+1)+v(i,j-1)-4*v(i,j));
            
            vt(i,j)=v(i,j)+dt*(-A+D);
            
        end
    end
    
end

% --------------------------------------------------------------------------------------
% Advection-diffusion equation for passive scalar
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
                
                % Passive scalar on the faces (linear interpolation)
                phie = 0.50*(phio(i+1,j)+phio(i,j));
                phiw = 0.50*(phio(i-1,j)+phio(i,j));
                phin = 0.50*(phio(i,j+1)+phio(i,j));
                phis = 0.50*(phio(i,j-1)+phio(i,j));
                
                % Gradients of passive scalar on the faces (centered)
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
