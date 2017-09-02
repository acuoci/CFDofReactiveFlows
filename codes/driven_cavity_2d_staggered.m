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
%	License                                                               %
%                                                                         %
%   Copyright(C) 2017 Alberto Cuoci                                       %
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
%  Code: 2D driven-cavity problem on a staggered grid                     %
%        The code is adapted and extended from Tryggvason, Computational  %
%        Fluid Dynamics http://www.nd.edu/~gtryggva/CFD-Course/           %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% ----------------------------------------------------------------------- %
% User data
% ----------------------------------------------------------------------- %

nx=48;              % number of (physical) cells along x
ny=nx;              % number of (physical) cells along y
L=1;                % length [m]
nu=0.1;             % kinematic viscosity [m2/s] 
tau=1;              % total time of simulation [-]

% Boundary conditions
un=1;       % north wall velocity [m/s]
us=0;       % south wall velocity [m/s]
ve=0;       % east wall velocity [m/s]
vw=0;       % west wall velocity [m/s]

% Parameters for SOR
max_iterations=1000;    % maximum number of iterations
beta=1.2;               % SOR coefficient
max_error=1e-4;         % error for convergence

% ----------------------------------------------------------------------- %
% Data processing
% ----------------------------------------------------------------------- %

% Grid step
h=L/nx;             % grid step [m]

% Time step
sigma = 0.5;                        % safety factor for time step (stability)
u2=(un^2+us^2+ve^2+vw^2);           % velocity measure [m2/s2]
dt_diff=h^2/4/nu;                   % time step (diffusion stability)
dt_conv=4*nu/u2;                    % time step (convection stability)
dt=sigma*min(dt_diff, dt_conv);     % time step (stability)
nsteps=tau/dt;                      % number of steps

fprintf('Time step: %f\n', dt);
fprintf(' - Diffusion:  %f\n', dt_diff);
fprintf(' - Convection: %f\n', dt_conv);

% Grid construction
x=zeros(nx+1,ny+1);         % grid coordinates (x axis)
y=zeros(nx+1,ny+1);         % grid coordinates (y axis)
for i=1:nx+1
    for j=1:ny+1
        x(i,j)=h*(i-1);
        y(i,j)=h*(j-1);
    end
end;

% ----------------------------------------------------------------------- %
% Memory allocation
% ----------------------------------------------------------------------- %

% Main fields (velocities and pressure)
u=zeros(nx+1,ny+2);
v=zeros(nx+2,ny+1);
p=zeros(nx+2,ny+2);

% Temporary velocity fields
ut=zeros(nx+1,ny+2);
vt=zeros(nx+2,ny+1);

% Temporary pressure field (convergence of SOR)
po=zeros(nx+2,ny+2);

% Fields used only for graphical post-processing purposes
uu=zeros(nx+1,ny+1);
vv=zeros(nx+1,ny+1);
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
    
    % Boundary conditions
    u(1:nx+1,1)=2*us-u(1:nx+1,2);
    u(1:nx+1,ny+2)=2*un-u(1:nx+1,ny+1);
    v(1,1:ny+1)=2*vw-v(2,1:ny+1);
    v(nx+2,1:ny+1)=2*ve-v(nx+1,1:ny+1);
    
    % Temporary u-velocity
    for i=2:nx
        for j=2:ny+1 
            ut(i,j)=u(i,j)+dt*(-(0.25/h)*((u(i+1,j)+u(i,j))^2-(u(i,j)+...
                    u(i-1,j))^2+(u(i,j+1)+u(i,j))*(v(i+1,j)+...
                    v(i,j))-(u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))+...
                    (nu/h^2)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j)));
        end
    end
    
    % Temporary v-velocity
    for i=2:nx+1
        for j=2:ny 
            vt(i,j)=v(i,j)+dt*(-(0.25/h)*((u(i,j+1)+u(i,j))*(v(i+1,j)+...
                    v(i,j))-(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j))+...
                    (v(i,j+1)+v(i,j))^2-(v(i,j)+v(i,j-1))^2)+...
                    (nu/h^2)*(v(i+1,j)+v(i-1,j)+v(i,j+1)+v(i,j-1)-4*v(i,j)));
        end
    end
    
    % Pressure equation (Poisson)
    for it=1:max_iterations
        
        po=p;
        for i=2:nx+1
            for j=2:ny+1
                p(i,j)=beta*c(i,j)*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1)-...
                        (h/dt)*(ut(i,j)-ut(i-1,j)+vt(i,j)-vt(i,j-1)))+...
                        (1-beta)*p(i,j);
            end
        end
        
        % Estimate the error
        epsilon=0.0; 
        for i=2:nx+1
            for j=2:ny+1
                epsilon=epsilon+abs(po(i,j)-p(i,j)); 
            end
        end
        epsilon = epsilon / (nx*ny);
        
        % Check the error
        if (epsilon <= max_error) % stop if converged
            break;
        end 
        
    end
    
    % Correct the velocity
    u(2:nx,2:ny+1)=ut(2:nx,2:ny+1)-(dt/h)*(p(3:nx+1,2:ny+1)-p(2:nx,2:ny+1));
    v(2:nx+1,2:ny)=vt(2:nx+1,2:ny)-(dt/h)*(p(2:nx+1,3:ny+1)-p(2:nx+1,2:ny));
    
    fprintf('Step: %d - Time: %f - Poisson iterations: %d\n', is, t, it);
 
    % Advance time
    t=t+dt;
 
    % Post-processing
    uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
    vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
    omega(1:nx+1,1:ny+1)=(  u(1:nx+1,2:ny+2)-u(1:nx+1,1:ny+1)-...
                            v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1) )/(2*h);
 
    % Contour: vorticity
    subplot(241);
    contour(x,y,omega,20);
    axis('square');
    xlabel('x [m]');ylabel('y [m]'); title('vorticity contour');
    hold off;
    
    % Contour: x-velocity
    subplot(242);
    contour(x,y,uu);
    axis('square');
    xlabel('x [m]');ylabel('y [m]'); title('x-velocity contour');
    hold off;
    
    % Contour: y-velocity
    subplot(246);
    contour(x,y,vv);
    axis('square');
    xlabel('x [m]');ylabel('y [m]'); title('y-velocity contour');
    
    % Plot: velocity components along the horizontal middle axis
    subplot(243);
    plot(x(:,round(ny/2)),uu(:, round(ny/2)));
    hold on;
    plot(x(:,round(ny/2)),vv(:, round(ny/2)));
    axis('square');
    xlabel('x [m]');ylabel('velocity [m/s]'); title('velocity along x-axis');
    legend('x-velocity', 'y-velocity');
    hold off;

    % Plot: velocity components along the vertical middle axis
    subplot(247);
    plot(y(round(nx/2),:),uu(round(nx/2),:));
    hold on;
    plot(y(round(nx/2),:),vv(round(nx/2),:));
    axis('square');
    xlabel('y [m]');ylabel('velocity [m/s]'); title('velocity along y-axis');
    legend('x-velocity', 'y-velocity');
    hold off;
    
    % Velocity vector field
    subplot(244);
    quiver(x,y,uu,vv);
    axis('square','equal');
    xlabel('x [m]');ylabel('y [m]'); title('velocity vector field');

    pause(0.01)
    
end

