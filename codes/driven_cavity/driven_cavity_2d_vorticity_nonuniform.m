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
%  Code: 2D driven-cavity problem in vorticity/streamline formulation     %
%        The code is adapted and extended from Tryggvason, Computational  %
%        Fluid Dynamics http://www.nd.edu/~gtryggva/CFD-Course/           %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% Basic setup
Ly_over_Lx = 1.;        % ratio between y and x lengths
nx=25;                  % number of grid points along x
ny=25;                  % number of grid points along y
deltax=2;               % stretching factor along x
deltay=2;               % stretching factor along y
Re=100;                 % Reynolds number [-]
tau=20;                 % total time of simulation [-]

% Data for reconstructing the velocity field
L=1;                     % length along x [m] (used as reference length)
nu=1e-3;                 % kinematic viscosity [m2/s] 
Uwall=nu*Re/L;           % wall velocity [m/s]
 
% Parameters for SOR
max_iterations=10000;   % maximum number of iterations
beta=1.5;               % SOR coefficient
max_error=0.0001;       % error for convergence

% Grid construction (dimensionless)
x = zeros(nx,1);
for i=1:nx
    x(i) = 0.5*(1+tanh(deltax*((i-1)/(nx-1)-0.5))/tanh(deltax/2));
end
y = zeros(ny,1);
for i=1:ny
    y(i) = 0.5*(1+tanh(deltay*((i-1)/(ny-1)-0.5))/tanh(deltay/2));
end
y = Ly_over_Lx*y;

% Time step
h2 = (x(2)-x(1))*(y(2)-y(1));       % minimum cell volume
sigma = 0.5;                        % safety factor for time step (stability)
dt_diff=h2*Re/4;                    % time step (diffusion stability)
dt_conv=4/Re;                       % time step (convection stability)
dt=sigma*min(dt_diff, dt_conv);     % time step (stability)
nsteps=tau/dt;                      % number of steps

fprintf('Time step: %f\n', dt);
fprintf(' - Diffusion:  %f\n', dt_diff);
fprintf(' - Convection: %f\n', dt_conv);

% Memory allocation
psi=zeros(nx,ny);       % streamline function
omega=zeros(nx,ny);     % vorticity
psio=zeros(nx,ny);      % streamline function at previous time
omegao=zeros(nx,ny);    % vorticity at previous time
u=zeros(nx,ny);         % reconstructed dimensionless x-velocity
v=zeros(nx,ny);         % reconstructed dimensionless y-velocity
U=zeros(nx,ny);         % reconstructed x-velocity [m/s]
V=zeros(nx,ny);         % reconstructed y-velocity [m/s]

% Mesh construction (only needed in graphical post-processing)
[X,Y] = meshgrid(x,y);  % mesh

% Allocation of vectors for non-uniform grid (see the Poisson equation function)
% Along the x direction
ae = zeros(nx,1);   ax = zeros(nx,1);   aw = zeros(nx,1);
for i=2:nx-1
    a = x(i)-x(i-1); b = x(i+1)-x(i-1); c = x(i+1)-x(i);
    ae(i) = 2/(b*c);
    ax(i) = 2/(a*c);
    aw(i) = 2/(a*b);
end
    
% Along the y direction
an = zeros(ny,1);   ay = zeros(ny,1);   as = zeros(ny,1);
for j=2:ny-1
    a = y(j)-y(j-1); b = y(j+1)-y(j-1); c = y(j+1)-y(j);
    an(j) = 2/(b*c);
    ay(j) = 2/(a*c);
    as(j) = 2/(a*b);
end

% Time loop
t = 0;
for istep=1:nsteps     
    
    % ------------------------------------------------------------------- %
    % Poisson equation (SOR)
    % ------------------------------------------------------------------- %
    [psi, iter] = Poisson2D( psi, x, y, omega, ...
                             beta, max_iterations, max_error, ...
                             ae, ax, aw, an, ay, as );
    
    % ------------------------------------------------------------------- %
    % Find vorticity on boundaries
    % ------------------------------------------------------------------- %
    omega(2:nx-1,1)=-2.0*psi(2:nx-1,2)/(y(2)-y(1))^2;               % south
    omega(2:nx-1,ny)=-2.0*psi(2:nx-1,ny-1)/(y(ny)-y(ny-1))^2 ...
                     -2.0/(y(ny)-y(ny-1))*1;                        % north
    omega(1,2:ny-1)=-2.0*psi(2,2:ny-1)/(x(2)-x(1))^2;               % east
    omega(nx,2:ny-1)=-2.0*psi(nx-1,2:ny-1)/(x(nx)-x(nx-1))^2;       % west
  
    % ------------------------------------------------------------------- %
    % Reconstruction of dimensionless velocity field
    % ------------------------------------------------------------------- %
    u(:,ny)=1;
    for i=2:nx-1 
         for j=2:ny-1
             u(i,j) =  (psi(i,j+1)-psi(i,j-1))/(y(j+1)-y(j-1));
             v(i,j) = -(psi(i+1,j)-psi(i-1,j))/(x(i+1)-x(i-1));
         end
    end
    
    % ------------------------------------------------------------------- %
    % Find new vorticity in interior points
    % ------------------------------------------------------------------- %
    [omega] = AdvectionDiffusion(omega, x, y, u, v, Re, dt);
   
    % Message on video
    if (mod(istep,50)==1)
        fprintf('Step: %d - Time: %f - Poisson iterations: %d\n', istep, t, iter);
    end
    
    % Advance time
    t=t+dt;

    % ------------------------------------------------------------------- %
    % Reconstruction of velocity field [m/s]
    % ------------------------------------------------------------------- %
    U = u*Uwall;
    V = v*Uwall;
    
    % ------------------------------------------------------------------- %
    % Graphics only
    % ------------------------------------------------------------------- %
    plot_2x4 = false;   % plotting the 2x4 plot
    
    if (plot_2x4 == true)
        
        subplot(241);
        contour(x,y,omega');
        axis('square'); title('omega'); xlabel('x'); ylabel('y');

        subplot(245);
        contour(x,y,psi');
        axis('square'); title('psi'); xlabel('x'); ylabel('y');

        subplot(242);
        contour(x,y,u');
        axis('square'); title('u'); xlabel('x'); ylabel('y');

        subplot(246);
        contour(x,y,v');
        axis('square'); title('v'); xlabel('x'); ylabel('y');

        subplot(243);
        plot(x,u(:, round(ny/2)));
        hold on;
        plot(x,v(:, round(ny/2)));
        axis('square'); legend('u', 'v');
        title('velocities along HA'); xlabel('x'); ylabel('velocities');
        hold off;
        
        subplot(247);
        plot(y,u(round(nx/2),:));
        hold on;
        plot(y,v(round(nx/2),:));
        axis('square'); legend('u', 'v');
        title('velocities along VA'); xlabel('y'); ylabel('velocities');
        hold off;

        subplot(244);
        quiver(x,y,u',v');
        axis('square', [0 1 0 Ly_over_Lx]);
        title('velocity vectors'); xlabel('x'); ylabel('y');
    
        pause(0.001);
        
    end
    
end

% ------------------------------------------------------------------- %
% Write final maps
% ------------------------------------------------------------------- %

subplot(231);
surface(x,y,u');
axis('square'); title('u'); xlabel('x'); ylabel('y');

subplot(234);
surface(x,y,v');
axis('square'); title('v'); xlabel('x'); ylabel('y');

subplot(232);
surface(x,y,omega');
axis('square'); title('omega'); xlabel('x'); ylabel('y');

subplot(235);
surface(x,y,psi');
axis('square'); title('psi'); xlabel('x'); ylabel('y');

subplot(233);
contour(x,y,psi', 30, 'b');
axis('square');
title('stream lines'); xlabel('x'); ylabel('y');

subplot(236);
quiver(x,y,u',v');
axis([0 1 0 Ly_over_Lx], 'square');
title('stream lines'); xlabel('x'); ylabel('y');

% ------------------------------------------------------------------- %
% Write velocity profiles along the centerlines for exp comparison
% ------------------------------------------------------------------- %
u_profile = u(round(nx/2),:);
fileVertical = fopen('vertical.out','w');
for i=1:ny 
    fprintf(fileVertical,'%f %f\n',y(i), u_profile(i));
end
fclose(fileVertical);

v_profile = v(:,round(ny/2));
fileHorizontal = fopen('horizontal.out','w');
for i=1:nx
    fprintf(fileHorizontal,'%f %f\n',x(i), v_profile(i));
end
fclose(fileHorizontal);

% ------------------------------------------------------------------- %
% Compare with exp data (available only for Re=100, 400, and 1000)
% ------------------------------------------------------------------- %
% Read experimental data from file
exp_u_along_y = dlmread('experimental_data/u_along_y.exp', '', 1, 0);
exp_v_along_x = dlmread('experimental_data/v_along_x.exp', '', 1, 0);

% Comparison with exp data
% Be careful: cols 1,2 for Re=100, 3,4 for Re=400, 5,6 for Re=1000
figure;
plot(exp_u_along_y(:,1), exp_u_along_y(:,2), 'o', y, u_profile, '-');
axis('square'); title('u along y (centerline)'); xlabel('y'); ylabel('u');

figure;
plot(exp_v_along_x(:,1), exp_v_along_x(:,2), 'o', x, v_profile, '-');
axis('square'); title('v along x (centerline)'); xlabel('x'); ylabel('v');


% --------------------------------------------------------------------------------------
% Poisson equation solver
% The second order derivative is discretized over a non uniform grid
% d2(psi)/dx2 = ae*psi(i+1)-ac*psi(i)+aw*psi(i-1)
%               ae = 2/(b*c), ax = 2/(a*c), aw = 2/(a*b)
%               a=x(i)-x(i-1), b=x(i+1)-x(i-1), c=x(i+1)-x(i)
% --------------------------------------------------------------------------------------
function [psi, iter] = Poisson2D(   psi, x, y, omega, beta, max_iterations, max_error, ...
                                    ae, ax, aw, an, ay, as)

    nx = length(x);
    ny = length(y);
    
    for iter=1:max_iterations
        
        % Update solution
        for i=2:nx-1
            for j=2:ny-1 
                psi(i,j)=( psi(i+1,j)*ae(i) + psi(i-1,j)*aw(i) + ...
                           psi(i,j+1)*an(j) + psi(i,j-1)*as(j) + ...
                           omega(i,j) ) / (ax(i)+ay(j))*beta + ...
                           (1.0-beta)*psi(i,j);

            end
        end
        
        % Estimate the error
        epsilon=0.0; 
        for i=2:nx-1
            for j=2:ny-1
                epsilon=epsilon+abs( psi(i+1,j)*ae(i) - psi(i,j)*ax(i) + psi(i-1,j)*aw(i) + ...
                                     psi(i,j+1)*an(j) - psi(i,j)*ay(j) + psi(i,j-1)*as(j) + ...
                                     omega(i,j) ); 
            end
        end
        epsilon = epsilon/(nx-2)/(ny-2);
        
        % Check the error
        if (epsilon <= max_error) % stop if converged
            break;
        end 
        
    end
    
end

% --------------------------------------------------------------------------------------
% Advection-diffusion equation: forward Euler + centered discretization
% --------------------------------------------------------------------------------------
function [omega] = AdvectionDiffusion(omega, x, y, u, v, Re, dt)

    nx = length(x);
    ny = length(y);
    
    omegao=omega;
    
    for i=2:nx-1 
    
        ax = x(i)-x(i-1); bx = x(i+1)-x(i-1); cx = x(i+1)-x(i);
        
        for j=2:ny-1
            ay = y(j)-y(j-1); by = y(j+1)-y(j-1); cy = y(j+1)-y(j);

            advection_x = -u(i,j)*(omegao(i+1,j)-omegao(i-1,j))/bx;
            advection_y = -v(i,j)*(omegao(i,j+1)-omegao(i,j-1))/by;
            
            diffusion_x = 1/Re*( ax*omegao(i+1,j)-bx*omegao(i,j)+...
                                 cx*omegao(i-1,j))/(0.5*ax*bx*cx);
            diffusion_y = 1/Re*( ay*omegao(i,j+1)-by*omegao(i,j)+...
                                 cy*omegao(i,j-1))/(0.5*ay*by*cy);

            omega(i,j)=omegao(i,j) + ...
                        dt*( advection_x + advection_y + ...
                             diffusion_x + diffusion_y );
        end
        
    end
    
end
