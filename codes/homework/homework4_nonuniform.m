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
% ----------------------------------------------------------------------- %

close all;
clear variables;

%% Input data
%--------------------------------------------------------------------------
% Basic setup
nx=25;                  % number of grid points along x
ny=nx;                  % number of grid points along y
Lx=1;                   % length along x [m]
Ly=Lx;                  % length along y [m]
deltax=2;               % stretching factor along x
deltay=2;               % stretching factor along y
Uwall=0.1;              % wall velocity [m/s]
nu=1e-3;                % kinematic viscosity [m2/s]
ttot=200;               % total time of simulation [s]

% Parameters for SOR
max_iterations=10000;   % maximum number of iterations
beta=1.5;               % SOR coefficient
max_error=0.0001;       % error for convergence


%% Pre-processing operations
%--------------------------------------------------------------------------
L=max(Lx,Ly);           % reference length [m]
V=Uwall;                % reference velocity [m/s]
Re=V*L/nu;              % Reynolds' number [-]
tref=L/V;               % reference time [s]
lengthx=Lx/L;           % dimensionless length along x [-] 
lengthy=Ly/L;           % dimensionless length along y [-] 
tautot=ttot/tref;       % total time of simulation (dimensionless) [-]
uwall = Uwall/V;        % dimensionless wall velocity [-]

% Grid
x = zeros(nx,1);
for i=1:nx
    x(i) = 0.5*(1+tanh(deltax*((i-1)/(nx-1)-0.5))/tanh(deltax/2));
end
x = lengthx*x;
y = zeros(ny,1);
for i=1:ny
    y(i) = 0.5*(1+tanh(deltay*((i-1)/(ny-1)-0.5))/tanh(deltay/2));
end
y = lengthy*y;


% Time step
hx_min = min( x(2:end)-x(1:end-1) );    % minimum grid step along x [-]
hy_min = min( y(2:end)-y(1:end-1) );    % minimum grid step along y [-]
h_min = min( hx_min, hy_min);           % minimum grid step [-]
sigma = 0.5;                            % safety factor for time step (stability)
dtau_diff=h_min^2*Re/4;                 % time step (diffusion stability) [-]
dtau_conv=4/Re;                         % time step (convection stability) [-]
dtau=sigma*min(dtau_diff, dtau_conv);   % time step (stability) [-]
nsteps=tautot/dtau;                     % number of steps

fprintf('Time step: %f\n', dtau);
fprintf(' - Diffusion:  %f\n', dtau_diff);
fprintf(' - Convection: %f\n', dtau_conv);

% Memory allocation
psi=zeros(nx,ny);       % streamline function [-]
omega=zeros(nx,ny);     % vorticity [-]
u=zeros(nx,ny);         % reconstructed dimensionless x-velocity [-]
v=zeros(nx,ny);         % reconstructed dimensionless y-velocity [-]

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

%% Numerical solution
%--------------------------------------------------------------------------
tau = 0;
for istep=1:nsteps     
    
    % ------------------------------------------------------------------- %
    % Poisson equation (SOR)
    % ------------------------------------------------------------------- %
    [psi,iter] = Poisson2D( psi,x,y,-omega, ...
                            beta,max_iterations,max_error, ...
                            ae, ax, aw, an, ay, as );
    
    % ------------------------------------------------------------------- %
    % Reconstruction of dimensionless velocity field
    % ------------------------------------------------------------------- % 
    [u,v] = ReconstructDimensionlessVelocity(u,v,psi,x,y,uwall);
    
    % ------------------------------------------------------------------- %
    % Find vorticity on boundaries
    % ------------------------------------------------------------------- %
    omega(2:nx-1,1) = -2.0*psi(2:nx-1,2)/(y(2)-y(1))^2;              % south
    omega(2:nx-1,ny)= -2.0*psi(2:nx-1,ny-1)/(y(ny)-y(ny-1))^2 ...
                      -2.0/(y(ny)-y(ny-1))*1;                        % north
    omega(1,2:ny-1) = -2.0*psi(2,2:ny-1)/(x(2)-x(1))^2;              % east
    omega(nx,2:ny-1)= -2.0*psi(nx-1,2:ny-1)/(x(nx)-x(nx-1))^2;       % west
    
    % ------------------------------------------------------------------- %
    % Advection-diffusion equation (new vorticity in interior points)
    % ------------------------------------------------------------------- %
    [omega] = AdvectionDiffusion2D(omega, x,y, u,v, Re, dtau);
    
    % ------------------------------------------------------------------- %
    % Advancing time
    % ------------------------------------------------------------------- %    
    if (mod(istep,25)==1)
        fprintf('Step: %d - Time: %f - Iterations: %d\n', ...
                istep, tau, iter);
    end
    
    tau=tau+dtau;
    
    % ------------------------------------------------------------------- %
    % On-the-fly graphical post-processing
    % ------------------------------------------------------------------- %
    if (mod(istep,25)==0)
        contour(x,y,psi', 30, 'b');
        axis('square');
        pause(0.01);
    end
    
end


%% Final post-processing operations
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
axis([0 lengthx 0 lengthy], 'square');
title('stream lines'); xlabel('x'); ylabel('y');


%% ------------------------------------------------------------------------
% Poisson equation solver
% The second order derivative is discretized over a non uniform grid
% d2(psi)/dx2 = ae*psi(i+1)-ac*psi(i)+aw*psi(i-1)
%               ae = 2/(b*c), ax = 2/(a*c), aw = 2/(a*b)
%               a=x(i)-x(i-1), b=x(i+1)-x(i-1), c=x(i+1)-x(i)
%  ------------------------------------------------------------------------
function [f,iter] = Poisson2D( f, x,y, S, beta,max_iterations,max_error, ...
                               ae,ax,aw,an,ay,as)

    nx = length(x);
    ny = length(y);
    
    for iter=1:max_iterations
        
        for i=2:nx-1
            for j=2:ny-1
                
                f(i,j) = beta*( ae(i)*f(i+1,j) + aw(i)*f(i-1,j) + ...
                                an(j)*f(i,j+1) + as(j)*f(i,j-1) + ...
                                -S(i,j) ) / (ax(i)+ay(j)) + ...
                         (1-beta)*f(i,j);
                
            end
        end
        
        res = 0;
        for i=2:nx-1
            for j=2:ny-1
                res = res+abs( f(i+1,j)*ae(i) - f(i,j)*ax(i) + f(i-1,j)*aw(i) + ...
                               f(i,j+1)*an(j) - f(i,j)*ay(j) + f(i,j-1)*as(j) + ...
                               -S(i,j) ); 
            end
        end
        res = res/(nx-2)/(ny-2);
        
        if (res <= max_error)
            break;
        end
        
    end

end


%% ------------------------------------------------------------------------
%  Reconstruction of velocity field (dimensionless)
%  ------------------------------------------------------------------------
function [u,v] = ReconstructDimensionlessVelocity(u,v,psi,x,y,uwall)

    nx = length(x);
    ny = length(y);
    
    u(:,ny) = uwall; 
    for i=2:nx-1
        for j=2:ny-1
            u(i,j) =  ( psi(i,j+1)-psi(i,j-1) )/(y(j+1)-y(j-1));
            v(i,j) = -( psi(i+1,j)-psi(i-1,j) )/(x(i+1)-x(i-1));
        end
    end

end


%% ------------------------------------------------------------------------
%  Advection-diffusion equation: forward Euler + centered discretization
%  ------------------------------------------------------------------------
function [f] = AdvectionDiffusion2D(f, x,y, u,v, Re, dtau)

    nx = length(x);
    ny = length(y);
    
    fo = f;
    for i=2:nx-1
        
        ax = x(i)-x(i-1); bx = x(i+1)-x(i-1); cx = x(i+1)-x(i);
        
        for j=2:ny-1
            
            ay = y(j)-y(j-1); by = y(j+1)-y(j-1); cy = y(j+1)-y(j);
            
            advection_x = -u(i,j)*(fo(i+1,j)-fo(i-1,j))/bx;
            advection_y = -v(i,j)*(fo(i,j+1)-fo(i,j-1))/by;
            
            diffusion_x = 1/Re*( ax*fo(i+1,j)-bx*fo(i,j)+...
                                 cx*fo(i-1,j))/(0.5*ax*bx*cx);
            diffusion_y = 1/Re*( ay*fo(i,j+1)-by*fo(i,j)+...
                                 cy*fo(i,j-1))/(0.5*ay*by*cy);

            f(i,j) = fo(i,j) + ...
                     dtau*( advection_x + advection_y + ...
                            diffusion_x + diffusion_y );
                         
        end
    end

end
