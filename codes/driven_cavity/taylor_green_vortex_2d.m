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
%  Code: Taylor-Green vortex in 2D
%  The Taylor-Green vortex is an exact closed form solution of 2D,        %
%  incompressible Navier-Stokes eqs. This 2D decaying vortex defined      %
%  in the square domain, 0-2pi, serves as a benchmark problem             %
%  for testing and validation of incompressible Navier-Stokes codes.      %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% ----------------------------------------------------------------------- %
% User data
% ----------------------------------------------------------------------- %

% Only even numbers of cells are acceptable
nx=60;      % number of (physical) cells along x
ny=nx;      % number of (physical) cells along y
L=2*pi;     % length [m]
nu=1;       % kinematic viscosity [m2/s]
tau=1;      % total time of simulation [s]

% Parameters for SOR
max_iterations=20000;   % maximum number of iterations
beta=1.5;               % SOR coefficient
max_error=1e-7;         % error for convergence

% Process the grid
h=L/nx;                 % grid step (uniform grid) [m]

% Initial Solution
t = 0;
[u,v] = AnalyticalSolutionVelocity(nx,ny, h, nu, t);
[p] = AnalyticalSolutionPressure(nx,ny, h, nu, t);

% Time step
umax = sqrt( max(max(u))^2 + ...
             max(max(v))^2 );       % maximum velocity [m/s]
sigma = 0.5;                        % safety factor for time step (stability)
dt_diff=h^2/4/nu;                   % time step (diffusion stability) [s]
dt_conv=4*nu/umax^2;                % time step (convection stability) [s]
dt=sigma*min(dt_diff, dt_conv);     % time step (stability) [s]
nsteps=tau/dt;                      % number of steps

fprintf('Time step: %f\n', dt);
fprintf(' - Diffusion:  %f\n', dt_diff);
fprintf(' - Convection: %f\n', dt_conv);

% Grid construction
x=0:h:L;                         % grid coordinates (x axis)
y=0:h:L;                         % grid coordinates (y axis)
[X,Y] = meshgrid(x,y);           % MATLAB grid

% Temporary velocity fields
ut=u;
vt=v;

% Coefficient for pressure equation
gamma=zeros(nx+2,ny+2)+1/4;
gamma(2,3:ny)=1/3;gamma(nx+1,3:ny)=1/3;gamma(3:nx,2)=1/3;gamma(3:nx,ny+1)=1/3;
gamma(2,2)=1/2;gamma(2,ny+1)=1/2;gamma(nx+1,2)=1/2;gamma(nx+1,ny+1)=1/2;

% ----------------------------------------------------------------------- %
% Solution over time
% ----------------------------------------------------------------------- %
t=0.0;
for is=1:1:nsteps
    
    % Boundary conditions
    u(1:nx+1,1)=u(1:nx+1,ny+1);         % south wall
    u(1:nx+1,ny+2)=u(1:nx+1,2);         % north wall
    v(1,1:ny+1)=v(nx+1,1:ny+1);         % west wall
    v(nx+2,1:ny+1)=v(2,1:ny+1);         % east wall
    
    % Advection-diffusion equation (predictor)
    [ut, vt] = AdvectionDiffusion2D( ut, vt, u, v, nx, ny, h, dt, nu);
    
    % Pressure equation (Poisson)
    [p, iter] = Poisson2D( p, ut, vt, gamma, nx, ny, h, dt, ...
                           beta, max_iterations, max_error );
    
    % Correction on the velocity
    u(2:nx,2:ny+1)=ut(2:nx,2:ny+1)-(dt/h)*(p(3:nx+1,2:ny+1)-p(2:nx,2:ny+1));
    v(2:nx+1,2:ny)=vt(2:nx+1,2:ny)-(dt/h)*(p(2:nx+1,3:ny+1)-p(2:nx+1,2:ny));
    
    % Print on the screen
    if (mod(is,100)==1)
        
        fprintf('Step: %d - Time: %f - Poisson iterations: %d\n', is, t, iter);
        
        % Reconstruct fields
        [uu, vv, pp] = FieldReconstruction(u,v,p, nx,ny);

        % Surface map: u-velocity
        subplot(221);
        surface(X,Y,uu');
        axis('square'); title('u'); xlabel('x'); ylabel('y');
        shading interp; colorbar;

        % Surface map: v-velocity
        subplot(222);
        surface(X,Y,vv');
        axis('square'); title('v'); xlabel('x'); ylabel('y');
        shading interp; colorbar;

        % Surface map: pressure
        subplot(223);
        surface(X,Y,pp');
        axis('square'); title('p'); xlabel('x'); ylabel('y');
        shading interp; colorbar;
        
        % Plot: velocity components along the horizontal middle axis
        subplot(224);
        plot(x,uu(:, round(ny/2)+1));
        hold on;
        plot(x,vv(:, round(ny/2)+1));
        axis('square');
        xlabel('x [m]');ylabel('velocity [m/s]'); title('velocity along x-axis');
        legend('x-velocity', 'y-velocity');
        hold off;

        pause(0.0001);
        
    end
    
    % Advance in time
    t=t+dt;
 
end

% ----------------------------------------------------------------------- %
% Final post-processing                                                   %
% ----------------------------------------------------------------------- %

% Field reconstruction
[uu, vv, pp] = FieldReconstruction(u,v,p, nx,ny);
                    
% Analytical Solution
[ua,va] = AnalyticalSolutionVelocity(nx,ny,h, nu, t);
[pa] = AnalyticalSolutionPressure(nx,ny,h,nu,t);
[uua, vva, ppa] = FieldReconstruction(ua,va,pa, nx,ny);

% Surface map: u-velocity
subplot(241);
surface(X,Y,uu'); shading interp; colorbar;
axis('square'); title('u (numerical)'); xlabel('x'); ylabel('y');

% Surface map: v-velocity
subplot(242);
surface(X,Y,vv'); shading interp; colorbar;
axis('square'); title('v (numerical)'); xlabel('x'); ylabel('y');

% Surface map: pressure
subplot(243);
surface(X,Y,pp'); shading interp; colorbar;
axis('square'); title('p (numerical)'); xlabel('x'); ylabel('y');

% Surface map: difference pressure
subplot(244);
surface(X,Y,(ppa-pp)'); shading interp; colorbar;
axis('square'); title('delta(p)'); xlabel('x'); ylabel('y');

% Surface map: u-velocity (analytical)
subplot(245);
surface(X,Y,uua'); shading interp; colorbar;
axis('square'); title('u (analytical)'); xlabel('x'); ylabel('y');

% Surface map: v-velocity (analytical)
subplot(246);
surface(X,Y,vva'); shading interp; colorbar;
axis('square'); title('v (analytical)'); xlabel('x'); ylabel('y');

% Surface map: pressure (analytical)
subplot(247);
surface(X,Y,ppa'); shading interp; colorbar;
axis('square'); title('p (analytical)'); xlabel('x'); ylabel('y');

% Streamlines
subplot(248);
sx = [0.25*L 0.75*L 0.25*L 0.75*L];
sy = [0.25*L 0.75*L 0.75*L 0.25*L];
streamline(X,Y,uu',vv',sx,sy)
axis([0 L 0 L], 'square');
title('streamlines'); xlabel('x'); ylabel('y');

 
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
% Field reconstruction (for graphical purposes only)
% --------------------------------------------------------------------------------------
function [uu, vv, pp] = FieldReconstruction(u,v,p, nx,ny)

    % u-velocity
    uu=zeros(nx+1,ny+1);
    uu(1:nx+1,1:ny+1)=0.50*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
    
    % v-velocity
    vv=zeros(nx+1,ny+1);
    vv(1:nx+1,1:ny+1)=0.50*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
    
    % Pressure
    pp=zeros(nx+1,ny+1);
    p(1:nx+2,1) = p(1:nx+2,ny+1);       % south wall
    p(1:nx+2,ny+2) = p(1:nx+2,2);       % north wall
    p(1,1:ny+2) = p(nx+1,1:ny+2);       % west wall
    p(nx+2, 1:ny+2) = p(2,1:ny+2);      % east wall
        
    pp(1:nx+1,1:ny+1)=0.25*(p(1:nx+1,1:ny+1)+p(1:nx+1,2:ny+2)+...
                            p(2:nx+2,1:ny+1)+p(2:nx+2,2:ny+2) );
                   
end

% --------------------------------------------------------------------------------------
% Analytical solution (velocity)
% --------------------------------------------------------------------------------------
function [u,v] = AnalyticalSolutionVelocity(nx,ny, h, nu, t)

    u=zeros(nx+1,ny+2);
    v=zeros(nx+2,ny+1);

    for i=1:nx+1
        for j=1:ny+2
            x = (i-1)*h;
            y = h/2+(j-2)*h;
            u(i,j) = -sin(x)*cos(y)*F(t,nu);
        end
    end
    
    for i=1:nx+2
        for j=1:ny+1
            x = h/2+(i-2)*h;
            y = (j-1)*h;
            v(i,j) = cos(x)*sin(y)*F(t,nu);
        end
    end

end

% --------------------------------------------------------------------------------------
% Analytical solution (pressure)
% --------------------------------------------------------------------------------------
function [p] = AnalyticalSolutionPressure(nx,ny, h, nu, t)

    p=zeros(nx+2,ny+2);

    for i=2:nx+1
        for j=2:ny+1
            x = h/2+(i-2)*h;
            y = h/2+(j-2)*h;
            p(i,j) = 0.25*(cos(2*x)+cos(2*y))*F(t,nu)^2;
        end
    end

end

% --------------------------------------------------------------------------------------
% Decaying function (for analytical solution)
% --------------------------------------------------------------------------------------
function f = F(t,nu)

    f = exp(-2.*nu*t);

end
