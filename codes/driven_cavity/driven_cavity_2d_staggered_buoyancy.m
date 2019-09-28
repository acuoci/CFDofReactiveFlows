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
%        with inclusion of a buoyancy term in the momentum equation       %
%        and a temperature equation as a passive scalar equation          %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% ----------------------------------------------------------------------- %
% User data
% ----------------------------------------------------------------------- %

% Only even numbers of cells are acceptable
nx=20;                % number of (physical) cells along x
ny=nx;                % number of (physical) cells along y
L=0.01;               % length [m]
nu=1.1e-5;            % kinematic viscosity of air [m2/s]
Beta=3.4e-3;          % volumetric expansion coeff. of air at 20C [1/K]
alpha=21e-6;          % thermal diffusion coefficient of air [m2/s]
T0=20;                % initial (and reference) value of temperature [C]
gx=0;                 % gravity vector (x component) [m/s2]
gy=-9.81;             % gravity vector (y component) []m/s2
tau=5;                % total time of simulation [s]
sigma = 0.01;         % safety factor for time step (stability)
maxsteps = 1000000;   % maximum number of time steps
Te=T0;                % east wall temperature (Dirichlet)
Tw=T0+10;             % west wall temperature (Dirichlet)

% Boundary conditions (velocities)
un=0;       % north wall velocity [m/s]
us=0;       % south wall velocity [m/s]
ve=0;       % east wall velocity [m/s]
vw=0;       % west wall velocity [m/s]

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

% Grashof number
g=sqrt(gx^2+gy^2);                        % gravity vector (module) [m2/s2]
Tmax=max([Te,Tw,T0]);                     % maximum temperature [C]
Tmin=min([Te,Tw,T0]);                     % minimum temperature [C]
Gr = g*Beta*(Tmax-Tmin)*L^3/nu^2;         % Grashof's number

fprintf('Grashof number: %f\n', Gr);

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
T=zeros(nx+2,ny+2);

% Temporary velocity fields
ut=zeros(nx+1,ny+2);
vt=zeros(nx+2,ny+1);

% Fields used only for graphical post-processing purposes
uu=zeros(nx+1,ny+1);
vv=zeros(nx+1,ny+1);
TT=zeros(nx+1,ny+1);
omega=zeros(nx+1,ny+1);

% Coefficient for pressure equation
c=zeros(nx+2,ny+2)+1/4;
c(2,3:ny)=1/3;c(nx+1,3:ny)=1/3;c(3:nx,2)=1/3;c(3:nx,ny+1)=1/3;
c(2,2)=1/2;c(2,ny+1)=1/2;c(nx+1,2)=1/2;c(nx+1,ny+1)=1/2;


% ----------------------------------------------------------------------- %
% Solution over time
% ----------------------------------------------------------------------- %
t=0.0;
for is=1:maxsteps

    % Choosing the time step
    umax = max(max(abs(u)));                    % maximum velocity [m/s]
    dt_conv=4*min(alpha,nu)/umax^2;             % time step (convection stability)
    dt_diff=h^2/4/max(nu,alpha);                % time step (diffusion stability)
    dt=sigma*min(dt_diff, dt_conv);             % time step (stability)

    % Boundary conditions (velocity)
    u(1:nx+1,1)=2*us-u(1:nx+1,2);
    u(1:nx+1,ny+2)=2*un-u(1:nx+1,ny+1);
    v(1,1:ny+1)=2*vw-v(2,1:ny+1);
    v(nx+2,1:ny+1)=2*ve-v(nx+1,1:ny+1);

    % Advection-diffusion equation (predictor)
    [ut, vt] = AdvectionDiffusion2D( ut, vt, u, v, nx, ny, h, dt, nu, ...
                                     gx, gy, Beta, T, T0);

    % Pressure equation (Poisson)
    [p, iter] = Poisson2D( p, ut, vt, c, nx, ny, h, dt, ...
                           beta, max_iterations, max_error);

    % Correct the velocity
    u(2:nx,2:ny+1)=ut(2:nx,2:ny+1)-(dt/h)*(p(3:nx+1,2:ny+1)-p(2:nx,2:ny+1));
    v(2:nx+1,2:ny)=vt(2:nx+1,2:ny)-(dt/h)*(p(2:nx+1,3:ny+1)-p(2:nx+1,2:ny));

    % Print time step on the screen
    if (mod(is,2000)==1)
        fprintf('Step: %d - Time: %f - Poisson: %d - Max v: %f - dt: %f\n', is, t, iter, umax, dt);
    end
    
    % ----------------------------------------------------------------------- %
    % Temperature advection-diffusion equation
    % ----------------------------------------------------------------------- %

    % Boundary conditions (temperature)
    T(2:nx+1,1)=T(2:nx+1,2);                % zero-gradient
    T(2:nx+1,ny+2)=T(2:nx+1,ny+1);          % zero-gradient
    T(1,2:ny+1)=2*Tw-T(2,2:ny+1);           % Dirichlet
    T(nx+2,2:ny+1)=2*Te-T(nx+1,2:ny+1);     % Dirichlet

    % Update temperature solution
    T = AdvectionDiffusion2DPassiveScalar( T, u, v, nx, ny, h, dt, alpha );
    
  
    % Advance time
    t=t+dt;
 
    % ----------------------------------------------------------------------- %
    % Update graphical output
    % ----------------------------------------------------------------------- %
    if (mod(is,2000)==1)
        
        % Post-processing (reconstructi on the corners of pressure cells)
        uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
        vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
        TT(1:nx+1,1:ny+1)=0.25*(T(1:nx+1,1:ny+1)+T(2:nx+2,1:ny+1) + ...
                                T(1:nx+1,2:ny+2)+T(2:nx+2,2:ny+2));
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
        
         % Surface map: velocity vectors
        subplot(234);
        quiver(X,Y,uu',vv');
        axis([0 L 0 L], 'square');
        title('velocity vector field'); xlabel('x'); ylabel('y');

        % Contour: temperature
        subplot(235);
        contour(X,Y,TT');
        axis('square'); xlabel('x [m]');ylabel('y [m]'); title('temperature contour');

        % Plot: temperature components along the horizontal middle axis
        subplot(236);
        plot(x,TT(:, round(ny/2)));
        hold on;
        plot(y,TT(round(ny/2),:));
        axis('square');
        xlabel('x or y [m]');ylabel('temperature'); title('temperature (centerlines)');
        legend('along x', 'along y'); ylim([min([Te,Tw]), max([Te,Tw])]);
        hold off;

        pause(0.001)
    
    end
    
    if (t>=tau) break; end

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
function [ut, vt] = AdvectionDiffusion2D( ut, vt, u, v, nx, ny, h, dt, nu, ...
                                          gx, gy, Beta, T, T0)
                            
    % Temporary u-velocity
    for i=2:nx
        for j=2:ny+1 
            
            % Temperature in the center of u-cell (linear interpolation)
            Tij = 0.50*(T(i+1,j)+T(i,j));
            
            ue2 = 0.25*( u(i+1,j)+u(i,j) )^2;
            uw2 = 0.25*( u(i,j)+u(i-1,j) )^2;
            unv = 0.25*( u(i,j+1)+u(i,j) )*( v(i+1,j)+v(i,j) );
            usv = 0.25*( u(i,j)+u(i,j-1) )*( v(i+1,j-1)+v(i,j-1) );
            
            A = (ue2-uw2+unv-usv)/h;
            D = (nu/h^2)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j));
            
            ut(i,j)=u(i,j)+dt*(-A+D) ... 
                    -dt*gx*Beta*(Tij-T0);
            
        end
    end
    
    % Temporary v-velocity
    for i=2:nx+1
        for j=2:ny 
            
            % Temperature in the center of v-cell (linear interpolation)
            Tij = 0.50*(T(i,j+1)+T(i,j));
            
            vn2 = 0.25*( v(i,j+1)+v(i,j) )^2;
            vs2 = 0.25*( v(i,j)+v(i,j-1) )^2;
            veu = 0.25*( u(i,j+1)+u(i,j) )*( v(i+1,j)+v(i,j) );
            vwu = 0.25*( u(i-1,j+1)+u(i-1,j) )*( v(i,j)+v(i-1,j) );
            A = (vn2 - vs2 + veu - vwu)/h;
            D = (nu/h^2)*(v(i+1,j)+v(i-1,j)+v(i,j+1)+v(i,j-1)-4*v(i,j));
            
            vt(i,j)=v(i,j)+dt*(-A+D) ... 
                    -dt*gy*Beta*(Tij-T0);
            
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
