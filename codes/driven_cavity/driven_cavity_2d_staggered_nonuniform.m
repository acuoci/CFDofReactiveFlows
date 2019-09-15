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
%  Code: 2D driven-cavity problem on a staggered non-uniform grid         %
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
nx=30;      % number of grid points along x
ny=30;      % number of grid points along y
deltax=2;   % stretching factor along x
deltay=2;   % stretching factor along y
Lx=1;       % length along the x axis [m]
Ly=1;       % length along the y axis [m]
nu=0.01;    % kinematic viscosity [m2/s] (if L=1 and un=1, then Re=1/nu)
tau=20;     % total time of simulation [s]

% Boundary conditions
un=1;       % north wall velocity [m/s]
us=0;       % south wall velocity [m/s]
ve=0;       % east wall velocity [m/s]
vw=0;       % west wall velocity [m/s]

% Parameters for SOR
max_iterations=10000;   % maximum number of iterations
beta=1.5;               % SOR coefficient
max_error=1e-4;         % error for convergence

% ----------------------------------------------------------------------- %
% Grid construction
% ----------------------------------------------------------------------- %
if (mod(nx,2)~=0 || mod(ny,2)~=0)
    error('Only even number of cells can be accepted (for graphical purposes only)');
end

% Grid points along x
x = zeros(nx+1,1);
for i=1:nx+1
    x(i) = 0.5*(1+tanh(deltax*((i-1)/(nx+1-1)-0.5))/tanh(deltax/2));
end
x = Lx*x;

% Grid points along y 
y = zeros(ny+1,1);
for i=1:ny+1
    y(i) = 0.5*(1+tanh(deltay*((i-1)/(ny+1-1)-0.5))/tanh(deltay/2));
end
y = Ly*y;

% Process the grid (centers and interpolation coefficients)
[xc, yc, interp_x, interp_y] = ProcessGrid(x,y);

% ----------------------------------------------------------------------- %
% Data processing
% ----------------------------------------------------------------------- %

% Time step
h2 = (x(2)-x(1))*(y(2)-y(1));       % minimum cell volume [m2]
sigma = 0.5;                        % safety factor for time step (stability)
dt_diff=h2/4/nu;                    % time step (diffusion stability) [s]
dt_conv=4*nu/un^2;                  % time step (convection stability) [s]
dt=sigma*min(dt_diff, dt_conv);     % time step (stability) [s]
nsteps=tau/dt;                      % number of steps
Re = un*Lx/nu;                      % Reynolds' number

fprintf('Time step: %f\n', dt);
fprintf(' - Diffusion:  %f\n', dt_diff);
fprintf(' - Convection: %f\n', dt_conv);
fprintf('Reynolds number: %f\n', Re);

% Grid construction
[X,Y] = meshgrid(x,y);  % MATLAB grid

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

% Fields used only for graphical post-processing purposes
uu=zeros(nx+1,ny+1);
vv=zeros(nx+1,ny+1);
pp=zeros(nx+1,ny+1);

% Allocation of vectors for non-uniform grid (see the Poisson equation function)
ae=zeros(nx+2,ny+2);    aw=zeros(nx+2,ny+2);    an=zeros(nx+2,ny+2);    as=zeros(nx+2,ny+2);
ap=zeros(nx+2,ny+2);    hx=zeros(nx+1,1);       hy=zeros(ny+1,1);

for i=2:nx+1
    
    hx(i) = x(i)-x(i-1);
    
    for j=2:ny+1

        hy(j) = y(j)-y(j-1);

        ae(i,j) = hy(j)/(xc(i+1)-xc(i));
        aw(i,j) = hy(j)/(xc(i)-xc(i-1));
        an(i,j) = hx(i)/(yc(j+1)-yc(j));
        as(i,j) = hx(i)/(yc(j)-yc(j-1));

        if (i>2)    ap(i,j) = ap(i,j)+aw(i,j); end
        if (i<nx+1) ap(i,j) = ap(i,j)+ae(i,j); end
        if (j>2)    ap(i,j) = ap(i,j)+as(i,j); end
        if (j<ny+1) ap(i,j) = ap(i,j)+an(i,j); end

    end
end

% ----------------------------------------------------------------------- %
% Solution over time
% ----------------------------------------------------------------------- %
t=0.0;
for is=1:nsteps
    
    % Boundary conditions
    u(1:nx+1,1)=2*us-u(1:nx+1,2);               % south wall
    u(1:nx+1,ny+2)=2*un-u(1:nx+1,ny+1);         % north wall
    v(1,1:ny+1)=2*vw-v(2,1:ny+1);               % west wall
    v(nx+2,1:ny+1)=2*ve-v(nx+1,1:ny+1);         % east wall
    
    % Advection-diffusion equation
    [ut, vt] = AdvectionDiffusion2D(   ut, vt, u, v, ...
                                       x, y, xc, yc, interp_x, interp_y, ...
                                       dt, nu);
    
    % Pressure equation (Poisson)
    [p, iter] = Poisson2D( p, ut, vt, ap, ae, aw, an, as, hx, hy, ...
                           dt, ...
                           beta, max_iterations, max_error);
    
    % Correct the velocity
    for i=2:nx
        for j=2:ny+1
            u(i,j)=ut(i,j) - dt*(p(i+1,j)-p(i,j))/(xc(i+1)-xc(i));
        end
    end           
    for i=2:nx+1
        for j=2:ny
            v(i,j)=vt(i,j) - dt*(p(i,j+1)-p(i,j))/(yc(j+1)-yc(j));
        end
    end       
    
    % Print on the screen
    if (mod(is,50)==1)
        fprintf('Step: %d - Time: %f - Poisson iterations: %d\n', is, t, iter);
    end
    
    % Advance time
    t=t+dt;
 
end

% ----------------------------------------------------------------------- %
% Final post-processing                                                   %
% ----------------------------------------------------------------------- %

% Field reconstruction (corners of pressure cells)
% The reconstruction is not accurate, since no interpolation factors
% are adopted
uu(1:nx+1,1:ny+1)=0.50*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
vv(1:nx+1,1:ny+1)=0.50*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
pp(1:nx+1,1:ny+1)=0.25*(p(1:nx+1,1:ny+1)+p(1:nx+1,2:ny+2)+...
                        p(2:nx+2,1:ny+1)+p(2:nx+2,2:ny+2));

for i=1:nx+1
    for j=1:ny+1
        uu(i,j)=( u(i,j+1)*(yc(j+1)-y(j)) + u(i,j)*(y(j)-yc(j)) ) / ( yc(j+1) - yc(j) );
        vv(i,j)=( v(i+1,j)*(xc(i+1)-x(i)) + v(i,j)*(x(i)-xc(i)) ) / ( xc(i+1) - xc(i) );
        pp(i,j)=( p(i+1,j+1)*(xc(i+1)-x(i))*(yc(j+1)-y(j)) + p(i+1,j)*(xc(i+1)-x(i))*(y(j)-yc(j)) + ...
                  p(i,j+1)*(x(i)-xc(i))*(yc(j+1)-y(j))     + p(i,j)*(x(i)-xc(i))*(y(j)-yc(j))  ) / ...
                ( xc(i+1) - xc(i) ) / ( yc(j+1) - yc(j) );
    end
end

% Surface map: u-velocity
subplot(231);
surface(X,Y,uu');
axis('square'); title('u'); xlabel('x'); ylabel('y');

% Surface map: v-velocity
subplot(234);
surface(X,Y,vv');
axis('square'); title('v'); xlabel('x'); ylabel('y');

% Surface map: pressure
% subplot(232);
% surface(X,Y,pp');
% axis('square'); title('pressure'); xlabel('x'); ylabel('y');

% Streamlines
subplot(233);
streamline(X,Y,uu',vv',x,y)
axis([0 Lx 0 Ly], 'square');
title('streamlines'); xlabel('x'); ylabel('y');

% Surface map: velocity vectors
subplot(236);
quiver(X,Y,uu',vv');
axis([0 1 0 1], 'square');
title('velocity vector field'); xlabel('x'); ylabel('y');

% Plot: velocity components along the horizontal middle axis
subplot(232);
plot(x,uu(:, round(ny/2)+1));
hold on;
plot(x,vv(:, round(ny/2)+1));
axis('square');
xlabel('x [m]');ylabel('velocity [m/s]'); title('velocity along x-axis');
legend('x-velocity', 'y-velocity');
hold off;

% Plot: velocity components along the horizontal middle axis
subplot(235);
plot(y,uu(round(nx/2)+1,:));
hold on;
plot(y,vv(round(nx/2)+1,:));
axis('square');
xlabel('y [m]');ylabel('velocity [m/s]'); title('velocity along y-axis');
legend('x-velocity', 'y-velocity');
hold off;


% ------------------------------------------------------------------- %
% Write velocity profiles along the centerlines for exp comparison
% ------------------------------------------------------------------- %
u_profile = uu(nx/2+1,:);
p_vertical_profile = pp(nx/2+1,:);
fileVertical = fopen('experimental_data/vertical.out','w');
for i=1:ny+1 
    fprintf(fileVertical,'%f %f %f\n', y(i), u_profile(i), p_vertical_profile(i) );
end
fclose(fileVertical);

v_profile = vv(:,ny/2+1);
p_horizontal_profile = pp(:,ny/2+1);
fileHorizontal = fopen('experimental_data/horizontal.out','w');
for i=1:nx+1
    fprintf(fileHorizontal,'%f %f %f\n',x(i), v_profile(i), p_horizontal_profile(i) );
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
axis('square'); title('u along y (centerline)'); xlabel('y'); ylabel('u'); xlim([0 Ly]);

figure;
plot(exp_v_along_x(:,1), exp_v_along_x(:,2), 'o', x, v_profile, '-');
axis('square'); title('v along x (centerline)'); xlabel('x'); ylabel('v'); xlim([0 Lx]);

figure;
plot(y(2:ny), p_vertical_profile(2:ny), '-');
axis('square'); title('p along y (centerline)'); xlabel('y'); ylabel('p'); xlim([0 Ly]);

figure;
plot(x(2:nx), p_horizontal_profile(2:nx), '-');
axis('square'); title('p along x (centerline)'); xlabel('x'); ylabel('p'); xlim([0 Lx]);


% --------------------------------------------------------------------------------------
% Poisson equation solver
% --------------------------------------------------------------------------------------
function [p, iter] = Poisson2D( p, ut, vt, ap, ae, aw, an, as, hx, hy, ...
                                dt, ...
                                beta, max_iterations, max_error)

    nx = length(hx)-1;
    ny = length(hy)-1;
    
    for iter=1:max_iterations
        
        for i=2:nx+1      
            for j=2:ny+1

                delta = ae(i,j)*p(i+1,j)+aw(i,j)*p(i-1,j) + ...
                        an(i,j)*p(i,j+1)+as(i,j)*p(i,j-1);
                S = (1/dt)*( (ut(i,j)-ut(i-1,j))*hy(j) + (vt(i,j)-vt(i,j-1))*hx(i) );
                
                p(i,j)=beta*( delta-S )/ap(i,j)+(1-beta)*p(i,j);
                
            end
        end
        
        % Estimate the error
        epsilon=0.0; 
        for i=2:nx+1
            for j=2:ny+1
                
                delta = ae(i,j)*p(i+1,j)+aw(i,j)*p(i-1,j) + ...
                        an(i,j)*p(i,j+1)+as(i,j)*p(i,j-1);
                S = (1/dt)*( (ut(i,j)-ut(i-1,j))*hy(j) + (vt(i,j)-vt(i,j-1))*hx(i) );              
                
                epsilon=epsilon+abs( p(i,j) - ( delta-S )/ap(i,j) );
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
function [ut, vt] = AdvectionDiffusion2D(   ut, vt, u, v, ...
                                            x, y, xc, yc, interp_x, interp_y, ...
                                            dt, nu)
                       
    nx = length(x)-1;
    ny = length(y)-1;
    
    % Temporary u-velocity
    for i=2:nx
        for j=2:ny+1 
            
            ue = (u(i,j)+u(i+1,j))/2;
            uw = (u(i,j)+u(i-1,j))/2;
            un = u(i,j)+(u(i,j+1)-u(i,j))*interp_y(j);
            us = u(i,j-1)+(u(i,j)-u(i,j-1))*interp_y(j-1);
            vn = v(i,j)+(v(i+1,j)-v(i,j))*interp_x(i);
            vs = v(i,j-1)+(v(i+1,j-1)-v(i,j-1))*interp_x(i);
            
            ue2 = ue^2 * (y(j)-y(j-1));
            uw2 = uw^2 * (y(j)-y(j-1));
            unv = un*vn * (xc(i+1)-xc(i));
            usv = us*vs * (xc(i+1)-xc(i));
            
            V = (xc(i+1)-xc(i)) * (y(j)-y(j-1));
            A = (ue2-uw2+unv-usv)/V;
            
            De = nu*(u(i+1,j)-u(i,j))/(x(i+1)-x(i)) * (y(j)-y(j-1));
            Dw = nu*(u(i,j)-u(i-1,j))/(x(i)-x(i-1)) * (y(j)-y(j-1));
            Dn = nu*(u(i,j+1)-u(i,j))/(yc(j+1)-yc(j)) * (xc(i+1)-xc(i));
            Ds = nu*(u(i,j)-u(i,j-1))/(yc(j)-yc(j-1)) * (xc(i+1)-xc(i));
            D = (De-Dw+Dn-Ds)/V;
            
            ut(i,j)=u(i,j)+dt*(-A+D);
            
        end
    end
    
    % Temporary v-velocity
    for i=2:nx+1
        for j=2:ny 
            
            vn = (v(i,j)+v(i,j+1))/2;
            vs = (v(i,j)+v(i,j-1))/2;
            ve = v(i,j)+(v(i+1,j)-v(i,j))*interp_x(i);
            vw = v(i-1,j)+(v(i,j)-v(i-1,j))*interp_x(i-1);
            ue = u(i,j)+(u(i,j+1)-u(i,j))*interp_y(j);
            uw = u(i-1,j)+(u(i-1,j+1)-u(i-1,j))*interp_y(j);
            
            vn2 = vn^2 * (x(i)-x(i-1));
            vs2 = vs^2 * (x(i)-x(i-1));
            veu = ve*ue * (yc(j+1)-yc(j));
            vwu = vw*uw * (yc(j+1)-yc(j));
            
            V = (x(i)-x(i-1)) * (yc(j+1)-yc(j));
            A = (vn2 - vs2 + veu - vwu)/V;
            
            De = nu*(v(i+1,j)-v(i,j))/(xc(i+1)-xc(i)) * (yc(j+1)-yc(j));
            Dw = nu*(v(i,j)-v(i-1,j))/(xc(i)-xc(i-1)) * (yc(j+1)-yc(j));
            Dn = nu*(v(i,j+1)-v(i,j))/(y(j+1)-y(j)) * (x(i)-x(i-1));
            Ds = nu*(v(i,j)-v(i,j-1))/(y(j)-y(j-1)) * (x(i)-x(i-1));
            D = (De-Dw+Dn-Ds)/V;
            
            vt(i,j)=v(i,j)+dt*(-A+D);
            
        end
    end
    
end


% --------------------------------------------------------------------------------------
% Process the grid
% --------------------------------------------------------------------------------------
function [xc, yc, interp_x, interp_y] = ProcessGrid(x,y)

    nx = length(x)-1;
    ny = length(y)-1;

    % Centers of cells along x
    xc = zeros(nx+2,1);
    for i=1:nx
        xc(i+1) = 0.5*(x(i)+x(i+1));
    end
    xc(1) = xc(2) - (xc(3)-xc(2));
    xc(ny+2) = xc(nx+1) + (xc(nx+1)-xc(nx));

    % Centers of cells along y
    yc = zeros(ny+2,1);
    for j=1:ny
        yc(j+1) = 0.5*(y(j)+y(j+1));
    end
    yc(1) = yc(2) - (yc(3)-yc(2));
    yc(ny+2) = yc(ny+1) + (yc(ny+1)-yc(ny));

    % Interpolation coefficients along x
    interp_x = zeros(nx+1,1);
    for i=1:nx+1
        interp_x(i) = (x(i)-xc(i))/(xc(i+1)-xc(i));
    end

    % Interpolation coefficients along y
    interp_y = zeros(ny+1,1);
    for j=1:ny+1
        interp_y(j) = (y(j)-yc(j))/(yc(j+1)-yc(j));
    end

end