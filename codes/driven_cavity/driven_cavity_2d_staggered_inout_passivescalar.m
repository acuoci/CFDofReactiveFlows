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
%        inlet and outlet sections and a passive scalar equation          %
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
nx=50;      % number of (physical) cells along x
ny=nx;      % number of (physical) cells along y
L=1;        % length [m]
nu=1e-2;    % kinematic viscosity [m2/s] (if L=1 and un=1, then Re=1/nu)
tau=10;     % total time of simulation [s]
alpha=0.01; % [PASSIVE] diffusion coefficient passive scalar [m2/s]
phi0=0;     % [PASSIVE] initial value of passive scalar

% Boundary conditions (velocity)
un=2;       % north wall velocity [m/s]
us=0;       % south wall velocity [m/s]
ve=0;       % east wall velocity [m/s]
vw=0;       % west wall velocity [m/s]
uin=1;      % [INOUT] inlet velocity on the west side [m/s]

% Boundary conditions for passive scalar
% We assume that the passive scalar enters the domain through the inlet
% section with a fixed, user-defined value
phi_inlet=1;  % [PASSIVE] inlet value

% Parameters for SOR
max_iterations=10000;   % maximum number of iterations
beta=1.5;               % SOR coefficient
max_error=1e-4;         % error for convergence

% [INOUT] Inlet section (west side)
nin_start = 1/2*(ny+2);         % first cell index 
nin_end = 3/4*(ny+2)+1;           % last cell index

% [INOUT] Outlet section (east side)
nout_start = 1/4*(ny+2);        % first cell index 
nout_end = 1/2*(ny+2)+1;          % last cell index

% ----------------------------------------------------------------------- %
% Data processing
% ----------------------------------------------------------------------- %
if (mod(nx,2)~=0 || mod(ny,2)~=0)
    error('Only even number of cells can be accepted (for graphical purposes only)');
end

% Process the grid
h=L/nx;                           % grid step (uniform grid) [m]

% [INOUT] Inlet/Outlet section areas
Ain = h*(nin_end-nin_start+1);      % inlet section area [m]
Aout = h*(nout_end-nout_start+1);   % outlet section area [m]

% [INOUT] Estimated max velocity
umax=max([abs(un),abs(us),abs(ve),abs(vw),...
            uin,uin*Ain/Aout]);     % maximum velocity [m/s]

% Time step
sigma = 0.50;                       % safety factor for time step (stability)
dt_diff=h^2/4/nu;                   % time step (diffusion stability) [s]
dt_conv=4*nu/umax^2;                % time step (convection stability) [s]
dt=sigma*min(dt_diff, dt_conv);     % time step (stability) [s]
nsteps=tau/dt;                      % number of steps
Re = umax*L/nu;                     % Reynolds' number

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
phi=zeros(nx+2,ny+2)+phi0;  % [PASSIVE] passive scalar field

% Temporary velocity fields
ut=zeros(nx+1,ny+2);
vt=zeros(nx+2,ny+1);

% Fields used only for graphical post-processing purposes
uu=zeros(nx+1,ny+1);
vv=zeros(nx+1,ny+1);
pp=zeros(nx+1,ny+1);
phiphi=zeros(nx+1,ny+1);    % [PASSIVE] passive scalar field

% Coefficient for pressure equation
gamma=zeros(nx+2,ny+2)+1/4;
gamma(2,3:ny)=1/3;gamma(nx+1,3:ny)=1/3;gamma(3:nx,2)=1/3;gamma(3:nx,ny+1)=1/3;
gamma(2,2)=1/2;gamma(2,ny+1)=1/2;gamma(nx+1,2)=1/2;gamma(nx+1,ny+1)=1/2;

% [INOUT] Correction of gamma coefficients for taking into account out sections
gamma(nx+1,nout_start:nout_end) = 1/4;

% [INOUT] Initial conditions
u(1,nin_start:nin_end) = uin;               % inlet section: fixed velocity [m/s]
u(nx+1,nout_start:nout_end) = uin;          % outlet section: fixed velocity [m/s]
u(2:nx,2:ny+1) = uin;                       % internal points: fixed velocity [m/s]
ut = u;                                     % temporary velocity [m/s]

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
    
    % [INOUT] Over-writing inlet conditions    
    u(1,nin_start:nin_end) = uin;               % fixed velocity [m/s]
    
    % [INOUT] Over-writing outlet conditions
    u(nx+1,nout_start:nout_end) = u(nx,nout_start:nout_end);    % zero-gradient     
    v(nx+2,nout_start:nout_end) = v(nx+1,nout_start:nout_end);  % zero-gradient    
    
    % Advection-diffusion equation (predictor)
    [ut, vt] = AdvectionDiffusion2D( ut, vt, u, v, nx, ny, h, dt, nu);
    
    % [INOUT] Update boundary conditions for temporary velocity
    ut(1,nin_start:nin_end) = u(1,nin_start:nin_end);            % fixed velocity [m/s]
    ut(nx+1,nout_start:nout_end) = u(nx+1,nout_start:nout_end);  % zero-gradient
    vt(nx+2,nout_start:nout_end) = v(nx+2,nout_start:nout_end);  % zero-gradient 
    
    % Pressure equation (Poisson)
    [p, iter] = Poisson2D( p, ut, vt, gamma, nx, ny, h, dt, ...
                           beta, max_iterations, max_error );
    
    % Correction on the velocity
    u(2:nx,2:ny+1)=ut(2:nx,2:ny+1)-(dt/h)*(p(3:nx+1,2:ny+1)-p(2:nx,2:ny+1));
    v(2:nx+1,2:ny)=vt(2:nx+1,2:ny)-(dt/h)*(p(2:nx+1,3:ny+1)-p(2:nx+1,2:ny));
    
    % [INOUT] Correction on outlet to ensure conservation of mass
    u(nx+1,nout_start:nout_end)=ut(nx+1,nout_start:nout_end) - ...
                                (dt/h)*(p(nx+2,nout_start:nout_end)-p(nx+1,nout_start:nout_end));
    
    % [INOUT] Because of numerical errors in the solution of the equations,
    %         the overall continuity equation, i.e. the conservation of
    %         mass cannot be guaranteed. It is better to correct the outlet
    %         velocity in order to force conservation of mass
    Qin = mean(u(1,nin_start:nin_end))*Ain;         % inlet flow rate [m2/s]
    Qout = mean(u(nx+1,nout_start:nout_end))*Aout;  % outlet flow rate [m2/s]    
    if (abs(Qout)>1.e-6)
        u(nx+1,nout_start:nout_end) = u(nx+1,nout_start:nout_end)*abs(Qin/Qout);
    end
    
    % ----------------------------------------------------------------------- %
    % [PASSIVE] Passive scalar advection-diffusion equation
    % ----------------------------------------------------------------------- %

    % [PASSIVE] Boundary conditions (zero-gradient)
    phi(2:nx+1,1)=phi(2:nx+1,2);
    phi(2:nx+1,ny+2)=phi(2:nx+1,ny+1);
    phi(1,2:ny+1)=phi(2,2:ny+1);
    phi(nx+2,2:ny+1)=phi(nx+1,2:ny+1);
    
    % [PASSIVE] Boundary conditions at inlet section (Dirichlet)
    phi(1,nin_start:nin_end)=2*phi_inlet-phi(2,nin_start:nin_end);

    % [PASSIVE] Update passive scalar solution
    phi = AdvectionDiffusion2DPassiveScalar( phi, u, v, nx, ny, h, dt, alpha );
    
    % [PASSIVE] Print on the screen
    if (mod(is,50)==1)
        fprintf( 'Step: %d - Time: %f - Poisson iterations: %d - dQ: %f%%\n', ...
                 is, t, iter, (Qout-Qin)/Qin*100. );
    end
    
    % Advance in time
    t=t+dt;
    
    % ----------------------------------------------------------------------- %
    % On the fly post-processing                                                   %
    % ----------------------------------------------------------------------- %
    if (mod(is,10)==1)
        
        % [PASSIVE] Adjust values on the corners (graphical purposes only)
        phi(1,1)=1/3*(phi(1,2)+phi(2,2)+phi(2,1));
        phi(nx+2,1)=1/3*(phi(nx+2,2)+phi(nx+1,1)+phi(nx+1,nx+1));
        phi(nx+2,1)=1/3*(phi(nx+2,2)+phi(nx+1,1)+phi(nx+1,2));
        phi(1,ny+2)=1/3*(phi(2,ny+1)+phi(1,ny+1)+phi(2,ny+1));
        phi(nx+2,ny+2)=1/3*(phi(nx+1,ny+2)+phi(nx+2,ny+1)+phi(nx+1,ny+1));


        % Field reconstruction
        phiphi(1:nx+1,1:ny+1)=0.25*(phi(1:nx+1,1:ny+1)+phi(2:nx+2,1:ny+1) + ...
                                    phi(1:nx+1,2:ny+2)+phi(2:nx+2,2:ny+2));

        % Surface: passive scalar
        surface(X,Y,phiphi');
        axis('square'); xlabel('x [m]');ylabel('y [m]'); title('passive scalar');
        
        pause(0.01);
        
    end
 
end

% ----------------------------------------------------------------------- %
% Final post-processing                                                   %
% ----------------------------------------------------------------------- %

% [PASSIVE] Adjust values on the corners (graphical purposes only)
phi(1,1)=1/3*(phi(1,2)+phi(2,2)+phi(2,1));
phi(nx+2,1)=1/3*(phi(nx+2,2)+phi(nx+1,1)+phi(nx+1,nx+1));
phi(nx+2,1)=1/3*(phi(nx+2,2)+phi(nx+1,1)+phi(nx+1,2));
phi(1,ny+2)=1/3*(phi(2,ny+1)+phi(1,ny+1)+phi(2,ny+1));
phi(nx+2,ny+2)=1/3*(phi(nx+1,ny+2)+phi(nx+2,ny+1)+phi(nx+1,ny+1));


% Field reconstruction
uu(1:nx+1,1:ny+1)=0.50*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
vv(1:nx+1,1:ny+1)=0.50*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
phiphi(1:nx+1,1:ny+1)=0.25*(phi(1:nx+1,1:ny+1)+phi(2:nx+2,1:ny+1) + ...
                            phi(1:nx+1,2:ny+2)+phi(2:nx+2,2:ny+2));
pp(1:nx+1,1:ny+1)=0.25*(p(1:nx+1,1:ny+1)+p(1:nx+1,2:ny+2)+...
                        p(2:nx+2,1:ny+1)+p(2:nx+2,2:ny+2));
                    
% Surface map: u-velocity
subplot(231);
surface(X,Y,uu');
axis('square'); title('u'); xlabel('x'); ylabel('y');

% Surface map: v-velocity
subplot(232);
surface(X,Y,vv');
axis('square'); title('v'); xlabel('x'); ylabel('y');

% Streamlines
subplot(233);
sy = 0:2*h:L;
sx = L*0.01*ones(length(sy),1);
streamline(X,Y,uu',vv',sx,sy)
axis([0 L 0 L], 'square');
title('streamlines'); xlabel('x'); ylabel('y');

% Surface map: velocity vectors
subplot(234);
quiver(X,Y,uu',vv');
axis([0 L 0 L], 'square');
title('velocity vector field'); xlabel('x'); ylabel('y');

% Surface: passive scalar
subplot(235);
surface(X,Y,phiphi');
axis('square'); xlabel('x [m]');ylabel('y [m]'); title('passive scalar');

% Plot: velocity components along the horizontal middle axis
subplot(236);
plot(x,phiphi(:, round(ny/2)));
hold on;
plot(y,phiphi(round(ny/2),:));
axis('square');
xlabel('x or y [m]');ylabel('passive scalar'); title('passive scalar (centerlines)');
legend('along x', 'along y');
hold off;


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
