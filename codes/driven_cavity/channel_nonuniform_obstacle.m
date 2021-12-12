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
%   Copyright(C) 2021 Alberto Cuoci                                       %
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
%  Code: channel problem on a staggered grid                              %
%        derived from 2d driven cavity + inlet and outlet sections        %
%        grid with different step sizes along the two directions          %
%        and internal obstacle                                            %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% ----------------------------------------------------------------------- %
% User data
% ----------------------------------------------------------------------- %

% Only even numbers of cells are acceptable
nx=200;     % number of (physical) cells along x
ny=50;      % number of (physical) cells along y 
Lx=5;       % length [m]
Ly=1;       % length [m]
nu=1e-2;    % kinematic viscosity [m2/s] (if L=1 and un=1, then Re=1/nu)
tau=20;     % total time of simulation [s]

% Boundary conditions
un=0;       % north wall velocity [m/s]
us=0;       % south wall velocity [m/s]
ve=0;       % east wall velocity [m/s]
vw=0;       % west wall velocity [m/s]
uin=1;      % [INOUT] inlet velocity on the west side [m/s]

% Parameters for Poisson solver
solver_type = 'SOR';    % Options: SOR | Direct | GMRES
max_iterations=50000;   % maximum number of iterations
beta=1.8;               % SOR coefficient
max_error=1e-4;         % error for convergence

% [INOUT] Inlet section (west side)
nin_start = 1/2*(ny+2);                 % first cell index 
nin_end = 3/4*(ny+2)+1;                 % last cell index

% [INOUT] Outlet section (east side)
nout_start = 1/4*(ny+2);                % first cell index 
nout_end = 1/2*(ny+2)+1;                % last cell index

% [OBST] Obstacle definition: rectangle with base xs:xe and height ys:ye
% Example: square obstacle in the center
xs=(nx+2)/2-12; xe=(nx+2)/2+12;
ys=(ny+2)/2-10; ye=(ny+2)/2+10;

% ----------------------------------------------------------------------- %
% Data processing
% ----------------------------------------------------------------------- %
if (mod(nx,2)~=0 || mod(ny,2)~=0)
    error('Only even number of cells can be accepted (for graphical purposes only)');
end
if (mod(ny+2,4)~=0)
    error('The total number of cells along y (i.e. ny+2) must be divisible by 4');
end

% Process the grid
hx=Lx/nx;                           % grid step (uniform grid) [m]
hy=Ly/ny;                           % grid step (uniform grid) [m]


% [INOUT] Inlet/Outlet section areas
Ain = hy*(nin_end-nin_start+1);      % inlet section area [m]
Aout = hy*(nout_end-nout_start+1);   % outlet section area [m]

% [INOUT] Estimated max velocity
umax=max([abs(un),abs(us),abs(ve),abs(vw),...
            uin,uin*Ain/Aout]);     % maximum velocity [m/s]

% Time step
sigma = 0.50;                       % safety factor for time step (stability)
dt_diff=min(hx,hy)^2/4/nu;                   % time step (diffusion stability) [s]
dt_conv=4*nu/umax^2;                % time step (convection stability) [s]
dt=sigma*min(dt_diff, dt_conv);     % time step (stability) [s]
nsteps=tau/dt;                      % number of steps
Re = umax*Ly/nu;                    % Reynolds' number

fprintf('Time step: %f\n', dt);
fprintf(' - Diffusion:  %f\n', dt_diff);
fprintf(' - Convection: %f\n', dt_conv);
fprintf('Reynolds number: %f\n', Re);

% Grid construction
x=0:hx:Lx;                         % grid coordinates (x axis)
y=0:hy:Ly;                         % grid coordinates (y axis)
[X,Y] = meshgrid(x,y);           % MATLAB grid

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

% Coefficient for pressure equation
gamma=zeros(nx+2,ny+2)+hx*hy/(2*hx^2+2*hy^2);       % internal cells
gamma(2,3:ny)=hx*hy/(2*hx^2+hy^2);                  % west cells
gamma(nx+1,3:ny)=hx*hy/(2*hx^2+hy^2);               % east cells
gamma(3:nx,2)=hx*hy/(hx^2+2*hy^2);                  % south cells
gamma(3:nx,ny+1)=hx*hy/(hx^2+2*hy^2);               % north cells
gamma(2,2)=hx*hy/(hx^2+hy^2); gamma(2,ny+1)=hx*hy/(hx^2+hy^2);          % corner cells
gamma(nx+1,2)=hx*hy/(hx^2+hy^2); gamma(nx+1,ny+1)=hx*hy/(hx^2+hy^2);    % corner cells

% [INOUT] Correction of gamma coefficients for taking into account in sections
% In inlet sections the velocity is prescribed, so there is no correction
% Cells belonging to inlet sections behave like cells belonging to a wall
gamma(nx+1,nout_start:nout_end) = hx*hy/(2*hx^2+hy^2);

% [INOUT] Correction of gamma coefficients for taking into account out sections
% In outlet section in principle we can have a correction, since we do not
% know the outlet velocity. So cells on the oultet section behave like
% internal cells
gamma(nx+1,nout_start:nout_end) = hx*hy/(2*hx^2+2*hy^2);

% [OBST] 
[flagu, flagv, flagp] = ObstaclePreProcessing(nx,ny, xs,xe, ys,ye);

% [OBST] Correction of gamma close to the obstacle
gamma(xs-1,ys:ye)=hx*hy/(2*hx^2+hy^2);
gamma(xe+1,ys:ye)=hx*hy/(2*hx^2+hy^2);
gamma(xs:xe,ys-1)=hx*hy/(hx^2+2*hy^2);
gamma(xs:xe,ye+1)=hx*hy/(hx^2+2*hy^2);

% [INOUT] Initial conditions (v component is automatically set equal to 0)
u(1,nin_start:nin_end) = uin;               % inlet section: fixed velocity [m/s]
u(nx+1,nout_start:nout_end) = uin;          % outlet section: fixed velocity [m/s]
for i=2:nx
        for j=2:ny+1 
            if (flagu(i,j)==0)
                u(i,j) = uin;
            end
        end
end
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

    % [OBST] Parallel velocities along the obstacle walls
    uwall = 0;
    u(xs-1:xe,ye)=2*uwall-u(xs-1:xe,ye+1);     % north
    u(xs-1:xe,ys)=2*uwall-u(xs-1:xe,ys-1);     % south
    v(xs,ys-1:ye)=2*uwall-v(xs-1,ys-1:ye);     % west
    v(xe,ys-1:ye)=2*uwall-v(xe+1,ys-1:ye);     % east    
    
    % [INOUT] Over-writing inlet conditions    
    u(1,nin_start:nin_end) = uin;               % fixed velocity [m/s]
    
    % [INOUT] Over-writing outlet conditions
    u(nx+1,nout_start:nout_end) = u(nx,nout_start:nout_end);    % zero-gradient     
    v(nx+2,nout_start:nout_end) = v(nx+1,nout_start:nout_end);  % zero-gradient    
    
    % [OBST] Advection-diffusion equation (predictor)
    [ut, vt] = AdvectionDiffusion2D( ut, vt, u, v, flagu, flagv, nx, ny, hx, hy, dt, nu);
    
    % [INOUT] Update boundary conditions for temporary velocity
    ut(1,nin_start:nin_end) = u(1,nin_start:nin_end);            % fixed velocity [m/s]
    ut(nx+1,nout_start:nout_end) = u(nx+1,nout_start:nout_end);  % zero-gradient
    vt(nx+2,nout_start:nout_end) = v(nx+2,nout_start:nout_end);  % zero-gradient 
    
    % [OBST] Pressure equation (Poisson)
    [p, iter] = Poisson2D( p, flagp, ut, vt, gamma, nx, ny, hx, hy, dt, ...
                            beta, max_iterations, max_error, solver_type);

    % [OBST] Correct the velocity (only outside the obstacle)
    for i=2:nx
        for j=2:ny+1
            if (flagu(i,j)==0)
                u(i,j)=ut(i,j)-(dt/hx)*(p(i+1,j)-p(i,j));
            end
        end
    end
    for i=2:nx+1
        for j=2:ny
            if (flagv(i,j)==0)
                v(i,j)=vt(i,j)-(dt/hy)*(p(i,j+1)-p(i,j));
            end
        end
    end
    
    % [INOUT] Correction on outlet to ensure conservation of mass
    u(nx+1,nout_start:nout_end)=ut(nx+1,nout_start:nout_end) - ...
                                (dt/hx)*(p(nx+2,nout_start:nout_end)-p(nx+1,nout_start:nout_end));
    
    % [INOUT] Because of numerical errors in the solution of the equations,
    %         the overall continuity equation, i.e. the conservation of
    %         mass cannot be guaranteed. It is better to correct the outlet
    %         velocity in order to force conservation of mass
    Qin = mean(u(1,nin_start:nin_end))*Ain;         % inlet flow rate [m2/s]
    Qout = mean(u(nx+1,nout_start:nout_end))*Aout;  % outlet flow rate [m2/s]    
    if (abs(Qout)>1.e-6)
        u(nx+1,nout_start:nout_end) = u(nx+1,nout_start:nout_end)*abs(Qin/Qout);
    end
    
    % Print on the screen
    if (mod(is,50)==1)
        fprintf( 'Step: %d - Time: %f - Poisson iterations: %d - dQ: %f%%\n', ...
                 is, t, iter, (Qout-Qin)/Qin*100. );
    end
    
    % Advance in time
    t=t+dt;
 
end

% ----------------------------------------------------------------------- %
% Final post-processing                                                   %
% ----------------------------------------------------------------------- %

% [OBST] Field reconstruction
[uu, vv, pp] = ReconstructFields(u, v, p, nx, ny, flagu, flagv, flagp);

% [OBST] Obstacle detailes
position = [hx*(xs-2),hy*(ys-2),hx*(xe-xs+1),hy*(ye-ys+1)];
color = [0.85 0.85 0.85];

% Surface map: u-velocity
subplot(231);
surface(X,Y,uu');
axis('square'); title('u'); xlabel('x'); ylabel('y'); shading interp;
rectangle( 'Position', position, 'FaceColor', color);
      
% Surface map: v-velocity
subplot(234);
surface(X,Y,vv');
axis('square'); title('v'); xlabel('x'); ylabel('y'); shading interp;
rectangle( 'Position', position, 'FaceColor', color);
      
% Surface map: pressure
% subplot(232);
% surface(X,Y,pp');
% axis('square'); title('pressure'); xlabel('x'); ylabel('y');
% rectangle( 'Position', position, 'FaceColor', color);
      
% Streamlines
subplot(233);
sx = [0:Lx/10:Lx 0:Lx/100:Lx/10];
sy = [0:Ly/10:Ly Ly/2:Ly/20:Ly];
streamline(X,Y,uu',vv',sx,sy)
axis([0 Lx 0 Ly], 'square');
title('streamlines'); xlabel('x'); ylabel('y');
rectangle( 'Position', position, 'FaceColor', color);
      
% Surface map: velocity vectors
subplot(236);
quiver(X,Y,uu',vv');
axis([0 Lx 0 Ly], 'square');
title('velocity vector field'); xlabel('x'); ylabel('y');
rectangle( 'Position', position, 'FaceColor', color);

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


% --------------------------------------------------------------------------------------
% [OBST] Poisson equation solver
% --------------------------------------------------------------------------------------
function [p, iter] = Poisson2D( p, flagp, ut, vt, gamma, nx, ny, hx, hy, dt, ...
                                beta, max_iterations, max_error, solver_type)

    % SOR solver
    if (strcmp(solver_type,'SOR'))

        % Main loop
        for iter=1:max_iterations
            
            for i=2:nx+1
                for j=2:ny+1
    
                    if (flagp(i,j)==0)
    
                        delta = hy/hx*(p(i+1,j)+p(i-1,j))+hx/hy*(p(i,j+1)+p(i,j-1));
                        S = 1/dt*( hy*(ut(i,j)-ut(i-1,j)) + hx*(vt(i,j)-vt(i,j-1)));
                        p(i,j)=beta*gamma(i,j)*( delta-S )+(1-beta)*p(i,j);
    
                    end
                end
            end
            
            % Estimate the error
            epsilon=0.0;
            count_cells = 0;
            for i=2:nx+1
                for j=2:ny+1
    
                    if (flagp(i,j)==0)
    
                        delta = hy/hx*(p(i+1,j)+p(i-1,j))+hx/hy*(p(i,j+1)+p(i,j-1));
                        S = 1/dt*( hy*(ut(i,j)-ut(i-1,j)) + hx*(vt(i,j)-vt(i,j-1)));           
                        epsilon=epsilon+abs( p(i,j) - gamma(i,j)*( delta-S ) );

                        count_cells = count_cells+1;
    
                    end
                end
            end
            epsilon = epsilon / count_cells;
            
            % Check the error
            if (epsilon <= max_error) % stop if converged
                break;
            end 
            
        end

    else    % Direct or GMRES solvers

        ne = (nx+2)*(ny+2);
        b = zeros(ne,1);

        % Fill main diagonal
        counter = 1;
        for i=1:ne
            I(counter) = i; J(counter) = i; V(counter) = 1.; counter = counter+1;
        end
    
        % Fill equations
        for i=2:nx+1
            for j=2:ny+1
                
                if (flagp(i,j)==0)
    
                        k = (nx+2)*(j-1) + i;
    
                        I(counter) = k; J(counter) = k+1; V(counter) = -gamma(i,j)*hy/hx; counter = counter+1;
                        I(counter) = k; J(counter) = k-1; V(counter) = -gamma(i,j)*hy/hx; counter = counter+1;
                        I(counter) = k; J(counter) = k+(nx+2); V(counter) = -gamma(i,j)*hx/hy; counter = counter+1;
                        I(counter) = k; J(counter) = k-(nx+2); V(counter) = -gamma(i,j)*hx/hy; counter = counter+1;

                        b(k) = -gamma(i,j)*(1/dt)*(hy*(ut(i,j)-ut(i-1,j))+hx*(vt(i,j)-vt(i,j-1)));
    
                end
            
            end
        end
    
        M = sparse(I,J,V,ne,ne);
    
        if (strcmp(solver_type,'Direct'))
    
            p = M\b;
            iter = 0;
    
        elseif (strcmp(solver_type,'GMRES'))
    
            tol = 1e-6;
            maxit = 10;
            [p,~,~,iter] = gmres(M,b,[],tol,maxit,[],[],p(:));
            iter=iter(2);
    
        end
    
        p = reshape(p,[nx+2 ny+2]);
    
    end

end


% --------------------------------------------------------------------------------------
% [OBST] Advection-diffusion equation
% --------------------------------------------------------------------------------------
function [ut, vt] = AdvectionDiffusion2D( ut, vt, u, v, flagu, flagv, nx, ny, hx, hy, dt, nu)
                       
    % Temporary u-velocity
    for i=2:nx
        for j=2:ny+1 
            
            if (flagu(i,j)==0)
    
                ue = (u(i,j)+u(i+1,j))/2;
                uw = (u(i,j)+u(i-1,j))/2;
                un = u(i,j)+(u(i,j+1)-u(i,j))*0.50;
                us = u(i,j-1)+(u(i,j)-u(i,j-1))*0.50;
                vn = v(i,j)+(v(i+1,j)-v(i,j))*0.50;
                vs = v(i,j-1)+(v(i+1,j-1)-v(i,j-1))*0.50;
                
                ue2 = ue^2 * hy;
                uw2 = uw^2 * hy;
                unv = un*vn * hx;
                usv = us*vs * hx;
                
                V = hx * hy;
                A = (ue2-uw2+unv-usv)/V;
                
                De = nu*(u(i+1,j)-u(i,j))/hx*hy;
                Dw = nu*(u(i,j)-u(i-1,j))/hx*hy;
                Dn = nu*(u(i,j+1)-u(i,j))/hy*hx;
                Ds = nu*(u(i,j)-u(i,j-1))/hy*hx;
                D = (De-Dw+Dn-Ds)/V;
                
                ut(i,j)=u(i,j)+dt*(-A+D);

            end
            
        end
    end
    
    % Temporary v-velocity
    for i=2:nx+1
        for j=2:ny 
            
            if (flagv(i,j)==0)

                vn = (v(i,j)+v(i,j+1))/2;
                vs = (v(i,j)+v(i,j-1))/2;
                ve = v(i,j)+(v(i+1,j)-v(i,j))*0.50;
                vw = v(i-1,j)+(v(i,j)-v(i-1,j))*0.50;
                ue = u(i,j)+(u(i,j+1)-u(i,j))*0.50;
                uw = u(i-1,j)+(u(i-1,j+1)-u(i-1,j))*0.50;
                
                vn2 = vn^2 * hx;
                vs2 = vs^2 * hx;
                veu = ve*ue * hy;
                vwu = vw*uw * hy;
                
                V = hx * hy;
                A = (vn2 - vs2 + veu - vwu)/V;
                
                De = nu*(v(i+1,j)-v(i,j))/hx*hy;
                Dw = nu*(v(i,j)-v(i-1,j))/hx*hy;
                Dn = nu*(v(i,j+1)-v(i,j))/hy*hx;
                Ds = nu*(v(i,j)-v(i,j-1))/hy*hx;
                D = (De-Dw+Dn-Ds)/V;
                
                vt(i,j)=v(i,j)+dt*(-A+D);

            end
            
        end
    end
    
end


% --------------------------------------------------------------------------------------
% [OBST] Identification of cells corresponding to the obstacle
% --------------------------------------------------------------------------------------
function [flagu, flagv, flagp] = ObstaclePreProcessing(nx,ny, xs,xe, ys,ye)

    flagu = zeros(nx+1,ny+2);   % u-cells corresponding to the obstacle
    for i=xs-1:xe
        for j=ys:ye
            flagu(i,j)=1;
        end
    end
    
    flagv = zeros(nx+2,ny+1);   % v-cells corresponding to the obstacle
    for i=xs:xe
        for j=ys-1:ye
            flagv(i,j)=1;
        end
    end
    
    flagp = zeros(nx+2,ny+2);   % p-cells corresponding to the obstacle
    for i=xs:xe
        for j=ys:ye
            flagp(i,j)=1;
        end
    end
    
end


% --------------------------------------------------------------------------------------
% [OBST] Reconstruct fields (graphical purposes only)
% --------------------------------------------------------------------------------------
function [uu, vv, pp] = ReconstructFields(u, v, p, nx, ny, flagu, flagv, flagp)

    % u-field reconstruction
    uu=zeros(nx+1,ny+1);
    for i=1:nx+1
        for j=1:ny+1
            if (flagu(i,j)==0)
                uu(i,j)=0.50*(u(i,j)+u(i,j+1));
            end
        end
    end

    % v-field reconstruction
    vv=zeros(nx+1,ny+1);
    for i=1:nx+1
        for j=1:ny+1
            if (flagv(i,j)==0)
                vv(i,j)=0.50*(v(i,j)+v(i+1,j));
            end
        end
    end

    % p-field reconstruction
    pp=zeros(nx+1,ny+1);
    for i=1:nx+1
        for j=1:ny+1
            if (flagp(i,j)==0)
                pp(i,j)=0.25*(p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1));
            end
        end
    end    

end