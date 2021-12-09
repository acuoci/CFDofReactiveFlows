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
%  Code: 2D driven-cavity problem on a staggered grid                     %
%        Transport equation for species A (with reaction A+B->C)          %
%        Second order reaction: r=kappa*CA*CB                             %
%        Implicit Backward Euler method + linearization of source term +  %
%        coupling of species                                              %
%        The solution of 9-diagonal system of equations is carried out    %
%        using a solver for banded systems                                %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% ----------------------------------------------------------------------- %
% User data
% ----------------------------------------------------------------------- %

% Only even numbers of cells are acceptable
nx=24;      % number of (physical) cells along x
ny=nx;      % number of (physical) cells along y
L=1;        % length [m]
nu=0.01;    % kinematic viscosity [m2/s] (if L=1 and un=1, then Re=1/nu)
tau=10;     % total time of simulation [s]

% Boundary conditions
un=1;       % north wall velocity [m/s]
us=0;       % south wall velocity [m/s]
ve=0;       % east wall velocity [m/s]
vw=0;       % west wall velocity [m/s]

% Parameters for SOR
max_iterations=10000;   % maximum number of iterations
beta=1.3;               % SOR coefficient
max_error=1e-5;         % error for convergence

% Data for transport equations of species
Gamma=0.01;             % diffusion coefficient [m2/s] 
CAin = 1;               % inlet concentration of A [kmol/m3]
CBin = 1;               % inlet concentration of B [kmol/m3]
kappa = 100;            % kinetic constant [m3/kmol/s]

% ----------------------------------------------------------------------- %
% Data processing
% ----------------------------------------------------------------------- %

% Grid step
h=L/nx;                                 % grid step (uniform grid) [m]

% Time step
sigma = 0.5;                            % safety factor for time step (stability)
dt_diff_ns=h^2/4/nu;                    % time step (diffusion stability) [s]
dt_conv_ns=4*nu/un^2;                   % time step (convection stability) [s]
dt_ns=min(dt_diff_ns, dt_conv_ns);      % time step (stability) [s]
dt_diff_sp=h^2/4/Gamma;                 % time step (species diffusion stability) [s]
dt_conv_sp=4*Gamma/un^2;                % time step (species convection stability) [s]
dt_sp=min(dt_conv_sp, dt_diff_sp);      % time step (stability) [s]
dt_reac=1/kappa/(CAin+CBin);            % time step (stability) [s]

dt=sigma*dt_ns;                         % [IMPL] time step (stability based on momentum equations only) [s]
nsteps=tau/dt;                          % number of steps
Re = un*L/nu;                           % Reynolds' number

% Summary on the screen
fprintf('Time step: %f\n', dt);
fprintf(' - Diffusion (NS):       %f\n', dt_diff_ns);
fprintf(' - Convection (NS):      %f\n', dt_conv_ns);
fprintf(' - Diffusion (Species):  %f\n', dt_diff_sp);
fprintf(' - Convection (Species): %f\n', dt_conv_sp);
fprintf(' - Reaction:             %f\n', dt_reac);
fprintf('Reynolds number: %f\n', Re);

% Grid construction
x=0:h:1;                % grid coordinates (x axis)
y=0:h:1;                % grid coordinates (y axis)
[X,Y] = meshgrid(x,y);  % MATLAB grid

% ----------------------------------------------------------------------- %
% Memory allocation
% ----------------------------------------------------------------------- %

% Main fields (velocities and pressure)
u=zeros(nx+1,ny+2);
v=zeros(nx+2,ny+1);
p=zeros(nx+2,ny+2);
CA=zeros(nx+2,ny+2);
CB=zeros(nx+2,ny+2);
CC=zeros(nx+2,ny+2);

% Temporary velocity fields
ut=zeros(nx+1,ny+2);
vt=zeros(nx+2,ny+1);

% Temporary pressure field (convergence of SOR)
po=zeros(nx+2,ny+2);

% Fields used only for graphical post-processing purposes
uu=zeros(nx+1,ny+1);
vv=zeros(nx+1,ny+1);

% Coefficient for pressure equation
gamma=zeros(nx+2,ny+2)+1/4;
gamma(2,3:ny)=1/3;gamma(nx+1,3:ny)=1/3;gamma(3:nx,2)=1/3;gamma(3:nx,ny+1)=1/3;
gamma(2,2)=1/2;gamma(2,ny+1)=1/2;gamma(nx+1,2)=1/2;gamma(nx+1,ny+1)=1/2;

% [IMPL] Total number of equations for a single scalar
ne=(nx+2)*(ny+2);

% ----------------------------------------------------------------------- %
% Solution over time
% ----------------------------------------------------------------------- %
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;

% Monitor CPU time (in s)
start_time = cputime;

t=0.0;
for is=1:nsteps
    
    % ------------------------------------------------------------------- %
    % 1. Projection algorithm                                             %
    % ------------------------------------------------------------------- %

    % Boundary conditions
    u(1:nx+1,1)=2*us-u(1:nx+1,2);               % south wall
    u(1:nx+1,ny+2)=2*un-u(1:nx+1,ny+1);         % north wall
    v(1,1:ny+1)=2*vw-v(2,1:ny+1);               % west wall
    v(nx+2,1:ny+1)=2*ve-v(nx+1,1:ny+1);         % east wall
    
    % Temporary velocity
    [ut, vt] = TemporaryVelocity(ut,vt, u,v, dt,h, nu, nx,ny);
    
    % Pressure equation (Poisson)
    [p, it] = PoissonSolver(p, ut,vt, dt,h, nx,ny, gamma, beta, max_iterations, max_error);
    
    % Correct the velocity
    u(2:nx,2:ny+1)=ut(2:nx,2:ny+1)-(dt/h)*(p(3:nx+1,2:ny+1)-p(2:nx,2:ny+1));
    v(2:nx+1,2:ny)=vt(2:nx+1,2:ny)-(dt/h)*(p(2:nx+1,3:ny+1)-p(2:nx+1,2:ny));
      

    % ------------------------------------------------------------------- %
    % 2. Transport of species                                             %
    % ------------------------------------------------------------------- %

    % Impermeable walls
    CA = ImpermeableWallsBoundaryConditions(CA, nx,ny);
    CB = ImpermeableWallsBoundaryConditions(CB, nx,ny);   
    CC = ImpermeableWallsBoundaryConditions(CC, nx,ny); 
    
    % Porous inlet (west)
    CA(1,ny/2:3/4*ny)=2*CAin-CA(2,ny/2:3/4*ny);
    CB(1,ny/2:3/4*ny)=2*0-CB(2,ny/2:3/4*ny);
    CC(1,ny/2:3/4*ny)=2*0-CC(2,ny/2:3/4*ny);

    % Porous inlet (east)
    CA(nx+2,ny/4:ny/2)=2*0-CA(nx+1,ny/4:ny/2);
    CB(nx+2,ny/4:ny/2)=2*CBin-CB(nx+1,ny/4:ny/2);
    CC(nx+2,ny/4:ny/2)=2*0-CC(nx+1,ny/4:ny/2);
    

    % [IMPL] Concentration of species

    tStart = cputime;   % track cpu time

    % 1. Storing current values
    CAo = CA;
    CBo = CB;
    CCo = CC;
    
    % 2. Setting the matrix M (common to every species)
    M = eye(ne,ne);
    for i=2:nx+1
            for j=2:ny+1
                
                k = (nx+2)*(j-1)+i;
                
                ue = u(i,j);    uw = u(i-1,j);
                vn = v(i,j);    vs = v(i,j-1);
                
                M(k,k)   = 1 + 4*Gamma*dt/h/h;
                M(k,k+1) = dt*( ue/2/h-Gamma/h^2);
                M(k,k-1) = dt*(-uw/2/h-Gamma/h^2);
                M(k,k+nx+2) = dt*( vn/2/h-Gamma/h^2);
                M(k,k-nx-2) = dt*(-vs/2/h-Gamma/h^2);                
                
            end
    end

    % 3. Setting the matrices for single species
    MA = M; MB = M; MC = M;
    for i=2:nx+1
            for j=2:ny+1
                
                k = (nx+2)*(j-1)+i;
                
                MA(k,k) = MA(k,k) + dt*kappa*CBo(i,j);
                MB(k,k) = MB(k,k) + dt*kappa*CAo(i,j);
            end
    end

    % 4. Assembling fully-coupled matrix
    MM = zeros(3*ne, 3*ne);
    MM(1:ne,1:ne) = MA;
    MM(ne+1:2*ne,ne+1:2*ne) = MB;
    MM(2*ne+1:3*ne,2*ne+1:3*ne) = MC;

    % 5. Adding coupling terms
    for i=2:nx+1
            for j=2:ny+1

                k = (nx+2)*(j-1)+i;
                MM(k,k+ne) = MM(k,k+ne) + dt*kappa*CAo(i,j);            % B on A
                MM(k+ne,k) = MM(k+ne,k) + dt*kappa*CBo(i,j);            % A on B
                MM(k+2*ne,k) = MM(k+2*ne,k) - dt*kappa*CBo(i,j);        % A on C
                MM(k+2*ne,k+ne) = MM(k+2*ne,k+ne) - dt*kappa*CAo(i,j);  % B on C
            end
    end

    % 6. Setting the RHS vectors
    bA=CAo(:);
    bB=CBo(:);
    bC=CCo(:);
    
    % 7. Setting the RHS vectors for single species
    for i=2:nx+1
            for j=2:ny+1
                k = (nx+2)*(j-1)+i;
                bA(k) = bA(k) + dt*kappa*CAo(i,j)*CBo(i,j);
                bB(k) = bB(k) + dt*kappa*CAo(i,j)*CBo(i,j);
                bC(k) = bC(k) - dt*kappa*CAo(i,j)*CBo(i,j);
            end
    end
    
    % 8. Assembling fully-coupled vector
    bb = zeros(3*ne,1);
    bb(1:ne) = bA;
    bb(1+ne:2*ne) = bB;
    bb(1+2*ne:3*ne) = bC;
    
    % 9. Solve the linear system
    C = decomposition(MM, 'banded')\bb;
    CA = C(1:ne);
    CB = C(1+ne:2*ne);
    CC = C(1+2*ne:3*ne);
    
    % 10. Reshape the solution
    CA = reshape(CA, [nx+2 ny+2]);
    CB = reshape(CB, [nx+2 ny+2]);
    CC = reshape(CC, [nx+2 ny+2]);    

    tEnd = cputime;   % track cpu time

    % Print info on the screen
    if (mod(is,25)==1)
        fprintf('Step: %d - Time: %f - Poisson iterations: %d - CPU time (species): %f s\n', is, t, it, tEnd-tStart);
    end

    % ------------------------------------------------------------------- %
    % 3. On-the-fly post-processing                                       %
    % ------------------------------------------------------------------- %
    if (mod(is,25) == 1)    % update figures every 25 steps
        
        % Reconstruction
        uu(1:nx+1,1:ny+1)=0.50*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
        vv(1:nx+1,1:ny+1)=0.50*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
        ca = ReconstructScalarField(CA, nx,ny);
        cb = ReconstructScalarField(CB, nx,ny);
        cc = ReconstructScalarField(CC, nx,ny);

        % Surface maps
        subplot(221);
        surface(X,Y,ca'+cb', 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
        axis('square'); title('CA+CB [kmol/m3]'); xlabel('x'); ylabel('y');
        colorbar; shading interp;
        subplot(222);
        surface(X,Y,cc', 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
        axis('square'); title('CC [kmol/m3]'); xlabel('x'); ylabel('y');
        colorbar; caxis([1e-11 max(max(cc))+1e-10]); shading interp;
        subplot(223);
        surface(X,Y,kappa*(ca.*cb)', 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
        axis('square'); title('reaction rate [kmol/m3/s]'); xlabel('x'); ylabel('y');
        colorbar; caxis([1e-10 max(max(kappa*(ca.*cb)))+1e-9]); shading interp;
        subplot(224);
        surface(X,Y,kappa*(ca'+cb'), 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
        axis('square'); title('keff [1/s]'); xlabel('x'); ylabel('y');
        colorbar; shading interp;
        
        pause(0.01);
    
    end
    
    % Advance time
    t=t+dt;
 
    pause(0.01);
    
end
hold off;

% Monitor CPU time
end_time = cputime;
fprintf('Elapsed time: %f s\n', end_time-start_time);

% ----------------------------------------------------------------------- %
% Final post-processing                                                   %
% ----------------------------------------------------------------------- %

% Field reconstruction
uu(1:nx+1,1:ny+1)=0.50*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
vv(1:nx+1,1:ny+1)=0.50*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
pp = ReconstructScalarField(p, nx,ny);
ca = ReconstructScalarField(CA, nx,ny);
cb = ReconstructScalarField(CB, nx,ny);
cc = ReconstructScalarField(CC, nx,ny);
      
% Surface map: CA
subplot(231);
surface(X,Y,ca');
axis('square'); title('CA [kmol/m3]'); xlabel('x'); ylabel('y');

% Surface map: CB
subplot(232);
surface(X,Y,cb');
axis('square'); title('CB [kmol/m3]'); xlabel('x'); ylabel('y');

% Surface map: CC
subplot(233);
surface(X,Y,cc');
axis('square'); title('CC [kmol/m3]'); xlabel('x'); ylabel('y');

% Plot: velocity components along the horizontal middle axis
subplot(234);
plot(x, 0.5*(ca(:,round(ny/2))+ca(:,round(ny/2)+1)) );
hold on;
plot(x, 0.5*(cb(:,round(ny/2))+cb(:,round(ny/2)+1)) );
hold on;
plot(x, 0.5*(cc(:,round(ny/2))+cc(:,round(ny/2)+1)) );
axis('square');
xlabel('x [m]');ylabel('concentrations [kmol/m3]'); title('concentrations along x-axis');
legend('CA', 'CB', 'CC');
hold off;

% Plot: velocity components along the horizontal middle axis
subplot(235);
plot(y, 0.5*(ca(round(nx/2),:)+ca(round(nx/2)+1,:)) );
hold on;
plot(y, 0.5*(cb(round(nx/2),:)+cb(round(nx/2)+1,:)) );
hold on;
plot(y, 0.5*(cc(round(nx/2),:)+cc(round(nx/2)+1,:)) );
axis('square');
xlabel('y [m]');ylabel('concentrations [kmol/m3]'); title('concentrations along y-axis');
legend('CA', 'CB', 'CC');
hold off;

% Surface map: r
subplot(236);
surface(X,Y,kappa*ca'*cb');
axis('square'); title('reaction rate [kmol/m3/s]'); xlabel('x'); ylabel('y');


% ----------------------------------------------------------------------- %
%                            FUNCTIONS                                    %
% ----------------------------------------------------------------------- %

% Temporary velocity calculation (Predictor step in Projection algorithm)
function [ut, vt] = TemporaryVelocity(ut,vt, u,v, dt,h, nu, nx,ny)

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


% Pressure equation (Poisson)
function [p, it] = PoissonSolver(p, ut,vt, dt,h, nx,ny, gamma, beta, max_iterations, max_error)

    for it=1:max_iterations
        
        po=p;
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
                epsilon=epsilon+abs(po(i,j)-p(i,j)); 
            end
        end
        epsilon = epsilon / (nx*ny);
        
        % Check the error
        if (epsilon <= max_error) % stop if converged
            break;
        end 
        
    end

end


% Reconstruct scalar field on a grid for graphical purposes
function [phir] = ReconstructScalarField(phi, nx,ny)
    
    phir(1:nx+1,1:ny+1)=0.25*(  phi(1:nx+1,1:ny+1)+phi(1:nx+1,2:ny+2)+...
                                phi(2:nx+2,1:ny+1)+phi(2:nx+2,2:ny+2));
end


% Set the boundary conditions on impermeable walls for scalars
function [phi] = ImpermeableWallsBoundaryConditions(phi, nx,ny)

    phi(2:nx+1,1)=phi(2:nx+1,2);          
    phi(2:nx+1,ny+2)=phi(2:nx+1,ny+1);    
    phi(1,2:ny+1)=phi(2,2:ny+1);          
    phi(nx+2,2:ny+1)=phi(nx+1,2:ny+1); 

end
