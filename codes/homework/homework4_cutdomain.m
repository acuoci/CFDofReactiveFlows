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
Uwnorth=0.1;            % north wall velocity [m/s]
Uwsouth=0.1;            % south wall velocity [m/s]
nu=1e-3;                % kinematic viscosity [m2/s]
ttot=100;               % total time of simulation [s]

% Parameters for SOR
max_iterations=10000;   % maximum number of iterations
beta=1.5;               % SOR coefficient
max_error=0.0001;       % error for convergence

% Active domain
active = ones(nx,ny);
ibc = ceil(nx/2);
jbc = ceil(ny/2);
for i=ibc:nx
    for j=1:jbc
        active(i,j) = 0;
    end
end

%% Pre-processing operations
%--------------------------------------------------------------------------
V=max(abs(Uwnorth),abs(Uwsouth)); % reference velocity [m/s]
L=max(Lx,Ly);           % reference length [m]
Re=V*L/nu;              % Reynolds' number [-]
tref=L/V;               % reference time [s]
lengthx=Lx/L;           % dimensionless length along x [-] 
lengthy=Ly/L;           % dimensionless length along y [-] 
tautot=ttot/tref;       % total time of simulation (dimensionless) [-]
uwnorth = Uwnorth/V;    % dimensionless north wall velocity [-]
uwsouth = Uwsouth/V;    % dimensionless south wall velocity [-]

% Grid step
hx=lengthx/(nx-1);      % grid step along x [-]
hy=lengthy/(ny-1);      % grid step along y [-]

% Time step
sigma = 0.5;                            % safety factor for time step (stability)
dtau_diff=min(hx,hy)^2*Re/4;            % time step (diffusion stability) [-]
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
x=0:hx:lengthx;         % grid coordinates (x axis) [-]
y=0:hy:lengthy;         % grid coordinates (y axis) [-]
[X,Y] = meshgrid(x,y);  % mesh

% Rectangle position
position = [lengthx/2 0 lengthx/2, lengthy/2];
color = [0.85 0.85 0.85];

%% Numerical solution
%--------------------------------------------------------------------------
tau = 0;
for istep=1:nsteps     
    
    % ------------------------------------------------------------------- %
    % Poisson equation (SOR)
    % ------------------------------------------------------------------- %
    [psi,iter] = Poisson2D( psi,nx,ny,hx,hy,-omega, ...
                            beta,max_iterations,max_error, active );
    
    % ------------------------------------------------------------------- %
    % Reconstruction of dimensionless velocity field
    % ------------------------------------------------------------------- % 
    [u,v] = ReconstructDimensionlessVelocity(u,v,psi,nx,ny,hx,hy,...
                                             uwnorth,uwsouth,active);
    
    % ------------------------------------------------------------------- %
    % Find vorticity on boundaries
    % ------------------------------------------------------------------- %
    omega(1,:)  = (psi(1,:)-psi(2,:))*2/hx^2 ;                      % west
    omega(:,ny) = (psi(:,ny)-psi(:,ny-1))*2/hy^2 -2/hy*uwnorth;     % north
    
    % South wall
    omega(1:ibc,1) = (psi(1:ibc,1)-psi(1:ibc,2))*2/hy^2 + 2/hy*uwsouth;                                
    omega(ibc:nx,jbc) = (psi(ibc:nx,jbc)-psi(ibc:nx,jbc+1))*2/hy^2;                    
    
    % East wall
    omega(nx,jbc:ny) = (psi(nx,jbc:ny)-psi(nx-1,jbc:ny))*2/hx^2 ;
    omega(ibc,1:jbc) = (psi(ibc,1:jbc)-psi(ibc-1,1:jbc))*2/hx^2 ;
    
    % ------------------------------------------------------------------- %
    % Advection-diffusion equation (new vorticity in interior points)
    % ------------------------------------------------------------------- %
    [omega] = AdvectionDiffusion2D(omega, u,v, Re, nx,ny, hx,hy, dtau, active);
    
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
        contour(x,y,psi', 40, 'b');
        axis('square');
        rectangle( 'Position', position, 'FaceColor', color);
        pause(0.01);
    end
    
end


%% Final post-processing operations
% ------------------------------------------------------------------- %

subplot(231);
surface(x,y,u');
axis('square'); title('u'); xlabel('x'); ylabel('y');
rectangle( 'Position', position, 'FaceColor', color);

subplot(234);
surface(x,y,v');
axis('square'); title('v'); xlabel('x'); ylabel('y');
rectangle( 'Position', position, 'FaceColor', color);

subplot(232);
surface(x,y,omega');
axis('square'); title('omega'); xlabel('x'); ylabel('y');
rectangle( 'Position', position, 'FaceColor', color);

subplot(235);
surface(x,y,psi');
axis('square'); title('psi'); xlabel('x'); ylabel('y');
rectangle( 'Position', position, 'FaceColor', color);

subplot(233);
contour(x,y,psi', 30, 'b');
axis('square');
title('stream lines'); xlabel('x'); ylabel('y');
rectangle( 'Position', position, 'FaceColor', color);

subplot(236);
quiver(x,y,u',v');
axis([0 lengthx 0 lengthy], 'square');
title('stream lines'); xlabel('x'); ylabel('y');
rectangle( 'Position', position, 'FaceColor', color);


%% ------------------------------------------------------------------------
% Poisson equation solver
%  ------------------------------------------------------------------------
function [f,iter] = Poisson2D(f,nx,ny,hx,hy,S,beta,...
                              max_iterations,max_error, active)

    B = (hx^2*hy^2)/2/(hx^2+hy^2);
    Ae = B/hx^2; Aw = Ae;
    An = B/hy^2; As = An;
    
    for iter=1:max_iterations
        
        for i=2:nx-1
            for j=2:ny-1
                
                if (active(i,j) == 1)
                    f(i,j) = beta*( Ae*f(i+1,j) + Aw*f(i-1,j) + ...
                                    An*f(i,j+1) + As*f(i,j-1) + ...
                                    -B*S(i,j) ) + ...
                             (1-beta)*f(i,j);
                end
                
            end
        end
        
        res = 0;
        for i=2:nx-1
            for j=2:ny-1
                if (active(i,j) == 1)
                    res = res + abs( (f(i+1,j)-2*f(i,j)+f(i-1,j))/hx^2 + ...
                                     (f(i,j+1)-2*f(i,j)+f(i,j-1))/hy^2 + ...
                                     -S(i,j) ) ;
                end
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
function [u,v] = ReconstructDimensionlessVelocity(u,v,psi,nx,ny,hx,hy,...
                 uwnorth, uwsouth, active)

    u(:,ny) = uwnorth;
    u(:,1)  = uwsouth; 
    for i=2:nx-1
        for j=2:ny-1
            if (active(i,j) == 1)
                u(i,j) =  ( psi(i,j+1)-psi(i,j-1) )/(2*hy);
                v(i,j) = -( psi(i+1,j)-psi(i-1,j) )/(2*hx);
            end
        end
    end

end


%% ------------------------------------------------------------------------
%  Advection-diffusion equation: forward Euler + centered discretization
%  ------------------------------------------------------------------------
function [f] = AdvectionDiffusion2D(f, u,v, Re, nx,ny, hx,hy, dtau, active)

    fo = f;
    for i=2:nx-1
        for j=2:ny-1
            
            if (active(i,j) == 1)
                
                A = u(i,j)*(fo(i+1,j)-fo(i-1,j))/(2*hx) + ...
                    v(i,j)*(fo(i,j+1)-fo(i,j-1))/(2*hy) ;

                D = 1/Re * ( (fo(i+1,j)-2*fo(i,j)+fo(i-1,j))/hx^2 + ...
                             (fo(i,j+1)-2*fo(i,j)+fo(i,j-1))/hy^2 ) ;

                f(i,j) = fo(i,j) + (-A + D)*dtau;
                
            end
        end
    end

end
