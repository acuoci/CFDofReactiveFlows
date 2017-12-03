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
%  Code: 1D advection-diffusion-reaction                                  %
%  Method: segregation + implicit Euler with linearization of source term %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% Data
%-------------------------------------------------------------------------%
nv=3;               % number of variables (species)
np=51;              % number of grid points                
nsteps=200;         % number of time steps
length=1.0;         % domain length [m]
h=length/(np-1);    % grid step [m]
dt=0.0025;          % time step [s]
v=1;                % velocity [m/s]
D=0.05;             % diffusion coefficient [m2/s]
kappa = 1;          % kinetic constant [kmol,m,s]
alpha = 0.8;        % order od reaction of species A
beta = 1.1;         % order od reaction of species B
C0 = [1 0 0];       % initial concentration of A, B and C [kmol/m3]
Cin = [1 1 0];      % inlet concentration of A, B and C [kmol/m3]

% Memory allocation
%-------------------------------------------------------------------------%
C=zeros(np,nv);     % concentrations of species [kmol/m3]

d=zeros(np,nv);     % tridiagonal system: central diagonal
u=zeros(np,nv);     % tridiagonal system: upper diagonal
l=zeros(np,nv);     % tridiagonal system: lower diagonal
b=zeros(np,nv);     % tridiagonal system: rhs vector

% Initial conditions
%-------------------------------------------------------------------------%
for k=1:nv
    C(:,k) = C0(k);
end

% Video setup
%-------------------------------------------------------------------------%
video = VideoWriter('advection_diffusion_reaction_1d_linearized.mp4', 'MPEG-4');
open(video);

% Loop over time
%-------------------------------------------------------------------------%
t=0.; 
for m=1:nsteps
    
    % Constant coefficients 
    Ae =  v*dt/2/h - D*dt/h^2;
    Aw = -v*dt/2/h - D*dt/h^2;
    
    % West BC (inlet conditions)
    for k=1:nv
        d(1,k) = 1;
        b(1,k) = Cin(k);
    end
    
    % Internal points
    for i=2:np-1 
    
        % Reaction term
        r  = kappa*C(i,1)^alpha*C(i,2)^beta;
        R(1) = -r;
        R(2) = -r;
        R(3) =  r;
        J(1) = -kappa*alpha*C(i,1)^(alpha-1)*C(i,2)^beta;
        J(2) = -kappa*beta*C(i,1)^alpha*C(i,2)^(beta-1);
        J(3) = 0;
    
        % Upper/Lower diagonals
        u(i,:) = Ae;
        l(i,:) = Aw;
        
        % Central diagonal
        for k=1:nv
            d(i,k) = 1+2*D*dt/h^2-dt*J(k);
        end

        % RHS vector
        for k=1:nv
            b(i,k) = C(i,k)*(1-dt*J(k)) +R(k)*dt;
        end
        
    end 
    
    % East BC
    d(np,:) = 1;
    l(np,:) = -1;
    b(np,:) = 0;
        
    % Triadiagonal solutions
    for k=1:nv
        C(:,k) = tridiagonal(d(:,k),l(:,k),u(:,k),b(:,k));
    end
    
    % Plots
	hold off; 
    plot(0:h:length, C(:,1), 'linewidth',2); hold on;
    plot(0:h:length, C(:,2), 'linewidth',2); hold on;
    plot(0:h:length, C(:,3), 'linewidth',2); hold on;
    legend('A', 'B', 'C');                  
    xlabel('spatial coordinate [m]');
    ylabel('concentrations');    
    axis([0 length 0 1]);
    frame = getframe(gcf);
    writeVideo(video,frame);
    
    % Advance time
	t=t+dt; 
    
end

close(video);

% ------------------------------------------------------------------- %
% Write final concentration profiles
% ------------------------------------------------------------------- %
outFile = fopen('linearized.txt','w');
for i=1:np 
    fprintf(outFile,'%f %f %f %f \n', h*(i-1), C(i,1), C(i,2), C(i,3));
end
fclose(outFile);
