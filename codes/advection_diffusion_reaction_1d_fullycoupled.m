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
%  Method: fully-coupled + method of lines                                %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% Global variables
%-------------------------------------------------------------------------%
global nv;
global np;
global Cin;
global kappa;
global alpha;
global beta;
global v;
global D;
global h;

% Data
%-------------------------------------------------------------------------%
nv=3;               % number of variables
np=51;             % number of grid points                
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
ne=nv*np;           % total number of equations
C=zeros(np,nv);     % concentrations of species [kmol/m3]
Y=zeros(ne,1);      % concentration vector [kmol/m3]

% Initial conditions
%-------------------------------------------------------------------------%
for i=1:np
    for k=1:nv
        j = (i-1)*nv+k;
        Y(j)=C0(k);
    end
end

% Video setup
%-------------------------------------------------------------------------%
video = VideoWriter('advection_diffusion_reaction_1d_fullycoupled.mp4', 'MPEG-4');
open(video);


% Fully-coupled method
%-------------------------------------------------------------------------%

% Mass matrix (1=differential, 0=algebraic)
Mdense = eye(ne, ne);
for k=1:nv
    Mdense(k,k)=0;                  % inlet
    Mdense(ne-nv+k,ne-nv+k)=0;      % outlet
end
M = sparse(Mdense);

% DAE system solution
options = odeset('Mass',M, 'RelTol',1e-6, 'AbsTol',1e-8);
[t,Ymatrix] = ode15s(@advection_diffusion_reaction,0:dt:nsteps*dt,Y,options);

for m=1:nsteps
    
    % Reconstruct current solution
    Y = Ymatrix(m,:);
    for k=1:nv
        for i=1:np
            j = (i-1)*nv + k;
            C(i,k) = Y(j);
        end
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

end
    
close(video);

% ------------------------------------------------------------------- %
% Write final concentration profiles
% ------------------------------------------------------------------- %
outFile = fopen('fullycoupled.txt','w');
for i=1:np 
    fprintf(outFile,'%f %f %f %f \n', h*(i-1), C(i,1), C(i,2), C(i,3));
end
fclose(outFile);


% ------------------------------------------------------------------- %
% DAE system
% ------------------------------------------------------------------- %
function dC = advection_diffusion_reaction(t,C)

    global nv;
    global np;
    global Cin;
    global kappa;
    global alpha;
    global beta;
    global v;
    global D;
    global h;
    
    dC = zeros(np*nv,1);

    % West BC (inlet)
    for k=1:nv
        j = k;
        dC(j) = C(j) - Cin(j);
    end

    % Internal points
    for i=2:np-1
        
        j = (i-1)*nv;
        
        r = kappa*(C(j+1)^alpha)*(C(j+2)^beta);
        R(1) = -r;
        R(2) = -r;
        R(3) =  r;
        
        for k=1:nv
            
            j = (i-1)*nv + k;

            dC(j)  = -v*(C(j+nv)-C(j-nv))/2/h + ...
                      D*(C(j+nv)-2*C(j)+C(j-nv))/h^2 + R(k);
        end
    end
    
    % East BC (outlet)
    for k=1:nv
        j = (np-1)*nv + k;
        dC(j) = C(j) - C(j-nv);
    end

end
