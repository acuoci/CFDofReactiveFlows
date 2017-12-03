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
%  Method: operator splitting (explicit transport)                        %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% Global variables
%-------------------------------------------------------------------------%
global nv;
global kappa;
global alpha;
global beta;

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
C=zeros(np,nv);     % current concentrations [kmol/m3]
Ct=zeros(nv,1);     % temporary concentrations [kmol/m3]

% Initial conditions
%-------------------------------------------------------------------------%
for k=1:nv
    C(:,k) = C0(k);
end

% Video setup
%-------------------------------------------------------------------------%
video = VideoWriter('advection_diffusion_reaction_1d_splitting.mp4', 'MPEG-4');
open(video);

% Loop over time
%-------------------------------------------------------------------------%
Co = C; % previous (old) concentrations [kmol/m3]

t=0.; 
for m=1:nsteps
    
    % West BC (inlet)
    for k=1:nv
        C(1,k) = Cin(k);
    end
    
    % Internal points
    for i=2:np-1 
        
        % Transport (explicit Euler method)
        for k=1:nv
            Ct(k) = Co(i,k) + dt*( -v*(Co(i+1,k)-Co(i-1,k))/2/h + ...
                                    D*(Co(i+1,k)-2*Co(i,k)+Co(i-1,k))/h^2);
        end
        
        % Reaction
        [dummy, Cf] = ode15s(@reaction,[t t+dt],Ct);
        for k=1:nv
            C(i,k) = Cf(end,k);
        end
        
    end 
    
    % East BC (outlet)
    for k=1:nv
        C(np,k) = C(np-1,k);
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
    Co = C;     % store old value
	t=t+dt; 
    
end

close(video);

% ------------------------------------------------------------------- %
% Write final concentration profiles
% ------------------------------------------------------------------- %
outFile = fopen('splitting.txt','w');
for i=1:np 
    fprintf(outFile,'%f %f %f %f \n', h*(i-1), C(i,1), C(i,2), C(i,3));
end
fclose(outFile);

% ------------------------------------------------------------------- %
% Stiff ODE system
% ------------------------------------------------------------------- %
function dC = reaction(t,C)

    global nv;
    global kappa;
    global alpha;
    global beta;
    
    dC = zeros(nv,1);
        
    r = kappa*(C(1)^alpha)*(C(2)^beta);
    R(1) = -r;
    R(2) = -r;
    R(3) =  r;
        
    for k=1:nv
        dC(k)  = R(k);
    end

end
