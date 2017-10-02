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
%  Code: 1D advection-diffusion by the FTCS scheme                        %
%        The code is adapted and extended from Tryggvason, Computational  %
%        Fluid Dynamics http://www.nd.edu/~gtryggva/CFD-Course/           %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% Data
%-------------------------------------------------------------------------%
np=21;              % number of grid points                
nstep=100;          % number of time steps
length=2.0;         % domain length [m]
h=length/(np-1);    % grid step [m]
dt=0.05;            % time step [s]
u=1;                % velocity [m/s]
D=0.05;             % diffusion coefficient [m2/s]
A=0.5;              % amplitude of initial solution
k=1;                % wave number [1/m]

% Memory allocation
%-------------------------------------------------------------------------%
f=zeros(np,1);      % temporary numerical solution
y=zeros(np,1);      % current numerical solution
a=zeros(np,1);      % exact solution

% Video setup
%-------------------------------------------------------------------------%
v = VideoWriter('advection_diffusion_1d.mp4', 'MPEG-4');
open(v);

% Initial solution
%-------------------------------------------------------------------------%
for i=1:np
	y(i)=A*sin(2*pi*k*h*(i-1)); 
end

% Loop over time
%-------------------------------------------------------------------------%
t=0.; 
for m=1:nstep 
	
    % Exact solution
    for i=1:np 
		a(i) = A*exp(-4*pi*pi*k*k*D*t)*sin(2*pi*k*(h*(i-1)-u*t)); 
    end
    
    % Integrals (post processing only)
    a2_int = 0.;
    y2_int = 0.;
    for i=1:np-1
         a2_int = a2_int + h/2*(a(i)^2+a(i+1)^2);
         y2_int = y2_int + h/2*(y(i)^2+y(i+1)^2);
    end
         
    % Plots
    message = sprintf('time=%d\na^2(int)=%d\ny^2(int)=%d', t, a2_int, y2_int);
	hold off; plot([0:h:length],y,'linewidth',2); axis([0 length -1, 1]); % plot num. 
	hold on; plot([0:h:length],a,'r','linewidth',2);                      % plot exact
    hold on; legend('numerical', 'exact');                  % legend
    xlabel('spatial coordinate [m]');
    ylabel('solution');    
    time = annotation('textbox',[0.15 0.8 0.1 0.1],'String',message,'EdgeColor','none');
    frame = getframe(gcf);
    writeVideo(v,frame);
    delete(time);
    
    % Temporary solution 
	f=y; 
    
    % Numerical solution (internal points)
    for i=2:np-1 
		y(i) = f(i)-(u*dt/2/h)*(f(i+1)-f(i-1))+...    % advection
			   D*(dt/h^2)*(f(i+1)-2*f(i)+f(i-1));     % diffusion
    end 
    
    % Numerical solution (periodic boundary conditions)
	y(np) = f(np)-(u*dt/2/h)*(f(2)-f(np-1))+...
            D*(dt/h^2)*(f(2)-2*f(np)+f(np-1)); 
	y(1)  = y(np);
    
    % Calculate error
    E = 0;
    for i=1:np 
        E = E + (y(i)-a(i))^2;
    end
    E = h*sqrt(E);
    
    % Advance time
	t=t+dt; 
    fprintf('time=%d E=%e\n', t, E);
    
end

close(v);
