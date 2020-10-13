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
%-------------------------------------------------------------------------%
%   Exact solutions from K. Masatsuka, "I do like CFD, Vol. 1", 2018      %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% List of problems
problem_1 = true;
problem_2 = true;
problem_3 = true;
problem_4 = true;

% User-defined data
%-------------------------------------------------------------------------%
D=0.01;             % diffusion coefficient [m2/s]
u=0.1;              % velocity [m/s]
L=1.0;              % domain length [m]
np=21;              % number of grid points                
nstep=200;          % number of time steps
dt=0.1;             % time step [s]


%% Pre-processing of user-defined data
%-------------------------------------------------------------------------%
% Grid step calculation
h=L/(np-1);         % grid step [m]
x = 0:h:L;          % grid coordinates [m]

% Memory allocation
fo=zeros(np,1);     % temporary numerical solution
f=zeros(np,1);      % current numerical solution
a=zeros(np,1);      % exact solution

% Check the stability conditions on time step
Di = D*dt/h^2;                      % Diffusion number
Co = u*dt/h;                        % Courant number
dt_max = min(h^2/D*1/2, h/u*1);     % Maximum allowed time step
fprintf('Di=%f, Co=%f, dt=%f, dt(max)=%f\n', Di, Co, dt, dt_max);

%% Problem 1
%-------------------------------------------------------------------------%
if problem_1 == true
    
    % IC: f(x,0) = sin(n*pi*x)
    % BC: f(0,t) = f(1,t) = 0
    % Analytical solution: f = exp(-D*n^2*pi^2*t)*sin(n*pi*x)

    % Additional data
    n = 2;  % arbitrary, integer number 

    % Initial solution
    f=sin(n*pi*x);

    % Advancing in time
    t = 0.;
    for m=1:nstep

        % Update the analytical solution
        a = exp(-D*n^2*pi^2*t)*sin(n*pi*x); 

        % Graphical output
        GraphicalComparison(x,f,a,-1,1); 

        % Dirichlet boundary conditions
        f(1)  = 0;
        f(np) = 0; 
        
        % Forward Euler + Centered Discretization (FECD)
        fo = f;
        f = FECD(f,fo, 0, D, h, dt);

        % Update the error between numerical and analytical solution
        E = ErrorEstimation(h,f,a);

        % New time step
        t=t+dt;

        % Print the current time (every 25 steps)
        if (mod(m,25)==1), fprintf('time=%d E=%e\n', t, E); end
    end
	
	fprintf('Press enter to continue...\n');
	pause;
    
end

%% Problem 2
%-------------------------------------------------------------------------%
if problem_2 == true
    
    % IC: f(x,0) = sqrt(A/B)*exp(-(x-1)^2/(4*D*B))
    % BC: f(0,t) = sqrt(A/(t+B))*exp(-1/(4*D*(t+B)))
    % BC: f(1,t) = sqrt(A/(t+B))
    % Analytical solution: f = sqrt(A/(t+B))*exp(-(x-1)^2/(4*D*(t+B)))

    % Additional data
    A=2; B=3;   % arbitrary constants 

    % Initial solution
    f=sqrt(A/B)*exp(-(x-1).^2/(4*D*B));

    % Advancing in time
    t = 0.;
    for m=1:nstep

        % Update the analytical solution
        a = sqrt(A/(t+B))*exp(-(x-1).^2/(4*D*(t+B))); 

        % Graphical output
        GraphicalComparison(x,f,a,0,1); 

        % Dirichlet boundary conditions
        f(1)  = sqrt(A/(t+B))*exp(-1/(4*D*(t+B)));
        f(np) = sqrt(A/(t+B));
        
        % Forward Euler + Centered Discretization (FECD)
        fo = f;
        f = FECD(f,fo, 0, D, h, dt); 

        % Update the error between numerical and analytical solution
        E = ErrorEstimation(h,f,a);

        % New time step
        t=t+dt;

        % Print the current time (every 25 steps)
        if (mod(m,25)==1), fprintf('time=%d E=%e\n', t, E); end
    end
    
    fprintf('Press enter to continue...\n');
	pause;
    
end

%% Problem 3
%-------------------------------------------------------------------------%
if problem_3 == true

    % IC: f(x,0) = x
    % BC: f(0,t) = 0, f(1,t) = 1
    % Analytical solution (steady-state): f = (1-exp(xRe))/(1-exp(Re))
    % where Re=u/D
    
    % Initial solution
    f=x;
    
    % Advancing in time
    t = 0.;
    for m=1:nstep

        % Graphical output
        GraphicalNumericalSolution(x,f,0,1); 
        
        % Dirichlet boundary conditions
        f(1)  = 0;
        f(np) = 1;
        
        % Forward Euler + Centered Discretization (FECD)
        fo = f;
        f = FECD(f,fo, u, D, h, dt); 

        % New time step
        t=t+dt;

    end
    
    % Steady-state analytical solution
    a = (1-exp(x*u/D))/(1-exp(u/D));
    
    % Graphical output
    GraphicalComparison(x,f,a,0,1); 
    
    fprintf('Press enter to continue...\n');
	pause;
    
end


%% Problem 4
%-------------------------------------------------------------------------%
if problem_4 == true

    % Source term: C*pi*(u*cos(pi*x)+pi*D*sin(pi*x))
    % IC: f(x,0) = x
    % BC: f(0,t) = 0, f(1,t) = 1
    % Analytical solution (steady-state): 
    % f = (1-exp(xRe))/(1-exp(Re)) + C*sin(pi*x)
    
    % Additional data
    C = 1.5;
    
    % Initial solution
    f=x;
    
    % Advancing in time
    t = 0.;
    for m=1:nstep

        % Graphical output
        GraphicalNumericalSolution(x,f,0,2); 
        
        % Dirichlet boundary conditions
        f(1)  = 0;
        f(np) = 1;
        
        % Forward Euler + Centered Discretization + Source term
        fo = f;
        for i=2:np-1 
            f(i) = fo(i) - 1/2*(u*dt/h)*(fo(i+1)-fo(i-1)) ...
                         + (D*dt/h^2)*(fo(i+1)-2*fo(i)+fo(i-1)) ...
                         + dt*C*pi*(u*cos(pi*x(i))+pi*D*sin(pi*x(i))) ;
        end 

        % New time step
        t=t+dt;

    end
    
    % Steady-state analytical solution
    a = (1-exp(x*u/D))/(1-exp(u/D)) + C*sin(pi*x);
    
    % Graphical output
    GraphicalComparison(x,f,a,0,2); 
    
    fprintf('Press enter to continue...\n');
	pause;
end

%-------------------------------------------------------------------------%
%% Function: Forward-Euler + Centered Difference
%-------------------------------------------------------------------------%
function f = FECD(f, fo, u, D, h, dt)  

    Co = u*dt/h;
    Di = D*dt/h^2;
    np = length(fo);
    for i=2:np-1 
		f(i) = fo(i) - 1/2*Co*(fo(i+1)-fo(i-1)) ...
                     + Di*(fo(i+1)-2*fo(i)+fo(i-1));
    end 
    
end

%-------------------------------------------------------------------------%
%% Function: error estimation
%-------------------------------------------------------------------------%
function E = ErrorEstimation(h,f,a)

    E = 0;
    for i=1:length(f) 
        E = E + (f(i)-a(i))^2;
    end
    E = h*sqrt(E);
    
end

%-------------------------------------------------------------------------%
%% Function: graphical comparison
%-------------------------------------------------------------------------%
function GraphicalComparison(x,f,a,ymin,ymax)

	hold off; plot(x,f,'linewidth',2); axis([0 max(x) ymin, ymax]); % plot num. 
	hold on; plot(x,a,'r--','linewidth',2);                         % plot exact
    hold on; legend('numerical', 'exact');                          % legend
    xlabel('spatial coordinate [m]');
    ylabel('solution');
    pause(0.01);
end

%-------------------------------------------------------------------------%
%% Function: graphical comparison
%-------------------------------------------------------------------------%
function GraphicalNumericalSolution(x,f,ymin,ymax)

	hold off; plot(x,f,'linewidth',2); axis([0 max(x) ymin, ymax]); % plot num. 
    xlabel('spatial coordinate [m]');
    ylabel('solution');
    pause(0.01);
    
end
