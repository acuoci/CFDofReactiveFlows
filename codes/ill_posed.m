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
%  Code: example of ill-posed problem                                     %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% Data
%-------------------------------------------------------------------------%
N=1000;             % number of grid points
h = 2*pi/(N-1);     % grid spacing
x=0:h:2*pi;         % grid
dt = 0.0005;        % time step
nsteps = 400;       % number of steps
alpha = 1;          % diffusion coefficient

klow = 1/pi;        % low frequency
kmedium = 5/pi;     % medium frequency
khigh = 40/pi;      % high frequency

wmedium = 0;        % weight medium frequency
whigh = 0.02;        % weight high frequency

% Calculations
%-------------------------------------------------------------------------%
alow = sin(2*pi*klow*x);                % low frequency wave
amedium = sin(2*pi*kmedium*x);          % medium frequency wave
ahigh = sin(2*pi*khigh*x);              % high frequency wave
a = alow+wmedium*amedium+whigh*ahigh;   % global wave
plot(x,a, x,alow, x,amedium, x,ahigh);  % plot contributions

% Video setup
%-------------------------------------------------------------------------%
v = VideoWriter('ill_posed.mp4', 'MPEG-4');
open(v);

% Time loop
%-------------------------------------------------------------------------%
f = a;          % current solution
fo = a;         % solution at n-1
foo = a;        % solution at n-2

t = 0;
for k=1:nsteps

    % Advancing solution
    for i=2:N-1
        f(i) = -dt^2/h^2*alpha*(fo(i+1)-2*fo(i)+fo(i-1)) + 2*fo(i) - foo(i);
    end
    f(1) = -dt^2/h^2*alpha*(fo(2)-2*fo(1)+fo(N-1)) + 2*fo(1) - foo(1);
    f(N) = f(1);
    
    % Plot solution
    hold off;
    plot(x,f);
    axis([0 6.28 -5 5]);
    frame = getframe(gcf);
    writeVideo(v,frame);

    % Store previous solution
    foo = fo;
    fo = f;
    
    t = t+dt;
    
end
    
close(v);
