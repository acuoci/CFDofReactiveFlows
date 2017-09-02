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
%  Code: 2D advection-diffusion by the FTCS scheme                        %
%        The code is adapted and extended from Tryggvason, Computational  %
%        Fluid Dynamics http://www.nd.edu/~gtryggva/CFD-Course/           %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;

% two-dimensional unsteady diffusion by the FTCS scheme
%------------------------------------------------------------
nx=33;                  % number of grid points along x
ny=33;                  % number of grid points along y
nstep=150;              % number of time steps
lengthx=2.0;            % domain length along x [m]
lengthy=2.0;            % domain length along y [m]
hx=lengthx/(nx-1);      % grid step along x [m]
hy=lengthy/(ny-1);      % grid step along y [m] 
D=0.025;                % diffusion coefficient [m2/s]

dt=1.0*0.125*hx*hy/D;   % time step [s]

u=-1;                   % velocity along x [m/s]
v= 0;                   % velocity along y [m/s]

f=zeros(nx,ny);     % current numerical solution
fo=zeros(nx,ny);    % previous numerical solution

% grid axes
xaxis = 0:hx:lengthx;
yaxis = 0:hy:lengthy;

% Boundary conditions
f(nx, ny*1/3:ny*2/3) = 1;

% prepare video
videompg4 = VideoWriter('advection_diffusion_2d.mp4', 'MPEG-4');
open(videompg4);

t = 0;
for l=1:nstep
    
    % Graphics only
	hold off; 
    mesh(xaxis, yaxis, f'); 
    axis([0 lengthx 0 lengthy 0 1.25]);
    xlabel('x'); ylabel('y'); zlabel('f');
    message = sprintf('time=%d\n', t);
    time = annotation('textbox',[0.15 0.8 0.15 0.15],'String',message,'EdgeColor','none');
    frame = getframe(gcf);
    writeVideo(videompg4,frame);
    delete(time);
    
    % Forward Euler
    fo=f;
    for i=2:nx-1
        for j=2:ny-1
            f(i,j) = fo(i,j)...
                    -(0.5*dt*u/hx)*(fo(i+1,j)-fo(i-1,j))...
                    -(0.5*dt*v/hy)*(fo(i,j+1)-fo(i,j-1))...
                    +(D*dt/hx^2)*(fo(i+1,j)-2*fo(i,j)+fo(i-1,j))...
                    +(D*dt/hy^2)*(fo(i,j+1)-2*fo(i,j)+fo(i,j-1));
        end
    end
    
    % Boundary conditions: south and north sides (Neumann)
    for i=1:nx
        f(i,1)=f(i,2);
        f(i,ny)=f(i,ny-1);
    end
    
    % Boundary conditions: west side (Neumann)
    for j=1:ny
        f(1,j)=f(2,j);
    end
    
    % Advance time
    t=t+dt;
    
end

close(videompg4);
