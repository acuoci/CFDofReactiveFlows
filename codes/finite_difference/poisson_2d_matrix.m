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
%                                                                         %
%  Code: 2D Poisson equation solved by explicit construction of the       %
%        sparse matrix associated to the linear system                    %
%                                                                         %
% ----------------------------------------------------------------------- %
close all;
clear variables;


%% User-defined data
%-------------------------------------------------------------------------%
lengthx=2.0;            % domain length along x [m]
lengthy=2.0;            % domain length along y [m]
nx=102;                  % number of grid points along x
ny=102;                  % number of grid points along y


%% Pre-processing of user-defined data
%-------------------------------------------------------------------------%
% Grid step calculation
hx=lengthx/(nx-1);      % grid step along x [m]
hy=lengthy/(ny-1);      % grid step along y [m]

% Memory allocation
ne = nx*ny;             % number of equations
S=zeros(nx,ny);         % source term


%% Construct sparsity matrix
%-------------------------------------------------------------------------%
As = -1/hy^2;
Ae = -1/hx^2;
Aw = -1/hx^2;
An = -1/hy^2;
Ap = -(As+Aw+Ae+An);

b = zeros(ne,1);
M=sparse(ne,ne);

% Internal points
for i=2:nx-1
    for j=2:ny-1
        k = (j-1)*nx+i;
        M(k,k) = Ap; 
        M(k,k-1)=Aw; M(k,k+1)=Ae; 
        M(k,k-nx)=As; M(k,k+nx)=An;
        b(k) = -S(i,j);
    end
end

% Boundaries (south & north)
for i=1:nx
    
    M(i,i) = 1.;
    b(i) = 0.;
    
    M(nx*(ny-1)+i,nx*(ny-1)+i) = 1.;
    b(nx*(ny-1)+i) = 0.;
    
end

% Boundaries (west & east)
for j=1:ny
    
    M((j-1)*nx+1,(j-1)*nx+1) = 1.;
    b((j-1)*nx+1) = 0;
    if (j>=ny*1/3 && j<=ny*2/3) 
        b((j-1)*nx+1) = 1.;
    end
    
    M(j*nx,j*nx) = 1.;
    b(j*nx) = 0.;
    
end


%% Solve linear system of equations
%-------------------------------------------------------------------------%
% Opt 1: Automatically choose best algorithm for the given pattern
f = M\b;   
% Opt 2: Generalized Minimum Residual Method
%        Be careful, no preconditioner used! Very dangerous!
%[f, flag, relres, iter, res] = gmres(M,b,[],1e-5,1000);


% Graphical output
%-------------------------------------------------------------------------%
figure;
xaxis = 0:hx:lengthx;
yaxis = 0:hy:lengthy;
f = reshape(f,nx,ny);
pcolor(xaxis, yaxis, f');
colorbar; shading interp; colormap(jet);
hcb=colorbar; hcb.Title.String = "f value";
xlabel('x-axis [m]'); ylabel('y-axis [m]');
