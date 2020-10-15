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
%   Copyright(C) 2019 Alberto Cuoci                                       %
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

% Input data
%-------------------------------------------------------------------------%
lengthx=1.0;            % domain length along x [m]
lengthy=1.0;            % domain length along y [m]
nx=41;                  % number of grid points along x
ny=41;                  % number of grid points along y
max_iterations=5000;    % max number of iterations
beta=1.;                % SOR coefficient (1 means Gauss-Siedler)
threshold = 0.0001;      % residual threshold

% Pre-processing of user-defined data
%-------------------------------------------------------------------------%
% Grid step calculation
hx=lengthx/(nx-1);      % grid step along x [m]
hy=lengthy/(ny-1);      % grid step along y [m]
x=0:hx:lengthx;         % x coordinates [m]
y=0:hy:lengthy;         % y coordinates [m]

% Memory allocation
f=zeros(nx,ny);         % numerical solution
fan=zeros(nx,ny);       % analytical solution
S=zeros(nx,ny);         % source term


%% Problem 1a: homogeneous source term
%-------------------------------------------------------------------------%

% Arbitrary constants
A = 2;
kappa = 1;

% Boundary conditions
f(1, :)  = exp(kappa*x(1))*sin(kappa*y) + A/4*(x(1)*x(1)+y.*y);
f(nx, :) = exp(kappa*x(nx))*sin(kappa*y) + A/4*(x(nx)*x(nx)+y.*y);
f(:, 1)  = exp(kappa*x)*sin(kappa*y(1)) + A/4*(x.*x+y(1)*y(1));
f(:, ny) = exp(kappa*x)*sin(kappa*y(ny)) + A/4*(x.*x+y(ny)*y(ny));

% Source term (homogeneous)
S(:,:) = A;

% Analytical solution
for i=1:nx
        for j=1:ny
            fan(i,j) = exp(kappa*x(i))*sin(kappa*y(j)) + A/4*(x(i)^2+y(j)^2);
        end
end

% Numerical solution
[f, tot_res] = SOR(f, S, nx, ny, hx, hy, beta, max_iterations, threshold);

% Graphical output
GraphicalOutput(f, fan, tot_res, x, y);


%% Problem 1b: non-homogeneous source term
%-------------------------------------------------------------------------%

% Boundary conditions
f(1, :)  = exp(x(1)*y);
f(nx, :) = exp(x(nx)*y);
f(:, 1)  = exp(x*y(1));
f(:, ny) = exp(x*y(ny));

% Source term (homogeneous)
for i=1:nx
        for j=1:ny
            S(i,j) = (x(i)^2+y(j)^2)*exp(x(i)*y(j));
        end
end

% Analytical solution
for i=1:nx
        for j=1:ny
            fan(i,j) = exp(x(i)*y(j));
        end
end

% Numerical solution
[f, tot_res] = SOR(f, S, nx, ny, hx, hy, beta, max_iterations, threshold);

% Graphical output
GraphicalOutput(f, fan, tot_res, x, y);


%% SOR Algorithm
%-------------------------------------------------------------------------%
function [f, tot_res] = SOR( f, S, nx, ny, hx, hy, ...
                             beta, max_iterations, threshold )

    for l=1:max_iterations

        for i=2:nx-1
            for j=2:ny-1
                f(i,j)= beta/(2*(hx^2+hy^2))*...
                        (   (f(i+1,j)+f(i-1,j))*hy^2 ...
                          + (f(i,j+1)+f(i,j-1))*hx^2 ...
                          - hx^2*hy^2*S(i,j) ...
                        ) + ...
                        (1.0-beta)*f(i,j);
            end
        end

        % Residual
        res=0;
        for i=2:nx-1
            for j=2:ny-1
                res=res+abs( (f(i+1,j)-2*f(i,j)+f(i-1,j))/hx^2 + ...
                             (f(i,j+1)-2*f(i,j)+f(i,j-1))/hy^2 - S(i,j) );
            end
        end
        tot_res(l) = res/((nx-2)*(ny-2));
        fprintf('Iteration: %d - Residual: %e\n', l, tot_res(l));

        if (tot_res(l) < threshold)
            break;
        end
    end
    
end


%% Graphical output
%-------------------------------------------------------------------------%
function GraphicalOutput(f, fan, tot_res, x, y)

    % Solution
    figure;
    pcolor(x, y, f');
    colorbar; shading interp; colormap(jet);
    hcb=colorbar; hcb.Title.String = "f value";
    title('numerical solution');
    xlabel('x-axis [m]'); ylabel('y-axis [m]');

    % Residual
    figure;
    semilogy(1:length(tot_res), tot_res);
    title('residual');
    xlabel('iteration'); ylabel('residual');

    % Error
    abserr = abs(f-fan);
    figure;
    pcolor(x, y, abserr');
    colorbar; shading interp; colormap(jet); 
    title('absolute error');
    xlabel('x-axis [m]'); ylabel('y-axis [m]');

end
