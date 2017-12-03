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

linearized = importdata('linearized.txt');
fullycoupled = importdata('fullycoupled.txt');
splitting = importdata('splitting.txt');

x = linearized(:,1);

figure;
plot(x, fullycoupled(:,2), x, linearized(:,2), x, splitting(:,2), ...
     'linewidth',2);
legend('fully-coupled', 'linearized', 'operator-splitting');
xlabel('spatial coordinate [m]');
ylabel('concentrations');
title('concentration of A');

figure;
plot(x, fullycoupled(:,3), x, linearized(:,3), x, splitting(:,3), ...
     'linewidth',2);
legend('fully-coupled', 'linearized', 'operator-splitting');
xlabel('spatial coordinate [m]');
ylabel('concentrations');
title('concentration of B');

figure;
plot(x, fullycoupled(:,4), x, linearized(:,4), x, splitting(:,4), ...
     'linewidth',2);
legend('fully-coupled', 'linearized', 'operator-splitting');
xlabel('spatial coordinate [m]');
ylabel('concentrations');
title('concentration of C');