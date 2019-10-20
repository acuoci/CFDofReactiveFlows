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
%  Code: solver for tridiagonal linear systems                            %
%                                                                         %
%  [ d(1) u(1)                                ] [  x(1)  ]   [  b(1)  ]   %
%  [ l(2) d(2) u(2)                           ] [  x(2)  ]   [  b(2)  ]   %
%  [      l(3) d(3) u(3)                      ] [        ]   [        ]   %
%  [           ...  ...  ...                  ] [  ...   ] = [  ...   ]   %
%  [                     l(n-1) d(n-1) u(n-1) ] [ x(n-1) ]   [ b(n-1) ]   %
%  [                              l(n)   d(n) ] [  x(n)  ]   [  b(n)  ]   %
% ----------------------------------------------------------------------- %

function x = tridiagonal( d, l, u, b )

    n = length(b);
    v = zeros(n,1);   
    x = v;
    w = d(1);
    x(1) = b(1)/w;
    for i=2:n
        v(i-1) = u(i-1)/w;
        w = d(i) - l(i)*v(i-1);
        x(i) = ( b(i) - l(i)*x(i-1) )/w;
    end
    for j=n-1:-1:1
       x(j) = x(j) - v(j)*x(j+1);
    end
    
end
