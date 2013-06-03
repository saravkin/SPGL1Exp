function [I, wts]=MSN_interp2D(C, f, C_o, ~, ~, s)
%This function performs a 2-D MSN interpolation.
%The following are the inputs:
%   C - Co-Ordinates of the input samples
%   f - Interpolating Samples
%   C_o - The output coordinates
%   s - the sobolev norm being minimized.
%   Mx - the order of polynomial along the x directio.
%   My - the order of the polynomial along the y direction.
%The following are the outputs:
%   I - The interpolated set of samples.
% Authors:  Karthik Jayaraman Raghuram, jrk@ece.ucsb.edu
% v1.0 Created 13-Oct-10
%
% Copyright (C) 2009  Karthik Jayaraman Raghuram (jrk@ece.ucsb.edu)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
maxgridsize = 1e6;
tol = 1/maxgridsize;
x = C(:,1);
y = C(:,2);
xprime = C_o(:,1);
yprime = C_o(:,2);

if(size(x,1) == 1), x = x'; end
if(size(y,1) == 1), y = y'; end
if(size(f,1) == 1), f = f'; end

x = 2*(x-ax)./(bx-ax) - 1;
xprime = 2*(xprime-ax)./(bx-ax) - 1;

y = 2*(y-ay)./(by-ay) - 1;
yprime = 2*(yprime-ay)./(by-ay) - 1;

 %Compute the mesh norm.
xc = acos(x);
yc = acos(y);

[Mx My] = mesh_norm2(xc, yc, tol);
M = Mx*My;
V = vandermonde2(x, y, Mx, My);
D_inv = zeros(M,1);
for l=0:Mx-1
    D_inv((l*My + 1):(l+1)*My,1) = (1+l^2+(0:My-1).^2).^(-s/2);
end

% D_inv = spdiags(D_inv, 0, M, M);
[~,ps] = sort(D_inv,'descend');
D_inv = spdiags(D_inv, 0, M, M);
[q r v ~] = coda_WLS(D_inv(ps,:), V', 2);
Vo = vandermonde2(xprime, yprime, Mx, My);
rhs = D_inv(ps,:)*Vo';
wts = v*(r\q'*rhs);
I = wts'*f;