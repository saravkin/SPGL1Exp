function [y CM]=MSN_interp_1D(x, a, b, f, xprime, M, Opt)
%This function is the Minimum Sobolev Norm interpolant, that accepts
%inputs: 1. x - the input nodes (Nx1)
%        2. a,b - Boundary points; x is in (a,b)
%        3. f - values of the function at the nodes. (Nx1)
%        4. xprime - nodes to interpolate. (Lx1)
%        5. M   - Order of the polynomial used in interpolation.
%        6. Opt - Controls the Sobolev Norm being minimized.
%outputs: 1. y - interpolated values (at xprime). (Lx1)

% Authors:  Karthik Jayaraman Raghuram, jrk@ece.ucsb.edu
% v1.0 Created 21-Okt-09 
%
% hss_lu: 2D msn interpolation 
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

%Default the Sobolev Norm minimized to 2
if nargin == 6
    Opt = 2;
end

%Default the Order of Polynomial to be used to 2*number of nodes.
if nargin == 5
    Opt = 2;
    M = 2*length(x);
end

%Normalize the input nodes to (-1,1)
x = 2*(x-a)/(b-a) - 1;

%Construct the Chebyshev Vandermonde
VNM = cos(acos(x)*[0:M-1]); %NxM matrix
D = diag([1:M].^(-Opt));

%Calculate the MSN Coeffs
%Note: This solution can be done much more effectively as the Matrix is
%Diagonal Plus Semi-Separable.
[q r] = qr(D*VNM'); %VNM = r'q';
% cond(D*VNM')
CM = D*q*(r'\f);%Mx1

%Normalize the output nodes to (-1,1)
xprime = 2*(xprime - a)/(b - a) - 1;

%Construct the output Chebyshev Vandermonde
VPM = cos( acos(xprime) * [0:M-1] );%PxM matrix

%Reconstruct the output samples
y = VPM*CM;%Px1

