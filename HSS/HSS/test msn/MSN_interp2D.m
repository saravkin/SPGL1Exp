function [I, Coeff, hss2]=MSN_interp2D(C, f, C_o, Mx, My, s)
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


%Note that we do not need to make any assumptions about the ordering of the
%coordinates.
%Make the data Double just in case.
f = double(f);

if nargin == 5
    s = 2;
end
%Step1 - Calculate the Input Vandermonde matrix at the input coordinates.
%At each of the input coordinates (a row), we need to calculate the product
%of cosines.

N = size(C, 1);%We have these many coordinates.
% display('Number of samples in the region:');
% display(N);

VNM2 = zeros(N, Mx*My);
x = C(:,1);
y = C(:,2);

nd=0;
if nd==1
    [xnd ynd fnd] = NestedDissection(x, y,f, -1, 1, -1, 1);
else
    xnd = x;
    ynd = y;
    fnd = f;
end

p=1;Ds = sparse(zeros(Mx*My));
for l=0:Mx-1,
    for m=0:My-1,
        VNM2(:,p) = (cos(l*acos(xnd)).*cos(m*acos(ynd)));
        Ds(p,p) = (1+l^2+m^2)^(-s);
        p = p+1;
    end
end

% [xd yd fndy] = NestedDissection(y, x,[1:length(x)]', -1, 1, -1, 1);
% 
% VX = cos(acos(x)*[0:Mx-1]);Dx = diag([1:Mx].^(-2*s));
% VY = cos(acos(y)*[0:My-1]);Dy = diag([1:My].^(-2*s));
% KX = VX*Dx*Dx'*VX';
% KY = VY*Dy*Dy'*VY';
% 
% [hssXND] = hss_buildA(KX,1e-10,2);
% [hssYND] = hss_buildA(KY,1e-10,2);
% 
% VX = cos(acos(x)*[0:Mx-1]);Dx = diag([1:Mx].^(-2*s));
% VY = cos(acos(y)*[0:My-1]);Dy = diag([1:My].^(-2*s));
% KX = VX*Dx*Dx'*VX';
% KY = VY*Dy*Dy'*VY';
% 
% [hssX] = hss_buildA(KX,1e-10,2);
% [hssY] = hss_buildA(KY,1e-10,2);
% 
% S1 = KX .* KY;
% 
% [hss1] = hss_buildA(S1,1e-10,2);

% Ds_vec = diag(Ds);
% [X, I] = sort(Ds_vec,1,'descend');
% VNM2 = VNM2(:, I);
% Ds = Ds(:,I);

VNM2 = VNM2*Ds;
% S2 = VNM2*VNM2';
% [hss2] = hss_buildA(S2,1e-7*norm(S2),5);
% figure, imagesc(S);
hss2 = [];
%Step2 - Calculate the MSN Coefficients for the inut 
[q r] = qr(VNM2',0);
Coeff = Ds*q*(r'\fnd);
%Coeff = pinv(VNM2)*f;

% [l u] = lu(VNM2');
% [q r] = qr(l);
% Coeff = Ds*q*(r'\((u'\f)));

%clear VNM2;
%Step3 - Calculate the output Vandermonde matrix at the output set of
%coordinates.
P = size(C_o, 1);
VPM2 = zeros(P, Mx*My);

x = C_o(:,1);
y = C_o(:,2);
%[x y ind] = NestedDissection(x, y, -1, 1, -1, 1);

p=1;
for l=0:Mx-1,
    for m=0:My-1,
        VPM2(:,p) = (cos(l*acos(x)).*cos(m*acos(y)));
        p = p+1;
    end
end

%VPM2 = VPM2(:,I);

%Step 4 - Calculate the output sample values using the coefficients and the
%output vandermonde.
I = VPM2*Coeff;