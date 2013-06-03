function hss = hss_buildSepK(P,Q,R,S,x,y,a,b,maxleafsize,minThetaDist)
%hss = hss_buildSepK(P,Q,R,S,x,y,a,b)
% Build a HSS structure containing a separable kernel.
% the x and the y values have to be sorted! x(1) <= x(2) <= ... <= x(end)
% and P,Q,R and S have to be in the order of x and y.
% The separable kernel has the following properties:
% f(x,y) = R*S'  if x >= y
% f(x,y) = P*Q'  if y >= x
%
%
% inputs:
%       1. P,Q,R,S - Describes the kernel matrix where,
%          f(x,y) = R*S'  if x >= y,
%          f(x,y) = P*Q'  if y >= x.
%       3. x,y - x has to be sorted s.t. x(1) <= x(2) <= ... <= x(end), and
%       equal for y.
%       2. a,b - Boundary points; y is in (a,b).
% optional inputs:
%       4. maxleafsize - the max size of the leaf nodes in the HSS
%       representation.
%          (default = 150);
%       5. minThetaDist - stop cutting if minThetaDist > abs(a-b). Needed
%       to prevent infinite loops for to many indistinguishable x, y
%       in machine accuracy
%          (default = 100*eps)
%
% output:
%       1. hss -  the HSS representation of the matrix
%

% Authors:  Stefan Pauli, stefan.pauli@alumni.ethz.ch
%           Karthik Jayaraman Raghuram, jrk@ece.ucsb.edu
% v1.0 Created 21-Okt-09
%
% hss_buildSepK: Build a HSS structure containing a separable kernel
% Copyright (C) 2009  Stefan Pauli (stefan.pauli@alumni.ethz.ch)
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

if nargin<10
    minThetaDist = 100*eps;
    if nargin<9
        maxleafsize = 150;
    end
end

% since the rank must be constant:
if  ~(size(P)==size(R))
    error('R and U must have the same size');
elseif ~(size(Q)==size(S))
    error('Q and V must have the same size');
else

    % We are not constructing the optimal HSS matrix, but still our code is
    % liner. The not optimality comes form the fact that we are forced to use
    % [P,U] as our U in the HSS representation, and the same for [Q,V] and
    % therefore the HSS matrix is uncompressed. The optimal HSS matrix
    % could be fund by compressing the HSS matrix. This can be done in
    % linear time.
    U = [P,R];
    V = [Q,S];
    m = size(U,1);
    n = size(V,1);
    rank = size(V,2);
    hss = hss_buildSepK_iter(m, n, U, V,a,b, rank, y, x,maxleafsize,minThetaDist);
end
end




function [hss] = hss_buildSepK_iter(m, n, U, V,a,b, rank, thetai, theta,maxleafsize, minThetaDist, k)
% Construct a HSS matrix:
% We are not constructing the optimal HSS matrix, but still our code is
% liner. The not optimality comes form the fact that we are forced to use
% [P,U] as our U in the HSS representation, and the same for [Q,V] and
% therefore the HSS matrix is uncompressed. The optimal HSS matrix
% could be fund by compressing the HSS matrix. This can be done in
% linear time.
%
% inputs:
%       1. m,n - Size of the matrix.
%       2. U,V - U and V of the HSS representation.
%       3. a,b - Boundary points; x is in (a,b).
%       4. rank - rank of the off-diagonal elements.
%       5. thetai,theta -  transformed x,xint (thetai = acos(x) , theta =
%       acos(xint)).
%       6. maxleafsize - the max size of the leaf nodes in the HSS
%       representation.
%       9. minThetaDist - stop cutting if minThetaDist > abs(a-b). Needed
%       to prevent infinite loops for to many indistinguishable x, y
%       in machine accuracy
% internal inputs:
%       10. k -  do not set k! (represents the depth of the actual node in
%       the tree)
%          (default = 1)
%
% output:
%       1. hss -  the HSS representation of the matrix
%
if nargin<12
    k=1;
end

if k==1;    % the top node of the tree representation
    hss = struct('topnode',1);
else        % not the top node of the tree representation
    hss = struct('topnode',0);
end

% the depth of the children nodes in the tree is one more
k=k+1;



if (((n>maxleafsize)||(m>maxleafsize)) && minThetaDist<abs(a-b)) % parent node
    hss.leafnode=0;
    % in order to divide the matrix in to 4 parts, compute where to divide
    % it.
    cut = (a-b)/2+b;

    % calculate the size of the children according to the divided matrix
    hss.nl = sum(thetai>=cut);
    hss.nr = n - hss.nl;

    hss.ml = sum(theta>=cut);
    hss.mr = m - hss.ml;

    %recursion part
    [hss.hssl] = hss_buildSepK_iter  (hss.ml, hss.nl, U(1:hss.ml,:), V(1:hss.nl,:), a, cut, rank, thetai(1:hss.nl), theta(1:hss.ml), maxleafsize, minThetaDist, k);
    [hss.hssr] = hss_buildSepK_iter  (hss.mr, hss.nr, U(hss.ml+1:end,:), V(hss.nl+1:end,:), cut, b, rank, thetai(hss.nl+1:end), theta(hss.ml+1:end), maxleafsize, minThetaDist, k);

    %
    Bu = spdiags([zeros(rank/2,1) ; ones(rank/2,1)], 0, rank, rank);
    Bl = spdiags([ones(rank/2,1) ; zeros(rank/2,1)], 0, rank, rank);

    % coppy the rigth U, V part of the children
    W = speye(rank);

    % Debug:
    %     Bu = full(Bu);
    %     Bl = full(Bl);
    %     W = full(W);

    % attach it to the structure
    hss.Rl=W;
    hss.Rr=W;
    hss.Bu=Bu;
    hss.Bl=Bl;
    hss.Wl=W;
    hss.Wr=W;

else    % leav node
    hss.leafnode= 1;
    hss.U=U;
    hss.V=V;

    % if theta >= thetai: A = U*V'
    % if thetai >= theta: A = P*Q'
    % and U = [P,U];  V = [Q,V];
    A = U(:,1:end/2)*V(:,1:end/2)'; % = P*Q';
    thetai=thetai(:);
    theta=theta(:);
    i = find((theta*ones(1,length(thetai)) - ones(length(theta),1)*thetai') > 0);
    B = U(:,end/2+1:end)*V(:,end/2+1:end)'; % = U*V';
    A(i) = B(i);
    hss.D = A;
end

end