function [y, n, time] = msn_interpol_1D(x,a,b,f,xev,s,sorted,maxleafsize,run,I, J, S,solverIterations)
%y = msn_interpol_1D(x,a,b,f,xev)
% This is the minimum Sobolev norm (MSN) interpolation in 1D. With the
% interpolations points x and there function values f = f(x)
% msn_interpol_1D calculates a approximative function value at the
% evaluation points xev. The algorithm of the interpolation uses linear
% time, to do the interpolation in bought variable length(x) and
% length(xev). But the program does not jet interpolate in linear time.
% Bought x and xev can be randomly distributed and must be within the
% boundary [a,b]. 
% Theoretical framework:
% * Fast Interpolation Schemes with Sobolev Type Norm Minimization, SIAM
%   Annual Meeting July 2009. Karthik Jayaraman Raghuram, Shivkumar
%   Chandrasekaran, Ming Gu, Hrushikesh Mhaskar 
% * Fast Interpolation Schemes with Sobolev Type Norm Minimization, SIAM
%   Annual Meeting July 2009. Karthik Jayaraman Raghuram, Shivkumar
%   Chandrasekaran, Ming Gu, Hrushikesh Mhaskar 
% * The Minimum Sobolev Norm Interpolation Scheme and its applications in
%   Image Processing, (To appear in) The Proceedings of IS&T/SPIE 2010,
%   Shivkumar Chandrasekaran, Karthik Jayaraman Raghuram, Hrushikesh
%   Mhaskar  
%
% inputs:
%         1. x - vector of interpolation points. For a O(n) algorithm, x has to
%         be sorted s.t. x(1) <= x(2) <= ... <= x(end).
%         2. a,b - Boundary points; x is in (a,b).
%         3. f - vector of values of the function at the interpolation points.
%         4. xev - vector of evaluation points. For a O(n) algorithm,
%         xev has to be sorted s.t. xev(1) <= xev(2) <= ... <=
%         xev(end).
% optional inputs:
%         5. s - sobolev parameter, determine the sobolev norm to be
%         minimized
%            (default = 2)
%         6. sorted - true if x and xev are already sorted s.t.
%         x(1) <= x(2) <= ... <= x (end) and equal for xev.
%         if sorted = false -> sort them in the program by calling sort(x)
%         -> takes longer than O(n)
%            (default = false)
%         7. maxleafsize - the max size of the leaf nodes in the HSS
%         representation.
%            (default = 150)
%         8. run - if true -> run the function
%                  if false -> precompute the memory usage for I,J and S.
%                  n = ( length of I, J and S).
%            (default = true)
%         9. I,J,S - preallocated Vector of the size (nx1).
%            (default: allocated during the process)
%         10. solverIterations - tries to improve the accuracy of the
%         solution by a number of iterative refinements. Did not
%         improve the accuracy of the solution!
%            (defult = 0)
% outputs:
%         1. y - evaluated values (at xev).
%         3. n - minimal length of I,J,S.
%         4. time - struct of different execution times measured.
%
%
% Example, a simple but slow way
%   f = @(x) 1./(1+100*x.^2);
%   a = -1;
%   b = 1;
%   x = sort(rand(1,100)*(b-a)+a);
%   xev = linspace(a,b);
%   y = msn_interpol_1D(x,a,b,f(x),xev);
%   plot(linspace(a,b), f(linspace(a,b)),'-',x,f(x),'o', xev,y,'--')
%
%
% Example, a bit more complicate but faster way. Its important to say to
% the interpolation, that your x and xev are sorted.
% If you do more than one Interpolation, its useful to preallocate I,J and
% S.
%   f = @(x) 1./(1+100*x.^2);
%   a = -1;
%   b = 1;
%   x = sort(rand(1,100)*(b-a)+a);
%   xev = linspace(a,b);
%   maxleafsize = 30;
%   [y,n] = msn_interpol_1D(x,a,b,f(x),xev,2,true,maxleafsize,false);
%   I = zeros(n,1);
%   J = zeros(n,1);
%   S = zeros(n,1);
%   y = msn_interpol_1D(x,a,b,f(x),xev,2,true,maxleafsize,true,I,J,S);
%   plot(linspace(a,b), f(linspace(a,b)),'-',x,f(x),'o', xev,y,'--')
%

% Authors:  Stefan Pauli, stefan.pauli@alumni.ethz.ch
%           Karthik Jayaraman Raghuram, jrk@ece.ucsb.edu
% v1.0 Created 21-Okt-09
%
% msn_interpol_1D: This is the minimum Sobolev norm (MSN) interpolation in 1D.
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


% measure the time of the whole interpolation
time = struct();
time.tTot = cputime;




% Default case for s: s = 2
if nargin <6
    s = 2;
end

% Default case for sorted: sorted = false
if nargin <7
    sorted = false;
end


% Default case for maxleafsize: maxleafsize = 150
if nargin <8
    maxleafsize = 150;
end

% Default case for run: run = true
if nargin < 9
    run = true;
end

% transform the inputs in to a column vector.
f=f(:);
x=x(:);
xev=xev(:);

% sort x and xev if needed
if (sorted ==  false)
    [x ind] = sort(x);
    f = f(ind); % interchange the f values equally to the x values
    [xev sortIntInd] =  sort(xev);
end

% Default case for solverIterations: solverIterations = 0
if nargin < 13
    solverIterations = 0;

    % Default case for I, J and S:
    %   run = true: allocate the space needed by those variables.
    %   run = false: I, J and S nead no space.
    if nargin < 12
        if run==true
            % Compute the needed size of I, J, and S and allocate them
            [y, n] = msn_interpol_1D(x,a,b,f,xev,s,true,maxleafsize,false);
            I = zeros(n,1);
            J = zeros(n,1);
            S = zeros(n,1);
        else
            I=[];
            J=[];
            S=[];
        end
    end
end

if s~=2 && s~=1
    error('just s = 2 and s = 1 is supported by this programm');
end

% normalize the x to [-1,1]
xN=2*(x-a)/(b-a)-1;

% normalize the xreq to [-1,1]
xevN=2*(xev-a)/(b-a)-1;

% make the transformati0n
thetai = acos(xN);
theta =  acos(xevN);

% Do some test's
if (length(x) ~= length(f))
    error('length of x must be the same as length of f')
end

%% the solver part
tic %tempTime = cputime;

% kernel = D + triu(U*V',1) + tril(P*Q',-1), where D = diagonal matrix with
% diag(D) = diag(U*V') = diag(P*Q'): compute U,V,P and Q
if s==1
    V = [ones(size(thetai)), thetai.^2];
    Q = [1 + pi^2/6 - pi/2*thetai + pi/4*thetai.^2, -1/4*ones(size(thetai))];
    P = V;
    U = Q;
elseif s==2
    V = [ones(size(thetai)), thetai.^2, thetai.^4];
    Q = [1 + pi^4/90 - pi^2/12*thetai.^2 + pi/12*thetai.^3 - 1/48*thetai.^4, ...
        -pi^2/12 + pi/4*thetai - 1/8*thetai.^2,    -1/48*ones(size(thetai))];
    P = V;
    U = Q;
else
    error('the s you choose is not supported.');
end

% put the kernel in a HSS matrix
hss = hss_buildSepK(P,Q,U,V,thetai,thetai,acos(-1),acos(1),maxleafsize);
tempTime = toc; %cputime-tempTime;

% solve the system kernel*ftemp = f
if (solverIterations==0)
    [ftemp,time.solve ,n] = hss_solve(hss,f,run,I, J, S);
    time.solve.buildhss = tempTime;
end

if (solverIterations~=0)
    [L, U, b, indx, time.solve ,n] = hss_lu(hss,f,run,I, J, S);
    time.solve.buildhss = tempTime;

    % iterations for improofing the accuracy of the solution
    % was a tried but did not help
    if (run==true)
        ftempLong = U\(L\b);
        ftemp = ftempLong(indx,:);
        for i = 1:solverIterations
            r = hss_mul(hss,ftemp) - f;   % residual = A*xtilde - b;
            rlong=sparse(length(b),1);
            rlong(indx) = r;
            ftempLong = ftempLong - U\(L\rlong);      % x = xtilde - A/r
            ftemp = ftempLong(indx,:);
        end
    end
end


%% the multiplication part
if run==true % not needed to calculate the space needed for I,J,S
    tic % time.multiply = cputime;

    % if theta >= thetai: kernel = U*V'
    % if thetai >= theta: kernel = P*Q'
    % compute U,V,P and Q:
    if s==1
        th=thetai;
        V = [ones(size(th)), th.^2];
        Q = [1 + pi^2/6 - pi/2*th + pi/4*th.^2, -1/4*ones(size(th))];

        th=theta;
        P = [ones(size(th)), th.^2];
        U = [1 + pi^2/6 - pi/2*th + pi/4*th.^2, -1/4*ones(size(th))];
    elseif s==2

        th=thetai;
        V = [ones(size(th)), th.^2, th.^4];
        Q = [1 + pi^4/90 - pi^2/12*th.^2 + pi/12*th.^3 - 1/48*th.^4, ...
            -pi^2/12 + pi/4*th - 1/8*th.^2,    -1/48*ones(size(th))];

        th=theta;
        P = [ones(size(th)), th.^2, th.^4];
        U = [1 + pi^4/90 - pi^2/12*th.^2 + pi/12*th.^3 - 1/48*th.^4, ...
            -pi^2/12 + pi/4*th - 1/8*th.^2,    -1/48*ones(size(th))];
    else
        error('the s you choose is not supported.');
    end
    % put the kernle in a HSS matrix
    hss = hss_buildSepK(P,Q,U,V,theta,thetai,acos(-1),acos(1),maxleafsize);

    % y = kernel*ftemp
    y = hss_mul(hss,ftemp);
    time.multiply = toc; %cputime-time.multiply;
else
    % create a void y in case run = false (Just compute the calculate the
    % space needed for I,J,S)
    y=[];
end

if sorted == false
    % y mast be reordert according to the unsorted input values
    y(sortIntInd) = y;
end

% measure the time of the whole interpolation
time.tTot = cputime-time.tTot;
end