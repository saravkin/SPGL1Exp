function msn2DTest()
% function to make a lot of plots to illustrate the capability of the msn
% solver in 2D. Use publish('msn2DTest.m') to store all the figures.

% Authors:  Stefan Pauli, stefan.pauli@alumni.ethz.ch
%           Karthik Jayaraman Raghuram, jrk@ece.ucsb.edu
% v1.0 Created 21-Okt-09
%
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

%% prepare some functions to test the interpolation capabilities.
% runga function
f0 = @(x,y) 1./((1+100*(x.^2 + y.^2)));
% Principles of Mathematical Analysis, Third Edition, Walter Rudin: p.239
% ex.6
f1 = @(x,y) (x.*y)./(x.^2 + y.^2);
% Principles of Mathematical Analysis, Third Edition, Walter Rudin: p.240
% ex.14
f2 = @(x,y) (x.^3)./(x.^2 + y.^2);
% Principles of Mathematical Analysis, Third Edition, Walter Rudin: p.240
% ex.15
f3 = @(x,y) x.^2 + y.^2 - 2.*x.^2.*y - (4*x.^6.*y.^2)./(x.^4 + y.^2).^2;


%% Adjustment
% choose the used function
f = f3;

% the interpolation boundary x in [a,b]
a = -1;
b = 1;

% Determine how many interpolation points are used in the calculations,
maxInterpolPoints = 10^6;

% Determine how many interpolation points are used in the calculations
% withs the old slow methode
maxInterpolPointsOld = 4*10^2;

% Determine how the parameter maxleafsize used in the calculations
maxleafsize =  100;

% Calculate Npoints different values for the plot
Npoints = 10;

% Calculate Npoints different values with the old method for the plot
NpointsOld = 4;

% Use sDefault as the sobolev parameter, determine the sobolev norm to be
% minimized
sDefault = 2;

%% convergence studies of the error

% allocate all the static values
sorted = true;
evaluationPoints = maxInterpolPoints*4;
% choose equally spaced grid points
N = round(exp([0:Npoints-1]*log(maxInterpolPoints)/(Npoints-1)));
N = floor(sqrt(N));    % points per side
N(find(mod(N,2)==1)) = N(find(mod(N,2)==1))+1; % allow just even number of points, avoid x=0
N = eliminateDuplicates(N);
% generate the evaluation points
xev = sort(rand(sqrt(evaluationPoints),1))*(b-a)+a;
[xm,ym]=meshgrid(xev);
indNotY = find(-10^-3< xm.^2+ym.^2 & xm.^2+ym.^2 < 10^-3);
% calcualte the exact y = f(x,y)
yex =  f(xm,ym);
% exclude the values close to the singularity
yex(indNotY) = 0;
clear err

% allocate I,J,S
run = false;
x = linspace(a,b,N(end));
[xm,ym]=meshgrid(x);
fval =  f(xm,ym);
s = 4;
[yslkr, n] = msn_interpol_2D(x,a,b,fval,xev,s,sorted,maxleafsize,run);
I = zeros(n,1);
J = zeros(n,1);
S = zeros(n,1);

% calculate all the values
run = true;
for s = 2:2:4
    for i = 1:length(N)
        i;
        x = linspace(a,b,N(i));
        [xm,ym]=meshgrid(x);
        fval =  f(xm,ym);
        y = msn_interpol_2D(x,a,b,fval,xev,s,sorted,maxleafsize,run,I, J, S);      
        y(indNotY) = 0; % exclude the values close to the singularity
        %figure
        %contour(y-yex)
        err(i,s/2) = max(max( abs(y - yex)));

    end
end
clear I J S y yex

% make the plot
N = N.^2;
a1 = regress(log(err(:,1)),[ones(length(err),1),log(N)']);
a2 = regress(log(err(:,2)),[ones(length(err),1),log(N)']);
loglog(N,err(:,1),N, exp(a1(2)*log(N)+ a1(1)),N,err(:,2),N, exp(a2(2)*log(N)+ a2(1)))
title(['loglog plot of the error for maxleafsize = ',num2str(maxleafsize), ' for f = ',func2str(f)])
legend('s = 2; O(N^-^1)',['measured: O(N^',num2str(a1(2),2),')'],'s = 4; O(N^-^2)',['measured: O(N^',num2str(a2(2),2),')'])
xlabel('number of interpolation points')
ylabel('maximal error at the evaluation points')



%% convergence in time for bigger x

% allocate all the static values
s=sDefault;
sorted = true;
evaluationPoints = 100;
% generate the evaluation points
xev = sort(rand(evaluationPoints,1))*(b-a)+a;
clear size Nnodes tTot tMultiply tComposeS tBackslash tBuildhss
% choose equally spaced grid points
N = round(linspace(1,maxInterpolPoints,Npoints));
N = floor(sqrt(N));    % points per dimension
NOld = round(linspace(1,maxInterpolPointsOld,NpointsOld));
NOld = floor(sqrt(NOld));    % points per dimension

% allocate I,J,S
run=false;
x = linspace(a,b,N(end));
[xm,ym]=meshgrid(x);
fval =  f(xm,ym);
[yslkr, n] = msn_interpol_2D(x,a,b,fval,xev,s,sorted,maxleafsize,run);
n=1.4*n;
I = zeros(n,1);
J = zeros(n,1);
S = zeros(n,1);

% calculate all the values
run=true;
for i = 1:Npoints
    i;
    x = linspace(a,b,N(i));
    [xm,ym]=meshgrid(x);
    fval =  f(xm,ym);
    [y, n, t] = msn_interpol_2D(x,a,b,fval,xev,s,sorted,maxleafsize,run,I, J, S);

    tTot(i) = t.tTot;

    tMultiply(i) = t.multiply;
    tComposeS(i) = t.solve.composeS;
    tLU(i) = t.solve.LU;
    tApplyLU(i) = t.solve.applyLU;
    tBuildhss(i) = t.solve.buildhss;  
end
clear I J S y
for i = 1:NpointsOld
    i;
    x = linspace(a,b,NOld(i));
    [xm,ym]=meshgrid(x);
    fval =  f(xm,ym);
    [xmev,ymev]=meshgrid(xev);
    tic
    [yref]=MSN_interp2D([xm(:),ym(:)], fval(:), [xmev(:),ymev(:)], 2*NOld(i),2*NOld(i),s);
    tOld(i) = toc;
end

% make the plot
N=N.^2;
NOld=NOld.^2;
plot (N,tBuildhss,N,tComposeS,N,tLU,N,tApplyLU,N,tMultiply,N,tTot,NOld,tOld)
title('time')
title(['time for ',num2str(evaluationPoints,2),' evaluation points, maxleafsize = ',num2str(maxleafsize),'and s = ',num2str(s) ])
legend('build HSS matrix for solver part', 'compose sparse matrix S','LU decompostion','apply LU','every thing in the multiply part','total time','total time with old methode','Location','North')
xlabel('number of interpolation points')
ylabel('time [s]')

%% find the optimal maxleafsize
% allocate all the static values
s = sDefault;
findMaxleafsize =  maxleafsize*3;
evaluationPoints = 10^6;
% generate the evaluation points
xev = sort(rand(evaluationPoints,1))*(b-a)+a;
% generate the evaluation points
N = round(linspace(10,findMaxleafsize,Npoints));
findMaxleafsize = N;
sorted = true;

% allocate I,J,S
run=false;
x = linspace(a,b,floor(sqrt(maxInterpolPoints)))';
[xm,ym]=meshgrid(x);
fval =  f(xm,ym);
[yslkr, n] = msn_interpol_2D(x,a,b,fval,xev,s,sorted,maxleafsize,run);
I = zeros(n,1);
J = zeros(n,1);
S = zeros(n,1);

% calculate all the values
fval =  f(xm,ym);
run=true;
for i = 1:Npoints
    i;

    [y, n, t] = msn_interpol_2D(x,a,b,fval,xev,s,sorted,maxleafsize,run,I, J, S);
    tTot(i) = t.tTot;

    tMultiply(i) = t.multiply;
    tComposeS(i) = t.solve.composeS;
    tLU(i) = t.solve.LU;
    tApplyLU(i) = t.solve.applyLU;
    tBuildhss(i) = t.solve.buildhss;  
end
clear I J S y

% make the plot
plot (N,tBuildhss,N,tComposeS,N,tLU,N,tApplyLU,N,tMultiply,N,tTot)
title({'find the optimal maxleafsize',[num2str(maxInterpolPoints,2),' interpolation and ',num2str(evaluationPoints,2),' evaluation points','and s = ',num2str(s)]})
legend('build HSS matrix for solver part', 'compose sparse matrix S','LU decompostion','apply LU','every thing in the multiply part','total time','Location','North')
xlabel('maxleafsize')
ylabel('time [s]')


%% Some small Examples
% Example msn_interpol_2D.m
  f = @(x,y) 1./((1+100*(x.^2 + y.^2)));
  a = -2;
  b = 2;
  x = rand(1,500)*(b-a)+a;
  [xm,ym]=meshgrid(x);
  xev = linspace(a,b,100);
  [xmref,ymref]=meshgrid(xev);
  yref = f(xmref,ymref);
  y = msn_interpol_2D(x,a,b,f(xm,ym),xev);
  max(max(abs(y-yref)))

end



function x = eliminateDuplicates(x)
i=1;
while i <length(x)
   ind = find ( x==x(i));
   x(ind(2:end))=[];
   i=i+1;
end
end

function t = factor(t)
    t = t(2:end) ./ t(1:end-1);
end
