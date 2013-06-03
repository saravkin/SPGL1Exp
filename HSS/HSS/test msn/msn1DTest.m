function msn1DTest
% function to make a lot of plots to illustrate the capability of the msn
% solver in 1D. Use publish('msn1DTest.m') to store all the figures.

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
f0 = @(x) 1./(1+100*x.^2);
% x = abs(x) if x<0,  x = sqrt(x) if x>=0
f1 = @(x)[-x(find(x<0));sqrt(x(find(x>=0)))];

%% Adjustment
% choose the used function
f = f0;

% the interpolation boundary x in [a,b]
a = -1;
b = 1;

% Determine how many interpolation points are used in the calculations
maxInterpolPoints=10^3;

% Determine how many interpolation points are used in the calculations
% withs the old slow methode
maxInterpolPointsOld = 2.5*10^2;

% Determine how the parameter maxleafsize used in the calculations
maxleafsize =  150;

% Calculate Npoints different values for the plot
Npoints = 20;

% Calculate Npoints different values with the Old Methode for the plot
NpointsOld = 6;

% Use sDefault as the sobolev parameter, determine the sobolev norm to be
% minimized
sDefault = 2;

% Try to improve the accuracy of the interpolation by maxrefinement
% iterative refinements in the 'test the iterative refinement' plot
maxrefinement = 10;


%% convergence studies of the error

% allocate all the static values
evaluationPoints = maxInterpolPoints*2;
% choose equally spaced grid points
N = round(exp([0:Npoints-1]*log(maxInterpolPoints)/(Npoints-1)));
N = eliminateDuplicates(N);
% generate the evaluation points
xev = sort(rand(evaluationPoints,1))*(b-a)+a;
% calcualte the exact y = f(x)
yex= f(xev);
clear err

% allocate I,J,S
x = linspace(a,b,maxInterpolPoints)';
[y, n] = msn_interpol_1D(x,a,b,f(x),xev,2,true,maxleafsize,false);
I = zeros(n,1);
J = zeros(n,1);
S = zeros(n,1);

% calculate all the values
for s = 1:2
    for i = 1:length(N)
        i;
        x = linspace(a,b,N(i))';
        y = msn_interpol_1D(x,a,b,f(x),xev,s,true,maxleafsize,true,I, J, S);
        err(i,s) = max( abs(y - yex));
    end
end
clear I J S y yex

% make the plot
a1 = regress(log(err(:,1)),[ones(length(err),1),log(N)']);
a2 = regress(log(err(:,2)),[ones(length(err),1),log(N)']);
loglog(N,err(:,1),N, exp(a1(2)*log(N)+ a1(1)),N,err(:,2),N, exp(a2(2)*log(N)+ a2(1)))
title({['loglog plot of the error for the function f = ',func2str(f)],['maxleafsize = ',num2str(maxleafsize)]})
legend('s = 1; O(N^-^1)',['fit: O(N^',num2str(a1(2),2),')'],'s = 2; O(N^-^2)',['fit: O(N^',num2str(a2(2),2),')'])
xlabel('number of  interpolation points')
ylabel('maximal reconstruction error')


%% convergence in time for a big number of interpolation points
% allocate all the static values
s=sDefault;
evaluationPoints = 1000;
% choose equally spaced grid points
N = round(linspace(1,maxInterpolPoints,Npoints));
NOld = round(linspace(1,maxInterpolPointsOld,NpointsOld));
% generate the evaluation points
xev = sort(rand(evaluationPoints,1))*(b-a)+a;
clear size Nnodes tTot tMultiply tComposeS tBackslash tBuildhss



% allocate I,J,S
run=false;
x = linspace(a,b,N(end))';
[y,n] = msn_interpol_1D(x,a,b,f(x),xev,s,true,maxleafsize,run);
n=1.4*n;
I = zeros(n,1);
J = zeros(n,1);
S = zeros(n,1);

% calculate all the values with the new methode
run=true;
for i = 1:Npoints
    i;
    x = linspace(a,b,N(i))';

    [y, n, t] = msn_interpol_1D(x,a,b,f(x),xev,s,true,maxleafsize,run,I, J, S);
    tTot(i) = t.tTot;

    tMultiply(i) = t.multiply;
    tComposeS(i) = t.solve.composeS;
    tBackslash(i) = t.solve.backslash;
    tBuildhss(i) = t.solve.buildhss;
end
clear I J S y
% calculate all the values with the old methode
for i = 1:NpointsOld
    i;
    x = linspace(a,b,NOld(i))';

    tic
    [yref]=MSN_interp_1D(x, a, b, f(x), xev);
    tOld(i) = toc;
end

% make the plot
plot (N,tBuildhss,N,tComposeS,N,tBackslash,N,tMultiply,N,tTot,NOld,tOld)
title({['time for ',num2str(evaluationPoints,2),' evaluation points'], ['maxleafsize = ',num2str(maxleafsize),' and s = ',num2str(s)]})
legend('build HSS matrix for solver part', 'compose sparse matrix S','solve Sx=f for x','every thing in the multiply part','total time','total time with old method','Location','North')
xlabel('number of  interpolation points')
ylabel('time [s]')
clear size Nnodes tTot tMultiply tComposeS tBackslash tBuildhss

%% find the optimal maxleafsize
% allocate all the static values
s=sDefault;
findMaxleafsize =  maxleafsize*3;
evaluationPoints = 1000;
% generate the evaluation points
xev = sort(rand(evaluationPoints,1))*(b-a)+a;
% generate the evaluation points
N = round(linspace(30,findMaxleafsize,Npoints));
findMaxleafsize = N;
x = linspace(a,b,maxInterpolPoints)';

% allocate I,J,S
run=false;
[y,n] = msn_interpol_1D(x,a,b,f(x),xev,s,true,findMaxleafsize(end),run);
I = zeros(n,1);
J = zeros(n,1);
S = zeros(n,1);

% calculate all the values
run=true;
for i = 1:Npoints
    i;

    [y, n, t] = msn_interpol_1D(x,a,b,f(x),xev,s,true,findMaxleafsize(i),run,I, J, S);
    tTot(i) = t.tTot;

    tMultiply(i) = t.multiply;
    tComposeS(i) = t.solve.composeS;
    tBackslash(i) = t.solve.backslash;
    tBuildhss(i) = t.solve.buildhss;
end
clear I J S y

% make the plot
plot (N,tBuildhss,N,tComposeS,N,tBackslash,N,tMultiply,N,tTot)
title({'find the optimal maxleafsize',[num2str(maxInterpolPoints,2),' interpolation and ',num2str(evaluationPoints,2),' evaluation points','and s = ',num2str(s)]})
legend('build HSS matrix for solver part', 'compose sparse matrix S','solve Sx=f for x','every thing in the multiply part','total time','total time with old method','Location','North')
xlabel('maxleafsize')
ylabel('time [s]')



%% test the iterative refinement
% allocate all the static values
evaluationPoints = maxInterpolPoints*2;
s=sDefault;
refinement = [0:maxrefinement];
% generate the evaluation points
xev = sort(rand(evaluationPoints,1))*(b-a)+a;
% calcualte the exact y = f(x)
yex= f(xev);
clear err

% allocate I,J,S
x = linspace(a,b,maxInterpolPoints)';
[y, n] = msn_interpol_1D(x,a,b,f(x),xev,2,true,maxleafsize,false);
I = zeros(n,1);
J = zeros(n,1);
S = zeros(n,1);

% generate the interpolation points
x = linspace(a,b,maxInterpolPoints)';

% calculate all the values
for i = 1:length(refinement)
    i;

    y = msn_interpol_1D(x,a,b,f(x),xev,s,true,maxleafsize,true,I, J, S,refinement(i));
    err(i) = max( abs(y - yex));

end
clear I J S y yex

% make the plot
plot(err);
title({'loglog plot of the error ',[num2str(maxInterpolPoints,2),' interpolation points, ',num2str(evaluationPoints,2),' evaluation points',' and s = ',num2str(s)]})
xlabel('number of iterative refinements')
ylabel('maximal reconstruction error')

%% Some small Examples

% Example hss_buildA.m, hss_mul.m
m = 10;
n = 11;
x = rand(n,1);
A = rand(m,n);
hss = hss_buildA(A);
b = hss_mul(hss,x);
max(abs(A*x-b))

if max(abs(A*x-b)) > 10^-5
    error('Error to big')
end

% Example hss_solve.m
m = 10;
n = 11;
b = rand(m,1);
A = rand(m,n);
hss = hss_buildA(A);
x = hss_solve(hss,b);
max(abs(A*x-b))

if max(abs(A*x-b)) > 10^-5
    error('Error to big')
end

% Example msn_interpol_1D.m
f = @(x) 1./(1+100*x.^2);
a = -1;
b = 1;
x = sort(rand(1,100)*(b-a)+a);
xev = linspace(a,b);
y = msn_interpol_1D(x,a,b,f(x),xev);
plot(linspace(a,b), f(linspace(a,b)),'-',x,f(x),'o', xev,y,'--')

if max(abs(f(xev)-y')) > .5
    error('Error to big')
end

% Example msn_interpol_1D.m
f = @(x) 1./(1+100*x.^2);
a = -1;
b = 1;
x = sort(rand(1,100)*(b-a)+a);
xev = linspace(a,b);
maxleafsize = 30;
[y,n] = msn_interpol_1D(x,a,b,f(x),xev,2,true,maxleafsize,false);
I = zeros(n,1);
J = zeros(n,1);
S = zeros(n,1);
y = msn_interpol_1D(x,a,b,f(x),xev,2,true,maxleafsize,true,I,J,S);
plot(linspace(a,b), f(linspace(a,b)),'-',x,f(x),'o', xev,y,'--')

if max(abs(f(xev)-y')) > .5
    error('Error to big')
end

% Example hss_lu.m
m = 10;
n = 11;
b = rand(m,1);
A = rand(m,n);
hss = hss_buildA(A);
[L, U, indb, indx]= hss_lu(hss,b);
bHSS = zeros(size(indb));
bHSS(indb) = b;
x = U\(L\bHSS);
x = x(indx);
max(abs(A*x-b))

if max(abs(A*x-b)) > 10^-5
    error('Error to big')
end

% Example hss_toSparse.m
m = 10;
n = 11;
b = rand(m,1);
A = rand(m,n);
hss = hss_buildA(A);
[S, indb, indx] = hss_toSparse(hss,b);
bHSS =  sparse(size(indb,1),size(b,2));
bHSS(indb,:)=b;
x=S\bHSS;
x = x(indx,:);
max(abs(A*x-b))

if max(abs(A*x-b)) > 10^-5
    error('Error to big')
end

% test Example
m = 10;
n = 11;
b = rand(m,2);
A = rand(m,n);
hss = hss_buildA(A);
x = hss_solve(hss,b);
max(abs(A*x(:,1)-b(:,1)))

if max(abs(A*x-b)) > 10^-5
    error('Error to big')
end

% test Example
f = @(x) 1./(1+100*x.^2);
a = -1;
b = 1;
x = (rand(1,100)*(b-a)+a);
xev = (rand(1,100)*(b-a)+a);
y = msn_interpol_1D(x,a,b,f(x),xev);
plot(linspace(a,b), f(linspace(a,b)),'-',x,f(x),'o', xev,y,'o')

if max(abs(f(xev)-y')) > .5
    error('Error to big')
end

% test Example
f = @(x) 1./(1+100*x.^2);
a = -1;
b = 1;
x = sort(rand(1,100)*(b-a)+a);
xev = linspace(a,b);
maxleafsize = 30;
[y,n] = msn_interpol_1D(x,a,b,f(x),xev,1,true,maxleafsize,false);
I = zeros(n,1);
J = zeros(n,1);
S = zeros(n,1);
y = msn_interpol_1D(x,a,b,f(x),xev,1,true,maxleafsize,true,I,J,S);
plot(linspace(a,b), f(linspace(a,b)),'-',x,f(x),'o', xev,y,'--')

if max(abs(f(xev)-y')) > .5
    error('Error to big')
end



end


function x = eliminateDuplicates(x)
i=1;
while i <length(x)
    ind = find ( x==x(i));
    x(ind(2:end))=[];
    i=i+1;
end
end
