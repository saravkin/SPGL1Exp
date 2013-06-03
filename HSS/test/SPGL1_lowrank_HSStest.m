% addpath(genpath('/Dropbox/spgl1Latest'));
% addpath(genpath('/Dropbox/MaxNorm/GulfofSuez350/code/Function'));
addpath(genpath('/users/slic/rkumar/spgl1Latest'));
%% fetch the data from cluster
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
clear all;
close all;
nt = 1024;
nr = 355;
ns = 355;
outlier =0;
D=ReadSuFast('/scratch/slic/rkumar/suezdatatest/SuezShots125-355shots.su');
D = reshape(D,nt,nr,ns);  
data = fft(D);
clear D
data = permute(reshape(data,nt,nr,ns),[3 2 1]);
freq = 0; % for low frequency else high frequency
a=270;
D = squeeze((data(1:354,1:354,a)));
nm = size(D);
nR = nm(2);
nC = nm(1);% length of frequency axis

%% Define Restriction Operator
inds=randperm(nR);
ind=inds(1:floor(nR/2));
R1 = opRestriction(nR,ind);
R2 = opKron(R1,opDirac(nR));
opSR=opSR2MH(354);
info=opSR([],0);
SR=opFunction(info{1},info{2},opSR);
Bfun = @(x)SR*R2'*R2*vec(x);
b = Bfun(D);
b=reshape(b,nR,2*nC);
B=b(:,nR+1:end);
A=B(1:nR/2,1:nR/2);
C=B(nR/2+1:end,1:nR/2);
F=B(1:nR/2,nR/2+1:end);
E=B(nR/2+1:end,nR/2+1:end);
%% SPGL1 for A
params.afunT = @(x)reshape(x,nR/2,nC/2);
A=vec(A);
params.Ind = find(A==0);
params.afun = @(x)afun(x,params.Ind);
params.numr = nR/2;
params.numc = nC/2;
params.nr = 30; 
params.funForward = @NLfunForward;
% options for trace norm
opts = spgSetParms('optTol',1e-5, ...
                   'bpTol', 1e-5,...
                    'decTol',1e-4,...
                   'project', @TraceNorm_project, ...
                   'primal_norm', @TraceNorm_primal, ...
                   'dual_norm', @TraceNorm_dual, ...
                   'proxy', 1, ...
                   'ignorePErr', 1, ...
                   'iterations', 200,...
                   'weights', []);

% Run the main algorithm

LInit   = (1e-2*randn(params.numr,params.nr)+1i*1e-2*randn(params.numr,params.nr));
RInit   = (1e-2*randn(params.numc,params.nr)+1i*1e-2*randn(params.numc,params.nr));
xinit  = [vec(LInit);vec(RInit)]; % Initial guess
tau = norm(xinit,1);
sigma = 2;
opts.funPenalty = @funLS;
[xLSA,r,g,info] = spgl1(@NLfunForward,A,tau,sigma,xinit,opts,params);
e = params.numr*params.nr;
L1 = xLSA(1:e);
R1 = xLSA(e+1:end);
L1 = reshape(L1,params.numr,params.nr);
R1 = reshape(R1,params.numc,params.nr);
xlsa = L1*R1';
%% SPGL1 for A
params.afunT = @(x)reshape(x,nR/2,nC/2);
C=vec(C);
params.Ind = find(C==0);
params.afun = @(x)afun(x,params.Ind);
params.numr = nR/2;
params.numc = nC/2;
params.nr = 30; 
params.funForward = @NLfunForward;
% options for trace norm
opts = spgSetParms('optTol',1e-5, ...
                   'bpTol', 1e-5,...
                    'decTol',1e-4,...
                   'project', @TraceNorm_project, ...
                   'primal_norm', @TraceNorm_primal, ...
                   'dual_norm', @TraceNorm_dual, ...
                   'proxy', 1, ...
                   'ignorePErr', 1, ...
                   'iterations', 200,...
                   'weights', []);

% Run the main algorithm

LInit   = (1e-2*randn(params.numr,params.nr)+1i*1e-2*randn(params.numr,params.nr));
RInit   = (1e-2*randn(params.numc,params.nr)+1i*1e-2*randn(params.numc,params.nr));
xinit  = [vec(LInit);vec(RInit)]; % Initial guess
tau = norm(xinit,1);
sigma = 2;
opts.funPenalty = @funLS;
[xLSC,r,g,info] = spgl1(@NLfunForward,C,tau,sigma,xinit,opts,params);
e = params.numr*params.nr;
L1 = xLSC(1:e);
R1 = xLSC(e+1:end);
L1 = reshape(L1,params.numr,params.nr);
R1 = reshape(R1,params.numc,params.nr);
xlsc = L1*R1';
%% SPGL1 for C
params.afunT = @(x)reshape(x,nR/2,nC/2);
F=vec(F);
params.Ind = find(F==0);
params.afun = @(x)afun(x,params.Ind);
params.numr = nR/2;
params.numc = nC/2;
params.nr = 30; 
params.funForward = @NLfunForward;
% options for trace norm
opts = spgSetParms('optTol',1e-5, ...
                   'bpTol', 1e-5,...
                    'decTol',1e-4,...
                   'project', @TraceNorm_project, ...
                   'primal_norm', @TraceNorm_primal, ...
                   'dual_norm', @TraceNorm_dual, ...
                   'proxy', 1, ...
                   'ignorePErr', 1, ...
                   'iterations', 200,...
                   'weights', []);

% Run the main algorithm

LInit   = (1e-2*randn(params.numr,params.nr)+1i*1e-2*randn(params.numr,params.nr));
RInit   = (1e-2*randn(params.numc,params.nr)+1i*1e-2*randn(params.numc,params.nr));
xinit  = [vec(LInit);vec(RInit)]; % Initial guess
tau = norm(xinit,1);
sigma = 2;
opts.funPenalty = @funLS;
[xLSF,r,g,info] = spgl1(@NLfunForward,F,tau,sigma,xinit,opts,params);
e = params.numr*params.nr;
L1 = xLSF(1:e);
R1 = xLSF(e+1:end);
L1 = reshape(L1,params.numr,params.nr);
R1 = reshape(R1,params.numc,params.nr);
xlsf = L1*R1';
%% SPGL1 for A
params.afunT = @(x)reshape(x,nR/2,nC/2);
E=vec(E);
params.Ind = find(E==0);
params.afun = @(x)afun(x,params.Ind);
params.numr = nR/2;
params.numc = nC/2;
params.nr = 30; 
params.funForward = @NLfunForward;
% options for trace norm
opts = spgSetParms('optTol',1e-5, ...
                   'bpTol', 1e-5,...
                    'decTol',1e-4,...
                   'project', @TraceNorm_project, ...
                   'primal_norm', @TraceNorm_primal, ...
                   'dual_norm', @TraceNorm_dual, ...
                   'proxy', 1, ...
                   'ignorePErr', 1, ...
                   'iterations', 200,...
                   'weights', []);

% Run the main algorithm

LInit   = (1e-2*randn(params.numr,params.nr)+1i*1e-2*randn(params.numr,params.nr));
RInit   = (1e-2*randn(params.numc,params.nr)+1i*1e-2*randn(params.numc,params.nr));
xinit  = [vec(LInit);vec(RInit)]; % Initial guess
tau = norm(xinit,1);
sigma = 2;
opts.funPenalty = @funLS;
[xLSE,r,g,info] = spgl1(@NLfunForward,E,tau,sigma,xinit,opts,params);
e = params.numr*params.nr;
L1 = xLSE(1:e);
R1 = xLSE(e+1:end);
L1 = reshape(L1,params.numr,params.nr);
R1 = reshape(R1,params.numc,params.nr);
xlse = L1*R1';
%%

XX=[xlsa xlsf;xlsc xlse]; 
Xmap=zeros(nR,2*nC);
Xmap(:,nC+1:end)=XX;
Xmap(:,1:nC)=XX(:,end:-1:1);
% imagesc(abs(Xmap))
xls2 = reshape(SR'*vec(Xmap),nR,nC);
SNR1 = -20*log10(norm(D-xls2,'fro')/norm(D,'fro')) 


