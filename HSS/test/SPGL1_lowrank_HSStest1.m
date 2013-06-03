% addpath(genpath('/Dropbox/spgl1Latest'));
% addpath(genpath('/Dropbox/MaxNorm/GulfofSuez350/code/Function'));
clear all;
addpath(genpath('/users/slic/rkumar/spgl1Latest'));
%% fetch the data from cluster
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
nt = 1024;
nr = 355;
ns = 355;
D=ReadSuFast('/scratch/slic/rkumar/suezdatatest/SuezShots125-355shots.su');
D = reshape(D,nt,nr,ns);  
data = fft(D);
clear D
data = permute(reshape(data,nt,nr,ns),[3 2 1]);
freq = 0; % for low frequency else high frequency
a=250;
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
params.nr = 30; 

initMult=1e-5;
%%
B=b;
B=vec(B);
params.afunT = @(x)reshape(x,nR,nC*2);
params.Ind = find(B==0);
params.afun = @(x)afun(x,params.Ind);
params.numr = nR;
params.numc = 2*nC;
params.nr = 30;
params.funForward = @NLfunForward;
opts = spgSetParms('optTol',1e-5, ...
                   'bpTol', 1e-5,...
                    'decTol',1e-4,...
                   'project', @TraceNorm_project, ...
                   'primal_norm', @TraceNorm_primal, ...
                   'dual_norm', @TraceNorm_dual, ...
                   'proxy', 1, ...
                   'ignorePErr', 1, ...
                   'iterations', 300,...
                   'weights', []);
LInit   = randn(params.numr,params.nr)+1i*randn(params.numr,params.nr);
RInit   = randn(params.numc,params.nr)+1i*randn(params.numc,params.nr);
xinit  = 1e-3*[vec(LInit);vec(RInit)]; % Initial guess
tau = norm(xinit,1);
sigmafact=0.001;
sigma=sigmafact*norm(B,2);
opts.funPenalty = @funLS;
[xLSA,r,g,info] = spgl1(@NLfunForward,B,tau,sigma,xinit,opts,params);
e = params.numr*params.nr;
L1 = xLSA(1:e);
R1 = xLSA(e+1:end);
L1 = reshape(L1,params.numr,params.nr);
R1 = reshape(R1,params.numc,params.nr);
xlsa = L1*R1';
Dtest=SR*vec(D);
Dtest=reshape(Dtest,nR,2*nC);
Ao=Dtest(:,nC-19:nC+20);
A1=xlsa(:,nC-19:nC+20);
SNR1 = -20*log10(norm(Ao-A1,'fro')/norm(Ao,'fro')) 
%% SPGL1 for A
A=b(:,nC-19:nC+20);
params.nr = 20;
params.numr = nR;
params.numc = size(A,2);
params.afunT = @(x)reshape(x,nR,size(A,2));
A=vec(A);
params.Ind = find(A==0);
params.afun = @(x)afun(x,params.Ind);

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
                   'iterations', 300,...
                   'weights', []);

% Run the main algorithm

LInit   = randn(params.numr,params.nr)+1i*randn(params.numr,params.nr);
RInit   = randn(params.numc,params.nr)+1i*randn(params.numc,params.nr);
xinit  = 1e-3*[vec(LInit);vec(RInit)]; % Initial guess
tau = norm(xinit,1);
sigmafact=0.001;
sigma=sigmafact*norm(A,2);
opts.funPenalty = @funLS;
[xLSA,r,g,info] = spgl1(@NLfunForward,A,tau,sigma,xinit,opts,params);
e = params.numr*params.nr;
L1 = xLSA(1:e);
R1 = xLSA(e+1:end);
L1 = reshape(L1,params.numr,params.nr);
R1 = reshape(R1,params.numc,params.nr);
xlsa = L1*R1';
Dtest=SR*vec(D);
Dtest=reshape(Dtest,nR,2*nC);
Ao=Dtest(:,nC-19:nC+20);
SNR1 = -20*log10(norm(Ao-xlsa,'fro')/norm(Ao,'fro')) 
norm(vec(A1-xlsa));
%% SPGL1 for B
B=b(:,1:nC-15);
params.numr = nR;
params.numc = size(B,2);
params.afunT = @(x)reshape(x,nR,size(B,2));
B=vec(B);
params.Ind = find(B==0);
params.afun = @(x)afun(x,params.Ind);
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
                   'iterations', 300,...
                   'weights', []);

% Run the main algorithm

LInit   = randn(params.numr,params.nr)+1i*randn(params.numr,params.nr);
RInit   = randn(params.numc,params.nr)+1i*randn(params.numc,params.nr);
xinit  = 1e-3*[vec(LInit);vec(RInit)]; % Initial guess
tau = norm(xinit,1);
sigmafact=0.001;
sigma=sigmafact*norm(B,2);
opts.funPenalty = @funLS;
[xLSB,r,g,info] = spgl1(@NLfunForward,B,tau,sigma,xinit,opts,params);
e = params.numr*params.nr;
L1 = xLSB(1:e);
R1 = xLSB(e+1:end);
L1 = reshape(L1,params.numr,params.nr);
R1 = reshape(R1,params.numc,params.nr);
xlsb = L1*R1';
Dtest=SR*vec(D);
Dtest=reshape(Dtest,nR,2*nC);
Ab=Dtest(:,1:nC-15);
SNR1 = -20*log10(norm(Ab-xlsb,'fro')/norm(Ab,'fro')) 
%% SPGL1 for C
C=b(:,nC+16:end);
params.numr = nR;
params.numc = size(C,2);
params.afunT = @(x)reshape(x,nR,size(C,2));
C=vec(C);
params.Ind = find(C==0);
params.afun = @(x)afun(x,params.Ind);
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
                   'iterations', 300,...
                   'weights', []);

% Run the main algorithm

LInit   = randn(params.numr,params.nr)+1i*randn(params.numr,params.nr);
RInit   = randn(params.numc,params.nr)+1i*randn(params.numc,params.nr);
xinit  = 1e-3*[vec(LInit);vec(RInit)]; % Initial guess
tau = norm(xinit,1);
sigmafact=0.001;
sigma=sigmafact*norm(C,2);
opts.funPenalty = @funLS;
[xLSC,r,g,info] = spgl1(@NLfunForward,C,tau,sigma,xinit,opts,params);
e = params.numr*params.nr;
L1 = xLSC(1:e);
R1 = xLSC(e+1:end);
L1 = reshape(L1,params.numr,params.nr);
R1 = reshape(R1,params.numc,params.nr);
xlsc = L1*R1';
Dtest=SR*vec(D);
Dtest=reshape(Dtest,nR,2*nC);
Ac=Dtest(:,nC+16:end);
SNR1 = -20*log10(norm(Ac-xlsc,'fro')/norm(Ac,'fro')) 
%%
XX=[xlsb xlsa xlsc]; 
xls2 = reshape(SR'*vec(XX),nR,nC);
SNR1 = -20*log10(norm(D-xls2,'fro')/norm(D,'fro')) 


