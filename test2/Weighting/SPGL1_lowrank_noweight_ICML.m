% addpath(genpath('/Dropbox/spgl1Latest'));
% addpath(genpath('/Dropbox/MaxNorm/GulfofSuez350/code/Function'));
clear all;
close all;
% addpath(genpath('/users/slic/rkumar/spgl1Latest'));

addpath(genpath('/Volumes/Users/rkumar/Dropbox/Research/Low-Rank'));
addpath(genpath('/Volumes/Users/rkumar/Dropbox/Research/Low-Rank/spgl1Latest'));
%% fetch the data from cluster
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
nt = 1024;
nr = 355;
ns = 355;
% outlier =0;
% D=ReadSuFast('/Users/rkumar/research/Data/SuezShots125.355shots.su');
% D = reshape(D,nt,nr,ns);  
% data = fft(D);
% clear D
% save('SuezShots125-355shots-fft.mat','data');
load SuezShots125-355shots-fft.mat;
a=40;
D = squeeze((data(a,1:354,1:354)));
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
%% Function Handle
params.afunT = @(x)reshape(x,nR,nC*2);
params.Ind = find(b==0);
params.afun = @(x)afun(x,params.Ind);
params.numr = nR;
params.numc = 2*nC;
params.nr = 30; 
params.funForward = @NLfunForward;
%% options for trace norm
opts = spgSetParms('optTol',1e-5, ...
                   'bpTol', 1e-5,...
                    'decTol',1e-4,...
                   'project', @TraceNorm_project, ...
                   'primal_norm', @TraceNorm_primal, ...
                   'dual_norm', @TraceNorm_dual, ...
                   'proxy', 1, ...
                   'ignorePErr', 1, ...
                   'iterations', 150,...
                   'weights', []);

%% Run the main algorithm
LInit   = randn(params.numr,params.nr)+1i*randn(params.numr,params.nr);
RInit   = randn(params.numc,params.nr)+1i*randn(params.numc,params.nr);
xinit  = 1e-3*[vec(LInit);vec(RInit)]; % Initial guess
tau = norm(xinit,1);
a=40:70;
for i=1:length(a)
    i
    D = squeeze((data(a(i),1:354,1:354)));
    b = Bfun(D);
    params.Ind = find(b==0);
    params.afun = @(x)afun(x,params.Ind);
    sigmafact=0.01;
    sigma=sigmafact*norm(b,2);
    opts.funPenalty = @funLS;
    [xLS2,r,g,info] = spgl1(@NLfunForward,b,tau,sigma,xinit,opts,params);
    e = params.numr*params.nr;
    L1 = xLS2(1:e);
    R1 = xLS2(e+1:end);
    L1 = reshape(L1,params.numr,params.nr);
    R1 = reshape(R1,params.numc,params.nr);
    xls2 = reshape(SR'*vec(L1*R1'),nR,nC);
    SNR1(i) = -20*log10(norm(D-xls2,'fro')/norm(D,'fro')) 
end
save('icmlnoweighted.mat','SNR1');
