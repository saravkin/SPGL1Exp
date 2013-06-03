%% Test script for the implementation of General SPG

% open matlab from qmatlabX because of large data file used in this example

% update here the location of spgl2 folder
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
% %for server
D=ReadSuFast('/scratch/slic/rkumar/suezdatatest/SuezShots125-355shots.su');
D = reshape(D,nt,nr,ns);  
data = fft(D);
clear D
data = permute(reshape(data,nt,nr,ns),[3 2 1]);
freq = 0; % for low frequency else high frequency

if freq== 0
    D = squeeze((data(1:354,1:354,40)));
else
    D = squeeze((data(1:354,1:354,280)));
end
clear data;
% load Data_Suez_freq10.mat;
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
%params.F = @(x)R2*SR'*vec(x);
params.Ind = find(b==0);
params.afun = @(x)afun(x,params.Ind);
params.numr = nR;
params.numc = 2*nC;
params.nr = 20; 
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
                   'iterations', 100,...
                   'weights', []);

%% Run the main algorithm

LInit   = (1e-2*randn(params.numr,params.nr)+1i*1e-2*randn(params.numr,params.nr));
RInit   = (1e-2*randn(params.numc,params.nr)+1i*1e-2*randn(params.numc,params.nr));
xinit  = [vec(LInit);vec(RInit)]; % Initial guess
tau = norm(xinit,1);
sigma = 6;
opts.funPenalty = @funLS;
tic;
[xLS,r,g,info] = spgl1(@NLfunForward,b,tau,sigma,xinit,opts,params);
toc;
e = params.numr*params.nr;
L1 = xLS(1:e);
R1 = xLS(e+1:end);
L1 = reshape(L1,params.numr,params.nr);
R1 = reshape(R1,params.numc,params.nr);
xls = reshape(SR'*vec(L1*R1'),nR,nC);
SNR = -20*log10(norm(D-xls,'fro')/norm(D,'fro'))

%%
q = params.nr; % low rank approximation
w = sqrt(0.3); % set weight value
xls=L1*R1';
[U E V] = svd(xls);
n1 = size(U,1);
n2 = size(V,1);
weigh=[w;q;n1;n2;vec(U);vec(V)];
params.weight=weigh;
%%
xinit  = [vec(LInit);vec(RInit)]; % Initial guess
tau = norm(xinit,1);
clear opts;
opts = spgSetParms('optTol',1e-5, ...
                   'bpTol', 1e-5,...
                    'decTol',1e-4,...
                   'project', @TraceNorm_projectweight, ...
                   'primal_norm', @TraceNorm_primalweight, ...
                   'dual_norm', @TraceNorm_dualweight, ...
                   'proxy', 1, ...
                   'ignorePErr', 1, ...
                   'iterations', 100,...
                   'weights', []);
opts.funPenalty = @funLS;
tic;
[xLS2,r1,g1,info1] = spgl1(@NLfunForward,b,tau,sigma,xinit,opts,params);
toc
e = params.numr*params.nr;
L2 = xLS2(1:e);
R2 = xLS2(e+1:end);
L2 = reshape(L2,params.numr,params.nr);
R2 = reshape(R2,params.numc,params.nr);
xls2 = reshape(SR'*vec(L2*R2'),nR,nC);
SNR1 = -20*log10(norm(D-xls2,'fro')/norm(D,'fro')) 
