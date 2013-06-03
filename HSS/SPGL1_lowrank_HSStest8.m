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
D = squeeze((data(1:352,1:352,a)));
nm = size(D);
nR = nm(2);
nC = nm(1);% length of frequency axis
opSR=opSR2MH(nR);
info=opSR([],0);
SR=opFunction(info{1},info{2},opSR);
%% Define Restriction Operator
inds=randperm(nR);
ind=inds(1:floor(nR/2));
R1 = opRestriction(nR,ind);
R2 = opKron(R1,opDirac(nR));
Bfun = @(x)R2'*R2*vec(x);
b = Bfun(D);
b=reshape(b,nR,nC);
A=SR*vec(b);
A=reshape(A,nR,2*nC);
params.afunT = @(x)reshape(x,size(A,1),size(A,2));
params.numr = size(A,1);
params.numc = size(A,2);
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
                   'iterations', 150,...
                   'weights', []);

% Run the main algorithm
ra =[5:5:20];
for i=1:length(ra)
  params.nr = ra(i);  
LInit   = randn(params.numr,params.nr)+1i*randn(params.numr,params.nr);
RInit   = randn(params.numc,params.nr)+1i*randn(params.numc,params.nr);
xinit  = [vec(LInit);vec(RInit)]; % Initial guess
tau = norm(xinit,1);
sigmafact=0.001;
sigma=sigmafact*norm(A,2);
opts.funPenalty = @funLS;
tic;
[xLRA,r,g,info] = spgl1(@NLfunForward,A,tau,sigma,xinit,opts,params);
toc;
e = params.numr*params.nr;
L1 = xLRA(1:e);
R1 = xLRA(e+1:end);
L1 = reshape(L1,params.numr,params.nr);
R1 = reshape(R1,params.numc,params.nr);
xlra = SR'*vec(L1*R1');
xlra=reshape(xlra,nR,nC);
SNR(i) = -20*log10(norm(D-xlra,'fro')/norm(D,'fro')) 
end
