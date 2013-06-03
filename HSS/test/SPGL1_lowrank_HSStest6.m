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
Column=[30 40 50 60 80 100 200 300 400 500 600 708];
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
    C=Column(i);
    C=floor(C/2);
A=b(:,nC-C+1:nC+C);
params.afunT = @(x)reshape(x,nR,size(A,2));
params.numr = nR;
params.numc = size(A,2);
A=vec(A);
params.Ind = find(A==0);
params.afun = @(x)afun(x,params.Ind);
params.funForward = @NLfunForward;
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
Ao=Dtest(:,nC-C+1:nC+C);
Al=vec(Ao);
indm=find(Al==0);
if ~isempty(indm)
    xlsa=vec(xlsa);
    xlsa(indm)=0;
    xlsa=reshape(xlsa,params.numr,params.numc);
end
s1=svd(xlsa);
s2=svd(Ao);
SVD1(1:length(s1),i)=s1;
SVD2(1:length(s2),i)=s2;
SNR1(i) = -20*log10(norm(Ao-xlsa,'fro')/norm(Ao,'fro')) 

