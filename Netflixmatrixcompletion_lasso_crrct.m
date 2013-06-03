clear all;
addpath(genpath('/users/slic/rkumar/spgl1Latest'));
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
% read the rating data
D=dlmread('ratings.dat');
userid=6040;
movie=3952;
%Define the intial incomplete matrix
Dinit=zeros(userid,movie);
l=length(D(:,1));
% fill the matrix with the given data
for Z=1:l
k=D(Z,1);
m=D(Z,3);
Dinit(k,m)=D(Z,5);
end
%%
Dtest=vec(Dinit);
ind=find(Dtest); % keep track of index where data is zero
ind1=randperm(length(ind));
%% data Initilization
ind2=ind1(1:floor(length(ind1)/2));
ind3=ind(ind2);
Dtest(ind3)=0;
b=vec(Dtest);
% Function Handle
params.afunT = @(x)reshape(x,userid,movie);
params.Ind = find(b==0);
params.afun = @(x)afun(x,params.Ind);
params.numr = userid;
params.numc = movie;
params.funForward = @NLfunForward;
opts = spgSetParms('optTol',1e-6, ...
'bpTol', 1e-6,...
'decTol',1e-4,...
'project', @TraceNorm_project, ...
'primal_norm', @TraceNorm_primal, ...
'dual_norm', @TraceNorm_dual, ...
'proxy', 1, ...
'ignorePErr', 1, ...
'iterations', 500);
rank=[5 10 30 50];
TAU=[1e4 1e5 1e8];

for k=1:length(TAU)
tau=TAU(k);
for i=1:length(rank)
params.nr=rank(i);
Linit=rand(userid,params.nr);
Rinit=rand(movie,params.nr);
xinit  = 0.001*[vec(Linit);vec(Rinit)]; % Initial guess
sigma=0.2*norm(b,2);
opts.funPenalty = @funLS;
[xLS,r1,g1,info] = spgl1(@NLfunForward,b,tau,sigma,xinit,opts,params);
e = params.numr*params.nr;
L1 = xLS(1:e);
R1 = xLS(e+1:end);
L1 = reshape(L1,params.numr,params.nr);
R1 = reshape(R1,params.numc,params.nr);
xls = L1*R1';
xls=vec(xls);
Drec=xls(ind3);
Dorig=vec(Dinit);
Dorigtest=Dorig(ind3);
SNR(k,i) = -20*log10(norm(Dorigtest-Drec,'fro')/norm(Dorigtest,'fro'));
RMSE(k,i)=sqrt(sum((Dorigtest-Drec).^2)/length(ind3));
Tau(k,i)= info.tau;
iter(k,i)=info.iter;
end
end

save('Netflix_lasso_crrct.mat','SNR','RMSE','Tau','iter');


