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
for i=1:l
k=D(i,1);
m=D(i,3);
Dinit(k,m)=D(i,5);
end
%%
Dtest=vec(Dinit);
ind=find(Dtest); % keep track of index where data is zero
ind1=randperm(length(ind));
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
rank=[5 10 30 50];
initMult = [0.1 0.01 .001];
options.itermax = 300;
options.tol = 1e-6;
for j=1:length(initMult);
IM=initMult(j);
for i=1:length(rank)
params.nr=rank(i);
Linit=rand(userid,params.nr);
Rinit=rand(movie,params.nr);
xinit  = IM*[vec(Linit);vec(Rinit)]; % Initial guess
funObj = @(x)funCompLBFGS(x, b, @NLfunForward, @funLS, params);
xLS= mylbfgs(funObj,xinit,options);
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
SNR(j,i) = -20*log10(norm(Dorigtest-Drec,'fro')/norm(Dorigtest,'fro'));
RMSE(j,i)=sqrt(sum((Dorigtest-Drec).^2)/length(ind3));
end
end

save('Netflix_lbfgs_crrct.mat','SNR','RMSE');