clear all;
%% Test of primal and dual for distributed
n=400;
A=rand(n);
params.numr=n;
params.numc=n/2;
spmd
dist = codistributor1d(2,[params.numc,params.numc],[params.numr,params.numr]);
AA=codistributed(A,dist);
%B=vec(getLocalPart(AA));
end
p = TraceNorm_primal(AA,1,params)
test=0.5*(norm(AA(:))^2)
params.N=n;
params.M=n;
d = TraceNorm_dual(AA,1, params)
max(svd(A))

%% test of projection for distributed
L=rand(10,5);
R=rand(10,5);
B=.001;
X=[L R];
params.numr=5;
params.nr=10;
spmd
    dist = codistributor1d(1,[5,5],[10,10]);
    Xtest=codistributed(X,dist);
end

[x] = TraceNorm_pproject(Xtest,1, B,params);

x=gather(x);
% L1=x(1:50,1);R1=x(51:end,:);
% L1=reshape(L1,10,5);R1=reshape(R1,10,5);

params.numr=10;
params.nr=5;
xs=[vec(L);vec(R)];
Xs=TraceNorm_project(xs,1, B,params);
L2=Xs(1:50,1);R2=Xs(51:end,:);
L2=reshape(L2,10,5);R2=reshape(R2,10,5);
xt2=[L2 R2];
x./xt2

%% Primal Norm Test
L=rand(10,5);
R=rand(10,5);
B=.001;
params.numr=10;
params.nr=5;
spmd
    dist = codistributor1d(1,[5,5],[10,5]);
    Ltest= codistributed(L,dist);
    Rtest= codistributed(R,dist);
    LLocal=getLocalPart(Ltest);
    LLocal=vec(LLocal);
    RLocal=getLocalPart(Rtest);
    RLocal=vec(RLocal);
    local_x = [LLocal;RLocal];
     new_codistr_vec = codistributor1d(1, [50 50], [100 1]);
     x = codistributed.build(local_x, new_codistr_vec,'noCommunication');
end
[X] = TraceNorm_primal(x,1,params)
0.5*norm([vec(L);vec(R)])^2

%% dual Norm Test
L=rand(10,5);
R=rand(10,5);
B=.001;
params.numr=10;
params.nr=5;
X=L*R';
spmd
    dist = codistributor1d(1,[5,5],[10,5]);
    Ltest= codistributed(L,dist);
    Rtest= codistributed(R,dist);
    LLocal=getLocalPart(Ltest);
    LLocal=vec(LLocal);
    RLocal=getLocalPart(Rtest);
    RLocal=vec(RLocal);
    local_x = [LLocal;RLocal];
     new_codistr_vec = codistributor1d(1, [50 50], [100 1]);
     x = codistributed.build(local_x, new_codistr_vec,'noCommunication');
end
[X] = TraceNorm_primal(x,1,params);
[X] = TraceNorm_pproject(x,1, B,params);
X=gather(X);

%%
params.numr=10;
params.nr=5;
xs=[vec(L);vec(R)];
Xs=TraceNorm_project(xs,1, B,params);
L2=Xs(1:50,1);R2=Xs(51:end,:);
L2=reshape(L2,10,5);R2=reshape(R2,10,5);
xt2=[vec(L2(1:5,:));vec(R2(1:5,:))];
x./xt2