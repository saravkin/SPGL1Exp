clear all;
L=randn(354,5);
R = randn(708,5);
params.numr=354;
params.numc=708;
params.nr=5;
w=1;
q=params.nr;
X=L*R';
[u,e,v]=svd(X);
n1 = size(u,1);
n2 = size(v,1);
weigh=[w;q;n1;n2;vec(u);vec(v)];
params.weight=weigh;

%% Projection
x=[vec(L);vec(R)];
B=50;
[x1] = TraceNorm_project(x,1, B,params);
e = params.numr*params.nr;
L1 = x1(1:e,:);
R1 = x1(e+1:end,:);
L1 = reshape(L1,params.numr,params.nr);
R1 = reshape(R1,params.numc,params.nr);
[xw] = TraceNorm_projectweight(x,1,B,params);
L2 = xw(1:e); %L = x(1:e,:);
R2 = xw(e+1:end); %R = x(e+1:end,:);
L2 = reshape(L2,params.numr,params.nr);
R2 = reshape(R2,params.numc,params.nr);
L1./L2
R1./R2
norm(x1-xw)
%% primal norm
x=[vec(L);vec(R)];
f=TraceNorm_primal(x,1,[]) 
p = TraceNorm_primalweight(x,1,params)
%% dual testing
X=L*R';
x=vec(X);
d = TraceNorm_dual(x,1, params)
d1 = TraceNorm_dualweight(x,1,params)




