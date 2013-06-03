function p = MaxNorm_primal(x,weights,params)

T1 = zeros(params.numr,1);
T2 = zeros(params.numc,1);
e = params.numr*params.nr;
L = x(1:e,:);
R = x(e+1:end,:);
L = reshape(L,params.numr,params.nr);
R = reshape(R,params.numc,params.nr);


L=L.*L;
R=R.*R;
Lmax=max(sqrt(sum(transp(L))));
Rmax=max(sqrt(sum(transp(R))));

p = max(Lmax , Rmax);