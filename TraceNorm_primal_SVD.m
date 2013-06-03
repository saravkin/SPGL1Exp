function xt = TraceNorm_primal_SVD(x,params)
e = params.numr*params.nr;
L = x(1:e,:);
R = x(e+1:end,:);
L = reshape(L,params.numr,params.nr);
R = reshape(R,params.numc,params.nr);
xt=L*R';
xt=sum(svd(xt));
end

