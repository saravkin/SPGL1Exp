function d = TraceNorm_pdual(x,weights, params)

% dual of trace norm is operator norm i.e maximum singular value
spmd
E = reshape(x,params.numr,params.numc); 
Dist=codistributor1d(2, [params.numc params.numc],[params.N,params.M]);
Eglobal = codistributed.build(E,Dist);
end
E=gather(Eglobal);
d=max(svd(E));

