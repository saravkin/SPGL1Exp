function d = TraceNorm_pdual(x,weights, params)

% dual of trace norm is operator norm i.e maximum singular value
spmd
E = reshape(x,params.N,params.M); 
dtest=max(svd(E));
end
d=dtest;


