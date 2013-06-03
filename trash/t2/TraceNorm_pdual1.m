function d = TraceNorm_pdual(x,weights, params)

% dual of trace norm is operator norm i.e maximum singular value

E = reshape(x,params.numr,params.numc);
spmd
    X_Local=norm(x);
    c=gcat(X_Local);
    normx_global=norm(c);
end
d = max(svd(E));

