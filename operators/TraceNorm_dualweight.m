function d = TraceNorm_dualweight(x,weights,params)
E = reshape(x,params.numr,params.numc);

% dual of trace norm is operator norm i.e maximum singular value
w=params.weight(1);
q=params.weight(2);
n=params.weight(3);
m=params.weight(4);
A=params.weight(5:end);
U=A(1:n*n);
V=A(n*n+1:end);
U=reshape(U,n,n);
V=reshape(V,m,m);
Q = U*diag([w^(-1)*ones(q,1); ones(n-q,1)])*U';
W = V*diag([w^(-1)*ones(q,1); ones(m-q,1)])*V';
X=Q*E*W;
% E2 = params.afunT(E1);
d = svds(X,1);
%d
end

