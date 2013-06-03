function p = TraceNorm_primalweight(x,weights,params)


e = params.numr*params.nr;
L = x(1:e); %L = x(1:e,:);
R = x(e+1:end); %R = x(e+1:end,:);
L = reshape(L,params.numr,params.nr);
R = reshape(R,params.numc,params.nr);

w=params.weight(1);q=params.weight(2);n=params.weight(3);m=params.weight(4);
A=params.weight(5:end);
U=A(1:n*n);
V=A(n*n+1:end);
U=reshape(U,n,n);
V=reshape(V,m,m);
Q = U*diag([w*ones(q,1); ones(n-q,1)])*U';
W = V*diag([w*ones(q,1); ones(m-q,1)])*V';

L1 = Q*L;
R1 = W*R;

E1 = [vec(L1);vec(R1)];
p = 0.5 * norm(E1)^2; 
%p
end