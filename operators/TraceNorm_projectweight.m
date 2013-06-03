function [x] = TraceNorm_projectweight(x,weights,B,params)
%%%% Force the rows of L and R to have norm at most B.
%
% (LOut, ROut) = arg min_(U,V) || [L;R]-[U;V] ||_F^2 s.t. || [U;V] ||_mr<B
%
% Where || A ||_mr is the maximum Euclidean norm of a row of A.
% x is a vector
e = params.numr*params.nr;
L = x(1:e); %L = x(1:e,:);
R = x(e+1:end); %R = x(e+1:end,:);
L = reshape(L,params.numr,params.nr);
R = reshape(R,params.numc,params.nr);
% dual of trace norm is operator norm i.e maximum singular value
w=params.weight(1);q=params.weight(2);n=params.weight(3);m=params.weight(4);
A=params.weight(5:end);
u=A(1:n*n);
v=A(n*n+1:end);
u=reshape(u,n,n);
v=reshape(v,m,m);
Q = u*diag([w*ones(q,1); ones(n-q,1)])*u';
W = v*diag([w*ones(q,1); ones(m-q,1)])*v';
%% checking to do projection or not
normLR = 0.5 * (norm(Q*L,'fro')^2 + norm(W*R,'fro')^2);

if(normLR > B)
% Projection for L and R
lam = 0;
tau=B;

fl = trace((L'*u)*diag([(w^2)*(lam*w^2 + 1)^-2*ones(q,1); (lam+1)^-2*ones(n-q,1)])*(u'*L));
fr = trace((R'*v)*diag([(w^2)*(lam*w^2 + 1)^-2*ones(q,1); (lam+1)^-2*ones(m-q,1)])*(v'*R));
f = 0.5*(fl+fr)-tau;
i=0;
while abs(f) >= 1e-6 
    i = i+1;
    fl = trace((L'*u)*diag([(w^2)*(lam*w^2 + 1)^-2*ones(q,1); (lam+1)^-2*ones(n-q,1)])*(u'*L));
    fr = trace((R'*v)*diag([(w^2)*(lam*w^2 + 1)^-2*ones(q,1); (lam+1)^-2*ones(m-q,1)])*(v'*R));
    f = 0.5*(fl+fr) - tau;
    gl = trace((L'*u)*diag([(-2*w^4)*(lam*w^2 + 1)^-3*ones(q,1); -2*(lam+1)^-3*ones(n-q,1)])*(u'*L));
    gr = trace((R'*v)*diag([(-2*w^4)*(lam*w^2 + 1)^-3*ones(q,1); -2*(lam+1)^-3*ones(m-q,1)])*(v'*R));
    g=(gl+gr);
    lam = lam - f/g;
    [i lam f g];
end

% X2 = (eye(m)+lam*Q'*Q)\L;
L1 = u*diag([(lam*w^2 + 1)^-1*ones(q,1); (lam+1)^-1*ones(n-q,1)])*u'*L;
R1 = v*diag([(lam*w^2 + 1)^-1*ones(q,1); (lam+1)^-1*ones(m-q,1)])*v'*R;
x=[vec(L1);vec(R1)];
end