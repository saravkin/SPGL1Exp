% CVX Test
m = 100;
Z = rand(m); 
% alpha = randn(m,1);alpha=alpha./norm(alpha);
% beta = randn(m,1);beta=beta./norm(beta);
A=ones(m,1);
cvx_begin sdp
    variables t(1) alpha(m) beta(m)
    minimize t
    subject to
        [diag(alpha) Z;Z' diag(beta)] > 0;
        ((A'*alpha)+(A'*beta)) <= t;
        alpha>0;
        beta>0;
cvx_end


%% Max Norm Matlab
U=sqrt(alpha./t);
V=sqrt(beta./t);
T=(diag(U.^-1))'*Z*diag(V.^-1);
Amax=svds(T,1)

%%
n=10;
C=zeros(n,1);
for i=1:n
U=rand(m,1);U=U./(sqrt(2)*norm(U));
V=rand(m,1);V=V./(sqrt(2)*norm(V));
T=(diag(U.^-1))'*Z*diag(V.^-1);
C(i)=svds(T,1);
end
B=2*svds(Z,1);

A=ones(m,1);
A'*diag(U.^2)*A









