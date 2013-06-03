m=10;
q=5;

L=1e2*rand(m,m);
R=1e2*rand(m,m);
% S=rand(m,m);
% %l=L(:);
% %s=S(:);
w=1;sqrt(0.3);
A=rand(m);
[u,r,e]=qr(A);
Q=u*diag([w*ones(q,1); ones(m-q,1)])*u';
(norm(Q*L,'fro'))
(norm(Q*R,'fro'))
%% SIngle variable
% tau=25;
% n=m;
% cvx_begin
% variables X(n,n)
% minimize(norm((X-L), 'fro'))
% subject to
% norm((Q*X),'fro') <=tau;
% cvx_end
% norm(Q*X,'fro')
% %
% lam = 0;
% g = trace(u*diag([(-2*w^4)*(lam*w^2 + 1)^-3*ones(q,1); -2*(lam+1)^-3*ones(m-q,1)])*u'*(L*L'));
% f = trace(u*diag([(w^2)*(lam*w^2 + 1)^-2*ones(q,1); (lam+1)^-2*ones(m-q,1)])*u'*(L*L')) - tau^2;
% i=0;
% while abs(f) >= 1e-6 
%     i = i+1
%     f = trace(u*diag([(w^2)*(lam*w^2 + 1)^-2*ones(q,1); (lam+1)^-2*ones(m-q,1)])*u'*(L*L')) - tau^2;
%     g = trace(u*diag([(-2*w^4)*(lam*w^2 + 1)^-3*ones(q,1); -2*(lam+1)^-3*ones(m-q,1)])*u'*(L*L'));
%     lam = lam - f/g;
%     [i lam f g];
% end
% 
% % X2 = (eye(m)+lam*Q'*Q)\L;
% X2 = u*diag([(lam*w^2 + 1)^-1*ones(q,1); (lam+1)^-1*ones(m-q,1)])*u'*L;
% X2./L
% X./L
 %% Double variable
tau=25;
n=m;
cvx_begin
variables X(n,n) Y(n,n)
% minimize(norm((X-L),'fro')^2+norm((Y-R), 'fro')^2)
minimize(sum_square(X(:)-L(:)) + sum_square(Y(:)-R(:)))
subject to
%norm(Q*X,'fro')^2+norm(Q*Y,'fro')^2 <=tau;
sum_square(vec(Q*X)) + sum_square(vec(Q*Y)) <=tau;
cvx_end
norm(Q*X,'fro')^2+norm(Q*Y,'fro')^2
%norm(Q*X,'fro')^2
%%
lam = 0;
fl = trace(u*diag([(w^2)*(lam*w^2 + 1)^-2*ones(q,1); (lam+1)^-2*ones(n-q,1)])*u'*(L*L'));
fr = trace(u*diag([(w^2)*(lam*w^2 + 1)^-2*ones(q,1); (lam+1)^-2*ones(m-q,1)])*u'*(R*R'));
f = (fl+fr)-tau;
i=0;
while abs(f) >= 1e-6 
    i  = i+1;
    fl = trace(u*diag([(w^2)*(lam*w^2 + 1)^-2*ones(q,1); (lam+1)^-2*ones(n-q,1)])*u'*(L*L'));
    fr = trace(u*diag([(w^2)*(lam*w^2 + 1)^-2*ones(q,1); (lam+1)^-2*ones(m-q,1)])*u'*(R*R'));
    f  = (fl+fr) - tau;
    gl = trace(u*diag([(-2*w^4)*(lam*w^2 + 1)^-3*ones(q,1); -2*(lam+1)^-3*ones(n-q,1)])*u'*(L*L'));
    gr = trace(u*diag([(-2*w^4)*(lam*w^2 + 1)^-3*ones(q,1); -2*(lam+1)^-3*ones(m-q,1)])*u'*(R*R'));
    g  =(gl+gr);
    lam = lam - f/g;
    [i lam f g];
end

% X2 = (eye(m)+lam*Q'*Q)\L;
L1 = u*diag([(lam*w^2 + 1)^-1*ones(q,1); (lam+1)^-1*ones(n-q,1)])*u'*L;
R1 = u*diag([(lam*w^2 + 1)^-1*ones(q,1); (lam+1)^-1*ones(m-q,1)])*u'*R;
L1./X
R1./Y

