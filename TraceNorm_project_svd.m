function [xout] = TraceNorm_project_svd(x,weights, B,params)
%%%% Force the rows of L and R to have norm at most B.
%
% (LOut, ROut) = arg min_(U,V) || [L;R]-[U;V] ||_F^2 s.t. || [U;V] ||_mr<B
%
% Where || A ||_mr is the maximum Euclidean norm of a row of A.
%
e = params.numr*params.nr;
L = x(1:e,:);
R = x(e+1:end,:);
L = reshape(L,params.numr,params.nr);
R = reshape(R,params.numc,params.nr);
xt=L*R';
xt=sum(svd(xt));
c=sqrt(B/(xt));
LOut = min(1,c)*L(:);
ROut = min(1,c)*R(:);
% xout = min(1,c)*x;
xout = [vec(LOut);vec(ROut)];
end