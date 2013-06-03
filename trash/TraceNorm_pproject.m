function [x] = TraceNorm_pproject(x,weights, B,params)
%%%% Force the rows of L and R to have norm at most B.
%
% (LOut, ROut) = arg min_(U,V) || [L;R]-[U;V] ||_F^2 s.t. || [U;V] ||_mr<B
%
% Where || A ||_mr is the maximum Euclidean norm of a row of A.
 
C=sqrt(B/(0.5*norm(x)^2));
x = min(1,C)*x;
end