function [X] = MaxNorm_project(x,weights,B,params)
%%%% Force the rows of L and R to have norm at most B.
%
% (LOut, ROut) = arg min_(U,V) || [L;R]-[U;V] ||_F^2 s.t. || [U;V] ||_mr<B
%
% Where || A ||_mr is the maximum Euclidean norm of a row of A.

% x = [vec(L);vec(R)];
%
e = params.numr*params.nr;
L = x(1:e,:);
R = x(e+1:end,:);
L = reshape(L,params.numr,params.nr);
R = reshape(R,params.numc,params.nr);
row_normsL = sqrt(sum(abs(L).^2,2));
row_normsR = sqrt(sum(abs(R).^2,2));

rescaleL = min(ones(params.numr,1), sqrt(B)./row_normsL);
rescaleR = min(ones(params.numc,1), sqrt(B)./row_normsR);

LOut = L.*(rescaleL*ones(1, size(L,2)));
ROut = R.*(rescaleR*ones(1, size(R,2)));
X = [vec(LOut);vec(ROut)];
end