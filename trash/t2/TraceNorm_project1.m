function [x] = TraceNorm_project(x,weights, B,params)
%%%% Force the rows of L and R to have norm at most B.
%
% (LOut, ROut) = arg min_(U,V) || [L;R]-[U;V] ||_F^2 s.t. || [U;V] ||_mr<B
%
% Where || A ||_mr is the maximum Euclidean norm of a row of A.
%

% x is a vector
normLR = norm(x)^2/2;

if(normLR > B)
  x = B*x/normLR;
end


% Let's leave this for now - we will probably need it as a template
% when we do the max norm. 

% e = params.numr*params.nr;
% L = x(1:e,:);
% R = x(e+1:end,:);
% L = reshape(L,params.numr,params.nr);
% R = reshape(R,params.numc,params.nr);
% 
% normLR = sqrt(norm(L(:))^2 + norm(R(:))^2);
% 
% if(normLR > B)
% 
%    LOut = B*L/normLR;
%    ROut = B*R/normLR;
% else
%    LOut = L;
%    ROut = R;
% end

%X = [vec(LOut);vec(ROut)];
end