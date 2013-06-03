function [x] = TraceNorm_projectweightnew(x,weights,B,params)
%%%% Force the rows of L and R to have norm at most B.
%
% (LOut, ROut) = arg min_(U,V) || [L;R]-[U;V] ||_F^2 s.t. || [U;V] ||_mr<B
%
% Where || A ||_mr is the maximum Euclidean norm of a row of A.
%

% x is a vector
e = params.numr*params.nr;
L = x(1:e); %L = x(1:e,:);
R = x(e+1:end); %R = x(e+1:end,:);
L = reshape(L,params.numr,params.nr);
R = reshape(R,params.numc,params.nr);
% dual of trace norm is operator norm i.e maximum singular value
w1=params.weight(1);w2=params.weight(2);
w3=params.weight(3);w4=params.weight(4);
W=params.weight(5:end);
QL1=W(w1*w2+w3*w4+1:w1*w2+w3*w4+w1*w2);
RL1=W(w1*w2+w3*w4+w1*w2+1:w1*w2+w3*w4+w1*w2+w3*w4);
QL1=reshape(QL1,w1,w2);
RL1=reshape(RL1,w3,w4);

L1 = QL1*L;
R1 = RL1*R;
numr = size(L1,1);
numc = size(R1,1);

row_normsL = sqrt(sum(L1.^2,2));
row_normsR = sqrt(sum(R1.^2,2));

rescaleL = min(ones(numr,1), sqrt(B)./row_normsL);
rescaleR = min(ones(numc,1), sqrt(B)./row_normsR);

LOut = (L1.*(rescaleL*ones(1, size(L1,2))));
ROut = (R1.*(rescaleR*ones(1, size(R1,2))));

x = [vec(LOut);vec(ROut)];
% E1 = [vec(L1);vec(R1)];
% normLR = norm(E1)^2/2;
% 
% if(normLR > B)
%   x = B*x/normLR;
% end

end