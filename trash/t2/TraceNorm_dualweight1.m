function d = TraceNorm_dualweight(x,weights,params)
E = reshape(x,params.numr,params.numc);

% dual of trace norm is operator norm i.e maximum singular value
w1=params.weight(1);w2=params.weight(2);
w3=params.weight(3);w4=params.weight(4);
W=params.weight(5:end);
QL=W(1:w1*w2);
RL=W(w1*w2+1:w1*w2+w3*w4);
QL=reshape(QL,w1,w2);
RL=reshape(RL,w3,w4);
X=QL*E*RL;
% E2 = params.afunT(E1);
d = svds(X,1);
end

