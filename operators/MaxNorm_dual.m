function D = MaxNorm_dual(x,weights, params)
% dual of max norm is one norm of Euclidean norm of a row of A.
% i.e  ||A||2,1

T = zeros(params.numr,1);
E = reshape(x,params.numr,params.numc);

for i=1:params.maxSample
U=rand(params.numr,1)+1i*rand(params.numr,1);U=U./norm(U);
V=rand(params.numc,1)+1i*rand(params.numc,1);V=V./norm(V);
T=(diag(U.^-1))'*E*diag(V.^-1);
d(i)=svds(T,1);
end
D=min(d);
d2 = svds(E, 1);

D = sqrt(D*d2);