function d = TraceNorm_dual(x,weights, params)

% dual of trace norm is operator norm i.e maximum singular value

E = reshape(x,params.numr,params.numc);
% e = params.numr*params.nr;
% L = x(1:e); %L = x(1:e,:);
% R = x(e+1:end); %R = x(e+1:end,:);
% L = reshape(L,params.numr,params.nr);
% R = reshape(R,params.numc,params.nr);
% E1 = vec(L*R');
% E2 = params.afunT(E1);
d = svds(E,1);

