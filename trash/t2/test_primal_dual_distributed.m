clear all;
L=rand(100,30);
R=rand(100,30);
params.numr=25;
params.nr=30;
X=L*R';
B=0.01;
spmd
    dist = codistributor1d(1,[25,25,25,25],[100,30]);
    Ltest= codistributed(L,dist);
    Rtest= codistributed(R,dist);
    LLocal=getLocalPart(Ltest);
    LLocal=vec(LLocal);
    RLocal=getLocalPart(Rtest);
    RLocal=vec(RLocal);
    local_x = [LLocal;RLocal];
     new_codistr_vec = codistributor1d(1, [1500 1500 1500 1500], [6000 1]);
     x = codistributed.build(local_x, new_codistr_vec,'noCommunication');
end
%projection test
Xtest=[vec(L);vec(R)];
[Xp] = TraceNorm_project(Xtest,1,B,params);
[Xpp] = TraceNorm_pproject(x,1, B,params);
norm(Xp)-norm(Xpp)
% primal norm test
[Xpp1] = TraceNorm_primal(x,1,params);
[Xp1] = TraceNorm_primal(Xtest,1,params);
norm(Xpp1)-norm(Xp1)