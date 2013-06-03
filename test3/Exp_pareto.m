%%  Test script to compare the residual norm versus Tau

%%
clear all;
clc;
%% Solving the factorized formulation in LASSO mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimize   ||A(LR') - b ||_2  subject to ||LR'||_*<=tau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
addpath(genpath('/users/slic/rkumar/solver/tuningspgl1/spgl1'));
addpath(genpath('/users/slic/rkumar/Self_Function'));
% Generate rank=70 matrix 
    sigmas = [1e2*rand(50,1); zeros((100-50), 1)];
    % set up U and V
    U = randn(100, 50);
    [Q R] = qr(U);
    U = Q;
    V = randn(100,50);
    [Q R] = qr(V);
    V = Q;
    % Form A
    D = U*diag(sigmas)*V';
    opts = spgSetParms('project', @TraceNorm_project, ...
        'primal_norm', @TraceNorm_primal, ...
        'dual_norm', @TraceNorm_dual, ...
        'proxy', 1, ...
        'verbosity',1,...
        'ignorePErr', 1, ...
        'weights', []);
    nm = size(D);
    nR = nm(2);
    nC = nm(1);
    [m,n] = size(D);
    idx = randperm(m*n);
    k   = min(m*n,max(1,round(m*n*60/100)));
    idx = sort(idx(1:k));
    op  = opRestriction(m*n,idx);
    Bfun = @(x)op'*op*vec(x);
    params.numr = nR;
    params.numc = nC;
    params.funForward = @NLfunForward;
    b = Bfun(D);
    params.afunT = @(x)reshape(x,nR,nC);
    params.Ind = find(b==0);
    params.afun = @(x)afun(x,params.Ind);
    rank = [1,2,5,10,20,30,50,70,80];
    tau = [1:2000]; 
    opts.funPenalty = @funLS;
for i = 1:length(rank)
    r = rank(i);   
    params.nr=r;
    LInit   = randn(params.numr,params.nr);
    RInit   = randn(params.numc,params.nr);
    xinit   = [vec(LInit);vec(RInit)]; % Initial guess
    for j = 1:length(tau)
        Tau = tau(j);
        % =====================================================================
        [xS,resi,grad,info] = spgl1(@NLfunForward,b,Tau,[],xinit,opts,params);
        % =====================================================================
        e = params.numr*params.nr;
        L1 = xS(1:e);
        R1 = xS(e+1:end);
        L1 = reshape(L1,params.numr,params.nr);
        R1 = reshape(R1,params.numc,params.nr);
        xs = L1*R1';
        tau_proxy(i,j) = 0.5*norm([L1; R1], 'fro')^2;
        tau_calcul(i,j) = sum(svd(xs));
        SNR_fact(i,j) = -20*log10(norm(D-xs,'fro')/norm(D,'fro'));
        Data_Fact{i,j}= info;
        value(i,j) = norm(resi,2);
        deri(i,j)  = norm(resi,2)*norm(grad,inf);
        i
        j
    end
    % can you compute A^T (A (LR') - b)? 
    
end

save('ExpPareto.mat','SNR_fact','Data_Fact','value','deri','tau_proxy','tau_calcul');
%%
load ExpPareto.mat
tau = [1:2000];
figure,
plot(tau,value(1,:),'*r',tau,value(2,:),'b',tau,value(3,:),'*k',tau,value(4,:),'*g',...
    tau,value(5,:),'r',tau,value(6,:),'y',tau,value(7,:),'cyan',tau,value(8,:),'g',...
    tau,value(9,:),'k');
legend('1','2','5','10','20','30','50','70','80');
xlabel('Tau');
ylabel('value=norm(residual,2)');

