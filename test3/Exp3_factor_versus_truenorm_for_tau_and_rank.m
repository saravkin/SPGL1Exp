%%  Test script to compare the classical versis factorized formulation
% For a given matrix with predefined rank, solve the factorized formulation
% for a given tau and then compute the tau from the final L and R. use this
% tau as a starting point in classical formulation and then compute
% derivative and value. also save intial and final tau in both the cases

%%
clear all;
clc;
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
%% Solving the factorized formulation in LASSO mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimize   ||A(LR') - b ||_2  subject to ||LR'||_*<=tau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% script1='/Volumes/Users/rkumar/Dropbox/Research/Low-Rank/';
% addpath(genpath([script1 'spgl1Latest']));
% Generate rank=70 matrix 
    sigmas = [1e2*rand(70,1); zeros((100-70), 1)];
    % set up U and V
    U = randn(100, 70);
    [Q R] = qr(U);
    U = Q;
    V = randn(100,70);
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
    % deleted options: 
%    'optTol',1e-8, ...
 %       'bpTol', 1e-8,...
 %       'decTol',1e-6,...
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
   % b = b + .2*randn(size(b));
    params.afunT = @(x)reshape(x,nR,nC);
    params.Ind = find(b==0);
    params.afun = @(x)afun(x,params.Ind);
    rank = [20,30,50,70];
    tau = [1e2 3e2 5e2 7e2 8e2 1e3 1.5e3]; 
    opts.funPenalty = @funLS;
  %  params.nu = 4;
    i = 1;
% for i = 1:length(rank)
    r = rank(i);   
    params.nr=r;
    LInit   = randn(params.numr,params.nr);
    RInit   = randn(params.numc,params.nr);
    xinit   = [vec(LInit);vec(RInit)]; % Initial guess
%     for j = 1:length(tau)
        j = 1;
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
        tau_proxy = 0.5*norm([L1; R1], 'fro')^2;
        
        tau_calcul(i,j) = sum(svd(xs)); 
        SNR_fact(i,j) = -20*log10(norm(D-xs,'fro')/norm(D,'fro'));
        Data_Fact{i,j}= info;
        value(i,j) = norm(resi,2);
        deri(i,j)  = norm(grad,inf);
%     end
%  
%     
% end

% save('Exp1_Factorized.mat','SNR_fact','Data_Fact','value','deri');
