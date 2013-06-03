%%  Test script to compare the classical versis factorized formulation
% For a given sigma compare the convergance of classical and
% factorized formulation. Compute the best solution from factorized and
% standard and calulating the norm of each and their residual. Also using
% the best approximation from factorized formulation, pluging into the
% classical and observe how fast classical solution converge.
clear all;
clc;
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part1 : Run SPGL1 for convergance using sigma=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r=70;
    sigmas = [1e2*rand(r,1); zeros((100-r), 1)];
    % set up U and V
    U = randn(100, r);
    [Q R] = qr(U);
    U = Q;
    V = randn(100,r);
    [Q R] = qr(V);
    V = Q;
    % Form A
    D = U*diag(sigmas)*V';
rank = [30];
 for i = 1:length(rank)
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Factorized formulation
    %%%%%%%%%%%%%%%%%%%%%%%%
    script1='/Volumes/Users/rkumar/Dropbox/Research/Low-Rank/';
    addpath(genpath([script1 'spgl1Latest']));
    opts = spgSetParms('optTol',1e-4, ...
    'bpTol', 1e-4,...
    'decTol',1e-6,...
    'project', @TraceNorm_project, ...
    'primal_norm', @TraceNorm_primal, ...
    'dual_norm', @TraceNorm_dual, ...
    'proxy', 1, ...
    'verbosity',1,...
    'ignorePErr', 1, ...
    'weights', []);
    r = rank(i); 
    nm = size(D);
    nR = nm(2);
    nC = nm(1);
    [m,n] = size(D);
    idx = randperm(m*n);
    k   = min(m*n,max(1,round(m*n*60/100)));
    idx = sort(idx(1:k));
    op  = opRestriction(m*n,idx);
    Bfun = @(x)op'*op*vec(x);
    b = Bfun(D);
    params.afunT = @(x)reshape(x,nR,nC);
    params.Ind = find(b==0);
    params.afun = @(x)afun(x,params.Ind);
    params.numr = nR;
    params.numc = nC;
    params.funForward = @NLfunForward;
    params.nr=r;
    LInit   = randn(params.numr,params.nr);
    RInit   = randn(params.numc,params.nr);
    xinit   = [vec(LInit);vec(RInit)]; % Initial guess
    sigma = 1e-8;
    tau = norm(xinit,1);
    opts.funPenalty = @funLS;
    % =====================================================================
    [xS,resi,grad,info] = spgl1(@NLfunForward,b,tau,sigma,xinit,opts,params);
    % =====================================================================
    e = params.numr*params.nr;
    L1 = xS(1:e);
    R1 = xS(e+1:end);
    L1 = reshape(L1,params.numr,params.nr);
    R1 = reshape(R1,params.numc,params.nr);
    xs = L1*R1';
    SNR_fact_sigma_zeros(i) = -20*log10(norm(D-xs,'fro')/norm(D,'fro'));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Insert best value from factorized in classical  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear opts b;
    clearvars op
    idx = randperm(m*n);
    k   = min(m*n,max(1,round(m*n*60/100)));
    idx = sort(idx(1:k));
    op  = opRestriction(m*n,idx);
    script1='/Volumes/Users/rkumar/Dropbox/Research/Low-Rank/Matrix_Pareto_Factorization/SLIM.Projects.MatrixParetoFact/scripts/Classical_NuclearNorm/';
    addpath(genpath([script1 'scripts-2']));
    addpath(genpath([script1 'scripts/dependencies']));
    run([script1 'TR201002-Figures/scripts/runme']);
    opts = spgSetParms('optTol',1e-4, ...
        'bpTol', 1e-4,...
        'decTol',1e-6,...
        'verbosity',1);
    opts.project     = @(x,weight,tau) NormNuc_project(m,n,x,tau);
    opts.primal_norm = @(x,weight)     NormNuc_primal(m,n,x);
    opts.dual_norm   = @(x,weight)     NormNuc_dual(m,n,x);
    b=D(:);
    sigma=1e-8;
%     xinit = vec(xs);
    % =====================================================================
    [XM,resi,grad,data] = spgl1(op, b(idx), [],sigma, [], opts);
    % =====================================================================
    XM=reshape(XM,nR,nC);
    SNR_Class_sigma_zeros_init_guess_from_factor(i) = -20*log10(norm(D-XM,'fro')/norm(D,'fro'));
    Data_Fact{i}= data;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Classic formulation
    %%%%%%%%%%%%%%%%%%%%%%%%
    clear opts b;
    clearvars op
    idx = randperm(m*n);
    k   = min(m*n,max(1,round(m*n*60/100)));
    idx = sort(idx(1:k));
    op  = opRestriction(m*n,idx);
    opts = spgSetParms('optTol',1e-4, ...
        'bpTol', 1e-4,...
        'decTol',1e-6,...
        'verbosity',1);
    opts.project     = @(x,weight,tau) NormNuc_project(m,n,x,tau);
    opts.primal_norm = @(x,weight)     NormNuc_primal(m,n,x);
    opts.dual_norm   = @(x,weight)     NormNuc_dual(m,n,x);
    b=D(:);
    sigma=1e-8;
    % =====================================================================
    [XM,resi,grad,data] = spgl1(op, b(idx), [],sigma, [], opts);
    % =====================================================================
    XM=reshape(XM,nR,nC);
    SNR_Class_sigma_zeros(i) = -20*log10(norm(D-XM,'fro')/norm(D,'fro'));
    Norm_fact(i) = norm(xs,'fro');
    Norm_Class(i) = norm(XM,'fro'); 
    clear opts b,
    clearvars op
end
save('Exp2_part1_rank30.mat','Norm_fact','Norm_Class','SNR_Class_sigma_zeros','SNR_Class_sigma_zeros_init_guess_from_factor',...
    'Data_Fact','SNR_fact_sigma_zeros');
