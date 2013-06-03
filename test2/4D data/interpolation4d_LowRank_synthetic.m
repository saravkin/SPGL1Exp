addpath(genpath('/users/slic/rkumar/spgl1Latest'));
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));

% dataDir = '/scratch/slic/shared/forRajiv/BG_raw_data/';
homeDir = '/scratch/slic/rkumar/BG5D/Interpolation/';

%% Sampled indices
% load([homeDir 'INDEX_YSorted.mat']);
% X=18875; % First coordinate position of source grid
% Y=29100; % First coordinate position of receiver grid
% 
% sampleX = 25; sampleY = 25;
% n=401^2;  % 401 is # of acquisition grid points in x and y direction
% nrecs=401;
% 
% INDX=(INDX-X)/sampleX+1; % Indicies corresponding to Source location X
% INDY=(INDY-Y)/sampleY+1; % Indicies corresponding to Source location Y
% 
% 
% I = [vec(INDX)'; vec(INDY)'];
% 

%% Data extraction
% files = dir([dataDir 'S*fftrealfreqindex75.bin']);
load /scratch/slic/shared/forRajiv/slice3d.mat;

% distributedMode = true;
% 
% recSubsamplingFactor = 10;
% 
% b = zeros(nrecs*nrecs,size(I,2));
% for k=1:size(I,2)
%     fid = fopen([dataDir files(k).name],'r','a');
%     b(:,k) = fread(fid,nrecs*nrecs,'single');
%     fclose(fid);
% end
% 
% b = b';
% 
% new_nrecs = ceil(nrecs/recSubsamplingFactor);
% c = zeros(size(I,2),new_nrecs*new_nrecs);
% for i=1:size(I,2)
%     t = reshape(b(i,:),nrecs,nrecs);
%     c(i,:) = vec(t(1:recSubsamplingFactor:end,1:recSubsamplingFactor:end));
% end
% 
% b = vec(c);
% nrecs = new_nrecs;
% b = b/norm(b);

%b is a vec of (src x, src y, rec x, recy)

%% Binning data
% I = round(I/10)+1;
  nsrcs = 41; 
  nrecs = 41; 
% idx_1d = sub2ind([nsrcs,nsrcs],I(1,:),I(2,:));
%% Set up operators
inds=randperm(nsrcs*nsrcs);
idx_1d=inds(1:floor(nsrcs*nsrcs/2));
S = opRestriction(nsrcs*nsrcs,idx_1d);
R = opKron(opDirac(nrecs*nrecs),S);

dims = [nsrcs, nrecs, nsrcs, nrecs];
P = opPermute(dims, [1 3 2 4]);

% A = R * P;
A = P*R'*R ;
rhs = A * vec(fslice);

nR = nrecs;
nC = nsrcs;% length of frequency axis

%% Function Handle
params.afunT = @(x)reshape(x,nR*nC,nC*nR);
params.Ind = find(rhs==0);
params.afun = @(x)afun(x,params.Ind);
params.numr = nR*nC;
params.numc = nR*nC;
params.nr = 5; 
params.funForward = @NLfunForward;

%% options for trace norm
opts = spgSetParms('optTol',1e-7, ...
                   'bpTol', 1e-7,...
                    'decTol',1e-4,...
                   'project', @TraceNorm_project, ...
                   'primal_norm', @TraceNorm_primal, ...
                   'dual_norm', @TraceNorm_dual, ...
                   'proxy', 1, ...
                   'ignorePErr', 1, ...
                   'iterations', 100,...
                   'weights', []);
LInit   = 1e-4*randn(params.numr,params.nr)+1i*1e-4*randn(params.numr,params.nr);
RInit   = 1e-4*randn(params.numc,params.nr)+1i*1e-4*randn(params.numr,params.nr);
xinit  = [vec(LInit);vec(RInit)]; % Initial guess
tau = norm(xinit,1);
sigma = 1e-4;
opts.funPenalty = @funLS;
tic;
[xLS,r,g,info] = spgl1(@NLfunForward,rhs,tau,sigma,xinit,opts,params);
toc;
e = params.numr*params.nr;
L1 = xLS(1:e);
R1 = xLS(e+1:end);
L1 = reshape(L1,params.numr,params.nr);
R1 = reshape(R1,params.numc,params.nr);
xls = vec(L1*R1');
Xrec=P'*xls;
Xrec=reshape(Xrec,41,41,41,41);


%% plot the data
Sx=10;
Indy=1:5:41;
for i=1:length(Indy)
  I=Indy(i);  
  Drec=squeeze(Xrec(Sx,I,:,:));
  Dorig=squeeze(fslice(Sx,I,:,:));
  SNR = -20*log10(norm(Dorig-Drec,'fro')/norm(Dorig,'fro'));
figure,
subplot(1,2,1);imagesc(squeeze(abs(Xrec(Sx,I,:,:))));
title(['After Interpolation for Sx= ' num2str(Sx) 'Sy=' num2str(I) 'SNR=' num2str(SNR)]);
subplot(1,2,2);imagesc(squeeze(abs(fslice(Sx,I,:,:))));
title(['Original data for Sx =' num2str(Sx) 'Sy=' num2str(I)]);

end

% save('Synthetic4D.mat');
