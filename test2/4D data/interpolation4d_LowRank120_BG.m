% clear all;
addpath(genpath('/users/slic/rkumar/spgl1Latest'));
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));

dataDir = '/scratch/slic/shared/forRajiv/BG_raw_data/';
homeDir = '/scratch/slic/rkumar/BG5D/Interpolation/';

%% Sampled indices
load([homeDir 'INDEX_YSorted.mat']);
X=18875; % First coordinate position of source grid
Y=29100; % First coordinate position of receiver grid

sampleX = 25; sampleY = 25;
n=401^2;  % 401 is # of acquisition grid points in x and y direction
nrecs=401;

INDX=(INDX-X)/sampleX+1; % Indicies corresponding to Source location X
INDY=(INDY-Y)/sampleY+1; % Indicies corresponding to Source location Y
% Indx=INDX(1:5:end);
% Indy=INDY(1:5:end);
I = [vec(INDX)'; vec(INDY)'];

%% Data extraction
files = dir([dataDir 'S*fftrealfreqindex75.bin']);
load([homeDir 'BG5D_Receiver_Header.mat']);

distributedMode = true;

recSubsamplingFactor = 4;

b = zeros(nrecs*nrecs,size(I,2));
for k=1:size(I,2)
    fid = fopen([dataDir files(k).name],'r','a');
    b(:,k) = fread(fid,nrecs*nrecs,'single');
    fclose(fid);
end
b = b';

new_nrecs = ceil(nrecs/recSubsamplingFactor);
c = zeros(size(I,2),new_nrecs*new_nrecs);
for i=1:size(I,2)
    t = reshape(b(i,:),nrecs,nrecs);
    c(i,:) = vec(t(1:recSubsamplingFactor:end,1:recSubsamplingFactor:end));
end

b = vec(c);
nrecs = new_nrecs;
% b = b/norm(b);
%b is a vec of (src x, src y, rec x, recy)
%% Binning data
I = round(I/4)+1;
nsrcs = 100; 
idx_1d = sub2ind([nsrcs,nsrcs],I(1,:),I(2,:));
%% Set up operators
S = opRestriction(nsrcs*nsrcs,idx_1d);
R = opKron(opDirac(nrecs*nrecs),S);
dims = [nsrcs, nrecs, nsrcs, nrecs];
P = opPermute(dims, [1 3 2 4]);
A = R * P;
rhs = A' * b;
nR = nrecs;
nC = nsrcs;% length of frequency axis
%% Function Handle
params.afunT = @(x)reshape(x,nR*nC,nC*nR);
params.Ind = find(rhs==0);
params.afun = @(x)afun(x,params.Ind);
params.numr = nR*nC;
params.numc = nR*nC;
params.nr = 150; 
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
                   'iterations', 300,...
                   'weights', []);
LInit   = 1e3*randn(params.numr,params.nr);
RInit   = 1e3*randn(params.numc,params.nr);
xinit  = [vec(LInit);vec(RInit)]; % Initial guess
tau = norm(xinit,1);
sigma = 1e-4;
opts.funPenalty = @funLS;
tic;
[xLS,r,g,info] = spgl1(@NLfunForward,rhs,tau,sigma,xinit,opts,params);
toc;
% save('/scratch/slic/rkumar/BG5D/Interpolation/Interp_BG5D_rank120.mat');
%% Plotting the results
% load Interp_BG5D_rank80.mat
% e = params.numr*params.nr;
% L1 = xLS(1:e);
% R1 = xLS(e+1:end);
% L1 = reshape(L1,params.numr,params.nr);
% R1 = reshape(R1,params.numc,params.nr);
% xls = vec(L1*R1');
% % Xrec=P*xls;
% % interp_data=reshape(Xrec,nsrcs,nsrcs,nrecs,nrecs);
% % interp_data=vec(interp_data);
% missingData = reshape(rhs,nsrcs,nrecs,nsrcs,nrecs);
% knownSrcs = [I(1,:);I(2,:)];     % known sources are (1,2) and (3,4)
% samplePoints = [31,64,63,22,30,73,64,13,81,39,39,64,64;36,25,66,18,40,18,51,3,81,17,18,50,51]; % want to display sources (1,1) and (1,2)
% upsampleRecs = 401;            % fourier upsampled grid to display figures on (set to 0 for no upsampling)
% isPermuted = true;           % true if data is in [src x rec x src y rec y] format, false if data is in [src x src y rec x rec y] format
% saveDir = '/scratch/slic/rkumar/BG5D/Interpolation/rank80iter300/';                % directory to save figures in (empty if no figures are to be saved)
% dispTitles = false;
% axisLabels = true;
% 
% interpolation4dFigures(missingData,xls,knownSrcs, ...
%                        'saveDir',saveDir,...
%                        'isPermuted',isPermuted,...
%                        'dispTitles',dispTitles, ...
%                        'axisLabels',axisLabels, ...
%                        'samplePoints',samplePoints,...
%                        'upsampleRecs',upsampleRecs);
