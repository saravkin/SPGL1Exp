% addpath(genpath('/Dropbox/spgl1Latest'));
% addpath(genpath('/Dropbox/MaxNorm/GulfofSuez350/code/Function'));
clear all;
addpath(genpath('/users/slic/rkumar/spgl1Latest'));
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
nt = 1024;
nr = 355;
ns = 355;
D=ReadSuFast('/scratch/slic/rkumar/suezdatatest/SuezShots125-355shots.su');
D = reshape(D,nt,nr,ns);  
data = fft(D);
clear D
data = permute(reshape(data,nt,nr,ns),[3 2 1]);
a=250;
D = squeeze((data(1:352,1:352,a)));
nm = size(D);
nR = nm(2);
nC = nm(1);% length of frequency axis
inds=randperm(nR);
ind=inds(1:floor(nR/2));
R1 = opRestriction(nR,ind);
R2 = opKron(R1,opDirac(nR));
Bfun = @(x)R2'*R2*vec(x);
b = Bfun(D);
b=reshape(b,nR,nC);

% [hss] = LowRank_HSS_test(b,2);


% options for trace norm
opts = spgSetParms('optTol',1e-5, ...
                   'bpTol', 1e-5,...
                    'decTol',1e-4,...
                   'project', @TraceNorm_project, ...
                   'primal_norm', @TraceNorm_primal, ...
                   'dual_norm', @TraceNorm_dual, ...
                   'proxy', 1, ...
                   'ignorePErr', 1, ...
                   'iterations', 150,...
                 'weights', []);
level=[1 2 3 4];
r=[5:5:20];
for j=1:length(level)
    kMax=level(j);
for i=1:length(r)
    rank=r(i);
sigmafact=1e-3;
initfact=1;
Z = zeros(size(b,1),size(b,2));
tic;
[hss,c] = LowRank_HSS(b,Z,kMax,opts,rank,sigmafact,initfact);
toc;
% XR = [hss.hssl.D hss.Bl;hss.Bu hss.hssr.D];
SNR(j,i) = -20*log10(norm(D-c,'fro')/norm(D,'fro'));
end
end

figure,
plot(r,SNR(1,:),'g',r,SNR(2,:),'b',r,SNR(3,:),'r',r,SNR(4,:),'k');


