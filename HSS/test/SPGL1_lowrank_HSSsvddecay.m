% addpath(genpath('/Dropbox/spgl1Latest'));
% addpath(genpath('/Dropbox/MaxNorm/GulfofSuez350/code/Function'));
clear all;
addpath(genpath('/users/slic/rkumar/spgl1Latest'));
%% fetch the data from cluster
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
nt = 1024;
nr = 355;
ns = 355;
D=ReadSuFast('/scratch/slic/rkumar/suezdatatest/SuezShots125-355shots.su');
D = reshape(D,nt,nr,ns);  
data = fft(D);
clear D
data = permute(reshape(data,nt,nr,ns),[3 2 1]);
freq = 0; % for low frequency else high frequency
a=250;
D = squeeze((data(1:354,1:354,a)));
nm = size(D);
nR = nm(2);
nC = nm(1);% length of frequency axis
%% Define Restriction Operator
inds=randperm(nR);
ind=inds(1:floor(nR/2));
R1 = opRestriction(nR,ind);
R2 = opKron(R1,opDirac(nR));
opSR=opSR2MH(354);
info=opSR([],0);
SR=opFunction(info{1},info{2},opSR);
Bfun = @(x)SR*R2'*R2*vec(x);
b = Bfun(D);
b=reshape(b,nR,2*nC);
Column=[30 40 50 60 80 100 200 300 400 500 600 708];
SVD1=zeros(354,length(Column));
SVD2=zeros(354,length(Column));
SVD3=zeros(354,length(Column));
SVD4=zeros(354,length(Column));
for i=1:length(Column)
    C=Column(i);
    C=floor(C/2);
A=b(:,nC-C+1:nC+C);
B=b(:,1:nC-C-1);
Dtest=SR*vec(D);
Dtest=reshape(Dtest,nR,2*nC);
Ao=Dtest(:,nC-C+1:nC+C);
A1=Dtest(:,1:nC-C-1);
s1=svd(A);
s2=svd(Ao);
s3=svd(B);
s4=svd(A1);
SVD1(1:length(s1),i)=s1;
SVD2(1:length(s2),i)=s2;
SVD3(1:length(s3),i)=s3;
SVD4(1:length(s4),i)=s4;
end

x=1:354;
figure,
semilogx(x,(SVD2(:,1)./max(SVD2(:,1))),'r',x,(SVD2(:,2)./max(SVD2(:,2))),'b',x,(SVD2(:,3)./max(SVD2(:,3))),'g',...
    x,(SVD2(:,4)./max(SVD2(:,4))),'k',x,(SVD2(:,5)./max(SVD2(:,5))),'m',x,(SVD2(:,6)./max(SVD2(:,6))),'*m',...
    x,(SVD2(:,7)./max(SVD2(:,7))),'*r',...
    x,(SVD2(:,8)./max(SVD2(:,8))),'*b',x,(SVD2(:,9)./max(SVD2(:,9))),'*g',...
    x,(SVD2(:,10)./max(SVD2(:,10))),'*b',x,(SVD2(:,11)./max(SVD2(:,11))),'*k');
title('SVD Decay of true data in mh domain');

figure,
semilogx(x,(SVD1(:,1)./max(SVD1(:,1))),'r',x,(SVD1(:,2)./max(SVD1(:,2))),'b',x,(SVD1(:,3)./max(SVD1(:,3))),'g',...
    x,(SVD1(:,4)./max(SVD1(:,4))),'k',x,(SVD1(:,5)./max(SVD1(:,5))),'m',x,(SVD1(:,6)./max(SVD1(:,6))),'*m',...
    x,(SVD1(:,7)./max(SVD1(:,7))),'*r',...
    x,(SVD1(:,8)./max(SVD1(:,8))),'*b',x,(SVD1(:,9)./max(SVD1(:,9))),'*g',...
    x,(SVD1(:,10)./max(SVD1(:,10))),'*b',x,(SVD1(:,11)./max(SVD1(:,11))),'*k');
title('SVD Decay of true data in mh domain with missing entries');

figure,
semilogx(x,(SVD3(:,1)./max(SVD3(:,1))),'r',x,(SVD3(:,2)./max(SVD3(:,2))),'b',x,(SVD3(:,3)./max(SVD3(:,3))),'g',...
    x,(SVD3(:,4)./max(SVD3(:,4))),'k',x,(SVD3(:,5)./max(SVD3(:,5))),'m',x,(SVD3(:,6)./max(SVD3(:,6))),'*m',...
    x,(SVD3(:,7)./max(SVD3(:,7))),'*r');
title('SVD Decay of true data in mh domain with missing entries');

figure,
semilogx(x,(SVD4(:,1)./max(SVD4(:,1))),'r',x,(SVD4(:,2)./max(SVD4(:,2))),'b',x,(SVD4(:,3)./max(SVD4(:,3))),'g',...
    x,(SVD4(:,4)./max(SVD4(:,4))),'k',x,(SVD4(:,5)./max(SVD4(:,5))),'m',x,(SVD4(:,6)./max(SVD4(:,6))),'*m',...
    x,(SVD4(:,7)./max(SVD4(:,7))),'*r');
title('SVD Decay of true data in mh domain');

figure,
semilogx(x,(SVD4(:,1)./max(SVD4(:,1))),'r',x,(SVD3(:,1)./max(SVD3(:,1))),'b',...
    x,(SVD1(:,12)./max(SVD1(:,12))),'g',x,(SVD2(:,1)./max(SVD2(:,1))),'k');


figure,
semilogx(x,(SVD1(:,5)),'r',x,(SVD2(:,5)),'b',...
    x,(SVD3(:,5)),'g',x,(SVD4(:,5)),'k');


