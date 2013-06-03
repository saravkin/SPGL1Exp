% x=[1 2 3 4 5];
x=rand(100,1);
B=1e-2;
a=B*x/(0.5*norm(x)^2);
0.5*norm(a)^2
%%
b=0.5*norm(x)^2;
scale=sqrt(B/b);
a=B*x*scale/norm(x);
0.5*norm(a)^2
%%
L=rand(10,5);
R=rand(10,5);
B=1e-6
normLR = norm(L(:)) + norm(R(:))
if(normLR > B)
   LOut = B*L/normLR;
   ROut = B*R/normLR;
else
   LOut = L;
   ROut = R;
end

X = [vec(LOut);vec(ROut)];
0.5*norm(X)^2
X = [vec(L);vec(R)];
%%
L=rand(10,1);
R=rand(10,1);
%% Projection in case of nuclear norm
L=LInit;
R=RInit;
B=tau;
numr = size(L,1);
numc = size(R,1);

row_normsL = sqrt(sum(L.^2,2));
row_normsR = sqrt(sum(R.^2,2));

rescaleL = min(ones(numr,1), sqrt(B)./row_normsL);
rescaleR = min(ones(numc,1), sqrt(B)./row_normsR);

LOut = (L.*(rescaleL*ones(1, size(L,2))));
ROut = (R.*(rescaleR*ones(1, size(R,2))));

X = [vec(LOut);vec(ROut)];
tfac=0.5*norm(X)^2;

Lfac=B/tfac;
LOut = (L.*(rescaleL*ones(1, size(L,2))))*sqrt(Lfac);
ROut = (R.*(rescaleR*ones(1, size(R,2))))*sqrt(Lfac);
X = [vec(LOut);vec(ROut)];
tfac=0.5*norm(X)^2

%%
L=rand(100,5);
R=rand(100,5);
x=[L(:);R(:)];
B=1e6
c=sqrt(B/(0.5*norm(x)^2));
L1=min(1,c)*L(:);
R1=min(1,c)*R(:);
X=[L1(:);R1(:)];
0.5*norm(X)^2
sqrt(c)
