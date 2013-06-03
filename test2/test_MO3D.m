addpath(genpath('/Volumes/Users/rkumar/Dropbox/Research/Low-Rank/spgl1Latest'));

D=rand(4,4,4,4);
sx=0:25:25*3;
gx=sx;

B  = MO3D(D,sx,gx,25,25,1);


C  = MO3D(B,sx,gx,25,25,-1);
%%
