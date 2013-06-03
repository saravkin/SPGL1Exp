function [y] = afun(x,Ind)
 x = vec(x);
 x(Ind)=0;
 y = x;
 end