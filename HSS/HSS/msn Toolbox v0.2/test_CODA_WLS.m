function test_CODA_WLS
%This script tests the coda_WLS module, that performs a complete orthogonal
%decomposition of the given ill-conditioned system. We shall use this to
%solve the MSN system and see it we get more accurate solutions.
% test1DMSN;

% Authors:  Karthik Jayaraman Raghuram, jrk@ece.ucsb.edu
% v1.0 Created 13-Oct-10
%
% Copyright (C) 2009  Karthik Jayaraman Raghuram (jrk@ece.ucsb.edu)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
testSystem1D;
end

function testSystem2D
s = [10 20 50];
N = 11;
Niter = length(s);
err_ref_ori = zeros(Niter,1);
err_new_ori = zeros(Niter,1);
err_ref_asc = zeros(Niter,1);
err_new_asc = zeros(Niter,1);
err_ref_desc = zeros(Niter,1);
err_new_desc = zeros(Niter,1);
err_ref_Ming = zeros(Niter,1);
err_new_Ming = zeros(Niter,1);
x = linspace(-1,1,N);
[X Y] = meshgrid(x);
x = reshape(X, numel(X), 1);
y = reshape(Y, numel(Y), 1);

Mx = 6*N;
My = 6*N;
V = zeros(N^2, Mx*My);
D = zeros(Mx*My,1);
p=1;
for m=0:Mx-1,
    for n=0:My-1,            
        V(:,p) = cos(acos(x)*m).*cos(acos(y)*n);
        D(p,1) = (1+m^2+n^2)^(-1);
        p = p+1;
    end
end
a0 = rand(N^2, 10);

%original ordering
P = eye(Mx*My);

for iter=1:Niter,
    iter
    Di = diag(D.^(s(iter)/2));
    f = P*Di*V'*a0;
    %Normal MSN Solution
    [q r] = qr(P*Di*V',0);
    %\[ \Rightarrow VDi = r^{T}q^{T} \]
    %\[ \Rightarrow a = q*(r^{t}\f) \]
    a = r\(q'*f);
    
    err_ref_ori(iter, 1) = max(abs(a-a0))/max(abs(a0));

    %New MSN Solution
    [q r v Av] = coda_WLS(P*Di, V', 2);
    anew = v*(r\(q'*f));

    err_new_ori(iter, 1) = max(max(abs(anew-a0)))/max(max(abs(a0)));    
end

%Sort the D (apply some permutation for numerical advantages)
[~,sp] = sort(D,'ascend');
Ds = D;
P = eye(Mx*My);
P = P(sp,:);
% V = V*P;
for iter=1:Niter,
    iter
    Di = diag(D.^(s(iter)/2));
    f = P*Di*V'*a0;
    %Normal MSN Solution
    [q r] = qr(P*Di*V',0);
    %\[ \Rightarrow VDi = r^{T}q^{T} \]
    %\[ \Rightarrow a = q*(r^{t}\f) \]
    a = r\(q'*f);
    
    err_ref_asc(iter, 1) = max(abs(a-a0))/max(abs(a0));

    %New MSN Solution
    [q r v Av] = coda_WLS(P*Di, V', 2);
    anew = v*(r\(q'*f));

    err_new_asc(iter, 1) = max(max(abs(anew-a0)))/max(max(abs(a0)));    
end


%Sort the D (apply some permutation for numerical advantages)
[~,sp] = sort(D,'descend');
Ds = D;
P = eye(Mx*My);
P = P(sp,:);
% V = V*P;
for iter=1:Niter,
    iter
    Di = diag(D.^(s(iter)/2));
    f = P*Di*V'*a0;
    %Normal MSN Solution
    [q r] = qr(P*Di*V',0);
    %\[ \Rightarrow VDi = r^{T}q^{T} \]
    %\[ \Rightarrow a = q*(r^{t}\f) \]
    a = r\(q'*f);
    
    err_ref_desc(iter, 1) = max(abs(a-a0))/max(abs(a0));

    %New MSN Solution
    [q r v Av] = coda_WLS(P*Di, V', 2);
    anew = v*(r\(q'*f));

    err_new_desc(iter, 1) = max(max(abs(anew-a0)))/max(max(abs(a0)));    
end

P = eye(Mx*My);
%Prof.Ming's code
for iter=1:Niter,
    iter
    Di = diag(D.^(s(iter)/2));
    f = P*Di*V'*a0;
    %Normal MSN Solution
    [q r] = qr(P*Di*V',0);
    %\[ \Rightarrow VDi = r^{T}q^{T} \]
    %\[ \Rightarrow a = q*(r^{t}\f) \]
    a = r\(q'*f);
    
    err_ref_Ming(iter, 1) = max(abs(a-a0))/max(abs(a0));

    %New MSN Solution
    [aMing out] = WLS2(P*Di*V', f);
    err_new_Ming(iter, 1) = max(max(abs(aMing-a0)))/max(max(abs(a0)));    
end

display('2D MSN LS System results');
display('Original_QR   CODA');
format short e;
[err_ref_ori err_new_ori err_ref_asc err_new_asc err_ref_desc err_new_desc err_ref_Ming err_new_Ming]
format compact;
figure, semilogy(s, [[err_ref_ori err_new_ori err_ref_asc err_new_asc err_ref_desc err_new_desc err_ref_Ming err_new_Ming]])
title('Comparison of Max relative error in solving 2D MSN LS System with 10 random RHS');
xlabel('s');ylabel('Max Relative Error');
legend('QR weights', 'CODA WEights', 'QR Ascending Weights', 'CODA Ascending Weights', 'QR Descending Wts', 'CODA Descending Wts', 'QR original', 'QR Ming');
end

function testSystem1D
s = [1:20];
N = 100;
Niter = length(s);
err_ref = zeros(Niter,1);
err_new = zeros(Niter,1);
x = linspace(-1,1,N)';

M = 2*N;
V = cos(acos(x)*[0:M-1]);
D = [1:M];
a0 = rand(N, 10);
for iter=1:Niter,
    iter
    Di = diag(D.^(-s(iter)));
    f = Di*V'*a0;
    %Normal MSN Solution
    [q r] = qr(Di*V',0);
    %\[ \Rightarrow VDi = r^{T}q^{T} \]
    %\[ \Rightarrow a = q*(r^{t}\f) \]
    a = r\(q'*f);
    
    err_ref(iter, 1) = max(abs(a-a0))/max(abs(a0));

    %New MSN Solution
    [q r v Av] = coda_WLS(Di, V', 2);
    anew = v*(r\(q'*f));

    err_new(iter, 1) = max(max(abs(anew-a0)))/max(max(abs(a0)));    
    
    aMing = WLS2(Di*V',f);
    errMing(iter, 1) = max(max(abs(aMing-a0)))/max(max(abs(a0)));    
    
end
display('1D MSN LS System results');
display('Original_QR   CODA');
format short e;
[err_ref err_new errMing]
format compact;
figure, semilogy(s, [err_ref err_new errMing])
end
