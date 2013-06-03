function B = sr2mht(A,mode)
% simple source-receiver to midpoint-offset conversion. 
%
% use:
%   sr2mht(A,mode)
%
% input:
%   A    - source-receiver slice of size n x n or midpoint-offset slice size 2*n-1 x 2*n-1
%   mode - 1: forward, -1: backward
%
% output:
%   B - midpoint-offset slice or source-receiver slice
%
% Tristan van Leeuwen, 2012.
 
switch mode
    case 1 % forward
        % initialize output matrix
        n  = size(A,1);
        nt = 2*n-1;
        B  = zeros(nt);
 
        % fill output matrix
        for k = 1-n:n-1
            tmp = [zeros(floor(abs(k)/2),1);diag(A,k);zeros(floor(abs(k)/2),1)];
            B((mod(k,2)+1):2:nt,k+n) = tmp(1:n-mod(k,2));
        end        
    case -1 % backward
        % initialize output matrix
        nt = size(A,1);
        n  = (nt+1)/2;
        B  = zeros(n);
        
        % fill output matrix
        for k = 1:n
            tmp    = diag(A,n-1+2*(1-k));
            B(k,:) = tmp((n-1)/2-abs(k-1-(n-1)/2) + [1:n]);
        end
end
