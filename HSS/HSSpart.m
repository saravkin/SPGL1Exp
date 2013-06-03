function [h] = HSSpart(A)

% This function will partitioning the matrix in 4 block 
[m,n] = size(A);
if floor(m/2)==0
m = [1:m];
n = [1:n]; 
    % Partition of a matrix in SR domain
%
%     nl    nr
% ----------------
%    |    |     |
% ml |    |     |
%    |    |     |
% ----------------
%    |    |     |
% mr |    |     |
%    |    |     |
% ----------------
%
hss = struct('ml',floor(m(end)/2),'mr',ceil(m(end)/2),'nl',floor(n(end)/2),'nr',ceil(n(end)/2));

ml = m(1):hss.ml;
mr = m(1)+hss.ml:m(end);
nl = n(1):hss.nl;
nr = n(1)+hss.nl:n(end);
h.ml = ml;
h.nl = nl;
h.mr = mr;
h.nr = nr;
% seperate out the diagonal and non-diagonal part fo matrix at various level
h.Ul = A(ml,nl);
h.Ll = A(mr,nl);
h.Ur = A(ml,nr);
h.Lr = A(mr,nr);
Ll = mcnnorm (h.Ll,opts,params,SR,R);
Ur = mcnnorm (h.Ur,opts);
else
    h = A;
end
end

