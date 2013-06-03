function [hss,C] = HSSlr(A,Z,kMax)
% Factorize matrix in HSS Structure
[hss,C] = hss_buildA_iter(A,Z,kMax);
end

function [hss,C] = hss_buildA_iter(A,C,kMax,k,m,n)

if nargin<4
    k=1;
    [m,n] = size(A);
    m = [1:m];
    n = [1:n];
else
    k=k+1;
end

if (k<kMax) % parent node
    hss = struct('ml',floor(length(m)/2),'mr',ceil(length(m)/2),'nl',floor(length(m)/2),'nr',ceil(length(m)/2));

    % indices of the current sub-block      
    %      nl nr
    %   ml
    %   mr
    ml = m(1):m(1)+hss.ml-1;
    mr = m(1)+hss.ml:m(end);
    nl = n(1):n(1)+hss.nl-1;
    nr = n(1)+hss.nl:n(end);

    % indices with out the indices of current sub-block
    B2=A(ml,nr);
    B3=A(mr,nl);
    hss.Bl=B2;
    hss.Bu=B3;
    C(ml,nr)=B2;
    C(mr,nl)=B3;
    %recursion part
    [hss.hssr,C] = hss_buildA_iter(A,C,kMax,k,mr,nr);
    [hss.hssl,C] = hss_buildA_iter(A,C,kMax,k,ml,nl);

else
    hss.D = A(m,n);
    C(m,n)=hss.D;
end

end
