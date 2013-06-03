function [hss] = LowRank_HSS_test(A,kMax)
% Factorize matrix in HSS Structure
[hss] = hss_buildA_iter(A,kMax);
end

function [hss] = hss_buildA_iter(A,kMax,k,m,n)

if nargin<3
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
    [hss.hssl] = hss_buildA_iter(A,kMax,k,ml,nl);
    [hss.hssr] = hss_buildA_iter(A,kMax,k,mr,nr);
else
    [hss.D] = A(m,n);
end

end
