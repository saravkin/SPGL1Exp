function [hss,C] = LowRank_HSS(A,Z,kMax,opts,rank,sigmafact,initfact)
% Factorize matrix in HSS Structure
[hss,C] = hss_buildA_iter(A,Z,kMax,opts,rank,sigmafact,initfact);
end

function [hss,C] = hss_buildA_iter(A,C,kMax,opts,rank,sigmafact,initfact,k,m,n)

if nargin<8
    k=0;
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
    [B] = lowrankinterp(B2,opts,rank,sigmafact,initfact);
    [M] = lowrankinterp(B3,opts,rank,sigmafact,initfact);
    hss.Bl=B;
    hss.Bu=M;
    C(ml,nr)=B;
    C(mr,nl)=M;
    %recursion part
    [hss.hssl,C] = hss_buildA_iter(A,C,kMax,opts,rank,sigmafact,initfact,k,ml,nl);
    [hss.hssr,C] = hss_buildA_iter(A,C,kMax,opts,rank,sigmafact,initfact,k,mr,nr);
else
    D1 = A(m,n);
    [hss.D] = lowrankinterp(D1,opts,rank,sigmafact,initfact);
    C(m,n) = hss.D;
end

end
