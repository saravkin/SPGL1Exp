function [hss] = HSS_Struct(A,kMax)
% Factorize matrix in HSS Structure
[hss] = hss_buildA(A,kMax);
end

function [hss] = hss_buildA(A,kMax,k,m,n)

if nargin<3
    k=0;
    [m,n] = size(A);
    m = [1:m];
    n = [1:n];
else
    k=k+1;
end

if (length(m)>2 && length(n)>2 && k<kMax) % parent node
    hss = struct('ml',floor((m(end)-m(1)+1)/2),'mr',ceil((m(end)-m(1)+1)/2),'nl',floor((n(end)-n(1)+1)/2),'nr',ceil((n(end)-n(1)+1)/2));
    ml = m(1):m(1)+hss.ml-1;
    mr = m(1)+hss.ml:m(end);
    nl = n(1):n(1)+hss.nl-1;
    nr = n(1)+hss.nl:n(end);

    % indices mapping of each HSS partitioning 
    indxml = transp(repmat(ml,length(nr),1));
    indxnr = repmat(nr,length(ml),1);
    indxmr = transp(repmat(mr,length(nl),1));
    indxnl = repmat(nl,length(mr),1);
    hss.Bl=sub2ind(size(A),vec(indxml),vec(indxnr));
    hss.Bu=sub2ind(size(A),vec(indxmr),vec(indxnl));
    %recursion part
    [hss.hssl] = hss_buildA(A,kMax,k,ml,nl);
    [hss.hssr] = hss_buildA(A,kMax,k,mr,nr);
else
    indxm = transp(repmat(m,length(n),1));
    indxn = repmat(n,length(m),1);
    hss.D = sub2ind(size(A),vec(indxm),vec(indxn));
end

end

function K= Block(kMax)

n = 1;
K = 2^(n+1);
for i =2:kMax
    t = n-(i-1);
    if t>0
       K = K+2^(n-(i-1)); 
    else
        K = K
    end
    
end
end



















