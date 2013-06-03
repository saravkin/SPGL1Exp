function [f g k sig2] = funSelfTuneST(r, params)

para.r = r;
n = length(r);

inVal = 5; % start search at 1.
[nu] = fminsearch(@(x)studentsDFfit(x, para), inVal);
%scale = 1;
%para.print = 0;
%[nu, dfmin, count] = NewtonAlg(@studentsDFfit, inVal, para, scale);

[val, g, h, eta] = studentsDFfit(nu, para);


f = eta*sum(log(1 + r.*conj(r)./nu));
f = f -n*gammaln(eta) + n*gammaln(eta - 0.5) + 0.5*n*log(nu) ;

g = 2*eta*r./(nu + r.*conj(r));

k = 2*eta - 1;
sig2 = nu/k;

end


function [y g h scale] = studentsDFfit(k, para)
% function [y] = studentsDFfit(k, para)

% Input:  k = df* scale^2 parameter
%         para.n  number of measurements
%         para.r  residual vector (length n)


r = para.r;
n = length(r);

r2 = r.*conj(r);                  % r_i^2
lr = log(1 + r2/k);         % log(1 + r_i^2/k)
slr = sum(lr);              % sum(log(1 + r_i^2/k)

% computing s(nu) quantity
r2dkpr2 = r2./(k + r2);     % r_i^2/(k + r_i^2)
sr2dkpr2 = sum(r2dkpr2);    % sum(r_i^2/(k + r_i^2))
s = n/(2*sr2dkpr2);         % n/(2*sum(r_i^2/(k + r_i^2)))

y = 0.5*n*log(k) + s*slr - n*gammaln(s) + n*gammaln(s - 0.5);

if(nargout > 2)
    % computing s'(nu)
    r2dkpr2sq = r2./(k + r2).^2;        % r_i^2/(k + r_i^2)^2
    sr2dkpr2sq = sum(r2dkpr2sq) ;        % sum(r_i^2/(k + r_i^2)^2)
    sprime = (2/n) * s^2 * sr2dkpr2sq;   % (2/n) * s^2 * sum(r_i^2/(k + r_i^2)^2)

    g = sprime*(-n*psi(0,s) + n*psi(0, s-0.5) + slr);
end

if(nargout > 3)
    % computing s''(nu)
    r2dkpr2cube = r2./(k + r2).^3;        % r_i^2/(k + r_i^2)^3
    sr2dkpr2cube = sum(r2dkpr2cube) ;        % sum(r_i^2/(k + r_i^2)^3)
    sdprime = (4/n)*s^2 *((2/n)*s*sr2dkpr2sq^2 -sr2dkpr2cube); % formula for s''

    h =  sdprime*(-n*psi(0,s) + n*psi(0, s-0.5) + slr) +  ...
        sprime^2*(-n*psi(1,s) + n*psi(1, s-0.5)) - n*sprime/(2*k*s);
end

scale = 0.5*n/sr2dkpr2;

end


% function [y scale] = studentsDFfit(k, para)
% % function [y] = studentsDFfit(k, para)
%
% % Input:  k = df* scale^2 parameter
% %         para.n  number of measurements
% %         para.r  residual vector (length n)
%
%
% r = para.r;
% n = length(r);
%
% r2 = r.*conj(r);                  % r_i^2
% lr = log(1 + r2/k);         % log(1 + r_i^2/k)
% slr = sum(lr);              % sum(log(1 + r_i^2/k)
% r2dkpr2 = r2./(k + r2);     % r_i^2/(k + r_i^2)
% sr2dkpr2 = sum(r2dkpr2);    % sum(r_i^2/(k + r_i^2))
%
% y = 0.5*n*log(k) + 0.5*n*slr/sr2dkpr2 - n*gammaln(0.5*n/sr2dkpr2) + n*gammaln(0.5*n/sr2dkpr2 - 0.5);
%
% scale = 0.5*n/sr2dkpr2;
%
% end
