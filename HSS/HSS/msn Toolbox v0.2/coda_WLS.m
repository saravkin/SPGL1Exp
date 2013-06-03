function [q r v A] = coda_WLS(D, A, eta)
%Computes the complete orthogonal decomposition corresponding to the
%weighted least squares problem min_{x}||D(Ax-b)||. This function factors
%the matrix DA=QRV, in a stable manner. We iteratively refine to get very
%high accuracy(tol) for the QR.
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
[Av v] = refine_v(D*A, eta);
[q r] = qr(Av, 0);
end

function [Ao vo] = refine_v(A, eta)
[~,s,v] = svd(A);
Ao = [];Av = A*v;vo = [];
if 0
    Ao = Av;
    vo = v;
    return;
end

while(1)
    [A1 A2 v1 v2 s term] = compute_refine(Av, v, eta, s);

    Ao = [Ao A1];
    vo = [vo v1];
    Av = A2;
    v = v2;
    if term == 1
        break;
    end    
end
Ao = [Ao A2];
vo = [vo v2];
end

function [A1 A2 v1 v2 s term] = compute_refine(Av, v, eta, s)
    tol = s(1,1)/eta;
    k = find(diag(s) < tol);
    if isempty(k)
        A1 = Av;
        A2 = [];
        v1 = v;
        v2 = [];
        term = 1;
        return;
    end
    k = k(1);
    A1 = Av(:,1:k-1);
    A2 = Av(:,k:end);
    v1 = v(:,1:k-1);
    v2 = v(:,k:end);
    [~,s,v] = svd(A2);
    A2 = A2*v;
    v2 = v2*v;
    if isempty(s)
        term = 1;
    else
        term = 0;
    end
end
