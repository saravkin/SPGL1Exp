


function [f g] = funST(r, params)
nu = params.nu;
f = 0.5*sum((nu).*log(1 + r.*conj(r)./nu));
g = r.*nu./(nu + r.*conj(r));

end
