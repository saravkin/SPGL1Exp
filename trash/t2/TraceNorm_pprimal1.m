function p = TraceNorm_pprimal(x, weights,params)
% is this primal norm should be on full vector x or for indivisual block
% will work??

spmd
    X_Local=norm(x);
    c=gcat(X_Local);
    normx_global=norm(c);
end
 
p = 0.5 * normx_global{1}^2; 

end