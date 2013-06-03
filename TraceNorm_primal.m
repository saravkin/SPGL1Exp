function p = TraceNorm_primal(x, weights,params)

p = 0.5 * norm(x.*weights)^2; 

end