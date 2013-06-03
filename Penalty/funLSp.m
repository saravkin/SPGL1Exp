function [f g] = funLSp(r, params)

spmd
    X_Local=norm(r,2);
    c=gcat(X_Local);
    normf_glob=norm(c,2);
    g = r./normf_glob;
end
   f=normf_glob{1};
end