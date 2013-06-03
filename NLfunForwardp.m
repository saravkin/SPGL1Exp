function [f1 f2] = NLfunForwardp(x,g,params,Ind)

spmd
    xlocal=getLocalPart(x);
    e = params.numr*params.nr;
    L = xlocal(1:e);
    R = xlocal(e+1:end);
    L = reshape(L,params.numr,params.nr);
    R = reshape(R,params.numc,params.nr);
    t=L*R';
    if isempty(g)
        t=vec(t);
        f1 = params.afun(t,Ind);
        f2 = 0;
    else
        glocal=getLocalPart(g);
        glocal=reshape(glocal,params.numr,params.numc);
        fp = params.afunT(g);
        f1 = [vec(glocal*R); vec(glocal'*L)];
        f2 = vec(fp);
    end
end
end