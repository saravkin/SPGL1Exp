function K= Block(kMax)

K = 2^(kMax+1);
for i =1:kMax-1
    n=i;
    t = kMax-i;
    if t>0
       K = K+2^(kMax-i); 
    end    
end
end


