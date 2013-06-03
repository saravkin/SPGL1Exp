function p = TraceNorm_primalweight(x,weights,params)

% p = 0.5 * norm(x.*weights)^2; 

% E = reshape(x,params.numr,params.numc);
e = params.numr*params.nr;
L = x(1:e); %L = x(1:e,:);
R = x(e+1:end); %R = x(e+1:end,:);
L = reshape(L,params.numr,params.nr);
R = reshape(R,params.numc,params.nr);

% dual of trace norm is operator norm i.e maximum singular value
w1=weights(1);w2=weights(2);
w3=weights(3);w4=weights(4);
w5=weights(5);w6=weights(6);
w7=weights(7);w8=weights(8);
W=weights(9:end);
QL=W(1:w1*w2);
RL=W(w1*w2+1:w1*w2+w3*w4);
QR=W(w1*w2+w3*w4+1:w1*w2+w3*w4+w5*w6);
RR=W(w1*w2+w3*w4+w5*w6+1:end);
QL=reshape(QL,w1,w2);
RL=reshape(RL,w3,w4);
QR=reshape(QR,w5,w6);
RR=reshape(RR,w7,w8);
L1 = QL*L*RL;
R1 = QR*R*RR;

E1 = [vec(L1);vec(R1)];
p = 0.5 * norm(E1)^2; 
end