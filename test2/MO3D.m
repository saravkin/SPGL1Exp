function [B,XM,XO] = MO3D(D,SX,GX,sampx,sampy,mode )

switch mode
    case 1 % forward
        n  = size(D,1);
        nt = 2*n-1;
        B  = zeros(nt,nt,nt,nt);
        XM=min(SX):(SX(2)-SX(1))/2:max(SX);
        XO=-max(SX):(SX(2)-SX(1)):max(SX);
        for i=1:length(SX)
            sx=SX(i);
            for j=1:length(GX);
                gx=GX(j);
                %   midpoint x  and offset x value
                xm = (sx + gx)/2;
                xo= (sx - gx);
                Ax=(sx/sampx)+1;
                Ay=(gx/sampy)+1;
                C=squeeze(D(Ax,:,Ay,:));
                Cmap=sr2mht(C,mode);
                % select subset of data and organize in data matrix
                Is   = find(XM==xm);
                Ij   = find(XO==xo);
                B(Is,:,Ij,:)=Cmap;
           end
        end
        
        
      
        
    case -1 % backward
        % initialize output matrix
        nt = size(A,1);
        n  = (nt+1)/2;
        B  = zeros(n,n,n,n);
        A=reshape(A,nt*nt,nt*nt);
        % fill output matrix
        for k = 1:n*n
            tmp    = diag(A,2*n+2*(1-k));
            B(k,:) = tmp((n-1)/2-abs(k-1-(n-1)/2) + [1:n]);
        end
%         nt = size(D,1);
%         n  = (nt+1)/2;
%         B  = zeros(n,n,n,n);
%         XM=min(SX):(SX(2)-SX(1))/2:max(SX);
%         XO=-max(SX):(SX(2)-SX(1)):max(SX);
%         for i=1:length(XM)
%             xm=XM(i);
%             for j=1:length(XO);
%                 xo=XO(j);
%                 %   midpoint x  and offset x value
%                 sx = abs((2*xm + xo));
%                 gx= abs((2*xm - xo));
%                 Ax=find(xm==XM);
%                 Ay=find(xo==XO);
%                 C=squeeze(D(Ax,:,Ay,:));
%                 Cmap=sr2mht(C,mode);
%                 % select subset of data and organize in data matrix
%                 Is   = find(sx==SX);
%                 Ij   = find(gx==GX);
%                 B(Is,:,Ij,:)=Cmap;
%             end
%         end
end
end

