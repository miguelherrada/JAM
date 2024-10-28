%% Expansion of xo
% We pass from a single unknown vector _x0_ to "readable" variables
% just for inspection/drawing purposes
%
nbi = 0;
for j=1:nbl
    bl = list_block{j};
    NVAR = eval(['NV' bl]);
    nunkn = eval(['nt' bl]); %number of gridpoints in block
    for i=1:NVAR
        inic = nbi+(i-1)*nunkn+1;
        ifin = nbi+i*nunkn;
        lv = eval(['list_var_' bl '{' num2str(i) '}']);
        order = [lv '= full(reshape(x0(' num2str(inic) ':' num2str(ifin) '),nr' bl ',nz' bl ',nx));'];
        evalc(order);
    end
    nbi = nbi + nunkn*NVAR; 
end
%% Parameters are set in a readeable way
%




 



 
 
 
% %     taurr(i,:)=filter6a(taurr(i,:),nz);
% %     tauzz(i,:)=filter6a(tauzz(i,:),nz);
% %    tautt(i,:)=filter6a(tautt(i,:),nz);
% %    tauzr(i,:)=filter6a(tautt(i,:),nz);
% end

%filteting




 
  
  
%           %% Streamfunction of Block A
%         for j=1:nzA
%             psiA(1,j)=0;
%             for i=2:nrA
%                 %2D
%                %psiA(i,j)=psiA(i-1,j)+0.5*fA(i,j)*(r0A(i)-r0A(i-1))*(wA(i,j)+wA(i-1,j));
%                
%                 %axis
%                 psiA(i,j)=psiA(i-1,j)+0.5*fA(i,j)*(r0A(i)-r0A(i-1))*(ra(i,j)*wA(i,j)+ra(i-1,j)*wA(i-1,j));
%                
%             end
%         end
%         factor=psiA(nrA,1)/psiA(nrA,nzA)
%         
%         wA(:,nzA)=wA(:,nzA)*factor;
        
        
        
  
  
  
  
  
  
  
  
  
  
      % We construct an unique vector of unknowns _x0_
    %
    

  
    
    order='[';
    for j=1:nbl
        bl = list_block{j};
        NVAR = eval(['NV' bl]);
        for i=1:NVAR
            lv = eval(['list_var_' bl '{' num2str(i) '}']);
            order = [order 'reshape(' lv ',nt' bl ',1);'];
        end
    end
    x0=eval([order ']']);
