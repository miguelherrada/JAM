
 

%time derivatives
xt = (3*x0-4*x0m + x0mm)/(2*dt);
bp=3/(2*dt); 
% Initializing the variables 
yfA= zeros(NVA,NDA,ntA);
for i=1:NVA
%pointer in whole domain 
l=JM1{i};
variable= x0(l);
variablet= xt(l);
yfA(i,1,:)=variable;
yfA(i,2,:) = ddr0A*variable;
yfA(i,3,:) = ddz0A*variable;
yfA(i,4,:) = ddrr0A*variable;
yfA(i,5,:) = ddzz0A*variable;
yfA(i,6,:) = ddrz0A*variable;
yfA(i,7,:)= variablet;
yfA(i,8,:) = ddx0A*variable;
yfA(i,9,:) = ddxx0A*variable;
yfA(i,10,:) = ddzx0A*variable;
yfA(i,11,:) = ddrx0A*variable;
end

[FAA,DFAA]=evaluatingMatricesA(nrA,nzA,zA,rA,xA,yfA,NVA,NDA,ndA,pa,nx);

%Righthandside
 for j=1:NVA %equations
 C=FAA(j,:);
 eval(['xa' num2str(j) '= -C;']) %radial

for k=1:NVA  %jacobians
km=(k-1)*NDA+1;
kp=k*NDA;   
C=squeeze(DFAA(j,km:kp,:));
B=spdiags(C(1,:)',0,ntA,ntA)      +spdiags(C(2,:)',0,ntA,ntA)*ddr0A   +spdiags(C(3,:)',0,ntA,ntA)*ddz0A +...
 +spdiags(C(4,:)',0,ntA,ntA)*ddrr0A  +spdiags(C(5,:)',0,ntA,ntA)*ddzz0A  +spdiags(C(6,:)',0,ntA,ntA)*ddrz0A+...
 +spdiags(C(7,:)',0,ntA,ntA)*bp +spdiags(C(8,:)',0,ntA,ntA)*ddx0A  +spdiags(C(9,:)',0,ntA,ntA)*ddxx0A  +spdiags(C(10,:)',0,ntA,ntA)*ddzx0A +spdiags(C(11,:)',0,ntA,ntA)*ddzx0A;

eval(['aa' num2str(j) num2str(k) '=sparse(B);']) 

end
end

 

 
 
%BLOCK A
xa=[];
aa=[];
for i=1:NVA
    ax=[];
for j=1:NVA
aa0=eval(cat(2,'aa',num2str(i),num2str(j)));
ax=[ax,aa0];
end
bb0=eval(cat(2,'xa',num2str(i)));
aa=[aa;ax];
xa=[xa;bb0'];
end



%MATRX a and vector a
a=[aa];
b=([xa]);
