w=wA;
u=uA;
v=vA;

 
 for i=1:nrA
 for j=1:nzA
 for k=1:nx
z0(i,j,k)=z0A(j);
t0(i,j,k)=xx(k);
 end
 end
 end

wc=- w.*sin(t0) - u.*cos(t0).*cos(z0) - v.*cos(t0).*sin(z0);


uc=- w.*sin(t0) - u.*cos(t0).*cos(z0) - v.*cos(t0).*sin(z0);
vc=  u.*cos(z0).*sin(t0) - w.*cos(t0) + v.*sin(t0).*sin(z0);