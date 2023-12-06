function[z,d,d2]=finites2thsparse(nz,H)
dz = H / (nz - 1);
z = 0:dz:H;
z=z;
dx = 0.5 / dz;
e=dx*ones(nz,1);
d= spdiags([-e e*0 e], -1:1, nz, nz);
%bordes
d(1,1)=-3*dx; %second order
d(1,2)=+4*dx;
d(1,3)=-dx;
d(nz,nz-2)=dx;
d(nz,nz-1)=-4*dx;
d(nz,nz)=3*dx;
dx=1/(dz*dz);
e=dx*ones(nz,1);
d2= spdiags([e -2*e e], -1:1, nz, nz);

d2(nz,nz)=2*dx;
d2(nz,nz-1)=-5*dx;
d2(nz,nz-2)=+4*dx;
d2(nz,nz-3)=-dx;
d2(1,1)=2*dx;
d2(1,2)=-5*dx;
d2(1,3)=+4*dx;
d2(1,4)=-dx;
