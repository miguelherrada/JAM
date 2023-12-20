function[z,d,d2]=finites2th(nz,H)
z=zeros(1,nz);
z=z;
dz=H/(nz-1);
for j=1:nz
z(j)=(j-1)*dz;
end
%use c,x,dx,dx2,dp,u,v,du
dx=0.5/dz;
% dx1=1/(12*dz);
% dx2=1/dz^2;
% Matrix of derivatives

d=zeros(nz,nz);


d(1,1)=-3*dx; %second order
d(1,2)=+4*dx;
d(1,3)=-dx;


for j=2:nz-1; %2forth order

d(j,j-1)=-dx;
d(j,j)=0;
d(j,j+1)=dx;

end
% d(nz-1,nz-3)=dx;
% d(nz-1,nz-2)=-4*dx;
% d(nz-1,nz-1)=3*dx;


d(nz,nz-2)=dx;
d(nz,nz-1)=-4*dx;
d(nz,nz)=3*dx;

dx=1/(dz*dz);
d2=zeros(nz,nz);
d2(1,1)=2*dx;
d2(1,2)=-5*dx;
d2(1,3)=+4*dx;
d2(1,4)=-dx;

for j=2:nz-1; %forth order

d2(j,j-1)=dx;
d2(j,j)=0-2*dx;
d2(j,j+1)=dx;

end
% 
%   d2(nz-1,nz-1)=2*dx;
%   d2(nz-1,nz-2)=-5*dx;
%   d2(nz-1,nz-3)=+4*dx;
%   d2(nz-1,nz-4)=-dx;


d2(nz,nz)=2*dx;
d2(nz,nz-1)=-5*dx;
d2(nz,nz-2)=+4*dx;
d2(nz,nz-3)=-dx;




d=sparse(d);
d2=sparse(d2);


