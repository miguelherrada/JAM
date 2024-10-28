function[z,d,d2]=finitas4ordentotal(nz,H)
z=zeros(1,nz);
dz=H/(nz-1);
for j=1:nz
z(j)=(j-1)*dz;
end
%use c,x,dx,dx2,dp,u,v,du
dx=0.5/dz;
dx1=1/(12*dz);
dx2=1/(12*dz^2);
% Matrix of derivatives

d=zeros(nz,nz);
d2=zeros(nz,nz);
% -(25*f0 - 48*f_1 + 36*f_2 - 16*f_3 + 3*f_4)/(12*dz)
d(1,1)=-25*dx1; %4thorder
d(1,2)=+48*dx1;
d(1,3)=-36*dx1;
d(1,4)=+16*dx1;
d(1,5)=-3*dx1;
%(35*f0 - 104*f_1 + 114*f_2 - 56*f_3 + 11*f_4)/(12*dz^2)
d2(1,1)=35*dx2; %4thorder
d2(1,2)=-104*dx2;
d2(1,3)=+114*dx2;
d2(1,4)=-56*dx2;
d2(1,5)=11*dx2;
%-(10*f0 + 3*f_1 - 18*f_2 + 6*f_3 - f_4)/(12*dz)
d(2,1)=-3*dx1; %4thorder
d(2,2)=-10*dx1;
d(2,3)=18*dx1;
d(2,4)=-6*dx1;
d(2,5)=+1*dx1;
%(11*f_1 - 20*f0 + 6*f_2 + 4*f_3 - f_4)/(12*dz^2)
d2(2,1)=+11*dx2; %4thorder
d2(2,2)=-20*dx2;
d2(2,3)=6*dx2;
d2(2,4)=4*dx2;
d2(2,5)=-1*dx2;


for j=3:nz-2; %forth order

d(j,j-2)=dx1;
d(j,j-1)=-8*dx1;
d(j,j)=0;
d(j,j+1)=8*dx1;
d(j,j+2)=-dx1;


d2(j,j-2)=-dx2;
d2(j,j-1)=16*dx2;
d2(j,j)=-30*dx2;
d2(j,j+1)=16*dx2;
d2(j,j+2)=-dx2;

end
 %(10*f0 + 3*f_1 - 18*f_2 + 6*f_3 - f_4)/(12*dz)
d(nz-1,nz-4)=-1*dx1; 
d(nz-1,nz-3)=+6*dx1; 
d(nz-1,nz-2)=-18*dx1;
d(nz-1,nz-1)=10*dx1;
d(nz-1,nz)=3*dx1;

% (25*f0 - 48*f_1 + 36*f_2 - 16*f_3 + 3*f_4)/(12*dz)
d(nz,nz)=25*dx1; %4thorder
d(nz,nz-1)=-48*dx1;
d(nz,nz-2)=+36*dx1;
d(nz,nz-3)=-16*dx1;
d(nz,nz-4)=+3*dx1;



%(11*f_1 - 20*f0 + 6*f_2 + 4*f_3 - f_4)/(12*dz^2)

d2(nz-1,nz)=+11*dx2; %4thorder
d2(nz-1,nz-1)=-20*dx2;
d2(nz-1,nz-2)=+6*dx2;
d2(nz-1,nz-3)=+4*dx2;
d2(nz-1,nz-4)=-1*dx2;




d2(nz,nz)=35*dx2; %4thorder
d2(nz,nz-1)=-104*dx2;
d2(nz,nz-2)=+114*dx2;
d2(nz,nz-3)=-56*dx2;
d2(nz,nz-4)=11*dx2;

d2=sparse(d2);
d=sparse(d);

