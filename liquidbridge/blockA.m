%% symbolic variables
syms s y t0 %mapped space(s,y) and time (t0)
syms w dwds dwdss dwdy dwdyy dwdt  dwdys dwdt0
syms u duds dudss dudy dudyy dudt  dudys dudt0
syms p dpds dpdss dpdy dpdyy dpdt  dpdys dpdt0
syms f dfds dfdss dfdy dfdyy dfdt  dfdys dfdt0

%Mapping to (r,z)
%radial position
  r=y.*f;
 %axial position
  z=s;
  %mapping factors
  fs=dfds;
  fss=dfdss;
  mi=1./f;
  mis=-fs./f.^2;
  miss=-fss./f.^2+2*fs.^2./f.^3;
  r1=r.*mis;
  r2=r1.^2;
  r3=r.*miss;
%Radial derivatives
dudr=dudy.*mi;
dpdr=dpdy.*mi;
dwdr=dwdy.*mi;
dudrr=dudyy.*mi.^2;
dwdrr=dwdyy.*mi.^2;
 %Axial derivatives      
 dwdz=dwds+r1.*dwdy;
 dudz=duds+r1.*dudy;
 dpdz=dpds+r1.*dpdy;
 dwdzz=dwdss+2*r1.*dwdys+r2.*dwdyy+r3.*dwdy;
 dudzz=dudss+2*r1.*dudys+r2.*dudyy+r3.*dudy; 
 %temporal   mapping 
 ft=dfdt0;
 tt=-ft.*r./(f.^2);
 dudt=dudt0+tt.*dudy;     
 dwdt=dwdt0+tt.*dwdy;   
 
 nv=4 %number of equations
%tangential and normal vector+curvature
n=1+fs.^2;
%tangent vector
t0z=1./n.^0.5;
t0r=fs./n.^0.5;     
%normal vector
n0z=t0r;
n0r=-t0z;     
%stress tensor
taurr=2*dudr;
tauzr=dwdr+dudz;
tauzz=2*dwdz;
%total stresses
%normal
tauni=n0r.*taurr.*n0r+n0z.*tauzz.*n0z+2*n0z.*tauzr.*n0r; 
%tangential
tauti=t0r.*taurr.*n0r+t0z.*tauzz.*n0z+t0z.*tauzr.*n0r+n0z.*tauzr.*t0r; 
 %Interface curvature
 k=(f.*fss-1-fs.^2)./(f.*(1+fs.^2).^1.5);



%% Equation for bulk+++++++++++++++++++
%1 Axial Momentum
FAAb(1)=-(dwdt+u.*dwdr+w.*dwdz+dpdz)+(dwdzz+dwdrr+dwdr./r)*Cmu;  
%2 Radial Momentum
FAAb(2)=-(dudt+u.*dudr+w.*dudz+dpdr)+(dudzz+dudrr+dudr./r-u./r.^2)*Cmu; 
%3 (p) continuity
FAAb(3)=dwdz+dudr+u./r;
%4 (f)
FAAb(4)=dfdyy; 
 
%% Boundary conditions
%left:wall
FAAl(1)=w;
FAAl(2)=u;
FAAl(3)=dpdz;
FAAl(4)=f-1;
%right:wall
FAAr(1)=w;
FAAr(2)=u;
FAAr(3)=dpdz;
FAAr(4)=f-1;
%top (interface)
%normal balance
FAAt(1)=p+Bo*z+k-tauni*Cmu;
%tangential balance
FAAt(2)=Cmu*tauti;
%pressure
FAAt(3)=dpdr;
%
FAAt(4)=ft+w.*fs-u;
%bottom(axis)
FAAd(1)=dwdr;
FAAd(2)=u;
FAAd(3)=dpdr;
FAAd(4)=dfdy;



%% Mounting vectors and jacobian tensors 
xs=    [w,dwds, dwdss, dwdy, dwdyy dwdys dwdt0]; %step 2; %variables
xs=[xs,[u,duds, dudss, dudy, dudyy dudys dudt0]];
xs=[xs,[p,dpds, dpdss, dpdy, dpdyy dpdys dpdt0]];
xs=[xs,[f,dfds, dfdss, dfdy, dfdyy dfdys dfdt0]];
%% Jacobians  (step 3)
for k=1:nv
dFAAb(k,:)=jacobian(FAAb(k),xs);
dFAAr(k,:)=jacobian(FAAr(k),xs);
dFAAl(k,:)=jacobian(FAAl(k),xs);
dFAAt(k,:)=jacobian(FAAt(k),xs);
dFAAd(k,:)=jacobian(FAAd(k),xs);
end
%% Saving Functions+Analitical jacobians (step 4)
matlabFunction(FAAb,dFAAb,'file',[path_jacobian 'equationFAAb.m'],'vars',{s,y,xs,pa});
matlabFunction(FAAr,dFAAr,'file',[path_jacobian 'equationFAAr.m'],'vars',{s,y,xs,pa});
matlabFunction(FAAl,dFAAl,'file',[path_jacobian 'equationFAAl.m'],'vars',{s,y,xs,pa});
matlabFunction(FAAt,dFAAt,'file',[path_jacobian 'equationFAAt.m'],'vars',{s,y,xs,pa});
matlabFunction(FAAd,dFAAd,'file',[path_jacobian 'equationFAAd.m'],'vars',{s,y,xs,pa});
