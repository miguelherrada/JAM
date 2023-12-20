

function Fshape=Fshape(x,nz,zz,dz,dz2,Bo,V,Lambda)%  Ecuación de Laplace-Young + conservación del volumen

 %x vector N+1, con N primeros la forma y x(N+1), C
 f=x(1:nz);
 pref=x(nz+1);
 z=zz;
 
 fz=f*dz';
 fzz=f*dz2';
  
  %Interface curvature
k=(f.*fzz-1-fz.^2)./(f.*(1+fz.^2).^1.5);  

 
Fshape=pref+Bo*z+k;
  
 % Laplace-Young
 Fshape(1)=f(1)-1; 
 Fshape(nz)=f(nz)-1; 
 %volumen
 V0=0.;
 for i=2:nz
V0=V0+0.5*(zz(i)-zz(i-1))*(f(i)^2+f(i-1)^2); 
 end
 V0=V0/(2*Lambda);
 
 Fshape(nz+1)=V-V0;


