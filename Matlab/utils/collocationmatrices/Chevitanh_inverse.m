%%  Collocation matrices (with a stretching function)
%z: position, dz: 1st derivative,dz2:2th derivative, Rin:init point Rout: end point alpha stretching parameter 
function[z,dz,dz2]=Chevitanh_inverse(nr1,Rin, Rout,alpha)
%Chebyshev collocation matrices 
[x,dx,dx2]=Chevigood(nr1,1,0.); %[0 1]  
% Transforming function z = Rin +(Rout-Rin)* Tanh(a*(1-x)/Tanh(a); 0 < x < 1;
% Accumulates points at Rout 
a = alpha;
z=Rout + (Rin-Rout)*tanh(a*(1-x))/tanh(a);
dzdx=-2.*a*(Rin-Rout)/tanh(a)./(1+cosh(2*a*(1-x)));
dz=diag(1./dzdx)*dx;
dz2=dz*dz;