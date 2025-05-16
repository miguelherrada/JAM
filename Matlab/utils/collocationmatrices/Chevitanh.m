%%  Collocation matrices (with a stretching function)
%z: position, dz: 1st derivative,dz2:2th derivative, Rin:init point Rout: end point alpha stretching parameter 
function[z,dz,dz2]=Chevitanh(nr1,Rin, Rout,alpha)
%Chebyshev collocation matrices 
[x,dx,dx2]=Chevigood(nr1,1,0.); %[0 1]  
% Transforming function z = Rin +(Rout-Rin)* Tanh(alphar*x)/Tanh(alphar); 0 < x < 1;
% Accumulates points at Rout 
z=Rin + (Rout-Rin)*tanh(alpha*x)/tanh(alpha);
dzdx=2.*alpha*(Rout-Rin)/tanh(alpha)./(1+cosh(2*alpha*x));
dz=diag(1./dzdx)*dx;
dz2=dz*dz;