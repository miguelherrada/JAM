%Program to solve the steady state case using Fsolve
%domain s [0,1]
%discretizatiinf of s using ns uniform points and 2th finite differences
ns=2000;
[s,ds,ds2]=finites2th(ns,1); %ds and dss collocation matrices (1th and 2 second derivatives) [STEP 5]
%important ds and dss are sparse matrices!

%Init u and g in the 
u=(1-2*s);
%uA=s;
g=s;
%properties
nu=0.1 % viscosity
alpha1=0.01; 
%alpha1=0;

% initial guess
x0=[u,g];

 options=optimset('Display','iter','MaxIter',200,'MaxFunEvals',1000000,'Algorithm','levenberg-marquardt');   
%options=optimset('Display','iter','MaxIter',200,'MaxFunEvals',1000000,'Algorithm','trust-region');   
options=optimset('Display','iter','MaxIter',200,'MaxFunEvals',1000000)
tic
x=fsolve(@(x)  Fshape(x,ns,ds,ds2,alpha1,nu),x0,options); 
 u=x(1:ns);
 g=x(ns+1:2*ns);
 toc

x=g(2:ns); %x position
dx=g(2:ns)-g(1:ns-1); %distance between points
namef="\nu="+nu
subplot(2,1,1); hold on  ;plot(x,dx,'.-'); xlabel('x'), ylabel('dx'); box('on')
subplot(2,1,2); hold on  ;plot(g,u,'.-','Displayname',namef),; xlabel('x'), ylabel('u')
 




%F=0
function Fshape=Fshape(x,ns,ds,ds2,alpha1,nu)%  Ecuación de Laplace-Young + conservación del volumen

 %x vector N+1, con N primeros la forma y x(N+1), C
 u=x(1:ns);
 g=x(ns+1:2*ns);
 
 
 duds=u*ds';
 dudss=u*ds2';
 dgds=g*ds';
 dgdss=g*ds2';
 
 %case 
dudx=duds./dgds;
dudxx=(dudss./dgds-duds.*dgdss./dgds.^2)./dgds;
dgdx=dgds./dgds;
dgdxx=(dgdss./dgds-dgds.*dgdss./dgds.^2)./dgds;
%streaching
M=1./(alpha1*dudx.^2+0.4);
Ms=M*ds';
n=(dgds.^2).^1.5;
%2
F1=u.*dudx-nu*dudxx; %equation for u
F2=dgds.*dgdss-n.*Ms; %equation for u
%BC
F1(1)=u(1)-1;
F1(ns)=u(ns)+1;
F2(1)=g(1);
F2(ns)=g(ns)-1;



 


 Fshape=[F1;F2];
 
 



a=Fshape;
 end
