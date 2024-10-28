%% Setting equations of Block B
%number of variables
NVAR = NVA;
%number of derivatice
NDER = NDA;


%velocities cartesians u (vx) w( vz)   v(vy) entrance: z=0 z0=0
%w=Poiseuille)
syms u(r0,z0,t0,time) w(r0,z0,t0,time)  v(r0,z0,t0,time) p(r0,z0,t0,time) %mapping functions

%Mapping variables r0= \eta z0=s, t0=\theta 
assume (r0,"real")
assume (z0,"real")
assume (t0,"real")

% Parametric equations for the volume of the tube (see equation (15))
%X:
F = (cos(z0)*(2*L - pi*r0*cos(t0)))/pi;
%Z:
H = (sin(z0)*(2*L - pi*r0*cos(t0)))/pi;
G = -r0*sin(t0);

u=u(r0,z0,t0,time);
w=w(r0,z0,t0,time);
v=v(r0,z0,t0,time);
p=p(r0,z0,t0,time);

%% Maping 3D Euclidean 
X=[F,G,H]; %x,y,z
%X=[r0*cos(t0),r0*sin(t0),z0];
%X=[r0,z0,t0];
xs=[r0,z0,t0];
JX=jacobian(X,xs); %Jacobian matrix
Jinv=inv(JX);  %inverse Jacobian matrix
Jinv=simplify(Jinv);

%% selecting vector basis: cartesinan e_x=[1,0,0] e_y=[0,1,0] e_z=[0,0,1]
%temporal mapping
syms dr0dtime0  dz0dtime0  dt0dtime0 real
x=X(1);
y=X(2);
z=X(3);
eqn1 = diff(x,time) +diff(x,z0)*dz0dtime0 +diff(x,r0)*dr0dtime0+diff(x,t0)*dt0dtime0==0;
eqn2 = diff(y,time) +diff(y,z0)*dz0dtime0 +diff(y,r0)*dr0dtime0+diff(y,t0)*dt0dtime0==0;
eqn3 = diff(z,time) +diff(z,z0)*dz0dtime0 +diff(z,r0)*dr0dtime0+diff(z,t0)*dt0dtime0==0;
 [Aeqn,Beqn]=equationsToMatrix([eqn1,eqn2,eqn3],[dr0dtime0,dz0dtime0,dt0dtime0]);
Xeqn=linsolve(Aeqn,Beqn);
dr0dtime0=simplify(Xeqn(1));
dz0dtime0=simplify(Xeqn(2));
dt0dtime0=simplify(Xeqn(3));

% Cartesian velocity (base bbase)
U3D=[u,v,w];
%Euclides derivatives

%first derivatives
dU3Ddx=Jinv(1,1)*diff(U3D,xs(1))+Jinv(2,1)*diff(U3D,xs(2))+Jinv(3,1)*diff(U3D,xs(3));
dU3Ddy=Jinv(1,2)*diff(U3D,xs(1))+Jinv(2,2)*diff(U3D,xs(2))+Jinv(3,2)*diff(U3D,xs(3));
dU3Ddz=Jinv(1,3)*diff(U3D,xs(1))+Jinv(2,3)*diff(U3D,xs(2))+Jinv(3,3)*diff(U3D,xs(3));
dpdx=Jinv(1,1)*diff(p,xs(1))+Jinv(2,1)*diff(p,xs(2))+Jinv(3,1)*diff(p,xs(3));
dpdy=Jinv(1,2)*diff(p,xs(1))+Jinv(2,2)*diff(p,xs(2))+Jinv(3,2)*diff(p,xs(3));
dpdz=Jinv(1,3)*diff(p,xs(1))+Jinv(2,3)*diff(p,xs(2))+Jinv(3,3)*diff(p,xs(3));


%second
dU3Ddxx=Jinv(1,1)*diff(dU3Ddx,xs(1))+Jinv(2,1)*diff(dU3Ddx,xs(2))+Jinv(3,1)*diff(dU3Ddx,xs(3));
dU3Ddyy=Jinv(1,2)*diff(dU3Ddy,xs(1))+Jinv(2,2)*diff(dU3Ddy,xs(2))+Jinv(3,2)*diff(dU3Ddy,xs(3));
dU3Ddzz=Jinv(1,3)*diff(dU3Ddz,xs(1))+Jinv(2,3)*diff(dU3Ddz,xs(2))+Jinv(3,3)*diff(dU3Ddz,xs(3));



%% velocites operators

%gradient U in cartesians
gradientU3D=[dU3Ddx.',dU3Ddy.',dU3Ddz.'];
%Ugradient
UgradientU3D=w*zeros(1,3);
for i=1:3
UgradientU3D(i)=simplify(U3D(1)*dU3Ddx(i)+U3D(2)*dU3Ddy(i)+U3D(3)*dU3Ddz(i));
end
%laplace operator
laplaceU3D=simplify(dU3Ddxx+dU3Ddyy+dU3Ddzz);
%gradient (p)
gradientp3D=[dpdx,dpdy,dpdz];
%gradientp3D=simplify(gradientp3D);
%divergence of the velocity filed
diverU3D=simplify(dU3Ddx(1)+dU3Ddy(2)+dU3Ddz(3));

%time derivative
dUdtime3D=diff(U3D,time)+diff(U3D,r0)*dr0dtime0+diff(U3D,z0)*dz0dtime0...
+diff(U3D,t0)*dt0dtime0;







%% Momentum equations in Euclidean space (equation (17) in the paper)
MX=(dUdtime3D+UgradientU3D)+gradientp3D-laplaceU3D/Re;

Eq1=MX(1); %x Momentum
Eq2=MX(2); %y Momentum
Eq3=MX(3); %z Momentum

%% Continuity equation (equatin (16) in the paper)
Eq4=diverU3D; 

dimension = 2;

%Replacing simbolic by real variables
for k=1:NVAR  %variable
    for i=1:NDER  %derivatives
        %%
        % We create an array of symbolic variables, _yfA_
        %real
        yFA(k,i)=sym([ 'yFA','v',num2str(k),'d',num2str(i)],'real');
       
    end
    tco = list_var_A {k};
    [bas nto] = size(tco);
    list_variables {k} = tco(1:nto-1);
end
list_derivatives = list_der;



for k=1:NVAR
    name=['yFAsym(k,1)=' list_varsymbolic{k} ';'];
  eval(name)  
    
for i=2:NDER    
name=['yFAsym(k,i)=diff(' list_varsymbolic{k},list_dersymbolic{i}     ';'];
eval(name);
end
end

%equations bulk
for k=1:NVAR
for i=NDER:-1:1    
Eq1=subs(Eq1,yFAsym(k,i),yFA(k,i));
Eq2=subs(Eq2,yFAsym(k,i),yFA(k,i));
Eq3=subs(Eq3,yFAsym(k,i),yFA(k,i));
Eq4=subs(Eq4,yFAsym(k,i),yFA(k,i));
end
end


%getting explicit expressions for the derivatives and functions
for k=1:NVAR
    eval([list_variables{k} list_derivatives{1} ' = yFA(k,1);']);
for i=2:NDER    
eval(['d' list_variables{k} 'd' list_derivatives{i} ' = yFA(k,i);']);
end
end


%creating analythical functions 
  
  %bulk 
  FAA(1)=Eq1;  
  FAA(2)=Eq2;
  FAA(3)=Eq3; 
  FAA(4)=Eq4;  

%% Boundary conditions

%wall top: eta=1
FAAt(1)=u;
FAAt(2)=v;
FAAt(3)=w;
FAAt(4)=dpdr0;

%entrance left: s=0
FAAl(1)=u;
FAAl(2)=v ;
FAAl(3)=w-2*(1-r0^2);
FAAl(4)=dpdz0;

%outflow rigth: s=1
FAAr(1)=dudz0;
FAAr(2)=dvdz0 ;
FAAr(3)=dwdz0;
FAAr(4)=p;


%% Mounting vectors and jacobian tensors
%

x=reshape(yFA',NVAR*NDER,1); %variables

%%
% Jacobians 
for k=1:NVAR
dFAA(k,:)=jacobian(FAA(k),x);
dFAAt(k,:)=jacobian(FAAt(k),x);
dFAAr(k,:)=jacobian(FAAr(k),x);
dFAAl(k,:)=jacobian(FAAl(k),x);
end

%saving functions
matlabFunction(FAA ,dFAA, 'file',[path_jacobian 'equationFAA.m' ],'vars',{z0,r0,t0,x,pa});
matlabFunction(FAAr,dFAAr,'file',[path_jacobian 'equationFAAr.m'],'vars',{z0,r0,t0,x,pa});
matlabFunction(FAAl,dFAAl,'file',[path_jacobian 'equationFAAl.m'],'vars',{z0,r0,t0,x,pa});
matlabFunction(FAAt,dFAAt,'file',[path_jacobian 'equationFAAt.m'],'vars',{z0,r0,t0,x,pa});
