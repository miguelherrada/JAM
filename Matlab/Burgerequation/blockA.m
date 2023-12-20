%% symbolic variables
syms s tau %space and time
syms u duds dudss dudtau g dgds dgdss dgdtau %u and g simbolics...

%1 mapping s to x+++++++++++++++++++++++++++++++++++++++++
dudx=duds/dgds;
dudxx=(dudss/dgds-duds*dgdss/dgds^2)/dgds;
dgdx=dgds/dgds;
dgdxx=(dgdss/dgds-dgds*dgdss/dgds^2)/dgds;
dudt=dudtau-dgdtau/dgds*duds;
dgdt=dgdtau-dgdtau/dgds*duds;

%step 1
%% Equation for bulk+++++++++++++++++++
%1 (u)
FAAb(1)=dudt+u*dudx-nu*dudxx;     
%x=g(s): we need to concetrate points around the gradients of u
%concentrating parameter M
M=1/(alpha1*dudx^2+0.4);%we are going to concetrate depending on the gradient of u
%you can choose your own M(s)
%M=1;
Ms=diff(M,duds)*dudss+diff(M,dgds)*dgdss;
%dl
n=(dgds^2)^1.5;
%2 (g)
FAAb(2)=dgds*dgdss-n*Ms; 
%FAAb(2)=dgdss; 
%% Boundary conditions
%left
FAAl(1)=u-1;
FAAl(2)=g;
%right
FAAr(1)=u+1;
FAAr(2)=g-1;

%% Mounting vectors and jacobian tensors (step2)
xs=[u, duds, dudss, dudtau, g, dgds, dgdss, dgdtau]; %step 2; %variables
%% Jacobians  (step 3)
for k=1:2
dFAAb(k,:)=jacobian(FAAb(k),xs);
dFAAr(k,:)=jacobian(FAAr(k),xs);
dFAAl(k,:)=jacobian(FAAl(k),xs);
end
%% Saving Functions+Analitical jacobians (step 4)
matlabFunction(FAAb,dFAAb,'file',[path_jacobian 'equationFAAb.m'],'vars',{s,xs,pa});
matlabFunction(FAAr,dFAAr,'file',[path_jacobian 'equationFAAr.m'],'vars',{s,xs,pa});
matlabFunction(FAAl,dFAAl,'file',[path_jacobian 'equationFAAl.m'],'vars',{s,xs,pa});
