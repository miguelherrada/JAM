%% Solver JAM (Miguel Ángel Herrada)
%Example for course in China

%% Initial setup 
% Lets clean before start and decide where to locate jacobians functions
% and auxiliary elements
clear;
restoredefaultpath;
path_jacobian= ['jacobians/'];
path(path,path_jacobian)
% number of parameters
Np = 2; %viscosity, streaching factor

lset=0 % 1 to compute the Analitical Jacobians
     if (lset==1)
 
%% Create symbolic functions by
for i=1:Np 
    pa(i,1)=sym(['pa_',num2str(i)],'real');
end
nu=pa(1); %the knimatic viscosity
alpha1=pa(2); % parameter to concetrate points 
%% Blocks evaluation
blockA

end                                           

% %% Initialization of the variables
%
%  _lstart_ control this initialization
% 
lstart=0;
if (lstart==0)
    %% Initializing parameters
%domain s [0,1]
%discretizatiinf of s using ns uniform points and 2th finite differences
ns=2000; %ns=5000

 % tic
 % [s,ds,dss]=finites2th(ns,1); %ds and dss collocation matrices (1th and 2 second derivatives) [STEP 5]
 % toc

tic
[s,ds,dss]=finites2thsparse(ns,1); %ds and dss collocation matrices (1th and 2 second derivatives) [STEP 5]
toc
% tic
% [s,ds,dss]=finitas2(ns,1);
% toc
%important ds and dss are sparse matrices!

%Init u and g 
u=(1-2*s);
g=s;
% u=s.^2
%  uss=u*dss'
%  plot(u,uss)
%  stop
%properties
nu=0.1 % viscosity
alpha1=0.01; % streching factor

pa=[nu;alpha1];
% We construct an unique vector of unknowns _x0 (guess solution)[STEP 7]
x0=[u';g']; %vector [2N,1]; 

%time step
dt=10^10; %we are looking for steady solutions
%solution in the previus time m (t-dt), and mm (t-2dt) 
x0m = x0;
x0mm = x0;


else
%reading from file...
load("solutions/solutionsteady.mat")
nu=pa(1); %the knimatic viscosity
alpha1=pa(2); % parameter to concetrate points 
x0m=x0;
x0mm=x0;
end



tic
for ll=1:1 %this a loop to decrease the diffusivity nu
%  if(ll>1)&&(ll<12)
      % nu=0.9*nu;
      % pa(1)=nu;
%  end
% if(ll>=12)
% nu=0.9*nu;
% end
%nu=0.95*nu


   %% Newton method  (step 10)
    error=1e9;
    iter=0;
    while (error > 1e-4 && iter < 300)
     %getting matrices a (numerical jacobian DF) and vector b (fuction F)
      matrixAB
    %inverting the Jacobian    
        dxa=a\b;
       
        %computing error
        error=norm(dxa)
        if(error > 10^9) %if it diverges stop
         stop
        end
        %%
        % Apply the correction
        x0=x0+dxa;
        
   %% end Newton method
    end
 %making variables readables...
 u=x0(1:ns)';
 g=x0(ns+1:2*ns)';


%plotting the results for different nu
x=g(2:ns); %x position
dx=g(2:ns)-g(1:ns-1); %distance between points
namef="\nu="+nu
 % subplot(2,1,1); hold on  ;plot(x,dx,'-'); xlabel('x'), ylabel('dx'); box('on')
 % subplot(2,1,2); hold on  ;plot(g,u,'-','Displayname',namef), xlabel('x'), ylabel('u')
      
end
toc
%saving file
lsave=0;
if(lsave==1)
save('solutions/solutionsteady','x0','pa','ns','s','ds','dss','dt')
end

%%computing eigenvalue problem
leigen=1;
if(leigen==1)
 n=30; %number of eigen values
 options.Display='off'
 omega=0; %to choose where to look at
 matrixABeigen
tic
 [V,d] = eigs(a,b,n,omega,options);
toc
%  tic
% [V,d] = eig(full(a),full(b));
% toc
 % [V,d] = eig(full(a),full(b));
omegar=real(diag(d));
omegai=imag(diag(d));


%figure; 
hold on; plot(omegar,omegai,'o');xlim([-0.001,0.001]);
end    


%stop

