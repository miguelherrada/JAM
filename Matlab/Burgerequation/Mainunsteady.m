%% JAM example: Burger equation (unsteady) 

%% Initial setup 
% Lets clean before start and decide where to locate jacobians functions
% and auxiliary elements
clear;
restoredefaultpath;
%Jacobians
path_jacobian= ['jacobians/'];
path(path,path_jacobian)
%Collocation matrices
path_jacobian1= ['../utils/collocationmatrices/'];
path(path,path_jacobian1)
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
ns=2001;
[s,ds,dss]=finites2thsparse(ns,1); %ds and dss collocation matrices (1th and 2 second derivatives) [STEP 5]
%important ds and dss are sparse matrices!

%Init u and g in the 

%Init u and g 
u=(1-2*s);
g=s;
%properties
nu=0.001 % viscosity (
alpha1=0.01; % streching factor

pa=[nu;alpha1];
% We construct an unique vector of unknowns _x0 (guess solution)[STEP 7]
x0=[u';g']; %vector [2N,1]; 
%time step
dt=0.01; %we are looking for steady solutions
%solution in the previus time m (t-dt), and mm (t-dt-dt1) 
x0m = x0;
x0mm = x0m;

else
%reading from file...
end


time(1)=0;

for ll=1:150 %this a loop to change the time
time(ll)=(ll-1)*dt;

   %% Newton method  (step 10)
    error=1e9;
    iter=0;
    while (error > 1e-2 && iter < 300)
     %getting matrices a (numerical jacobian DF) and vector b (fuction F)
     tic
     matrixAB
       %inverting the Jacobian     
        dxa=a\b;
      
        %computing error
        error=norm(dxa)
        if(error > 10^9) %if it diverges
         stop
        end
        %%
        % Apply the correction
        x0=x0+dxa;
      

   %% end Newton method
 
    end
  %interpolating solution in a fixed point x=0.25
  ux25(ll)=interp1(g,u,0.25)

  %interpolating solution in a fixed point x=0.25
  ux04(ll)=interp1(g,u,0.45)
  %updatting solutions
  x0mm=x0m;
  x0m=x0;

 %plotting the results for different nu
x=g(2:ns); %x position
dx=g(2:ns)-g(1:ns-1); %distance between points
namef="\nu="+nu
% subplot(2,1,1); hold on  ;plot(x,dx,'-'); xlabel('x'), ylabel('dx'); box('on')
% subplot(2,1,2); hold on  ;plot(g,u,'-','Displayname',namef); xlabel('x'), ylabel('u')
        end

%figure
hold on
nameP="ns="+ns+ ",dt="+dt +"\nu="+nu
plot(time, ux25,'*','Displayname',nameP); xlabel('t');ylabel('u_{25}')
plot(time, ux04,'o','Displayname',nameP); xlabel('t');ylabel('u_{25}')


%stop

