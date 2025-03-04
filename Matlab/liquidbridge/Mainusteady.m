%% Solver JAM (Miguel Ángel Herrada)
%Exercise:  Nonlinear oscillations in a liquid bridge:
%[M. A. Herrada; J. M. Montanero. 2016. A numerical method to study the dynamics of capillary fluid systems. Journal of Computational Physics, 306:137-147, 2016.

%% Initial setup 
% Lets clean before start and decide where to locate jacobians functions
% and auxiliary elements
clear;
restoredefaultpath;
path_jacobian= ['jacobians/'];
path(path,path_jacobian)
%Collocation matrices
path_jacobian1= ['../utils/collocationmatrices/'];
path(path,path_jacobian1)
lset=0 % 1 to compute the Analytical Jacobians
if (lset==1)
 
%% Create symbolic functions by
% number of parameters
Np = 4; 
for i=1:Np 
    pa(i,1)=sym(['pa_',num2str(i)],'real');
end
   Cmu=pa(1); %Capilar number, Cmu=mu/(rho*sigma*r_o)^0.5
    Bo=pa(2); %Bond number, Bo=rho*g*r_o^2/sigma
Lambda=pa(3); %slenderness parameter Lambda=L/(r_o)
     V=pa(4); %dimimensioless volume
%% Blocks evaluation
blockA
stop
end                                           

% %% Initialization of the variables

lstart=0; % start from rest
if (lstart==0)
%Parameters (figure 4 Herrada&Montanero 2016)
Cmu=0.0349;
Bo=0.478;
Lambda=1.72;
V=0.925; 

%% Initializing parameters
%domains: mapped domain: [y,s]= [0,1]x[0,2Lambda]. s is "\xi" and y is "\eta" in Herrada&Montanero (2016)
%          Real  domain:  z=s  and r=f(s,t). f is the shape (R)  and t is time
%discretization of y(\xi in the paper) using ny chebyshev points 
ny=9;%number of y points
[y,dy,dy2]=Chevigood(ny-1,1,0);
%discretization of s (\eta in the paper) using nx uniform points and 2th finite differences
ns=61;%number of  s points
%discretization of s using 4th order finite differences
[s,ds,ds2]=finites4thsparse(ns,2*Lambda); 
%Getting the full collocation matrices
matrixalg_A



%Init field: u (radial),w (axial) ,p (pressure), f (liquid bridge shape)
w=zeros(ny,ns);
u=zeros(ny,ns);
p=ones(ny,ns);
%Getting the initial shape (r=f_i(s)) by solving the Laplace-Young equation
f=ones(1,ns); 
pref=1;
x0=[f,pref];
options=optimset('Display','iter');   
x=fsolve(@(x) Fshape(x,ns,s,ds,ds2,Bo,V,Lambda),x0,options);
f=x(1:ns);
pref=x(ns+1);


%reshaping
f=repmat(f, [ny 1]);
p=pref*ones(ny,ns);
%total points
N=ns*ny;
% We construct an unique vector of unknowns _x0 (guess solution)[STEP 7]
x0=[];
x0=[x0;reshape(w,N,1)]; %vector [2N,1]; 
x0=[x0;reshape(u,N,1)];
x0=[x0;reshape(p,N,1)];
x0=[x0;reshape(f,N,1)];
%time step
dt=0.1
dt1=0.1
%solution in the previus time m (t-dt), and mm (t-2dt) 
x0m = x0;
x0mm = x0;
%introducing perturbation in Bo
Bo=Bo+0.085*Bo;
%time step


else
%reading from file...

end
%evaluating parameters
pa=[Cmu;Bo;Lambda;V];
%end of the simulation
tend=130
%number of time step
NT=tend/dt+1
for ll=1:NT %this a loop to decrease the diffusivity nu
  %% Newton method  (step 10)
    error=1e9;
    iter=0;
    while (error > 1e-4 && iter < 300)
     %getting matrices a (numerical jacobian DF) and vector b (fuction F)
      matrixAB
    %inverting the Jacobian    
        dxa=a\b;
       
        %computing error
        error=norm(dxa);
        if(error > 10^9) %if it diverges stop
         stop
        end
        %%
        % Apply the correction
        x0=x0+dxa;
 %making variables readables...
 w=reshape(x0(lw),ny,ns);
 u=reshape(x0(lu),ny,ns);
 p=reshape(x0(lp),ny,ns);
 f=reshape(x0(lf),ny,ns);
   
    end
    %% end Newton method


timef(ll)=(ll-1)*dt; %time
%computing the shape at z=Lambda
f0(ll)=f(ny, (ns-1)/4+1);
t=timef(ll)
f00=f0(ll)
%updatting solutions
x0mm=x0m;
x0m=x0;
end
hold on

%plotting figure 4 Herrada&Montanero (2016)
plot(timef,f0); xlabel('t'); ylabel('R(\Lambda/2,t)')
