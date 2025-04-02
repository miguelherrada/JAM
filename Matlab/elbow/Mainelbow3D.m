%% JAM example II of the article "Goodbye Christoffel Symbols: A Flexible and Efficient Approach for Solving Physical Problems in Curved Spaces" by Miguel A. Herrada

%% Initial Setup
clear all;% Clear workspace
restoredefaultpath;
path_jacobian= ['jacobians3D/'];% Restore default MATLAB paths
path(path,[pwd '/' path_jacobian]);% Path for Jacobian functions
%Collocation matrices
path_jacobian1= ['../utils/collocationmatrices/'];% Path for collocation matrices
path(path,path_jacobian1)

%% Problem Configuration

% Number of parameters:
Np = 2; %  Re and L
% List of blocks
list_block = {'A' };
%List of variables: Cartesian velocity filed V=[v_x,v_y,v_z] and pressure
list_var_A = {'uA', 'vA',  'wA' 'pA'};% uA is  v_x, vA is v_y, wA is v_z and p is the pressure      
%symbolic velocity and pressure field
list_varsymbolic = {'u', 'v', 'w', 'p'}; 
%Number of variables
NV1=length(list_var_A);
[basura nbl ] = size(list_block);
% The list of derivatives: r0 is \eta z0 is s  and t0 is \thetata time is time 
list_der = {'', 'r0', 'z0',                'rr0',     'zz0'   ,  'rz0',   'time0',  't0' ,'tt0'   ,'tz0'    ,  'tr0'};
% The list of symbolic derivatives: 
list_dersymbolic = {'', ',r0)', ',z0)',  ',r0,r0)', ',z0,z0)' , ',r0,z0)' ,',time)' ',t0)',',t0,t0)',',t0,z0)',',t0,r0)' };
%Number of derivatives
ND1=length(list_der);

for i=1:nbl
    bl = list_block{i};
    order = ['[bas NV' bl '] = size(list_var_' bl ');'];
    evalc(order);
    order = ['[bas ND' bl '] = size(list_der);' ];
    evalc(order);
end


 %% Compute symbolic functions
lset=0  %If lset=1 we compile the equations
if(lset==1)

%parameters
for i=1:Np 
pa(i,1)=sym(['pa_',num2str(i)],'real');
end
Re=pa(1); % Reynolds number
L=pa(2); %legnt of the pipe (timesx radius)
blockAcartesianvelocity 
end

%%  Domain Discretization 
%Control parameters
%Reynolds number
Re=0.0001;
%Dimensionless length of the tube
L=10;
pa=[Re;L];
   
%\eta discretization
nrA=8; %number Chebyshev spectral collocation points in\theta
epsilon=0.001; %epsilon value
%x
[r0A,dr0A,drr0A]=Chevigood(nrA-1,1,0.001); %eta
%s  discretization
nzA=101; %4th order finite differences
[z0A,dz0A,dzz0A]=finites4thsparse(nzA,pi/2); %s
%theta dervivative
nx=8; %Fourier collocation points
[xx,dx] = fourdif(nx,1); %theta
[xx,dx2] = fourdif(nx,2);   

%Getting geometry
plottingelbowinit

% initialization of variable of block
uA=0*ones(nrA,nzA,nx);
wA=0*ones(nrA,nzA,nx);
vA=0*ones(nrA,nzA,nx);
pA=0*ones(nrA,nzA,nx);

rA=repmat(r0A', [1 nzA nx]);
%full collocation matrices
[ddr0A,ddrr0A,ddz0A,ddzz0A,ddx0A,ddxx0A,ddrz0A,ddzx0A,ddrx0A]= matrices3D (nrA,nzA,nx,dr0A,drr0A,dz0A,dzz0A,dx,dx2);

ntA=nrA*nzA*nx; %total number of grid points block A
nt1=ntA;
% pointers to BC
  pointerstoBC;
       
  %% Newton method 
    %Inittial guess solution x0
    order='[';
    for j=1:nbl
        bl = list_block{j};
        NVAR = eval(['NV' bl]);
        for i=1:NVAR
            lv = eval(['list_var_' bl '{' num2str(i) '}']);
            order = [order 'reshape(' lv ',nt' bl ',1);'];
        end
    end
    x0=eval([order ']']);
    

    
%single time step
dt=10^(10);
%previous solutions in t-dt, t-2*dt  
x0m = x0;
x0mm = x0m;



     error=1e9;
    iter=0;
    while (error > 1e-2 && iter < 300)
        % Increase the iteration counter
        iter=iter+1;

         tic  
  %getting matrices a (numerical jacobian DF) and vector b (fuction F)
       matrixAB3Doptima1
  %inversion
       dxa=(a\b);
            toc
       
   %errors   
        error=max(abs(dxa))
        
        if(error > 10^9)
            stop
        end
       % Apply the correction
        x0=x0+dxa;
       
       
       
      
    end
%making x0 redeable
 xotoredeable

 %% Plotting figure 3 in the paper
plottingslide

   
    




