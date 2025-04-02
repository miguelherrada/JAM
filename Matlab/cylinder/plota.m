
%% Solver de Miguel Angel Herrada

%% Initial setup 
% Lets clean before start and decide where to locate jacobians functions
% and auxiliary elements
%

clear all;
restoredefaultpath;
path_jacobian= ['eigen/jacobians2D1/'];
path(path,[pwd '/' path_jacobian]);
path_jacobian1='C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Dropbox\subroutinesMatlab/';
path_jacobian1='C:\Users\usuario\Dropbox\subroutinesMatlab/';

%path_jacobian1='subroutinesMatlab/';
path(path,path_jacobian1)

%% Details of the block, variables and connections
%
%
% Fix the number of blocks (subdomains). In each block it has to be defined
% the number of variables and the number of symbolic variables. 
%
%%
% number of parameters

Np = 4;

% To construct the numerical scheme we need to use the list of variables of
% each block (*Important: In the same order as defined in the corresponding
% block*). 

list_block = {'B'};
list_var_B = {'wB', 'uB', 'pB', 'FB' 'GB', 'VTB','p0B','VB','HB' , 'CB'}; %P0B internal pressure VB mass conservaion
name_Boundary_B = { 'BB'  'LineaBA'   'VertexBA1' 'VertexBA2' 'LineaBt1', 'LineaBt2'  'LineaBl', 'LineaBr'};
list_Boundary_B =(1:1:length(name_Boundary_B));




[basura nbl ] = size(list_block);


% The list of derivatives. We have assumend that the number and distribution
% is the same in all the block. This information is relevant in _Matrix_
%
%list_der = {'', 'r0', 'z0', 'rr0', 'zz0', 'rz0',  'time'  };
list_der = {'', 'r0', 'z0', 'rr0', 'zz0', 'rz0', 't0'  };



for i=1:nbl
    bl = list_block{i};
    order = ['[bas NV' bl '] = size(list_var_' bl ');'];
    evalc(order);
    order = ['[bas ND' bl '] = size(list_der);' ];
    evalc(order);
end

namef=['soluciones/mixta/']
%namef=['soluciones/cilindrical/']

%namef=['soluciones/spherical/']
 aaa = dir(fullfile(namef, '*.mat'));
%aaa = dir(fullfile('C:\Users\herrada\Dropbox\Marco\Jens1\cases\FULL\case81D/', '*.mat'));




 [~,index] = sortrows({aaa.date}.'); 
 aaa = aaa(index);
 clear index
  %[~,index] = sortrows({aaa.name}.'); aaa = aaa(indeox); clear index
% Camino B: Subiendo pext manteniendo sigma=cte
jmax=length(aaa)



lplot=0  %shapes
lplot=3 %comparison with experiment
lplot=1 %streamlineas and voritcity contours
lplot=5 %curvature 
lplot=3
lplot=5 %length
lplot=3
if(lplot==4)

 
 mov = VideoWriter('We10Oh0,1a','MPEG-4')
 %profiles = VideoWriter.getProfiles()
mov.FrameRate=2
 open(mov)
end



%% List of type of plotting
%


%2,4
j0=2
for js=jmax:jmax
    name=[namef,aaa(js).name]
   % name= 'soluciones/water/mesh2/Radius0,0006Rout15nrB1nrB91nzB351nzB351eps0alphar2'
%
hold on
%filename='soluciones/lambda1e-7/caseCB_0,33Re0 Rout6Hout6rho0lambda1e-07VF1nrB5nrB21nzB1301alphar2,5alphaz3,5'
load(name);
x0=x0m;
Bond=pa(1); %Stone G*R*mu2  %G shear rate, R radius(a) mu external viscosity  
    Re=pa(2);%Reynolds
    Rout=pa(3); % Box radius 
    VF=pa(4); %dimensionless volume VF=1 is a semisphere or radius 1
    Ma=pa(5); %Marangini
    Pes=pa(6); %Peclet S
    Peb=pa(7); %PecletB
    Ci=pa(8);
    AL=pa(9); % ratio mol/m^3
   %alphaf=2
Cireal=0.002
  %Radius=Radius1;
   % parameterPEsciDs8 



 Radius=Rbubble/1000
% 
%sigma=2

 % muw=5*muw
  %teminal velocity
    % Re=rhow*sigma*Radius/muw^2
    % 
    % Bond=rhow*g*Radius^2/sigma
  

 

    matrixalg_B;



    %%
    % The grid can be refined/coarsened. If _interpol_ = 1 the existing solution
    % will be interpolated to the new grid.
       xotoredeable
       hold on
       gettingcilindricalsvelocites
figure
contourf(GB,FB,wc)
stop
          %computing legnth    
   Lbubble(js)=0;
   fz=FB(1,:)*dz0B';
   gz=GB(1,:)*dz0B';
   ds=(fz.^2+gz.^2).^0.5;
   for j=2:nzB
   Lbubble(js)=Lbubble(js)+0.5*(z0B(j)-z0B(j-1))*(ds(j)+ds(j-1));
   end   
%  plot(z0B,wB(nrB,:))
%     stop


hold on
%       figure
 
plottingmesh0
%caracterisc time
omega0=0.0043
t0=2*pi/omega0
t0new=t0/1/VTB(1,1)
for i=1:nrB
pB(i,:)=filter6a(pB(i,:),nzB)  ; 
end 
j0=(nzB-1)/2+1;
pB(:,j0)=pB(:,j0-1);
 

for i=1:nrB

    pB(i,:)=filter6a(pB(i,:),nzB)  ; 
    CB(i,:)=filter6a(CB(i,:),nzB)  ; 
end 

%plot(z0B,pB)



 Vol=0;
     f=FB(1,:);    
 fz=FB(1,:)*dz0B';
 gz=GB(1,:)*dz0B';
 fza=FB(1,:);
 n=real((fz.^2+gz.^2).^0.5);
     dVol=pi*(f.^2.*gz);   
     for j=2:nzB-1
     Vol=Vol+0.5*(z0B(j)-z0B(j-1))*(dVol(j)+dVol(j-1));
     end
     a=Vol/(4/3*pi)
     rho=0.001
Etrust=(1-rho)*Bond*Vol
 a=Etrust/(4*pi*VTB(1,1))
VT=VTB(1,1);

%Vc=(sigma/(rhow*Radius))^0.5;
Vc=sigma/muw;
%computing 
VT=VT*Vc;
VTf(js)=VT;
Xia(js)=max(2*FB(1,:))/abs((GB(1,1)-GB(1,nzB)));
We1(js)=rhow*VTf(js)^2*Radius/(sigma*Xia(js));
We2(js)=rhow*VTf(js)^2*Radius*abs(GB(1,1)-GB(1,nzB))/(sigma);


Radiusf(js)=Radius;

Oh(js)=Oh;


gettingcilindricalsvelocites
% stop
switch lplot;
    case 0
         name1='R= '+string(Radiusf(js)*1000)+' mm';
            name1="nrB="+nrB+"nrB="+nrB+"nzB="+nzB+"alphar="+alphar
     %  velocities(GB,FB,GB,FB,pB,pB)
       gettingcilindricalsvelocites
     figure
  contourf(FB,GB,pB)
  stop
     hold on

     plot(GB(1,:),wB(1,:),'ro-','DisplayName',name1)
     plot(GB(:,1),wB(:,1),'go-','DisplayName',name1)
     plot(GB(:,nzB),wB(:,nzB),'go-','DisplayName',name1)

   
%       plot(-FB(1,:),GB(1,:),'m-')
%      hold on
%      axis('equal')
%        title(name1)
%        xlim([-1.5,1.5])
%         ylim([-1.2,1.2])
%        xlabel('r_s/R')
%        ylabel('z_s/R')
% 
%        Xi=2*max(FB(1,:))/(GB(1,1)-GB(1,nzB));

    case 1

         matrixBBCvorticity
         a1=min(min(vorticityBa))
         a2=max(max(vorticityBa))
         v1=[a1:-a1/10:-0.0001]
contourf(FB,GB,vorticityBa,v1)

 %creating cartesing grid
x=[reshape(GB,nrB*nzB,1);reshape(GB,nrB*nzB,1);reshape(GB,nrB*nzB,1);reshape(GB,nrB*nzB,1)];
y=[reshape(FB,nrB*nzB,1);reshape(FB,nrB*nzB,1);reshape(-FB,nrB*nzB,1);reshape(-FB,nrB*nzB,1)];
w=[reshape(wB,nrB*nzB,1);reshape(wB,nrB*nzB,1);reshape(wB,nrB*nzB,1);reshape(wB,nrB*nzB,1)];
u=[reshape(uB,nrB*nzB,1);reshape(uB,nrB*nzB,1);reshape(-uB,nrB*nzB,1);reshape(-uB,nrB*nzB,1)];
p=[reshape(pB,nrB*nzB,1);reshape(pB,nrB*nzB,1);reshape(pB,nrB*nzB,1);reshape(pB,nrB*nzB,1)];
Fw = scatteredInterpolant(x,y,w);
Fu = scatteredInterpolant(x,y,u);
Fp = scatteredInterpolant(x,y,p);

x0=linspace(min(-3),max(3),600);
y0=linspace(min(-3),max(0),600);
[X,Y]=meshgrid(x0,y0);
W=Fw(X,Y);
U=Fu(X,Y);
P=Fp(X,Y);

plot(FB(nrB,:),GB(nrB,:),'m-','LineWidth',2)
plot(-FB(nrB,:),GB(nrB,:),'m-','LineWidth',2)
axis('equal')

% starty=[-5:0.125:5]
% startx=5*ones(size(starty))

starty=[-3:0.1:0]
startx=2.5*ones(size(starty))

% startx=9.99;
% starty=1;
options=[0.1,500000];
%contourf(Y,X,P)
h=streamline(Y',X',U',W',starty,startx,options)
Y0=Y(1:10:600,1:10:600);
X0=X(1:10:600,1:10:600);
U0=U(1:10:600,1:10:600);
W0=W(1:10:600,1:10:600);
scale_factor = 0.5;
q=quiver(Y0',X0',U0',W0','b-','ButoScale','on','ButoScaleFactor',1.5)
%q.ShowBrrowHead = 'off';
%q.Marker = '.';


set( h, 'Color', 'red' )
starty=[-1:0.1:1]
startx=0.5*ones(size(starty))


options=[0.1,5000000];
%contourf(Y,X,P)
starty=[-0.045:0.04:0]
startx=3*ones(size(starty))
h=streamline(Y',X',U',W',starty,startx,options)

set( h, 'Color', 'r' )


% h=streamline(Y',X',U',W',starty,startx,options)
% set( h, 'Color', 'red' )

xlim([-2,2])
ylim([-1.5,1.5])

if ( abs(Radius- 1.0000e-03)<1e-6)
    options=[0.001,200000];
%contourf(Y,X,P)
starty=[-0.2:0.05:0]
startx=-0.6*ones(size(starty))

% clear starx 
% clear stary
% startx(1)=-0.195;
% starty(1)-0.39
% startx(1)=-0.195;
% starty(1)-0.39
h=streamline(Y',X',U',W',starty,startx,options)
set( h, 'Color', 'g' )
end

 if ( abs(Radius- 9.5000e-04)<1e-6)
    options=[0.001,200000];
%contourf(Y,X,P)
starty=[-0.2:0.05:0]
startx=-0.4*ones(size(starty))

% clear starx 
% clear stary
% startx(1)=-0.195;
% starty(1)-0.39
% startx(1)=-0.195;
% starty(1)-0.39
h=streamline(Y',X',U',W',starty,startx,options)
set( h, 'Color', 'r' )
end    







 name1='R= '+string(Radiusf(js)*1000)+' mm';

%  xlabel('x/R')
%  ylabel('z/R')
% title(name1)
   


    
    
    
    case 2
        hold on
         plot(GB(1,:),wB(1,:),'c-')
%  
      plot(GB(:,1),wB(:,1),'g-')
       plot(GB(:,nzB),wB(:,nzB),'r-')
      vz=wB(nrB,:);
      vr=uB(nrB,:);
      fz=FB(nrB,:)*dz0B';
      hz=GB(nrB,:)*dz0B';
      n=(fz.^2+hz.^2).^0.5;
      vt=(fz.*vr+hz.*vr)./n;

     % plot(GB(nrB,:),vt)

       %contourf(GB,FB,wB)
   
          




    case 4
    axis tight
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

    h=figure
 namea=['t=', po2com(timef(length(timef)))]
 title(namea)
    hold on
% plottingmesh
% axis('equal')
 %velocities(GB,FB,GB,FB,wB,wB)
%  plottingmesh
%  axis('equal')
%  ylim([1,2.2])
%  xlim([0,1.2])
    hold on
       plot(FB(nrB,:),GB(nrB,:),'b-')
       plot(-FB(nrB,:),GB(nrB,:),'b-')
 j0=(nzB-1)/2+1
 s=1:j0;
      plot(0.5*FB(nrB,s),0.5*GB(nrB,s),'r-')
      hold on
      plot(-0.5*FB(nrB,s),0.5*GB(nrB,s),'r-')
 s=j0:nzB   

      plot(0.5*FB(nrB,s),0.5*GB(nrB,s),'b-')
      hold on
      plot(-0.5*FB(nrB,s),0.5*GB(nrB,s),'b-')

      axis('equal')   
      ylim([-7.5,7.5])
      xlim([-7.5,7.5])
   drawnow
   % F(j) = getframe;
    frame=getframe(h),
    writeVideo(mov,frame);
    close(h)
%  ontour(zBinterp 

      case 5

        % Initializing the variables 
yfB= zeros(NVB,NDB,ntB);
xt=0*x0;
for i=1:NVB
%pointer in whole domain i variable j block
l=JM1{i};
variable= x0(l);
variablet= xt(l);
yfB(i,1,:)=variable;
yfB(i,2,:) = ddr0B*variable;
yfB(i,3,:) = ddz0B*variable;
yfB(i,4,:) = ddrr0B*variable;
yfB(i,5,:) = ddzz0B*variable;
yfB(i,6,:) = ddrz0B*variable;
yfB(i,7,:)= variablet;
end
r1B=repmat(r0B', [1 nzB]);
z1B=repmat(z0B, [nrB 1]);
z1B=reshape(z1B,ntB,1);
r1B=reshape(r1B,ntB,1);
r1B=repmat(r0B', [1 nzB]);
z1B=repmat(z0B, [nrB 1]);
z1B=reshape(z1B,ntB,1);
r1B=reshape(r1B,ntB,1);
%curvature and tangential velocites
curvature=zeros(nrB,nzB);
vt=zeros(nrB,nzB);
for l=1:ntB  
xa=reshape(yfB(:,:,l)',NVB*NDB,1); 
  [iss,jss]=ind2sub([nrB,nzB],l);
  [curvature(iss,jss),vt(iss,jss)]= Curvature(z1B(l),r1B(l),xa,pa);
end


Curvaturemax(js)=max(curvature(nrB,2:nzB-1));
Wek(js)=rhow*VTf(js)^2*Radius/(sigma*Curvaturemax(js));


%         hold on
%     name1='R= '+string(Radius*1000)+' mm';
% subplot(4,1,1); hold on; plot(z0B,curvature(nrB,:),'DisplayName', name1); 
% totalpressure=pB(nrB,nzB)-pB(nrB,nzB)+0.5*Re*(VTB(1,1)^2-(wB.^2+uB.^2))+curvature;
% 
% 
% subplot(4,1,2);hold on; plot(z0B,-vt(nrB,:)); 
% matrixBBCvorticity
% subplot(4,1,3);hold on; plot(z0B,vorticityBa(1,:))


end
end


if(lplot==5)
    hold on
x=Radiusf*1000;
[a,I]=sort(x);
y=Curvaturemax;
%plot(x(I),Wek(I),'o-');
plot(x(I),y(I),'o-');
hold on
%plot(Xia,We2);
end

if(lplot==6)
    hold on
x=Radiusf*1000;
[a,I]=sort(x);

plot(x(I),Lbubble(I),'o-');
hold on
%plot(Xia,We2);
end


%F=[Radiusf'*1000*2,VTf'*100];
if(lplot==3)
    hold on
x0=load('resultados/desdethe/P5.txt')
name='experiments with pure water DUINEVELD (1995)'
%plot(x0(:,1),x0(:,2),'x','Displayname',name)

x0=load('resultados/desdethe/P6.txt')
name='Moore (1965)'
%plot(x0(:,1),x0(:,2),'-','Displayname',name)

[a,I]=sort(Radiusf);

x=Radiusf(I)*1000;
y=VTf(I)*100;
%y=Re(I);

plot(x,y,'x')

% x0=linspace(min(x),max(x),60);
% y0=interp1(x,y,x0,'spline');

%plot(x,y,'o-','Displayname','2D simlations')



xlabel('R(mm)') 
ylabel('Velocity(cm/s)')

end
close(mov)


