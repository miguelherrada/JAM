%% Fundamental subroutine of the 'JAM' method, which allows the equations to be expressed in Cartesian coordinates by non-singular mappings and generates the analytical Jacobians required by the Newton method.

%% Expressing the equations in 3D Euclidean space
%The variables in the 3D box [z0,r0,t0] and time  
syms r0 z0 t0 time 
% Define symbolic variables for all blocks and variables
for kk = 1:Nblock
    for j = 1:NVAR(kk)
        % Declare symbolic variable
        variable_name = sprintf('%s%d(r0,z0, time)', list_var{kk}{j}, kk);
        eval(['syms ' variable_name ';']);
        % Assign variable to its name
        assume (r0,"real")
        assume (z0,"real")
        assume (time,"real")
        assume (t0,"real")
        eval(sprintf('%s%d= %s;', list_var{kk}{j}, kk, variable_name));
    end
end

%Mometum equations 
EqMomentumy=cell(1,Nblock); %radial
EqMomentumx=cell(1,Nblock); %axial
%Divergence
Eqdiver=cell(1,Nblock); %axial
%Mappings
EqMapx=cell(1,Nblock); %axial
EqMapy=cell(1,Nblock); %axial


for kk = 1:Nblock 
%Auxiliar variable "0"
for j=1:NVAR(kk)
eval(sprintf('%s0= %s%d;', list_var{kk}{j}, list_var{kk}{j},kk));
end

% Maping 3D Euclidean 
X=[F0,G0,t0]; %y,x,z: The problem is 2D t0 (z) is not used
xs=[r0,z0,t0];% the mapped box
J=jacobian(X,xs);
Jinv=simplify(inv(J));

%Getting vector of the space bc{i}
Xc=[r0,z0,t0]; 
%1:r0-y
bc{1}=diff(Xc,r0);
%a:z0-x
bc{2}=diff(Xc,z0);
%c:t0-z
bc{3}=diff(Xc,t0);
% Selecting vector base ( bc)
bbase=bc; 

% Temporal mapping
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
% 2D Cartesian velocity (base bbase)
U3D=vy0*bbase{1}+vx0*bbase{2}+0*bbase{3};

% Derivatives in the 3D Euclidean space 
       for j = 1:NVAR(kk)
        var_name = sprintf('%s%d', list_var{kk}{j},kk);
        % First derivatives in X
        eval(sprintf('d%sdy = Jinv(1, 1) * diff(%s, xs(1)) + Jinv(2, 1) * diff(%s, xs(2))+ Jinv(3, 1) * diff(%s, xs(3));', var_name, var_name, var_name,var_name));
        eval(sprintf('d%sdx = Jinv(1, 2) * diff(%s, xs(1)) + Jinv(2, 2) * diff(%s, xs(2))+ Jinv(3, 2) * diff(%s, xs(3));', var_name, var_name, var_name,var_name));
        eval(sprintf('d%sdz = Jinv(1, 3) * diff(%s, xs(1)) + Jinv(2, 3) * diff(%s, xs(2))+ Jinv(3, 3) * diff(%s, xs(3));', var_name, var_name, var_name,var_name));
        % Second derivatives in X
        eval(sprintf('d%sdyy = Jinv(1, 1) * diff(d%sdy, xs(1)) + Jinv(2, 1) * diff(d%sdy, xs(2))+ Jinv(3, 1) * diff(d%sdy, xs(3));', var_name, var_name, var_name,var_name));
        eval(sprintf('d%sdxx = Jinv(1, 2) * diff(d%sdx, xs(1)) + Jinv(2, 2) * diff(d%sdx, xs(2))+ Jinv(3, 2) * diff(d%sdx, xs(3));', var_name, var_name, var_name,var_name));
        eval(sprintf('d%sdzz = Jinv(1, 3) * diff(d%sdz, xs(1)) + Jinv(2, 3) * diff(d%sdz, xs(2))+ Jinv(3, 3) * diff(d%sdz, xs(3));', var_name, var_name, var_name,var_name));
        % Derivatives in the mapped space
        eval(sprintf('d%sdr0 = diff(%s, r0);', var_name, var_name));
        eval(sprintf('d%sdz0 = diff(%s, z0);', var_name, var_name));
        eval(sprintf('d%sdrr0 = diff(%s, r0, r0);', var_name, var_name));
        eval(sprintf('d%sdzz0 = diff(%s, z0, z0);', var_name, var_name));
        eval(sprintf('d%sdrz0 = diff(%s, r0, z0);', var_name, var_name));
        eval(sprintf('d%sdtime0 = diff(%s, time);', var_name, var_name));
         %temporal derivatives
        eval(sprintf('d%sdtime3D = diff(%s,time)+diff(%s,r0)*dr0dtime0+diff(%s,z0)*dz0dtime0+diff(%s,t0)*dt0dtime0;', var_name, var_name, var_name,var_name,var_name));
    end

%Additional derivatives in X ...
%x and y derivatives of all
list_more={'U3D','F0','G0','p0'};
Nmore=length(list_more)
   for j = 1:Nmore
        var_name = sprintf('%s', list_more{j});
        %First derivatives in X
        eval(sprintf('d%sdy = Jinv(1, 1) * diff(%s, xs(1)) + Jinv(2, 1) * diff(%s, xs(2))+ Jinv(3, 1) * diff(%s, xs(3));', var_name, var_name, var_name,var_name));
        eval(sprintf('d%sdx = Jinv(1, 2) * diff(%s, xs(1)) + Jinv(2, 2) * diff(%s, xs(2))+ Jinv(3, 2) * diff(%s, xs(3));', var_name, var_name, var_name,var_name));
        eval(sprintf('d%sdz = Jinv(1, 3) * diff(%s, xs(1)) + Jinv(2, 3) * diff(%s, xs(2))+ Jinv(3, 3) * diff(%s, xs(3));', var_name, var_name, var_name,var_name));
        %Second derivatives in X
        eval(sprintf('d%sdyy = Jinv(1, 1) * diff(d%sdy, xs(1)) + Jinv(2, 1) * diff(d%sdy, xs(2))+ Jinv(3, 1) * diff(d%sdy, xs(3));', var_name, var_name, var_name,var_name));
        eval(sprintf('d%sdxx = Jinv(1, 2) * diff(d%sdx, xs(1)) + Jinv(2, 2) * diff(d%sdx, xs(2))+ Jinv(3, 2) * diff(d%sdx, xs(3));', var_name, var_name, var_name,var_name));
        eval(sprintf('d%sdzz = Jinv(1, 3) * diff(d%sdz, xs(1)) + Jinv(2, 3) * diff(d%sdz, xs(2))+ Jinv(3, 3) * diff(d%sdz, xs(3));', var_name, var_name, var_name,var_name));
        %Derivatives in the mapped space
        eval(sprintf('d%sdr0 = diff(%s, r0);', var_name, var_name));
        eval(sprintf('d%sdz0 = diff(%s, z0);', var_name, var_name));
        eval(sprintf('d%sdrr0 = diff(%s, r0, r0);', var_name, var_name));
        eval(sprintf('d%sdzz0 = diff(%s, z0, z0);', var_name, var_name));
        eval(sprintf('d%sdrz0 = diff(%s, r0, z0);', var_name, var_name));
        eval(sprintf('d%sdtime0 = diff(%s, time);', var_name, var_name));
        %Temporal derivatives
        eval(sprintf('d%sdtime3D = diff(%s,time)+diff(%s,r0)*dr0dtime0+diff(%s,z0)*dz0dtime0+diff(%s,t0)*dt0dtime0;', var_name, var_name, var_name,var_name,var_name));
   end

   

% Velocites operators
%gradient U in cartesians
gradientU3D=[dU3Ddy.',dU3Ddx.',dU3Ddz.'];
%Tensor D
DD=(gradientU3D+gradientU3D.')/2;
%tensor W
WW=(gradientU3D-gradientU3D.')/2;
%Ugradient
UgradientU3D=vx1*zeros(1,3);
for i=1:3
UgradientU3D(i)=simplify(U3D(1)*dU3Ddy(i)+U3D(2)*dU3Ddx(i)+U3D(3)*dU3Ddz(i));
end
%Newtonian stress 
TAUN3D=2*DD;
%total stress
TAUT3D=TAUN3D;
%Laplace operatoir
laplaceU3D=simplify(dU3Ddxx+dU3Ddyy+dU3Ddzz);
%gradient p
gradientp3D=[dp0dy,dp0dx,dp0dz];
%diverU
diverU3D=simplify(dU3Ddx(2)+dU3Ddy(1)+dU3Ddz(3));



% Momentum and continuity equations in 3D Euclidean space
EQMM3D=-(dU3Ddtime3D+UgradientU3D)-gradientp3D+laplaceU3D/Re;
%Proyecting in base bc+++++++++++++++++++++++++++++++++++++++++++
EQMM=0*EQMM3D;
for i=1:3
EQMM(i)=(EQMM3D*(bc{i}.'));
end
%2D momentum
EqMomentumy{kk}=EQMM(1); %axial momentum
EqMomentumx{kk}=EQMM(2); %radial mometum

 Eqdiver{kk}=diverU3D; %continuity

 % Mapping :Dimakopoulos, Tsamopoulos - 2003  quasi-elliptic transformation 
 g11=dG0dz0^2+dF0dz0^2;
 g22=dG0dr0^2+dF0dr0^2;
 g12=dG0dr0*dG0dz0+dF0dr0*dF0dz0;
 J=dG0dr0*dF0dz0-dG0dz0*dF0dr0;
 eps1=0.3;
 D1=eps1*((dF0dz0^2+dG0dz0^2)/(dF0dr0^2+dG0dr0^2))^0.5+(1-eps1);
 dD1dr0=diff(D1,r0);
 dD1dz0=diff(D1,z0);
 Q1=-(dD1dr0*dF0dz0-dD1dz0*dF0dr0)*J/D1;

 EqMapy{kk}=g22*dF0dzz0+g11*dF0drr0-2*g12*dF0drz0-Q1;
 EqMapx{kk}=g22*dG0dzz0+g11*dG0drr0-2*g12*dG0drz0;

end



%%  Real yF matrix to replace simbolic matrix yFsym
yFF=cell(1,Nblock); %real matrix
yFFsym=cell(1,Nblock);%symbolic matrix

for kk=1:Nblock
yFF{kk} = sym('yF', [NVAR(kk), NDA], 'real'); % Preallocating symbolic array

% Preallocate yFsym for symbolic derivatives
yFFsym{kk} = sym('yFsym', [NVAR(kk), NDA], 'real');
% Getting symbolic matrices with all
for l = 1:NVAR(kk)
  

    for k = 1:NDA  % Derivatives
        % Assign symbolic variable to yF
        yFF{kk}(l, k) = sym(['yb', num2str(kk), 'v', num2str(l), 'd', num2str(k)], 'real');

        % For the base case (k == 1)
        if k == 1
            name = sprintf('yFFsym{kk}(l, 1) = %s%d;', list_var{kk}{l},kk);
        else
            % For derivatives, use list_dersymbolic to construct the diff command
            name = sprintf('yFFsym{kk}(l, %d) = diff(%s%d%s;', k, list_var{kk}{l},kk, list_dersymbolic{k});
        end

       % Execute the command to assign yFsym
        eval(name);
    end
end
end
%full matrices
yFsym = vertcat(yFFsym{:});
yF = vertcat(yFF{:});

%% Define equations
for k = 1:Nequations
    eq_name = sprintf('FF%d', k);
    eval(sprintf('%s = sym(zeros(NVA, 1));', eq_name));
lv=0;
for j = 1:Nblock
    for i = 1:NVAR(j)
    variable_name = list_var{j}{i}; % Get the variable name
lv=lv+1;
eval(sprintf('%s(%d)=%s%d;', eq_name,lv,variable_name,j));
    end 
end
end
%% Populate equations
%pointers
jp0=0;
jp1=NVAR(1);
% Bulk 1
FF1(1)=EqMomentumy{1};
FF1(2)=EqMomentumx{1};
FF1(3)=Eqdiver{1};
FF1(4)=EqMapy{1};
FF1(5)=EqMapx{1};
% Bulk 2
FF2(1+jp1)=EqMomentumy{2};
FF2(2+jp1)=EqMomentumx{2};
FF2(3+jp1)=Eqdiver{2};
FF2(4+jp1)=EqMapy{2};
FF2(5+jp1)=EqMapx{2};
% left: entrance
FF3(1+jp0)=vy1;
FF3(2+jp0)=vx1-1;
FF3(3+jp0)=dp1dx;
FF3(4+jp0)=F1-Y1;
FF3(5+jp0)=G1-X1;
FF3(1+jp1)=vy2;
FF3(2+jp1)=vx2-1;
FF3(3+jp1)=dp2dx;
FF3(4+jp1)=F2-Y2;
FF3(5+jp1)=G2-X2;
% right: exit
FF4(1+jp0)=vy1;
FF4(2+jp0)=dvx1dx;
FF4(3+jp0)=p1;
FF4(4+jp0)=F1-Y1;
FF4(5+jp0)=G1-X1;
FF4(1+jp1)=vy2;
FF4(2+jp1)=dvx2dx;
FF4(3+jp1)=p2;
FF4(4+jp1)=F2-Y2;
FF4(5+jp1)=G2-X2;
% Midle plane 
FF5(1)=vy1-vy2;
FF5(2)=vx1-vx2;
FF5(3)=p1-p2;
FF5(4)=F1-Y1;
FF5(5)=G1-X1;
FF5(1+jp1)=dvy1dy-dvy2dy;
FF5(2+jp1)=dvx1dy-dvx2dy;
FF5(3+jp1)=dp1dy-dp2dy;
FF5(4+jp1)=F2-Y2;
FF5(5+jp1)=G2-X2;
% 6 bottom wall
FF6(1+jp0)=vy1;
FF6(2+jp0)=dvx1dy;
FF6(3+jp0)=dp1dy;
FF6(4+jp0)=F1-Y1;
FF6(5+jp0)=dG1dzz0;
% 7 top wall
FF7(1+jp1)=vy2;
FF7(2+jp1)=dvx2dy;
FF7(3+jp1)=dp2dy;
FF7(4+jp1)=F2-Y2;
FF7(5+jp1)=dG2dzz0;
% Cylinder 
FF8(1)=vy1;
FF8(2)=vx1;
FF8(3)=dp1dr0;
FF8(4)=F1-Y1;
FF8(5)=G1-X1;
FF8(1+jp1)=vy2;
FF8(2+jp1)=vx2;
FF8(3+jp1)=dp2dr0;
FF8(4+jp1)=F2-Y2;
FF8(5+jp1)=G2-X2;



%% Replace symbolic variables with real ones in equations
for le = 1:Nequations
    eq_name = sprintf('FF%d', le);
    eval(sprintf('Eq = %s;', eq_name));

     k = 1:NVA;
     i = NDA:-1:1;
     Eq = subs(Eq, yFsym(k, i), yF(k, i));
    %2D putting t0=0
     Eq = subs(Eq, t0, 0);
    eval(sprintf('%s = Eq;', eq_name));
 %   kn=subs(kn, t0, 0);
end

%% Generate MATLAB functions for equations and their Jacobians
x = reshape(yF', NVA * NDA, 1); % Variables and derivatives

for i = 1:Nequations
    eq_name = sprintf('FF%d', i);
    dFF_name = sprintf('dFF%d', i);
       eval(sprintf('%s = sym(zeros(NVA, NVA * NDA));', dFF_name));

     k = 1:NVA;
        eval(sprintf('%s(k, :) = jacobian(%s(k), x);', dFF_name, eq_name));
    

    % Save equations and Jacobians as MATLAB functions
    file_name = sprintf('%sequationF%d', path_jacobian, i);
    eval(sprintf('matlabFunction(%s, %s, ''File'', ''%s'', ''Vars'', {z0, r0, x, pa});', eq_name, dFF_name, file_name));
end


