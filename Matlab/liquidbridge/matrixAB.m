%% Time derivative 
 %second order backwar differences
 dta=1/(2*dt);
 bm= -4*dta; 
 bmm= 1*dta; 
 bp=3*dta;  %[STEP 6]
 %temporal deriviatives
 x0t = bp*x0 + bm*x0m + bmm*x0mm;
 lw=1:N; %pointer
 w=x0(lw);
 dwds=dds*x0(lw);
 dwdss=dds2*x0(lw);
 dwdy=ddy*x0(lw);
 dwdyy=ddy2*x0(lw);
 dwdys=ddsy*x0(lw);
 dwdt0=x0t(lw);
 lu=N+1:2*N; %pointer
 u=x0(lu);
 duds=dds*x0(lu);
 dudss=dds2*x0(lu);
 dudy=ddy*x0(lu);
 dudyy=ddy2*x0(lu);
 dudys=ddsy*x0(lu);
 dudt0=x0t(lu);
 lp=2*N+1:3*N; %pointer
 p=x0(lp);
 dpds=dds*x0(lp);
 dpdss=dds2*x0(lp);
 dpdy=ddy*x0(lp);
 dpdyy=ddy2*x0(lp);
 dpdys=ddsy*x0(lu);
 dpdt0=x0t(lu);
 lf=3*N+1:4*N;%pointer
 f=x0(lf);
 dfds=dds*x0(lf);
 dfdss=dds2*x0(lf);
 dfdy=ddy*x0(lf);
 dfdyy=ddy2*x0(lf);
 dfdys=ddsy*x0(lf);
 dfdt0=x0t(lf);

%number of variables
nv=4;
%number of symbolic derivatives
nd=7;

FAA=zeros(nv,N);
DFAA=zeros(nv,nv*nd,N);
xs=zeros(1,nv*nd);
%getting analitical espresions[STEP8]
for l=1:N
[i,j]=ind2sub([ny,ns],l);

xs(1:nd)       =[w(l),dwds(l), dwdss(l), dwdy(l), dwdyy(l), dwdys(l), dwdt0(l)]; %step 2; %variables
xs(nd+1:2*nd)  =[u(l),duds(l), dudss(l), dudy(l), dudyy(l), dudys(l), dudt0(l)];
xs(2*nd+1:3*nd)=[p(l),dpds(l), dpdss(l), dpdy(l), dpdyy(l), dpdys(l), dpdt0(l)];
xs(3*nd+1:4*nd)=[f(l),dfds(l), dfdss(l), dfdy(l), dfdyy(l), dfdys(l), dfdt0(l)];
%Bulk
if (( i>1) && (i<ny))&& (( j>1) && (j<ns))
[FAA(:,l),DFAA(:,:,l)]=equationFAAb(s(j),y(i),xs,pa);
end
%bottom
 if ( i==1 )
[FAA(:,l),DFAA(:,:,l)]=equationFAAd(s(j),y(i),xs,pa);     
 end
%top
if ( i==ny )
[FAA(:,l),DFAA(:,:,l)]=equationFAAt(s(j),y(i),xs,pa);  
end


%Left
 if ( j==1 )
[FAA(:,l),DFAA(:,:,l)]=equationFAAl(s(j),y(i),xs,pa);    
 end
%Right
if ( j==ns )
[FAA(:,l),DFAA(:,:,l)]=equationFAAr(s(j),y(i),xs,pa);    
end
   


end
% Getting numerical Jacobian (step 9)
ablock = cell(nv, nv);
bablock = cell(nv, 1);

for j = 1:nv
    C1 = FAA(j, :);
    bablock{j} = -C1';
 
    for k = 1:nv
        km = (k - 1) * nd + 1;
        kp = k * nd;
        C = squeeze(DFAA(j, km:kp, :));
        B = spdiags(C(1, :)', 0, N, N) + spdiags(C(2, :)', 0, N, N) * dds + spdiags(C(3, :)', 0, N, N) * dds2 +  + spdiags(C(4, :)', 0, N, N) * ddy + spdiags(C(5, :)', 0, N, N) * ddy2+  spdiags(C(6, :)', 0, N, N) * ddsy  +spdiags(C(7, :)', 0, N, N) * bp;
    
        ablock{j, k} = sparse(B);
       
    end
end

% Mouting BLOCK 
ba = cell(nv, 1);
aa = cell(nv, nv);

for i = 1:nv
    ba{i} = bablock{i};
    for j = 1:nv
        aa{i, j} = ablock{i, j};
    end
end

% Evaluate sparse matrices
a = cell2mat(aa);
b = cell2mat(ba);
