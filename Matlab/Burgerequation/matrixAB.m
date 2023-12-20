%% Time derivative 
 %second order backwar differencesd
 dta=1/(2*dt);
 bm= -4*dta; 
 bmm= 1*dta; 
 bp=3*dta;  %[STEP 6]
 %temporal deriviatives
 x0t = bp*x0 + bm*x0m + bmm*x0mm;
 u=x0(1:ns);
 duds=ds*x0(1:ns);
 dudss=dss*x0(1:ns);
 dudt=x0t(1:ns);
 g=x0(ns+1:2*ns);
 dgds=ds*x0(ns+1:2*ns);
 dgdss=dss*x0(ns+1:2*ns);
 dgdt=x0t(ns+1:2*ns);

%number of variables
nv=2;
%number of symbolic derivatives
nd=4;

FAA=zeros(nv,ns);
DFAA=zeros(nv,nv*nd,ns);
%getting analitical expresions[STEP8]
for i=1:ns
%[is,js]=ind2sub([nrA,nzA],l);
xs=[u(i),duds(i),dudss(i),dudt(i),g(i),dgds(i),dgdss(i),dgdt(i)];  %symbolic derivatives evaluation..
%Bulk
if (( i>1) && (i<ns))
[FAA(:,i),DFAA(:,:,i)]=equationFAAb(s(i),xs,pa);
end
%Left
 if ( i==1 )
[FAA(:,i),DFAA(:,:,i)]=equationFAAl(s(i),xs,pa);     
 end
%Right
if ( i==ns )
[FAA(:,i),DFAA(:,:,i)]=equationFAAr(s(i),xs,pa);    
end
   
end
% Getting numerical Jacobian (step 9)
ablock = cell(nv, nv);
xablock = cell(nv, 1);

for j = 1:nv
    C1 = FAA(j, :);
    xablock{j} = -C1';
 
    for k = 1:nv
        km = (k - 1) * nd + 1;
        kp = k * nd;
        C = squeeze(DFAA(j, km:kp, :));
        B = spdiags(C(1, :)', 0, ns, ns) + spdiags(C(2, :)', 0, ns, ns) * ds + spdiags(C(3, :)', 0, ns, ns) * dss + spdiags(C(4, :)', 0, ns, ns) * bp;
    
        ablock{j, k} = sparse(B);
       
    end
end

% Mouting BLOCK 
xa = cell(nv, 1);
aa = cell(nv, nv);

for i = 1:nv
    xa{i} = xablock{i};
    for j = 1:nv
        aa{i, j} = ablock{i, j};
    end
end

% Evaluate sparse matrices
a = cell2mat(aa);
b = cell2mat(xa);