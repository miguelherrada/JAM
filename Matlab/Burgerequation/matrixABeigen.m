 xt = 0*x0;
 %% temporal derivatives
bp=-1i;

 u=x0(1:ns);
 duds=ds*x0(1:ns);
 dudss=dss*x0(1:ns);
 dudt=x0t(1:ns);
 g=x0(ns+1:2*ns);
 dgds=ds*x0(ns+1:2*ns);
 dgdss=dss*x0(ns+1:2*ns);
 dgdt=x0t(ns+1:2*ns);



%number of variables
nv=2
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
% Getting numerical a and b matrices for the eigenvalue provlem
ablock = cell(nv, nv);
bblock = cell(nv, nv);

for j = 1:nv
   
    for k = 1:nv
        km = (k - 1) * nd + 1;
        kp = k * nd;
        C = squeeze(DFAA(j, km:kp, :));
        B = spdiags(C(1, :)', 0, ns, ns) + spdiags(C(2, :)', 0, ns, ns) * ds + spdiags(C(3, :)', 0, ns, ns) * dss ;
        B1=  spdiags(C(4, :)', 0, ns, ns) * bp;   
        ablock{j, k} = B;
        bblock{j, k} = -B1;
    end
end

% % Mouting BLOCK 
% xa = cell(nv, 1);
% aa = cell(nv, nv);
% bb = cell(nv, nv);
% 
% 
% for i = 1:nv
%     for j = 1:nv
%         aa{i, j} = ablock{i, j};
%         bb{i, j} = bblock{i, j};
%     end
% end

% Evaluate sparse matrices
a = cell2mat(ablock);
b = cell2mat(bblock);