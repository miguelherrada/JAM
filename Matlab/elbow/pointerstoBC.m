
%pointers to BC

JM1 = cell(NVA,1);

%creating pointers
nbi = 0;
nunkn = ntA; %number of gridpoints in block    
for i=1:NVA
 inic = nbi+(i-1)*nunkn+1;
        ifin = nbi+i*nunkn;    
JM1{i}=(inic:ifin); %columns 
end

%
ndA=zeros(nrA,nzA,nx);
zA=zeros(nrA,nzA,nx);
rA=zeros(nrA,nzA,nx);
for k=1:nx
for j=1:nzA
for i=1:nrA
rA(i,j,k)=r0A(i);
zA(i,j,k)=z0A(j);
xA(i,j,k)=xx(k);



if(j==1) % entrance
ndA(i,j,k)=3; %t
end 

if(j==nzA) % exit
ndA(i,j,k)=1; %t
end 


if(i==nrA) % wall
ndA(i,j,k)=2; %top
end

end
end
end






ndA=reshape(ndA,1,ntA);

rA=reshape(rA,1,ntA);

zA=reshape(zA,1,ntA);

xA=reshape(xA,1,ntA);


