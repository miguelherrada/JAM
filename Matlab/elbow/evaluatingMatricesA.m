function [FAA,DFAA] = evaluatingMatricesA(nrA,nzA,zA,rA,tA,yfA,NVA,NDA,ndA,pa,nx) 
ntA = nrA*nzA*nx;

FAA = zeros(NVA,ntA);
DFAA =zeros(NVA,NVA*NDA,ntA);

for l=1:ntA
xa=reshape(yfA(:,:,l)',NVA*NDA,1);
if(ndA(l)==0) %bulk
[FAA(:,l),DFAA(:,:,l)]=equationFAA(zA(l),rA(l),tA(l),xa,pa);
end
if(ndA(l)==1) %outflow
[FAA(:,l),DFAA(:,:,l)]=equationFAAr(zA(l),rA(l),tA(l),xa,pa);
end
if(ndA(l)==2)%wall
[FAA(:,l),DFAA(:,:,l)]=equationFAAt(zA(l),rA(l),tA(l),xa,pa);
end

if(ndA(l)==3)%entrances
[FAA(:,l),DFAA(:,:,l)]=equationFAAl(zA(l),rA(l),tA(l),xa,pa);
end



end





end
