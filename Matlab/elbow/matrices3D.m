
function [ddr,ddr2,ddz,ddz2,ddx,ddx2,ddrz,ddzx,ddrx]= matrices3D (nr,nz,nx,dr,dr2,dz,dz2,dx,dx2)


  ddr=kron(speye(nz,nz),dr);
  ddr=kron(speye(nx,nx),ddr);
  %2 r derivatives
  ddr2=kron(speye(nz,nz),dr2);
  ddr2=kron(speye(nx,nx),ddr2);
   %1 z derivatives
  ddz=kron(dz,speye(nr,nr));
  ddz=kron(speye(nx,nx),ddz);
  
  %2 z derivatives
  ddz2=kron(dz2,speye(nr,nr));
  ddz2=kron(speye(nx,nx),ddz2);
  %1 thetaderivatives
  ddx=kron(dx,speye(nz,nz));
  ddx=kron(ddx,speye(nr,nr));
  %2 theta derivatives
  ddx2=kron(dx2,speye(nz,nz));
  ddx2=kron(ddx2,speye(nr,nr));
  %mixing derivatives (r-z)
  ddrz=kron(dz,dr); 
  ddrz=kron(speye(nx,nx),ddrz); 
  %mixing derivatives (z-theta)
  ddzx=kron(dx,dz); 
  ddzx=kron(ddzx,speye(nr,nr)); 
  %mixing derivatives (r-theta)
  ddrx=kron(dx,speye(nz,nz)); 
  ddrx=kron(ddrx,dr);