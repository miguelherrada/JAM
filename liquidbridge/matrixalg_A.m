%% Matrix Algebra for Block A


ss=repmat(s, [ny 1]);
yy=repmat(y', [1 ns]);

ddy=kron(speye(ns,ns),dy);
dds=kron(ds,speye(ny,ny));
delta=kron(speye(ny,ny),speye(ns,ns));
ddy2=kron(speye(ns,ns),dy2);
dds2=kron(ds2,speye(ny,ny));
ddsy=kron(ds,dy);


