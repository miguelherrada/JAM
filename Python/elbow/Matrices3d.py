import scipy.sparse as sp

def matrices3D(nr, nz, nx, dr, dr2, dz, dz2, dx, dx2):
    """
    Computes 3D differentiation sparse matrices using Kronecker products.
    Translated from MATLAB matrices3D.m
    """
    ddr = sp.kron(sp.eye(nx, format='csr'), sp.kron(sp.eye(nz, format='csr'), dr, format='csr'), format='csr')
    ddr2 = sp.kron(sp.eye(nx, format='csr'), sp.kron(sp.eye(nz, format='csr'), dr2, format='csr'), format='csr')
    
    ddz = sp.kron(sp.eye(nx, format='csr'), sp.kron(dz, sp.eye(nr, format='csr'), format='csr'), format='csr')
    ddz2 = sp.kron(sp.eye(nx, format='csr'), sp.kron(dz2, sp.eye(nr, format='csr'), format='csr'), format='csr')
    
    ddx = sp.kron(sp.kron(dx, sp.eye(nz, format='csr'), format='csr'), sp.eye(nr, format='csr'), format='csr')
    ddx2 = sp.kron(sp.kron(dx2, sp.eye(nz, format='csr'), format='csr'), sp.eye(nr, format='csr'), format='csr')
    
   # ddrz = sp.kron(sp.eye(nx, format='csr'), sp.kron(dz, dr, format='csr'), format='csr')
   # ddzx = sp.kron(sp.kron(dx, dz, format='csr'), sp.eye(nr, format='csr'), format='csr')
   # ddrx = sp.kron(sp.kron(dx, sp.eye(nz, format='csr'), format='csr'), dr, format='csr')
    ddrz=ddr*ddz
    ddzx=ddz*ddx
    ddrx=ddr*ddx
    return ddr, ddr2, ddz, ddz2, ddx, ddx2, ddrz, ddzx, ddrx