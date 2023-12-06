import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse import diags
# second order central differences


def finitas2(nz, H):
    z = np.zeros(nz)
    dz = H / (nz - 1)
    
    for j in range(nz):
        z[j] = (j) * dz
    
    dx = 0.5 / dz
    dx1 = 1 / (12 * dz)
    dx2 = 1 / dz**2
    
    d = lil_matrix((nz, nz))
    
    d[0, 0] = -3 * dx  # second order
    d[0, 1] = 4 * dx
    d[0, 2] = -dx
    
    for j in range(1, nz - 1):  # forth order
        d[j, j - 1] = -dx
        d[j, j] = 0
        d[j, j + 1] = dx
    
    d[nz - 1, nz - 3] = dx
    d[nz - 1, nz - 2] = -4 * dx
    d[nz - 1, nz - 1] = 3 * dx
    
    dx = 1 / (dz * dz)
    d2 = lil_matrix((nz, nz))
    
    d2[0, 0] = 2 * dx
    d2[0, 1] = -5 * dx
    d2[0, 2] = 4 * dx
    d2[0, 3] = -dx
    
    for j in range(1, nz - 1):  # forth order
        d2[j, j - 1] = dx
        d2[j, j] = -2 * dx
        d2[j, j + 1] = dx
    
    d2[nz - 1, nz - 1] = 2 * dx
    d2[nz - 1, nz - 2] = -5 * dx
    d2[nz - 1, nz - 3] = 4 * dx
    d2[nz - 1, nz - 4] = -dx
    
    return z, d.tocsr(), d2.tocsr()
import numpy as np
from scipy.sparse import lil_matrix, diags

def finitas4(nz, H):
    z = np.zeros(nz)
    dz = H / (nz - 1)
    for j in range(nz):
        z[j] = j * dz

    dx = 0.5 / dz
    dx1 = 1 / (12 * dz)
    dx2 = 1 / (12 * dz ** 2)

    e1 = dx1 * np.ones(nz)
    e2 = dx2 * np.ones(nz)

    d = lil_matrix((nz, nz), dtype=np.float64)
    d2 = lil_matrix((nz, nz), dtype=np.float64)

    d[0, 0] = -25 * dx1
    d[0, 1] = 48 * dx1
    d[0, 2] = -36 * dx1
    d[0, 3] = 16 * dx1
    d[0, 4] = -3 * dx1

    d2[0, 0] = 35 * dx2
    d2[0, 1] = -104 * dx2
    d2[0, 2] = 114 * dx2
    d2[0, 3] = -56 * dx2
    d2[0, 4] = 11 * dx2

    d[1, 0] = -3 * dx1
    d[1, 1] = -10 * dx1
    d[1, 2] = 18 * dx1
    d[1, 3] = -6 * dx1
    d[1, 4] = 1 * dx1

    d2[1, 0] = 11 * dx2
    d2[1, 1] = -20 * dx2
    d2[1, 2] = 6 * dx2
    d2[1, 3] = 4 * dx2
    d2[1, 4] = -1 * dx2

    for j in range(2, nz - 2):
        d[j, j - 2] = dx1
        d[j, j - 1] = -8 * dx1
        d[j, j] = 0
        d[j, j + 1] = 8 * dx1
        d[j, j + 2] = -dx1

        d2[j, j - 2] = -dx2
        d2[j, j - 1] = 16 * dx2
        d2[j, j] = -30 * dx2
        d2[j, j + 1] = 16 * dx2
        d2[j, j + 2] = -dx2

    d[nz - 2, nz - 5] = -1 * dx1
    d[nz - 2, nz - 4] = 6 * dx1
    d[nz - 2, nz - 3] = -18 * dx1
    d[nz - 2, nz - 2] = 10 * dx1
    d[nz - 2, nz - 1] = 3 * dx1

    d[nz - 1, nz - 5] = -1 * dx1
    d[nz - 1, nz - 4] = 6 * dx1
    d[nz - 1, nz - 3] = -18 * dx1
    d[nz - 1, nz - 2] = 10 * dx1
    d[nz - 1, nz - 1] = 3 * dx1

    d2[nz - 2, nz - 5] = -1 * dx2
    d2[nz - 2, nz - 4] = 4 * dx2
    d2[nz - 2, nz - 3] = 6 * dx2
    d2[nz - 2, nz - 2] = -20 * dx2
    d2[nz - 2, nz - 1] = 11 * dx2

    d2[nz - 1, nz - 5] = -1 * dx2
    d2[nz - 1, nz - 4] = 4 * dx2
    d2[nz - 1, nz - 3] = 6 * dx2
    d2[nz - 1, nz - 2] = -20 * dx2
    d2[nz - 1, nz - 1] = 11 * dx2

   

    return z, d.tocsr(), d2.tocsr()
