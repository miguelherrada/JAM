import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse import diags
from scipy.sparse import csr_matrix
# second order central differences


def finites2th(nz, H):
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


def finites4th(nz, H):

    z = np.zeros(nz)
    dz = H / (nz - 1)
    for j in range(nz):
        z[j] = j * dz

    dx = 0.5 / dz
    dx1 = 1 / (12 * dz)
    dx2 = 1 / (12 * dz**2)

    d = np.zeros((nz, nz))
    d2 = np.zeros((nz, nz))

    d[0, :5] = np.array([-25, 48, -36, 16, -3]) * dx1
    d2[0, :5] = np.array([35, -104, 114, -56, 11]) * dx2

    d[1, :5] = np.array([-3, -10, 18, -6, 1]) * dx1
    d2[1, :5] = np.array([11, -20, 6, 4, -1]) * dx2

    for j in range(2, nz - 2):
        d[j, j-2:j+3] = np.array([1, -8, 0, 8, -1]) * dx1
        d2[j, j-2:j+3] = np.array([-1, 16, -30, 16, -1]) * dx2

    d[-2, -5:] = np.array([-1, 6, -18, 10, 3]) * dx1
    d2[-2, -5:] = np.array([-1, 4, 6, -20, 11]) * dx2

    d[-1, -5:] = np.array([3, -16, 36, -48, 25]) * dx1
    d2[-1, -5:] = np.array([11, -56, 114, -104, 35]) * dx2
   

    d = csr_matrix(d)
    d2 = csr_matrix(d2)

    return z, d, d2








def Chevigood(nr, Rout, rin):
    nr1 = nr + 1
    pi = np.arccos(-1.0)

    c = np.zeros(nr1)
    c[1:nr] = 1.0
    c[0] = 2.0
    c[nr] = 2.0

    x = np.zeros(nr1)
    xc = np.zeros(nr1)
    dx = np.zeros(nr1)
    dx2 = np.zeros(nr1)
    d = np.zeros((nr1, nr1))
    d2 = np.zeros((nr1, nr1))

    x[0] = 1.0
    xc[0] = rin
    dx[0] = -2.0 / Rout
    dx2[0] = 0.0
    for j in range(1, nr1):
        x[j] = np.cos((j * pi) / nr)
        xc[j] = rin + 0.5 * (1.0 - x[j]) * Rout
        dx[j] = -2.0 / Rout
        dx2[j] = 0

    for l in range(nr1):
        for j in range(nr1):
            if j == l and 1 <= j < nr:
                d[l, j] = -x[l] / (2.0 * (1 - (x[j] ** 2.0)))
            elif j != l:
                d[l, j] = c[l] * (-1.0) ** (l + j) / (c[j] * (x[l] - x[j]))

    d[0, 0] = (2.0 * (nr ** 2.0) + 1.0) / 6.0
    d[nr, nr] = -d[0, 0]

    for l in range(nr1):
        for j in range(nr1):
            d[l, j] = d[l, j] * dx[l]

    d2 = np.dot(d, d)

    return xc, d, d2
