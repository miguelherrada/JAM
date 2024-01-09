


import numpy as np
from scipy.sparse import kron, eye
#to compute the initial shape of the liquid bridge
def Fshape(x, nz, zz, dz, dz2, Bo, V, Lambda):
    # Extracting variables from x
    f = x[:nz]
    pref = x[nz]

    z = zz
    fz = dz@f
    fzz = dz2@f

    # Interface curvature
    k = (f * fzz - 1 - fz**2) / (f * (1 + fz**2)**1.5)
    c=pref + Bo * z[1:nz-1] + k[1:nz-1]
    # Laplace-Young equation + Volume conservation
    Fshape = np.concatenate([c, [f[0] - 1], [f[nz-1] - 1]])


    # Volume conservation
    V0 = 0.0
    for i in range(1, nz):
        V0 += 0.5 * (zz[i] - zz[i - 1]) * (f[i]**2 + f[i - 1]**2)
    V0 = V0 / (2 * Lambda)

    Fshape = np.concatenate([Fshape, [V - V0]])

    return Fshape



#to compute the full saptial positon [ss yy] and full colocation matrices
def matrix_A(s,ds,dss,ns,y,dy,dyy,ny):
    # repmat(s, [ny, 1])
 ss = np.tile(s, (ny, 1))
 # repmat(y', [1, ns])
 yy = np.tile(y[:, np.newaxis], (1, ns))
 # kron(speye(ns, ns), dy)
 ddy = kron(dy,eye(ns))
 # kron(ds, speye(ny, ny))
 dds = kron(eye(ny),ds)
 # kron(speye(ns, ns), dy2)
 ddyy = kron(dyy,eye(ns))
 # kron(ds2, speye(ny, ny))
 ddss = kron(eye(ny),dss)
 # kron(ds, dy)
 ddsy = kron(dy, ds)
 return ss, yy, dds, ddss, ddy,ddyy,ddsy

#to compute the pointers to the BC
def computingpointers(ny,ns):
  #pointers to Boundaries
 #left
 l_l = []
 for i in range(ny):
    l_l.append(np.ravel_multi_index((i, 0), (ny, ns)))
 #right
 l_r = []
 for i in range(ny):
    l_r.append(np.ravel_multi_index((i, ns-1), (ny, ns)))
#down
 l_d = []
 for i in range(ns):
    l_d.append(np.ravel_multi_index((0, i), (ny, ns)))
#top
 l_t = []
 for i in range(ns):
    l_t.append(np.ravel_multi_index((ny-1,i), (ny, ns)))    
 # bulk
 # Crear una matriz de índices
 l_b = np.array([], dtype=int)
 for i in range(1, ny-1):
    for j in range(1, ns-1):
        l_b = np.append(l_b, np.ravel_multi_index((i, j), (ny, ns), order='F'))

 l_b = []
 for i in range(2, ny):
    for j in range(2, ns):
        l_b.append((i-1)*ns + (j-1))

 # Convert to a NumPy array if needed
 l_b = np.array(l_b)

 #print(f"Array de índices:\n{l_b}")


 return l_l, l_r, l_d, l_t, l_b






