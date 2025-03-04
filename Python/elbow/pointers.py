import numpy as np

def pointers_to_BC(NVA, ntA, nrA, nzA, nx, r0A, z0A, xx):
    """
    Generates pointers for boundary conditions and reshapes coordinate arrays.
    Translated from MATLAB pointersToBC.m
    """
    # Creating pointers
    JM1 = [None] * NVA
    nbi = 0
    nunkn = ntA  # Number of grid points in block
    
    for i in range(NVA):
        inic = nbi + i * nunkn
        ifin = nbi + (i + 1) * nunkn
        JM1[i] = np.arange(inic, ifin)
    
    # Initialize arrays
    ndA = np.zeros((nrA, nzA, nx), dtype=int)
    rA = np.zeros((nrA, nzA, nx))
    zA = np.zeros((nrA, nzA, nx))               
    xA = np.zeros((nrA, nzA, nx))
    
    # Assign values and boundary conditions
    for k in range(nx):
        for j in range(nzA):
            for i in range(nrA):
               
               # rA[i,j,k]=r0A[i]    # r
               # zA[i,j,k]=z0A[j]    # z                     
               # xA[i,j,k]=xx[k]
                if j == 0:  # Entrance
                    ndA[i, j, k] = 3  # t
                if j == nzA - 1:  # Exit
                    ndA[i, j, k] = 1  # t
                if i == nrA - 1:  # Wall
                    ndA[i, j, k] = 2  # top
    
    # Reshape arrays
    ndA = ndA.reshape(1, ntA,order='F')
    rA = np.tile(r0A[:, np.newaxis, np.newaxis], (1, nzA, nx))
    zA = np.tile(z0A[np.newaxis, :, np.newaxis], (nrA, 1, nx))
    xA = np.tile(xx[np.newaxis, np.newaxis, :], (nrA, nzA, 1))
    
    zA = zA.reshape(ntA, order='F')  # Antes: order='F'
    rA = rA.reshape(ntA, order='F')
    xA = xA.reshape(ntA, order='F')
  
   # zA=np.squeeze(zA)
   # rA=np.squeeze(rA)
   # xA=np.squeeze(xA)
    
    return JM1, ndA, rA, zA, xA