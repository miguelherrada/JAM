import numpy as np
def initial_guess_x0(list_var_A, NVA, ntA):
    """
    Constructs the initial guess solution x0 by reshaping variables.
    Translated from MATLAB code.
    """
    x01_list = [np.reshape(list_var_A[i], ntA,order='C') for i in range(NVA)]
    #x01 = np.vstack(x01_list)
    x01 = np.concatenate(x01_list)
    return x01
def inverse_initial_guess_x0(x0, list_var_A, nrA, nzA, nxA):
    new_list = []
    start = 0
    for var in list_var_A:
        n_el = np.prod(var.shape)
        # Se reorganiza el segmento con orden Fortran
        new_var = np.reshape(x0[start:start+n_el], (nrA, nzA, nxA), order='F')
        new_list.append(new_var)
        start += n_el
    return new_list


def computingpointers(nrA, nzA, nxA):
    """
    Calculate pointers for boundary and internal points matching computingpointers1 output
    """
    # Left boundary (x = 0)
    l_l = [np.ravel_multi_index((i, 0, k), (nrA, nzA, nxA)) for i in range(nrA) for k in range(nxA)]
    
    # Right boundary (x = L)
    l_r = [np.ravel_multi_index((i, nzA - 1, k), (nrA, nzA, nxA)) for i in range(nrA) for k in range(nxA)]
    
   
    
    # Top boundary (y = 1)
    l_t = [np.ravel_multi_index((nrA - 1, j, k), (nrA, nzA, nxA)) for j in range(nzA) for k in range(nxA)]
    
    # Internal points - all points not in boundaries
    mask = np.ones((nrA, nzA, nxA), dtype=bool)
    mask[:, 0, :] = False    # Remove left boundary
    mask[:, -1, :] = False   # Remove right boundary
    mask[-1, :, :] = False   # Remove top boundary
    l_b = np.ravel_multi_index(np.where(mask), (nrA, nzA, nxA))
    
    return l_l, l_r, l_t, l_b

def computingpointersF(nrA, nzA, nxA):
    """
    Calculate pointers for boundary and internal points using Fortran (column-major) order.
    """
    # Left boundary (x = 0)
    l_l = [np.ravel_multi_index((i, 0, k), (nrA, nzA, nxA), order='F')
           for i in range(nrA) for k in range(nxA)]
    
    # Right boundary (x = L)
    l_r = [np.ravel_multi_index((i, nzA - 1, k), (nrA, nzA, nxA), order='F')
           for i in range(nrA) for k in range(nxA)]
    
    # Top boundary (y = 1)
    l_t = [np.ravel_multi_index((nrA - 1, j, k), (nrA, nzA, nxA), order='F')
           for j in range(nzA) for k in range(nxA)]
    
    # Internal points - all points not in boundaries
    mask = np.ones((nrA, nzA, nxA), dtype=bool)
    mask[:, 0, :] = False    # Remove left boundary
    mask[:, -1, :] = False   # Remove right boundary
    mask[-1, :, :] = False   # Remove top boundary
    l_b = np.ravel_multi_index(np.where(mask), (nrA, nzA, nxA), order='F')
    
    return l_l, l_r, l_t, l_b

def fix_result(res, target_length):
    # Si es escalar, crear un arreglo lleno de ese valor.
    if np.isscalar(res):
        return np.full(target_length, res)
    # Asegurarse de que sea un arreglo unidimensional
    arr = np.atleast_1d(res).flatten()
    if arr.size == 1:
        return np.full(target_length, arr.item())
    elif arr.size == target_length:
        return arr
    else:
        raise ValueError(f"Size mismatch: se esperaba {target_length} elementos, pero se obtuvo {arr.size}.")
