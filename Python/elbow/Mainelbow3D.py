import pickle
import os
import sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse import spdiags, csc_matrix, hstack, vstack
import matplotlib.pyplot as plt
import sympy as sym
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(os.path.join(parent_dir, 'utils/collocationmatrices'))


import collocation
from plot_elbow_geometry import generate_elbow_geometry, plot_elbow_geometry, plot_solution
from Matrices3d import matrices3D
from utils import initial_guess_x0, inverse_initial_guess_x0, computingpointersF, fix_result
from pointers import pointers_to_BC

# JAM Example II: Flow in a Bent Cylindrical Tube (MATLAB to Python translation)

# Parameters (Reynolds number and dimensionless tube length)
Np = 2
Re = 0.0001  # Reynolds number
L = 10       # Dimensionless tube length
pa = np.array([Re, L])

# Variables: uA (v_x), vA (v_y), wA (v_z), pA (pressure)
list_var_A = ['uA', 'vA', 'wA', 'pA']
NVA = len(list_var_A)

# Derivative list (for indexing purposes)
list_der = ['', 'r0', 'z0', 'rr0', 'zz0', 'rz0', 'time0', 't0', 'tt0', 'tz0', 'tr0']
NDA = len(list_der)

# Read symbolic Jacobians from pickle files
files = ['FAA', 'FAAl', 'FAAr', 'FAAt', 'DFAA', 'DFAAl', 'DFAAr', 'DFAAt']
loaded_data = {}
for file in files:
    with open(f'jacobians/{file}.pkl', 'rb') as f:
        loaded_data[file] = pickle.load(f)

# Define symbolic variables (all derivatives and additional parameters)
yFA_vars = sym.symbols(' '.join([f'yFAv{i}d{j}' for i in range(NVA) for j in range(NDA)]))
symbolic_vars = yFA_vars + sym.symbols('r0 z0 t0 Re L')

# Convert symbolic expressions to numerical functions using lambdify
FAA_func  = [sym.lambdify(symbolic_vars, expr, 'numpy') for expr in loaded_data['FAA']]
FAAl_func = [sym.lambdify(symbolic_vars, expr, 'numpy') for expr in loaded_data['FAAl']]
FAAr_func = [sym.lambdify(symbolic_vars, expr, 'numpy') for expr in loaded_data['FAAr']]
FAAt_func = [sym.lambdify(symbolic_vars, expr, 'numpy') for expr in loaded_data['FAAt']]
DFAA_func = [sym.lambdify(symbolic_vars, expr, 'numpy') for expr in loaded_data['DFAA']]
DFAAl_func = [sym.lambdify(symbolic_vars, expr, 'numpy') for expr in loaded_data['DFAAl']]
DFAAr_func = [sym.lambdify(symbolic_vars, expr, 'numpy') for expr in loaded_data['DFAAr']]
DFAAt_func = [sym.lambdify(symbolic_vars, expr, 'numpy') for expr in loaded_data['DFAAt']]

# Discretization in eta (Chebyshev collocation)
nrA = 4  # Number of Chebyshev collocation points in eta
epsilon = 0.001
r0A = np.linspace(0, 1, nrA)
r0A, dr0A, drr0A = collocation.Chevigood(nrA - 1, 1.0, 0.001)

# Discretization in s (4th order finite differences)
nzA = 101
z0A = np.linspace(0, np.pi/2, nzA)
z0A, dz0A, dzz0A = collocation.finites4th(nzA, np.pi/2)

# Discretization in theta (Fourier collocation)
nxA = 8  # Fourier collocation points
xxA, dxA = collocation.fourdif(nxA, 1)
xxA, dxxA = collocation.fourdif(nxA, 2)

# Compute boundary pointers (using Fortran order)
l_l, l_r, l_t, l_b = computingpointersF(nrA, nzA, nxA)

# Build full derivative matrices
ddr, ddr2, ddz, ddz2, ddx, ddx2, ddrz, ddzx, ddrx = matrices3D(
    nrA, nzA, nxA, dr0A, drr0A, dz0A, dzz0A, dxA, dxxA
)

# Generate elbow geometry (and optionally plot it)
X, Y, Z = generate_elbow_geometry(nzA, nrA, nxA, z0A, r0A, xxA, L)
# plot_elbow_geometry(X, Y, Z)

# Initialize solution variables
uA = np.zeros((nrA, nzA, nxA))
vA = np.zeros((nrA, nzA, nxA))
wA = np.zeros((nrA, nzA, nxA))
pA = np.zeros((nrA, nzA, nxA))

# Total number of grid points
ntA = nrA * nzA * nxA

# Compute boundary pointers for the full domain
JM1, ndA, rA, zA, xA = pointers_to_BC(NVA, ntA, nrA, nzA, nxA, r0A, z0A, xxA)

# Assemble initial guess vector (and previous time steps)
list_var_A_data = [uA, vA, wA, pA]
x0 = initial_guess_x0(list_var_A_data, NVA, ntA)
x0m = x0
x0mm = x0

# Newton's Method parameters
dt = 1e15  # Time step
error = 1e9
iter_count = 0
max_iter = 300

# List of derivative operators (for reference)
derivative_ops = [ddr, ddz, ddr2, ddz2, ddrz, ddx, ddx2, ddzx, ddrx]

# Newton's Iteration
while error > 1e-2 and iter_count < max_iter:
    iter_count += 1

    # Compute time derivative and scaling parameter
    xt = (3 * x0 - 4 * x0m + x0mm) / (2 * dt)
    bp = 3 / (2 * dt)

    # Initialize yfA and construct list xs for function evaluations
    yfA = np.zeros((NVA, NDA, ntA))
    xs = []

    for i in range(NVA):
        l = JM1[i]  # Global indices for variable i
        variable = x0[l]
        variablet = xt[l]
        yfA[i, 0, :] = variable
        yfA[i, 1, :] = ddr @ variable
        yfA[i, 2, :] = ddz @ variable
        yfA[i, 3, :] = ddr2 @ variable
        yfA[i, 4, :] = ddz2 @ variable
        yfA[i, 5, :] = ddrz @ variable
        yfA[i, 6, :] = variablet
        yfA[i, 7, :] = ddx @ variable
        yfA[i, 8, :] = ddx2 @ variable
        yfA[i, 9, :] = ddzx @ variable
        yfA[i, 10, :] = ddrx @ variable

        xs.append(variable)
        xs.append(ddr @ variable)
        xs.append(ddz @ variable)
        xs.append(ddr2 @ variable)
        xs.append(ddz2 @ variable)
        xs.append(ddrz @ variable)
        xs.append(variablet)
        xs.append(ddx @ variable)
        xs.append(ddx2 @ variable)
        xs.append(ddzx @ variable)
        xs.append(ddrx @ variable)

    # Append coordinate arrays
    xs.append(rA)
    xs.append(zA)
    xs.append(xA)

    # Extract subarrays for different boundaries
    xs_t = [param[l_t] for param in xs]  # Top boundary
    xs_l = [param[l_l] for param in xs]  # Left boundary
    xs_r = [param[l_r] for param in xs]  # Right boundary
    xs_b = [param[l_b] for param in xs]  # Bulk (internal points)

    # Evaluate analytical functions at boundaries
    result_b = [func(*xs_b, Re, L) for func in FAA_func]  # Bulk
    result_t = [func(*xs_t, Re, L) for func in FAAt_func]   # Top
    result_l = [func(*xs_l, Re, L) for func in FAAl_func]   # Left
    result_r = [func(*xs_r, Re, L) for func in FAAr_func]   # Right

    FAA = np.zeros((NVA, ntA))
    mask = np.zeros(ntA, dtype=bool)
    mask[l_b] = True
    mask[l_t] = True
    mask[l_r] = True
    mask[l_l] = True
    FAA[:, l_l] = np.array(result_l)  # Left (highest priority)
    FAA[:, l_r] = np.array(result_r)  # Right
    FAA[:, l_t] = np.array(result_t)  # Top
    FAA[:, l_b] = np.array(result_b)  # Bulk (lowest priority)

    # Evaluate analytical Jacobians at boundaries and fix dimensions
    result = [func(*xs_b, Re, L) for func in DFAA_func]
    result_b = [fix_result(x, len(l_b)) for x in result]
    result = [func(*xs_l, Re, L) for func in DFAAl_func]
    result_l = [fix_result(x, len(l_l)) for x in result]
    result = [func(*xs_r, Re, L) for func in DFAAr_func]
    result_r = [fix_result(x, len(l_r)) for x in result]
    result = [func(*xs_t, Re, L) for func in DFAAt_func]
    result_t = [fix_result(x, len(l_t)) for x in result]

    # Assemble the Jacobian array
    DFAA = np.zeros((NVA * NVA * NDA, ntA))
    DFAA[:, l_b] = np.array(result_b)
    DFAA[:, l_l] = np.array(result_l)
    DFAA[:, l_r] = np.array(result_r)
    DFAA[:, l_t] = np.array(result_t)
    DFAA = DFAA.reshape(NVA, NVA * NDA, ntA)

    # Build the numerical Jacobian matrix
    ablock = [[None] * NVA for _ in range(NVA)]
    xablock = [None] * NVA

    for j in range(NVA):
        C1 = FAA[j, :]
        xablock[j] = -C1
        for k in range(NVA):
            km = k * NDA
            kp = (k + 1) * NDA
            C = DFAA[j, range(km, kp), :]
            c0  = C[0, :].reshape(1, -1)
            c1  = C[1, :].reshape(1, -1)
            c2  = C[2, :].reshape(1, -1)
            c3  = C[3, :].reshape(1, -1)
            c4  = C[4, :].reshape(1, -1)
            c5  = C[5, :].reshape(1, -1)
            c6  = C[6, :].reshape(1, -1)
            c7  = C[7, :].reshape(1, -1)
            c8  = C[8, :].reshape(1, -1)
            c9  = C[9, :].reshape(1, -1)
            c10 = C[10, :].reshape(1, -1)

            B = (spdiags(c0, 0, ntA, ntA) +
                 spdiags(c1, 0, ntA, ntA) @ ddr +
                 spdiags(c2, 0, ntA, ntA) @ ddz +
                 spdiags(c3, 0, ntA, ntA) @ ddr2 +
                 spdiags(c4, 0, ntA, ntA) @ ddz2 +
                 spdiags(c5, 0, ntA, ntA) @ ddrz +
                 spdiags(c6, 0, ntA, ntA) * bp +
                 spdiags(c7, 0, ntA, ntA) @ ddx +
                 spdiags(c8, 0, ntA, ntA) @ ddx2 +
                 spdiags(c9, 0, ntA, ntA) @ ddzx +
                 spdiags(c10, 0, ntA, ntA) @ ddrx)
            ablock[j][k] = B

    b = np.concatenate(xablock)
    a = vstack([hstack([m for m in row if m is not None])
                for row in ablock if any(m is not None for m in row)])

    # Solve the linear system and compute error
    dxa = sp.linalg.spsolve(a, b)
    error = np.max(np.abs(dxa))
    print(error)
    if error > 1e9:
        print("Divergence detected, stopping Newton's method.")
        break

    # Update the solution
    x0 = x0 + dxa

print("Newton's iteration completed.")

# Convert the solution vector back to 3D arrays and plot the results
[uA, vA, wA, pA] = inverse_initial_guess_x0(x0, list_var_A_data, nrA, nzA, nxA)
plot_solution(X, Y, Z, uA, wA, pA)
