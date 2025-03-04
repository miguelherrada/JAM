import sympy as sp
import pickle
import os
from sympy import symbols, Matrix, diff, simplify, linsolve

def generate_blockA():
    # Define symbolic variables and parameters
    r0, z0, t0, time = symbols('r0 z0 t0 time', real=True)
    L, Re = symbols('L Re', real=True)
    
    # Define velocity components and pressure as functions
    u = sp.Function('u')(r0, z0, t0, time)
    v = sp.Function('v')(r0, z0, t0, time)
    w = sp.Function('w')(r0, z0, t0, time)
    p = sp.Function('p')(r0, z0, t0, time)
    
    # List of derivatives matching MATLAB's structure
    list_der = ['', 'r0', 'z0', 'rr0', 'zz0', 'rz0', 'time0', 't0', 'tt0', 'tz0', 'tr0']
    NDER = len(list_der)  # Dynamic NDER based on list_der
    NVAR = 4  # u, v, w, p

    # Parametric mapping equations (Equation 15)
    F = (sp.cos(z0)*(2*L - sp.pi*r0*sp.cos(t0)))/sp.pi
    G = -r0*sp.sin(t0)
    H = (sp.sin(z0)*(2*L - sp.pi*r0*sp.cos(t0)))/sp.pi
    X = Matrix([F, G, H])
    xs = Matrix([r0, z0, t0])
    
    # Jacobian matrix and inverse
    JX = X.jacobian(xs)
    Jinv = JX.inv().applyfunc(simplify)
    
    # Helper function to parse derivative strings
    def parse_derivative(der_str):
        if der_str == '': return []
        if der_str == 'time0': return [time]
        vars = []
        for c in der_str[:-1]:  # Remove trailing '0'
            if c == 'r': vars.append(r0)
            elif c == 'z': vars.append(z0)
            elif c == 't': vars.append(t0)
        return vars
    
    # Create substitution dictionary for all derivatives
    yFA = sp.Matrix([[symbols(f'yFAv{v}d{d}') for d in range(NDER)] for v in range(NVAR)])
    subs_dict = {}
    
    for var_idx, func in enumerate([u, v, w, p]):
        # Add function itself (no derivative)
        subs_dict[func] = yFA[var_idx, 0]
        
        # Add all derivatives from list_der
        for der_idx, der_str in enumerate(list_der[1:], 1):
            der_vars = parse_derivative(der_str)
            if der_vars:
                der_expr = func.diff(*der_vars)
                subs_dict[der_expr] = yFA[var_idx, der_idx]

    # Compute temporal derivatives (dr0dtime0, etc.)
    x, y, z = X
    dr0, dz0, dt0 = symbols('dr0dtime0 dz0dtime0 dt0dtime0')
    eq1 = diff(x, time) + diff(x, r0)*dr0 + diff(x, z0)*dz0 + diff(x, t0)*dt0
    eq2 = diff(y, time) + diff(y, r0)*dr0 + diff(y, z0)*dz0 + diff(y, t0)*dt0
    eq3 = diff(z, time) + diff(z, r0)*dr0 + diff(z, z0)*dz0 + diff(z, t0)*dt0
    sol = linsolve([eq1, eq2, eq3], [dr0, dz0, dt0])
    dr0dtime0, dz0dtime0, dt0dtime0 = list(sol)[0]

    # Compute spatial derivatives using inverse Jacobian
    def apply_jacobian(func, Jinv):
        df_dr = diff(func, r0)
        df_dz = diff(func, z0)
        df_dt = diff(func, t0)
        return (Jinv[0,0]*df_dr + Jinv[1,0]*df_dz + Jinv[2,0]*df_dt,
                Jinv[0,1]*df_dr + Jinv[1,1]*df_dz + Jinv[2,1]*df_dt,
                Jinv[0,2]*df_dr + Jinv[1,2]*df_dz + Jinv[2,2]*df_dt)

    # First derivatives of velocities
    du_dx, du_dy, du_dz = apply_jacobian(u, Jinv)
    dv_dx, dv_dy, dv_dz = apply_jacobian(v, Jinv)
    dw_dx, dw_dy, dw_dz = apply_jacobian(w, Jinv)

    # Second derivatives for Laplacian
    def compute_laplacian(deriv_x, deriv_y, deriv_z, Jinv):
        d2x = apply_jacobian(deriv_x, Jinv)[0]
        d2y = apply_jacobian(deriv_y, Jinv)[1]
        d2z = apply_jacobian(deriv_z, Jinv)[2]
        return (d2x + d2y + d2z)

    laplace_u = compute_laplacian(du_dx, du_dy, du_dz, Jinv)
    laplace_v = compute_laplacian(dv_dx, dv_dy, dv_dz, Jinv)
    laplace_w = compute_laplacian(dw_dx, dw_dy, dw_dz, Jinv)

    # Pressure gradient
    dp_dx, dp_dy, dp_dz = apply_jacobian(p, Jinv)

    # Material derivatives
    du_dt = diff(u, time) + dr0dtime0*diff(u, r0) + dz0dtime0*diff(u, z0) + dt0dtime0*diff(u, t0)
    dv_dt = diff(v, time) + dr0dtime0*diff(v, r0) + dz0dtime0*diff(v, z0) + dt0dtime0*diff(v, t0)
    dw_dt = diff(w, time) + dr0dtime0*diff(w, r0) + dz0dtime0*diff(w, z0) + dt0dtime0*diff(w, t0)

    # Convective terms
    u_convec = u*du_dx + v*du_dy + w*du_dz
    v_convec = u*dv_dx + v*dv_dy + w*dv_dz
    w_convec = u*dw_dx + v*dw_dy + w*dw_dz

    # Final equations (FAA)
    FAA = [
        (du_dt + u_convec + dp_dx - laplace_u/Re),
        (dv_dt + v_convec + dp_dy - laplace_v/Re),
        (dw_dt + w_convec + dp_dz - laplace_w/Re),
        (du_dx + dv_dy + dw_dz)
    ]

    # Substitute variables and derivatives
    FAA_subs = [eq.subs(subs_dict) for eq in FAA]

    # Boundary conditions
    FAAt = [u, v, w, diff(p, r0)]
    FAAl = [u, v, w - 2*(1 - r0**2), diff(p, z0)]
    FAAr = [diff(u, z0), diff(v, z0), diff(w, z0), p]

    # Substitute BCs and flatten for saving
    FAAt_subs = [cond.subs(subs_dict) for cond in FAAt]
    FAAl_subs = [cond.subs(subs_dict) for cond in FAAl]
    FAAr_subs = [cond.subs(subs_dict) for cond in FAAr]
 # Define symbolic variables for function evaluation
    yFA = sp.Matrix([[symbols(f'yFAv{v}d{d}') for d in range(NDER)] for v in range(NVAR)])

    # Compute Jacobian calculation correctly
    def compute_jacobian(equations, variables):
        return Matrix([[diff(eq, var) for var in variables] for eq in equations])

    # Compute all Jacobians
    c=yFA.reshape(NVAR * NDER, 1)

    #print(c)
    DFAA = compute_jacobian(FAA_subs, yFA.reshape(NVAR * NDER, 1))
    DFAAt = compute_jacobian(FAAt_subs, yFA.reshape(NVAR * NDER, 1)) 
    DFAAl = compute_jacobian(FAAl_subs, yFA.reshape(NVAR * NDER, 1))
    DFAAr = compute_jacobian(FAAr_subs, yFA.reshape(NVAR * NDER, 1) )   
    print(DFAA[0,1])
    #print(DFAA[0,0])

    # Save each component separately
    os.makedirs('jacobians', exist_ok=True)
    components = {
        'FAA': FAA_subs, 'DFAA': DFAA,
        'FAAt': FAAt_subs, 'DFAAt': DFAAt,
        'FAAl': FAAl_subs, 'DFAAl': DFAAl,
        'FAAr': FAAr_subs, 'DFAAr': DFAAr
    }

    for name, data in components.items():
        with open(f'jacobians/{name}.pkl', 'wb') as f:
            pickle.dump(data, f)

    print("All components saved successfully in 'jacobians/' directory")

if __name__ == "__main__":
    generate_blockA()
