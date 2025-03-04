import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def generate_elbow_geometry(nzA, nrA, nx, z0A, r0A, xx, L):
    """
    Generates the 3D geometry of the bent cylindrical tube.
    Translated from MATLAB plottingelbowinit.m
    """
    ns = nzA
    neta = nrA
    ntheta = nx
    
    # Define the radius of the center curve
    Rc = 2 * L / np.pi
    
    # Initialize coordinate arrays
    X = np.zeros((ns, neta, ntheta))
    Y = np.zeros((ns, neta, ntheta))
    Z = np.zeros((ns, neta, ntheta))
    
    # Compute the coordinates
    for i in range(ns):
        for j in range(neta):
            for k in range(ntheta):
                X[i, j, k] = (np.cos(z0A[i]) * (2 * L - np.pi * r0A[j] * np.cos(xx[k]))) / np.pi
                Y[i, j, k] = -r0A[j] * np.sin(xx[k])
                Z[i, j, k] = (np.sin(z0A[i]) * (2 * L - np.pi * r0A[j] * np.cos(xx[k]))) / np.pi
    
    # Permute dimensions for correct ordering
    X = np.transpose(X, (1, 0, 2))
    Y = np.transpose(Y, (1, 0, 2))
    Z = np.transpose(Z, (1, 0, 2))
    
    return X, Y, Z

def plot_elbow_geometry(X, Y, Z):
    """
    Plots the mesh of the cylindrical tube bent into a semicircle.
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot mesh lines along s-direction
    for i in range(X.shape[1]):
        for j in range(X.shape[0]):
            ax.plot(X[j, i, :], Y[j, i, :], Z[j, i, :], 'b')
    
    # Plot mesh lines along theta-direction
    for j in range(X.shape[0]):
        for k in range(X.shape[2]):
            ax.plot(X[j, :, k], Y[j, :, k], Z[j, :, k], 'r')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Mesh of the Cylindrical Tube Bent into a Semicircle')
    
    plt.show()

if __name__ == "__main__":
    # Example usage
    nzA, nrA, nx = 10, 10, 10
    z0A = np.linspace(0, np.pi / 2, nzA)
    r0A = np.linspace(0, 1, nrA)
    xx = np.linspace(0, 2 * np.pi, nx)
    L = 1.0
    
    X, Y, Z = generate_elbow_geometry(nzA, nrA, nx, z0A, r0A, xx, L)
    plot_elbow_geometry(X, Y, Z)
    
def plot_solution(X, Y, Z, uA, wA, pA):
    """
    Plots v_z, v_x, and p in three subplots, each combining two 'theta' slices
    (for example j=0 and j=3) using 20 contour levels and a shared color scale 
    for each variable.
    
    Parameters
    ----------
    X, Y, Z : 3D arrays of shape (nrA, nzA, nxA)
        Coordinates in the bent cylindrical tube.
    uA, wA, pA : 3D arrays of shape (nrA, nzA, nxA)
        Velocity components (v_x and v_z) and pressure, respectively.
    """

    fig, axs = plt.subplots(1, 3, figsize=(16, 5))

    # 1) v_z = wA
    min_vz = wA.min()
    max_vz = wA.max()
    levels_vz = np.linspace(min_vz, max_vz, 20)

    # Plot half 1 (e.g. j=0)
    j1 = 0
    x1 = X[:, :, j1]
    z1 = Z[:, :, j1]
    w1 = wA[:, :, j1]
    cs1 = axs[0].contourf(z1, x1, w1, levels=levels_vz, cmap='viridis')

    # Plot half 2 (e.g. j=3)
    j2 = 3
    x2 = X[:, :, j2]
    z2 = Z[:, :, j2]
    w2 = wA[:, :, j2]
    cs2 = axs[0].contourf(z2, x2, w2, levels=levels_vz, cmap='viridis')

    axs[0].set_xlabel('z')
    axs[0].set_ylabel('x')
    axs[0].set_title('v_z')
    axs[0].axis('equal')
    fig.colorbar(cs2, ax=axs[0])

    # 2) v_x = uA
    min_vx = uA.min()
    max_vx = uA.max()
    levels_vx = np.linspace(min_vx, max_vx, 20)

    j1 = 0
    x3 = X[:, :, j1]
    z3 = Z[:, :, j1]
    u3 = uA[:, :, j1]
    cs3 = axs[1].contourf(z3, x3, u3, levels=levels_vx, cmap='viridis')

    j2 = 3
    x4 = X[:, :, j2]
    z4 = Z[:, :, j2]
    u4 = uA[:, :, j2]
    cs4 = axs[1].contourf(z4, x4, u4, levels=levels_vx, cmap='viridis')

    axs[1].set_xlabel('z')
    axs[1].set_ylabel('x')
    axs[1].set_title('v_x')
    axs[1].axis('equal')
    fig.colorbar(cs4, ax=axs[1])

    # 3) p = pA
    min_p = pA.min()
    max_p = pA.max()
    levels_p = np.linspace(min_p, max_p, 20)

    j1 = 0
    x5 = X[:, :, j1]
    z5 = Z[:, :, j1]
    p3 = pA[:, :, j1]
    cs5 = axs[2].contourf(z5, x5, p3, levels=levels_p, cmap='viridis')

    j2 = 3
    x6 = X[:, :, j2]
    z6 = Z[:, :, j2]
    p4 = pA[:, :, j2]
    cs6 = axs[2].contourf(z6, x6, p4, levels=levels_p, cmap='viridis')

    axs[2].set_xlabel('z')
    axs[2].set_ylabel('x')
    axs[2].set_title('p')
    axs[2].axis('equal')
    fig.colorbar(cs6, ax=axs[2])

    plt.tight_layout()
    plt.show()

