import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt

def compute_topological_charge_profile_by_unitcell(positions: np.ndarray, spins: np.ndarray, 
                                                    lattice_constant_z: float, plot=True) -> np.ndarray:
    """
    Computes and plots the topological charge Q per unit cell in the z-direction,
    aggregating spins in all basis atoms within a single unit cell.

    Parameters:
        positions (np.ndarray): Array of shape (N, 3), spin positions (x, y, z)
        spins (np.ndarray): Array of shape (N, 3), spin vectors
        lattice_constant_z (float): Distance between unit cells along z
        plot (bool): Whether to show a plot of Q vs z unit cell index

    Returns:
        np.ndarray: Array of shape (n_layers, 2), [z_center, Q(z)]
    """
    z_min = positions[:, 2].min()
    z_indices = np.floor((positions[:, 2] - z_min) / lattice_constant_z).astype(int)
    n_cells = z_indices.max() + 1

    q_profile = []

    for layer_index in range(n_cells):
        mask = z_indices == layer_index
        pos_layer = positions[mask]
        spin_layer = spins[mask]

        Q = compute_topological_charge(pos_layer, spin_layer)
        z_center = z_min + (layer_index + 0.5) * lattice_constant_z
        q_profile.append((z_center, Q))

    q_profile = np.array(q_profile)

    if plot:
        plt.figure(figsize=(6, 4))
        plt.plot(q_profile[:, 0], q_profile[:, 1], marker='o')
        plt.xlabel('Layer depth z (center of unit cell)')
        plt.ylabel('Topological charge Q(z)')
        plt.title('Topological charge per unit cell layer')
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    return q_profile

def compute_topological_charge(positions: np.ndarray, spins: np.ndarray) -> float:
    xy = positions[:, :2]
    xs = np.unique(xy[:, 0])
    ys = np.unique(xy[:, 1])
    Nx, Ny = len(xs), len(ys)

    grid_spins = np.zeros((Nx, Ny, 3))
    for idx, (x, y) in enumerate(xy):
        i = np.where(xs == x)[0][0]
        j = np.where(ys == y)[0][0]
        grid_spins[i, j] = spins[idx]

    Q = 0.0
    for i in range(Nx - 1):
        for j in range(Ny - 1):
            S1 = grid_spins[i, j]
            S2 = grid_spins[i + 1, j]
            S3 = grid_spins[i + 1, j + 1]
            Q += triangle_contribution(S1, S2, S3)

            S1 = grid_spins[i, j]
            S2 = grid_spins[i + 1, j + 1]
            S3 = grid_spins[i, j + 1]
            Q += triangle_contribution(S1, S2, S3)

    return Q / (4 * np.pi)


def triangle_contribution(S1: np.ndarray, S2: np.ndarray, S3: np.ndarray) -> float:
    S1 = S1 / np.linalg.norm(S1)
    S2 = S2 / np.linalg.norm(S2)
    S3 = S3 / np.linalg.norm(S3)

    omega = np.arctan2(
        np.dot(S1, np.cross(S2, S3)),
        1.0 + np.dot(S1, S2) + np.dot(S2, S3) + np.dot(S3, S1)
    )
    return omega

