import numpy as np
from typing import List
from scipy.spatial import cKDTree


def compute_neighbor_array(positions: np.ndarray, cutoff: float) -> List[List[int]]:
    """
    Computes the neighbor list for each spin based on a cutoff distance.
    
    Parameters:
    positions (np.ndarray): Array of shape (N, 3) containing the positions of spins.
    cutoff (float): Cutoff distance to consider neighbors.
    
    Returns:
    List[List[int]]: A list where each element is a list of indices of neighboring spins.
    """
    tree = cKDTree(positions)
    neighbor_array = np.zeros(positions.shape[0])  # Initialize with empty lists
    
    #for i,pos in enumerate(positions):
    #    indices = tree.query_ball_point(pos, r=cutoff)
    #    #print(f"Spin {i} neighbors: {indices}")
    #    indices.remove(i)
    
    for i, pos in enumerate(positions):
        indices = tree.query_ball_point(pos, r=cutoff)
        indices.remove(i)  # Remove self from neighbors
        neighbor_array[i] = indices
    
    return neighbor_array


def compute_local_J(spins: np.ndarray, neighbor_array: List[List[int]], T: float):
    local_J = np.zeros(spins.shape[0])
    for i, Si in enumerate(spins):
        Jij_sum = sum(np.dot(Si, spins[j]) for j in neighbor_array[i])
        local_J[i] = -Jij_sum / (T + 1e-12)
    return local_J.reshape(lattice_shape[:2])  # 2D grid


def compute_local_K(spins: np.ndarray):
    Sz2 = spins[:, 2] ** 2
    return Sz2.reshape(lattice_shape[:2])



def compute_local_D(spins: np.ndarray, positions: np.ndarray, neighbor_array: List[List[int]]):
    local_D = np.zeros(len(spins))
    for i, Si in enumerate(spins):
        for j in neighbor_array[i]:
            rij = positions[j] - positions[i]
            rij /= np.linalg.norm(rij)
            chirality = np.dot(rij, np.cross(Si, spins[j]))
            local_D[i] += chirality
    return local_D.reshape(lattice_shape[:2])






