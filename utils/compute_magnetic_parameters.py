import numpy as np
from typing import List

def compute_local_J(self, spins: np.ndarray, neighbor_list: List[List[int]], T: float):
    local_J = np.zeros(len(spins))
    for i, Si in enumerate(spins):
        Jij_sum = sum(np.dot(Si, spins[j]) for j in neighbor_list[i])
        local_J[i] = -Jij_sum / (T + 1e-12)
    return local_J.reshape(self.lattice_shape[:2])  # 2D grid


def compute_local_K(self, spins: np.ndarray):
    Sz2 = spins[:, 2] ** 2
    return Sz2.reshape(self.lattice_shape[:2])



def compute_local_D(self, spins: np.ndarray, positions: np.ndarray, neighbor_list: List[List[int]]):
    local_D = np.zeros(len(spins))
    for i, Si in enumerate(spins):
        for j in neighbor_list[i]:
            rij = positions[j] - positions[i]
            rij /= np.linalg.norm(rij)
            chirality = np.dot(rij, np.cross(Si, spins[j]))
            local_D[i] += chirality
    return local_D.reshape(self.lattice_shape[:2])






