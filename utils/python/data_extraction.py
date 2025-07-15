import numpy as np


def load_energy(path="./output/energies.txt"):
    return np.loadtxt(path)  # Each row: [step, energy]

def load_magnetization(path="output/magnetization.txt"):
    return np.loadtxt(path)  # Each row: [step, Mx, My, Mz]

def load_ovf_ascii(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    start = next(i for i, l in enumerate(lines) if "# Begin: Data Text" in l) + 1
    data = np.loadtxt(lines[start:])
    return data#data.reshape(-1, 3)  # [Nx*Ny*Nz, 3]

