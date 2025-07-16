from utils.python.menu import menu
from utils.python.data_extraction import load_ovf_ascii
from temp_params import TemperatureDependentParameters
from utils.python.visualization import visualize_spin_vectors
from utils.python.compute_magnetic_parameters import compute_neighbor_array
#from utils.python.compute_topological_charge import compute_topological_charge_profile


from spirit import system, quantities
import numpy as np
import matplotlib.pyplot as plt
from spirit.quantities import get_topological_charge

import ctypes

from spirit import simulation, state

if __name__ == "__main__":

    usr_choice = menu()

    if usr_choice == 1:
        with state.State("config.cfg") as p_state:
            topologicalCharge = get_topological_charge(p_state)
            topologicalChargeDensity = quantities.get_topological_charge_density(p_state)
            simulation.start(p_state,simulation.METHOD_LLG,simulation.SOLVER_DEPONDT)
        
        print(f"Topological charge: {topologicalCharge}")
        print(f"Topological charge density: {topologicalChargeDensity}")
            
    elif usr_choice == 2:
        # Load simulation final state
        spins = load_ovf_ascii("positions.txt")
        positions = load_ovf_ascii("spins.txt")
        
        print(f"Spins shape: {spins.shape}")
        print(f"Positions shape: {positions.shape}")
        
        visualize_spin_vectors(positions, spins, './')
        
        # Initialize temperature-dependent parameters
        #TDP = TemperatureDependentParameters(spins, positions, lattice_shape=(10, 10, 10))  
        
        
    elif usr_choice == 3:
        # Load simulation final state
        #spins = load_ovf_ascii("output/2025-07-14_18-18-28_Image-00_Spins-final.ovf")
        #positions = load_ovf_ascii("output/2025-07-14_18-18-28_positions_final.txt")
        #
        #compute_neighbor_array(positions, cutoff=1.0)
        
        #compute_topological_charge_profile(positions, spins, plot=True)
        
        #visualize_spin_vectors(positions, spins, './')
        
        # Initialize temperature-dependent parameters
        #TDP = TemperatureDependentParameters(spins, positions, lattice_shape=(10, 10, 10))  
        #
        #TDP.calculate_temperature_dependent_parameters()

        N = 21
        positions = []
        spins = []

        for i in range(N):
            for j in range(N):
                x = i
                y = j
                z = 0
                dx = i - N//2
                dy = j - N//2
                r = np.sqrt(dx**2 + dy**2)
                phi = np.arctan2(dy, dx)
                # Skyrmion profile: center down, edge up
                theta = np.pi * r / (N//2)
                Sx = np.sin(theta) * np.cos(phi)
                Sy = np.sin(theta) * np.sin(phi)
                Sz = np.cos(theta)
                positions.append(f"{x} {y} {z}")
                spins.append(f"{Sx} {Sy} {Sz}")

        with open("positions.txt", "w") as f:
            f.write("\n".join(positions))
        with open("spins.txt", "w") as f:
            f.write("\n".join(spins))
        
    elif usr_choice == 4:
        
        lib = ctypes.CDLL('./utils/c++/libtopo.so')

        lib.compute_layered_topo_charge.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
        lib.compute_layered_topo_charge.restype = None

        pos_file = b"./positions.txt"
        spin_file = b"./spins.txt"

        lib.compute_layered_topo_charge(pos_file, spin_file)
        
