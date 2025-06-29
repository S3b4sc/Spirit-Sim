from utils.menu import menu
from utils.data_extraction import load_ovf_ascii
from temp_params import TemperatureDependentParameters
from utils.visualization import visualize_spin_vectors

from spirit import simulation, state

if __name__ == "__main__":

    usr_choice = menu()

    if usr_choice == 1:
        with state.State("config.cfg") as p_state:
            simulation.start(p_state,simulation.METHOD_LLG,simulation.SOLVER_DEPONDT)
            
    elif usr_choice == 2:
        # Load simulation final state
        spins = load_ovf_ascii("output/2025-06-29_09-06-31_Image-00_Spins-final.ovf")
        positions = load_ovf_ascii("output/2025-06-29_09-06-31_positions_final.txt")
        
        visualize_spin_vectors(positions, spins, './')
        
        # Initialize temperature-dependent parameters
        #TDP = TemperatureDependentParameters(spins, positions, lattice_shape=(10, 10, 10))  
        
        
        
        