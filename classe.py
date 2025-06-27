import os
import numpy as np
from typing import List, Dict
from spirit import state, simulation, configuration, quantities, system
import matplotlib.pyplot as plt

class TemperatureDependentParameters:
    def __init__(self, temperatures:List[float], magnetic_parameters:Dict[str,List[float]]) -> None:
        self.temperatures = temperatures
        self.magnetic_paremters = magnetic_parameters
        self.results = {
        "M": [],
        "E": [],
        "Eani": [],
        "Q": [],
        "T": temperatures
        }

        
    def update_input_cfg(self, template_path: str, output_path: str, temperature: float):
        with open(template_path, "r") as f:
            text = f.read()

        # Combine temperature and all magnetic parameters
        values = {"TEMP": temperature, **self.magnetic_parameters}

        for key, val in values.items():
            text = text.replace(f"__{key}__", str(val))

        with open(output_path, "w") as f:
            f.write(text)
        
        return output_path
    
    def compute_D(self):
        if "J_eff" not in self.results:
            self.compute_J()
        D_eff = []
        for i, T in enumerate(self.results["T"]):
            q = self._get_spiral_q_from_spin_config(T)  # you define this
            D = 2 * self.results["J_eff"][i] * q
            D_eff.append(D)
        self.results["D_eff"] = D_eff
        return D_eff
    
    def compute_J(self):
        J_eff = []
        for T, M in zip(self.results["T"], self.results["M"]):
            Mnorm = np.linalg.norm(M)
            J = Mnorm**2 / (T + 1e-12)
            J_eff.append(J)
        self.results["J_eff"] = J_eff
        return J_eff
    
    def compute_K(self):
        K_eff = []
        N = self.num_spins  # should be defined from geometry
        for T, M, Eani in zip(self.results["T"], self.results["M"], self.results["Eani"]):
            Mz = M[2]
            K = -Eani / (N * Mz**2 + 1e-12)  # avoid division by zero
            K_eff.append(K)
        self.results["K_eff"] = K_eff
        return K_eff
    
    def _initialize_simulation(self, T):
        
        p_state = state.setup("geometry.json", "hamiltonian.json")
        simulation.set_temperature(p_state, T, damping=0.5)
        return p_state
    
    def _sample(self, p_state, n_steps):
        
        M, E, Eani = [], [], []
        for _ in range(n_steps // 100):
            simulation.iterate(p_state, 100)
            M.append(quantities.magnetization(p_state))
            E.append(system.get_energy(p_state))
            Eani.append(quantities.energy_contributions(p_state)["anisotropy"])
        return np.mean(M, axis=0), np.mean(E), np.mean(Eani)
    
    def run(self, equil_steps=5000, sample_steps=10000, repeats=3):
        for T in self.temperatures:
            M_accum, E_accum, Eani_accum = [], [], []
            for _ in range(repeats):
                p_state = self._initialize_simulation(T)
                self._equilibrate(p_state, equil_steps)
                M, E, Eani = self._sample(p_state, sample_steps)
                M_accum.append(M)
                E_accum.append(E)
                Eani_accum.append(Eani)
            self.results["M"].append(np.mean(M_accum, axis=0))
            self.results["E"].append(np.mean(E_accum))
            self.results["Eani"].append(np.mean(Eani_accum))
            
            filename = f"spins_T{T:.1f}.npy"
            spins = configuration.get(p_state)
            np.save(f"spins/{filename}", spins)
            
            
    def plot_magnetization(self):
        
        T = self.results["T"]
        Mz = [np.abs(m[2]) for m in self.results["M"]]
        plt.plot(T, Mz, marker='o')
        plt.xlabel("Temperature (K)")
        plt.ylabel("|M_z|")
        plt.title("Magnetization vs Temperature")
        plt.show()
    