from classe import TemperatureDependentParameters
from spirit import simulation, state

if __name__ == "__main__":

    with state.State("config.cfg") as p_state:
        simulation.start(p_state,simulation.METHOD_LLG,simulation.SOLVER_DEPONDT)
    