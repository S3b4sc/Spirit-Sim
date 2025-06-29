import pyvista as pv
import numpy as np


def visualize_spin_vectors(positions: np.ndarray, spins: np.ndarray, route:str) -> None:
    """
    Visualizes the spin vectors in 3D using PyVista.

    Parameters:
    pos (np.ndarray): Array of positions of the spins.
    spins (np.ndarray): Array of spin vectors.
    """
    
    if positions.shape[0] != spins.shape[0]:
        raise ValueError("positions and spins must have same length.")

    mesh = pv.PolyData(positions)
    # Attach the spin vectors to the mesh
    mesh["vectors"] = spins
    mesh.set_active_vectors("vectors")

    # Use glyph filter to add arrows oriented/scaled to spins
    
    #glyphs = mesh.glyph(orient="vectors", scale="vectors", factor=0.5, geom=arrow)
    
    # Use S_z component for scalar coloring
    mesh["S_z"] = spins[:, 2]
    mesh.set_active_scalars("S_z")
    arrow = pv.Arrow(tip_length=0.5, tip_radius=0.7, shaft_radius=0.6)  # Personalize arrow size
    glyphs = mesh.glyph(orient="vectors", scale="vectors", factor=1)
    plotter = pv.Plotter()
    plotter.add_mesh(glyphs, scalars="S_z", cmap="inferno", show_scalar_bar=True)

    # Plot the glyphs
    plotter.set_background("c2c2c2")
    
    #plotter.add_mesh(glyphs, color="black", show_scalar_bar=False)
    plotter.export_html(route + "interactive_spin_plot.html")
    plotter.show()
