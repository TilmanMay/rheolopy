import matplotlib.pyplot as plt
from rheolopy.materials import materials, get_material_by_id
from rheolopy.core import compute_dsigma
import numpy as np

def main():
    """
    Demonstrates how to compute and plot a basic Yield Strength Envelope (YSE)
    for a single material (Quartzite) assuming a linear temperature gradient.
    """
    print("Loading material database...")
    mats = materials()
    
    # Let's use wet quartzite from the database
    quartzite = get_material_by_id(mats, "quartzite_kirby_wet")
    
    if quartzite is None:
        print("Material not found!")
        return

    print(f"Material loaded: {quartzite.type} ({quartzite.source})")

    # Define a depth array (0 to -30 km, negative downwards)
    z_max = 30000.0
    depths = np.linspace(0, -z_max, 100)
    
    # Assume a simple linear geotherm: 20 K/km, surface temp 273.15 K
    temperature_gradient = 20.0 / 1000.0  # K / m
    surface_temp = 273.15
    temperatures = surface_temp - depths * temperature_gradient
    
    # Define strain rate
    strain_rate = 1e-15  # 1/s
    
    print("Computing differential stresses...")
    # compute_dsigma returns dsigma (stress) and depths for both compression and extension
    dsigma, d_plot = compute_dsigma(
        background=quartzite,
        z=depths,
        T=temperatures,
        strain_rate=strain_rate
    )
    
    # Plotting
    plt.figure(figsize=(6, 8))
    
    # Plot Yield Strength Envelope
    # Convert stress to GPa, depth to km
    plt.plot(dsigma / 1e9, -d_plot / 1000, label="YSE (Kirby 1983 Wet Quartzite)", color="cornflowerblue", linewidth=2)
    
    plt.axvline(0, color="black", linewidth=1)
    
    # Formatting
    plt.gca().invert_yaxis()  # Depth increases downwards
    plt.xlabel("Differential Stress (GPa)")
    plt.ylabel("Depth (km)")
    plt.title("Yield Strength Envelope for a Single Material")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig("01_basic_yield_strength.png", dpi=300)
    print("Plot saved as 01_basic_yield_strength.png")
    
    # If running interactively, uncomment this to see the plot:
    # plt.show()

if __name__ == "__main__":
    main()
