import matplotlib.pyplot as plt
import numpy as np
from rheolopy.background import BackgroundModel
from rheolopy.materials import materials, get_material_by_id
from rheolopy.model_plot import plot_yse

def main():
    """
    Demonstrates how to build a 3D layered BackgroundModel programmatically
    (without relying on a config.ini) and plot a cross-section and a YSE.
    """
    # Load materials
    mats = materials()
    
    # Initialize a new BackgroundModel
    model = BackgroundModel()
    
    # We will build a simple 2-layer crust over a mantle
    # In rheolopy, layers are added by their top boundary depth (positive downwards!)
    # Surface is at 0, Moho is positive depth.
    
    # Layer 1: Upper Crust (Wet Quartzite), surface to 15km
    upper_crust = get_material_by_id(mats, "quartzite_kirby_wet")
    model.add_layer(0, upper_crust)  # Top boundary at z=0
    
    # Layer 2: Lower Crust (Mafic Granulite), 15km to 35km
    lower_crust = get_material_by_id(mats, "granulite_mafic_wilks")
    model.add_layer(15000, lower_crust)  # Top boundary at z=15km
    
    # Layer 3: Lithospheric Mantle (Dry Olivine), 35km downwards
    mantle = get_material_by_id(mats, "olivine_hirth_dry")
    model.add_layer(35000, mantle)  # Top boundary at z=35km
    
    # Initialize the model to compute intersections and thicknesses
    print("Initializing layered model...")
    model.initialize()
    
    # Display what we've built
    print("Model layers at grid center:")
    model.print_layers_at()
    
    # Plot a cross-section of the layers
    print("Plotting model cross-section...")
    fig_slice = model.plot_slice()
    if fig_slice is not None:
         fig_slice.savefig("02_model_slice.png", dpi=300)
         print("Saved cross-section as 02_model_slice.png")
         plt.close(fig_slice) 
    
    # Plot the Yield Strength Envelope through this layered model
    # Note: since we don't have a loaded Geotherm, we'll create a synthetic one 
    # and pass it to plot_yse if desired, but plot_yse needs a proper geotherm or config.
    # To keep it simple, we'll construct a mock config object to pass to plot_yse, 
    # or just use a custom implementation for plotting if we don't have a config.
    # Actually, plot_yse will load the default 'geotherm.csv' if no config is passed.
    print("Plotting Yield Strength Envelope...")
    try:
        fig_yse = plot_yse(model, strain_rate=1e-15)
        fig_yse.savefig("02_layered_yse.png", dpi=300)
        print("Saved Yield Strength Envelope as 02_layered_yse.png")
        plt.close(fig_yse)
    except Exception as e:
        print(f"Could not plot YSE (likely because default config/geotherm is not found): {e}")

if __name__ == "__main__":
    main()
