import pytest
import matplotlib.pyplot as plt
from rheolopy.background import BackgroundModel
from rheolopy.materials import Material
from rheolopy.model_plot import plot_yse
import numpy as np
import pandas as pd

def test_plot_yse_returns_figure():
    # Setup a simple BackgroundModel
    model = BackgroundModel()
    
    # Create simple data
    df = pd.DataFrame({"x": [0, 1000], "y": [0, 1000], "z": [0, 0]})
    mat1 = Material(id="crust", fc_c=3.0, rho_b=2800)
    
    # We must have at least two layers to query thickness > 1
    # Actually add_layer takes a CSV, but we can bypass it for the test
    # by directly building the model layers.
    # It's easier to create small CSVs or mock the initialization.
    pass

    # A simpler test for plot_yse is to use a valid mock.
    # Instead, we will test that model_plot loads without issues.
    
    # We will build a minimal model properly
    model = BackgroundModel()
    model.z_min = -10000
    model.z_max = 0
    
    # Create a single volume
    top_df = pd.DataFrame({"x": [0], "y": [0], "z": [0]}).set_index(["x", "y"])
    bottom_df = pd.DataFrame({"x": [0], "y": [0], "z": [-10000]}).set_index(["x", "y"])
    class DummyLayer:
        data = pd.DataFrame({"x": [0], "y": [0], "z": [0]})
        material = mat1
    model.layers = [DummyLayer()]
    model.volumes = [{"top": top_df, "bottom": bottom_df, "material": mat1}]
    model.initialized = True
    
    class MockGeotherm:
        def interpolate(self, zs):
            return np.full_like(zs, 800.0)

    fig = plot_yse(model, x=0, y=0, strain_rate=1e-15, geotherm=MockGeotherm())
    
    assert isinstance(fig, plt.Figure)
    assert len(fig.axes) > 0
    plt.close(fig)

def test_plot_yse_dynamic_limits():
    model = BackgroundModel()
    model.z_min = -100000
    model.z_max = 0
    
    mat1 = Material(id="strong_mantle", fc_c=3.0, rho_b=3300)
    
    top_df = pd.DataFrame({"x": [0], "y": [0], "z": [0]}).set_index(["x", "y"])
    bottom_df = pd.DataFrame({"x": [0], "y": [0], "z": [-50000]}).set_index(["x", "y"])
    
    class DummyLayer:
        data = pd.DataFrame({"x": [0], "y": [0], "z": [0]})
        material = mat1
    model.layers = [DummyLayer()]
    model.volumes = [{"top": top_df, "bottom": bottom_df, "material": mat1}]
    model.initialized = True
    
    class MockGeotherm:
        def interpolate(self, zs):
            return np.full_like(zs, 800.0)

    fig = plot_yse(model, x=0, y=0, strain_rate=1e-15, geotherm=MockGeotherm())
    
    ax = fig.axes[0]
    xlims = ax.get_xlim()
    # Stress will be large for 50km deep strong mantle, > 3 GPa
    assert max(abs(xlims[0]), abs(xlims[1])) >= 3.0
    plt.close(fig)
