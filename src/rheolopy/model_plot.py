import numpy as np
import logging
from cmcrameri import cm

from .io_util import load_config
from .core import compute_dsigma
from .geotherm import Geotherm

import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def plot_yse(model, x=None, y=None, strain_rate=None, ax=None, geotherm=None):
    """
    Plot the yield strength envelope (YSE) at a given (x, y) index.

    Includes background layers as colored bands and the geotherm overlay.
    If x or y is None, the function uses the grid midpoint.

    Parameters
    ----------
    model : BackgroundModel
        The geological model.
    x : int, optional
        Grid index for the x-coordinate.
    y : int, optional
        Grid index for the y-coordinate.
    strain_rate : float, optional
        Reference strain rate in 1/s. Defaults to config value.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, a new figure is created.
    geotherm : Geotherm, optional
        Geotherm instance to use for temperature. Defaults to config value.

    Returns
    -------
    fig : matplotlib.figure.Figure or None
        The created Figure if `ax` was None, else None.
    """
    config = getattr(model, "config", None) or load_config("config.ini")
    if strain_rate is None:
        strain_rate = config.getfloat("General", "strain_rate", fallback=1e-17)
    if geotherm is None:
        geotherm_obj = Geotherm(config=config)
        geotherm = geotherm_obj.as_array()
    else:
        geotherm_obj = geotherm  # if user passes a Geotherm instance

    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 8))
        created_fig = True

    # Get the (x, y) indices and values
    x_vals, y_vals = model._get_unique_xy()
    if x is None:
        x_idx = len(x_vals) // 2
    else:
        x_idx = x
    if y is None:
        y_idx = len(y_vals) // 2
    else:
        y_idx = y

    # Bounds check
    if not (0 <= x_idx < len(x_vals)):
        raise IndexError(f"x index {x_idx} out of bounds (0 to {len(x_vals)-1})")
    if not (0 <= y_idx < len(y_vals)):
        raise IndexError(f"y index {y_idx} out of bounds (0 to {len(y_vals)-1})")

    x_val = x_vals[x_idx]
    y_val = y_vals[y_idx]

    logger.info(f"Plotting YSE at (x={x_val}, y={y_val})")

    # Get top and bottom z at this (x, y)
    top_layer = model.layers[0]
    row = top_layer.data[
        (top_layer.data["x"] == x_val) & (top_layer.data["y"] == y_val)
    ]
    if row.empty:
        raise ValueError(
            f"No data found for (x={x_val}, y={y_val}). Check your grid indices or coordinates."
        )
    top_z = row["z"].values[0]
    bottom_z = model.z_min
    # Depth and temperature arrays specific to (x, y)
    zmin = bottom_z
    zs = np.linspace(0, zmin, 300)

    # Load geotherm and interpolate to this (x, y)
    T = geotherm_obj.interpolate(zs)

    # Compute YSE
    sigma_plot, z_plot = compute_dsigma(
        model, zs, T, strain_rate, x_idx=x_idx, y_idx=y_idx
    )
    z_plot_shifted = z_plot - top_z

    # Plot background layers as colored bands for this (x, y)
    layers = model.print_layers_at(x_idx=x_idx, y_idx=y_idx, tag=True)
    layer_materials = [layer.material for layer in model.layers]
    material_to_index = {mat: i for i, mat in enumerate(layer_materials)}
    n_layers = len(layer_materials)
    colormap = cm.lipari

    for z_top, z_bottom, material in layers:
        # Use the global index for color assignment
        idx = material_to_index.get(material, 0)
        color = colormap(idx / max(1, n_layers - 1))
        top_km = -z_top / 1000
        bottom_km = -z_bottom / 1000
        span = ax.axhspan(ymin=top_km, ymax=bottom_km, color=color, alpha=0.3, zorder=0)
        mid_km = (top_km + bottom_km) / 2
        ax.text(
            ax.get_xlim()[1] if ax.get_xlim()[1] != 1.0 else 1.1,
            y=mid_km,
            s=material.id if hasattr(material, "id") else str(material),
            fontsize=8,
            color="black",
            ha="right",
            va="center",
            alpha=0.6,
            zorder=1,
        )
    # Plot YSE
    ax.plot(
        sigma_plot / 1e9,
        -z_plot_shifted / 1000,
        label="YSE",
        color="cornflowerblue",
        zorder=2,
    )
    # Add vertical black line at x = 0 (0 GPa)
    ax.axvline(0, color="black", linestyle="-", linewidth=1.2, zorder=3)

    # Plot geotherm (temperature) as dashed orange line
    ax2 = ax.twiny()
    ax2.plot(
        T - 273.15,
        -(zs + top_z) / 1000,  # shift by top_z so depths are absolute
        linestyle="--",
        color="orange",
        label="Temperature",
    )
    ax2.set_xlabel("Temperature (°C)")
    ax2.set_xlim(0, np.nanmax(T - 273.15) * 1.05)
    ax2.get_yaxis().set_visible(False)

    # Labels and legend
    ax.set_xlim(-3, 3)
    ax.set_xlabel("Differential Stress (GPa)")
    ax.set_ylabel("Depth (km)")
    ax.set_title(f"Yield Strength Envelope at (x={x_val:.1f}, y={y_val:.1f})")
    ax.invert_yaxis()
    ax.legend(loc="lower left")

    if created_fig:
        return fig
