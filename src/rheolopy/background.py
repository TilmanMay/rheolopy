import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cmcrameri import cm

from .materials import get_material_by_id, materials
from .io_util import load_config


class Layer3D:
    """Represents a 3D layer with associated material and grid data."""

    def __init__(self, csv_path, material_id):
        self.material = material_id
        self.data = pd.read_csv(csv_path, delimiter=" ")  # expects columns x, y, z


class BackgroundModel:
    """
    Represents a geological background model composed of multiple 3D layers.
    Provides methods for initialization, querying, and visualization.
    """

    def __init__(self):
        self.layers = []
        self.grid_info = None
        self.z_min = -100_000  # Set default minimum depth
        self.z_max = 0  # Set default maximum depth
        self.initialized = False

    def _update_grid_info(self, layer):
        x_vals = np.unique(layer.data["x"].values)
        y_vals = np.unique(layer.data["y"].values)
        z_vals = layer.data["z"].values
        self.grid_info = (
            x_vals.min(),
            x_vals.max(),
            y_vals.min(),
            y_vals.max(),
            len(x_vals),
            len(y_vals),
        )
        # Set z_min and z_max from data, or default if not present
        self.z_min = z_vals.min() if len(z_vals) > 0 else 0
        self.z_max = z_vals.max() if len(z_vals) > 0 else 100_000

    # Add this to your BackgroundModel class
    def _get_unique_xy(self):
        x_vals = np.unique(self.layers[0].data["x"].values)
        y_vals = np.unique(self.layers[0].data["y"].values)
        return x_vals, y_vals

    def _xy_index_to_value(self, x_idx=None, y_idx=None):
        x_vals, y_vals = self._get_unique_xy()
        if x_idx is None:
            x_idx = len(x_vals) // 2
        if y_idx is None:
            y_idx = len(y_vals) // 2
        if not (0 <= x_idx < len(x_vals)) or not (0 <= y_idx < len(y_vals)):
            raise IndexError(f"x_idx or y_idx out of bounds: {x_idx}, {y_idx}")
        return x_vals[x_idx], y_vals[y_idx]

    def add_layer(self, csv_or_depth: str | float | int, material_id) -> None:
        """
        Add a layer to the model, either from a CSV file or a constant depth.

        Parameters
        ----------
        csv_or_depth : str, float, or int
            Path to a CSV file (expects columns x, y, z) or a constant depth value.
        material_id : str or Material
            The material or its id for this layer.
        """
        if self.initialized:
            raise RuntimeError("Cannot add layers after initialization")

        def make_constant_depth_layer(depth: float, material):
            if self.layers:
                # Use x and y from the first loaded layer for exact shape
                ref_layer = self.layers[0].data
                df = ref_layer.copy()
                df["z"] = depth
            elif self.grid_info is not None:
                x_min, x_max, y_min, y_max, nx, ny = self.grid_info
                x_vals = np.linspace(x_min, x_max, nx)
                y_vals = np.linspace(y_min, y_max, ny)
                xx, yy = np.meshgrid(x_vals, y_vals)
                zz = np.full_like(xx, depth, dtype=float)
                df = pd.DataFrame({"x": xx.ravel(), "y": yy.ravel(), "z": zz.ravel()})
            else:
                x_min, x_max, nx = 0, 100_000, 101
                y_min, y_max, ny = 0, 100_000, 101
                x_vals = np.linspace(x_min, x_max, nx)
                y_vals = np.linspace(y_min, y_max, ny)
                xx, yy = np.meshgrid(x_vals, y_vals)
                zz = np.full_like(xx, depth, dtype=float)
                df = pd.DataFrame({"x": xx.ravel(), "y": yy.ravel(), "z": zz.ravel()})
            layer = Layer3D.__new__(Layer3D)
            layer.material = material
            layer.data = df
            return layer

        if isinstance(csv_or_depth, (int, float)):
            layer = make_constant_depth_layer(csv_or_depth, material_id)
        else:
            layer = Layer3D(csv_or_depth, material_id)
            if self.grid_info is None:
                self._update_grid_info(layer)
            else:
                # Update z_min/z_max if this layer extends the model
                z_vals = layer.data["z"].values
                if len(z_vals) > 0:
                    self.z_min = min(self.z_min, z_vals.min())
                    self.z_max = max(self.z_max, z_vals.max())
        self.layers.append(layer)

    def initialize(self):
        # Standardize columns and sort
        for layer in self.layers:
            df = layer.data
            df = df.rename(
                columns={df.columns[0]: "x", df.columns[1]: "y", df.columns[2]: "z"}
            )
            df = df.sort_values(["y", "x"]).reset_index(drop=True)
            layer.data = df

        # Sort layers by mean z (descending: topmost first)
        self.layers.sort(
            key=lambda layer: (
                layer.data["z"].mean() if not layer.data.empty else -np.inf
            ),
            reverse=True,
        )

        # Find the deepest z in all layers
        all_z = np.concatenate(
            [layer.data["z"].values for layer in self.layers if not layer.data.empty]
        )
        deepest_z = all_z.min()
        self.z_min = deepest_z - 10_000  # 10 km below the deepest boundary

        # Crossing check: no lower layer can be above an upper layer at any (x, y)
        for i in range(len(self.layers) - 1):
            top = self.layers[i].data
            bottom = self.layers[i + 1].data
            merged = pd.merge(top, bottom, on=["x", "y"], suffixes=("_top", "_bottom"))
            crossings = merged[merged["z_bottom"] > merged["z_top"]]
            if not crossings.empty:
                raise ValueError(
                    f"Layer crossing detected between layer {i} and {i+1} at points:\n"
                    f"{crossings[['x', 'y', 'z_top', 'z_bottom']]}"
                )

        # Build volumes: each is (top_surface, bottom_surface, material)
        self.volumes = []
        for i in range(len(self.layers) - 1):
            top = self.layers[i]
            bottom = self.layers[i + 1]
            # Merge to compute thickness at each (x, y)
            merged = pd.merge(
                top.data, bottom.data, on=["x", "y"], suffixes=("_top", "_bottom")
            )
            merged["thickness"] = np.abs(merged["z_top"] - merged["z_bottom"])
            # Keep only points where thickness >= 1.0 m
            thick_mask = merged["thickness"] >= 1.0
            top_filtered = merged.loc[thick_mask, ["x", "y", "z_top"]].rename(
                columns={"z_top": "z"}
            )
            bottom_filtered = merged.loc[thick_mask, ["x", "y", "z_bottom"]].rename(
                columns={"z_bottom": "z"}
            )
            self.volumes.append(
                {
                    "top": top_filtered.reset_index(drop=True),
                    "bottom": bottom_filtered.reset_index(drop=True),
                    "material": top.material,
                }
            )
        # Add bottom-most volume: from last layer's bottom to z_min
        last_layer = self.layers[-1]
        bottom_df = last_layer.data.copy()
        bottom_df_zmin = bottom_df.copy()
        bottom_df_zmin["z"] = self.z_min
        self.volumes.append(
            {
                "top": bottom_df.reset_index(drop=True),
                "bottom": bottom_df_zmin.reset_index(drop=True),
                "material": last_layer.material,
            }
        )
        self.initialized = True

        for v in self.volumes:
            if not v["top"].index.names == ["x", "y"]:
                v["top"].set_index(["x", "y"], inplace=True)
            if not v["bottom"].index.names == ["x", "y"]:
                v["bottom"].set_index(["x", "y"], inplace=True)

    def get_material_at(
        self, x_idx=None, y_idx=None, z=None, all_depths=False, n_z=200
    ):
        """
        Returns the material at (x_idx, y_idx, z).
        If all_depths=True, returns (z_sampled, materials) where z_sampled is an array of 200 depths
        (from self.z_max to self.z_min) and materials is a list of materials at each depth.
        """
        x, y = self._xy_index_to_value(x_idx, y_idx)
        if all_depths:
            # Sample 200 z values from top (z_max) to bottom (z_min)
            z_sampled = np.linspace(self.z_max, self.z_min, n_z)
            materials = []
            for z_query in z_sampled:
                mat = None
                for i, v in enumerate(self.volumes):
                    try:
                        z_top = v["top"].loc[(x, y), "z"]
                        z_bottom = v["bottom"].loc[(x, y), "z"]
                    except KeyError:
                        try:
                            z_top = v["top"].loc[(round(x, 2), round(y, 2)), "z"]
                            z_bottom = v["bottom"].loc[(round(x, 2), round(y, 2)), "z"]
                        except KeyError:
                            continue
                    thickness = abs(z_top - z_bottom)
                    if thickness < 1.0:
                        continue
                    if i == len(self.volumes) - 1:
                        if z_top >= z_query >= z_bottom:
                            mat = v["material"]
                            break
                    else:
                        if z_top >= z_query > z_bottom:
                            mat = v["material"]
                            break
                materials.append(mat)
            return z_sampled, materials

        # Default: single depth
        if z is None:
            z = (self.z_min + self.z_max) / 2

        for i, v in enumerate(self.volumes):
            try:
                z_top = v["top"].loc[(x, y), "z"]
                z_bottom = v["bottom"].loc[(x, y), "z"]
            except KeyError:
                try:
                    z_top = v["top"].loc[(round(x, 2), round(y, 2)), "z"]
                    z_bottom = v["bottom"].loc[(round(x, 2), round(y, 2)), "z"]
                except KeyError:
                    continue
            thickness = abs(z_top - z_bottom)
            if thickness < 1.0:
                continue
            if i == len(self.volumes) - 1:
                if z_top >= z >= z_bottom:
                    return v["material"]
            else:
                if z_top >= z > z_bottom:
                    return v["material"]
        return None

    def plot_slice(self, y_index=None, x_index=None):
        """
        Plot a cross-section at the y_index-th unique y value (default) or x_index-th unique x value.
        Black lines for layer surfaces, colored fill between/below for materials.
        Returns the matplotlib Figure object.
        """
        if y_index is not None and x_index is not None:
            raise ValueError("Specify only one of y_index or x_index, not both.")

        x_vals, y_vals = self._get_unique_xy()

        if x_index is not None:
            if x_index < 0 or x_index >= len(x_vals):
                print(f"x_index {x_index} out of range (0 to {len(x_vals)-1})")
                return
            x_val = x_vals[x_index]
            # For each layer, get (y, z) at this x
            layer_yz = []
            for layer in self.layers:
                df = layer.data[layer.data["x"] == x_val]
                if not df.empty:
                    ys = df["y"].values
                    zs = df["z"].values
                    idx = np.argsort(ys)
                    ys = ys[idx]
                    zs = zs[idx]
                    layer_yz.append((ys, zs, layer.material))
            if not layer_yz:
                print(f"No data for x={x_val}")
                return
            all_y = np.unique(np.concatenate([ys for ys, _, _ in layer_yz]))
            y_min, y_max = np.min(all_y), np.max(all_y)
            z_min, z_max = self.z_min, self.z_max
            all_mats = [
                str(mat.id) if hasattr(mat, "id") else str(mat)
                for _, _, mat in layer_yz
            ]
            color_map = {
                mat: cm.lipari(i / max(1, len(all_mats) - 1))
                for i, mat in enumerate(all_mats)
            }
            fig, ax = plt.subplots(figsize=(10, 6))
            # Fill between layers (from top to bottom)
            for i in range(len(layer_yz) - 1):
                ys_top, zs_top, mat_top = layer_yz[i]
                ys_bot, zs_bot, _ = layer_yz[i + 1]
                if not np.array_equal(ys_top, ys_bot):
                    zs_bot_interp = np.interp(ys_top, ys_bot, zs_bot)
                else:
                    zs_bot_interp = zs_bot
                ax.fill_between(
                    ys_top,
                    zs_top,
                    zs_bot_interp,
                    color=color_map[
                        str(mat_top.id) if hasattr(mat_top, "id") else str(mat_top)
                    ],
                    alpha=0.8,
                )
            # Fill below the lowest layer to z_min
            ys_low, zs_low, mat_low = layer_yz[-1]
            ax.fill_between(
                ys_low,
                zs_low,
                z_min,
                color=color_map[
                    str(mat_low.id) if hasattr(mat_low, "id") else str(mat_low)
                ],
                alpha=0.8,
            )
            # Draw black lines for each layer
            for ys, zs, _ in layer_yz:
                ax.plot(ys, zs, color="black", linewidth=1.5)
            from matplotlib.patches import Patch

            handles = [Patch(color=color_map[mat], label=mat) for mat in all_mats]
            ax.legend(
                handles=handles,
                title="Material",
                bbox_to_anchor=(1.05, 1),
                loc="upper left",
            )
            ax.set_xlabel("y")
            ax.set_ylabel("z")
            ax.set_ylim(z_min, z_max)
            ax.set_xlim(y_min, y_max)
            ax.set_title(f"Cross-section at x index={x_index} (x={x_val})")
            plt.show()
            return fig

        # Default: y_index cross-section
        if y_index is None:
            y_index = len(y_vals) // 2

        if y_index < 0 or y_index >= len(y_vals):
            print(f"y_index {y_index} out of range (0 to {len(y_vals)-1})")
            return
        y_val = y_vals[y_index]

        # For each layer, get (x, z) at this y
        layer_xz = []
        for layer in self.layers:
            df = layer.data[layer.data["y"] == y_val]
            if not df.empty:
                xs = df["x"].values
                zs = df["z"].values
                idx = np.argsort(xs)
                xs = xs[idx]
                zs = zs[idx]
                layer_xz.append((xs, zs, layer.material))

        if not layer_xz:
            print(f"No data for y={y_val}")
            return

        all_x = np.unique(np.concatenate([xs for xs, _, _ in layer_xz]))
        x_min, x_max = np.min(all_x), np.max(all_x)
        z_min, z_max = self.z_min, self.z_max

        all_mats = [
            str(mat.id) if hasattr(mat, "id") else str(mat) for _, _, mat in layer_xz
        ]
        color_map = {
            mat: cm.lipari(i / max(1, len(all_mats) - 1))
            for i, mat in enumerate(all_mats)
        }

        fig, ax = plt.subplots(figsize=(10, 6))

        # Fill between layers (from top to bottom)
        for i in range(len(layer_xz) - 1):
            xs_top, zs_top, mat_top = layer_xz[i]
            xs_bot, zs_bot, _ = layer_xz[i + 1]
            if not np.array_equal(xs_top, xs_bot):
                zs_bot_interp = np.interp(xs_top, xs_bot, zs_bot)
            else:
                zs_bot_interp = zs_bot
            ax.fill_between(
                xs_top,
                zs_top,
                zs_bot_interp,
                color=color_map[
                    str(mat_top.id) if hasattr(mat_top, "id") else str(mat_top)
                ],
                alpha=0.8,
            )

        # Fill below the lowest layer to z_min
        xs_low, zs_low, mat_low = layer_xz[-1]
        ax.fill_between(
            xs_low,
            zs_low,
            z_min,
            color=color_map[
                str(mat_low.id) if hasattr(mat_low, "id") else str(mat_low)
            ],
            alpha=0.8,
        )

        # Draw black lines for each layer
        for xs, zs, _ in layer_xz:
            ax.plot(xs, zs, color="black", linewidth=1.5)

        from matplotlib.patches import Patch

        handles = [Patch(color=color_map[mat], label=mat) for mat in all_mats]
        ax.legend(
            handles=handles,
            title="Material",
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
        )
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_ylim(z_min, z_max)
        ax.set_xlim(x_min, x_max)
        ax.set_title(f"Cross-section at y index={y_index} (y={y_val})")
        plt.show()
        return fig

    def print_layers_at(self, x_idx=None, y_idx=None, as_string=False, tag=False):
        """
        Print or return the model's layers at a given (x_idx, y_idx) index.
        If x_idx or y_idx is None, use the model's grid midpoints.
        Skips layers thinner than 1 meter.
        If tag=True, return a list of (z_top, z_bottom, mat) for present layers.
        """
        x, y = self._xy_index_to_value(x_idx, y_idx)

        output_lines = []
        tag_list = []
        for v in self.volumes:
            try:
                z_top = v["top"].loc[(x, y), "z"]
                z_bottom = v["bottom"].loc[(x, y), "z"]
            except KeyError:
                continue
            thickness = abs(z_top - z_bottom)
            if thickness < 1.0:
                continue  # Skip layers thinner than 1 meter
            mat = v["material"]
            mat_id = mat.id if hasattr(mat, "id") else str(mat)
            output_lines.append((z_top, z_bottom, mat_id))
            tag_list.append((z_top, z_bottom, mat))

        # Sort by z_top descending (top to bottom)
        output_lines_sorted = sorted(output_lines, key=lambda tup: tup[0], reverse=True)
        tag_list_sorted = sorted(tag_list, key=lambda tup: tup[0], reverse=True)

        if tag:
            return tag_list_sorted

        if output_lines_sorted:
            result = (
                f"Layers at x_idx={x_idx}, y_idx={y_idx} (x={x:.2f}, y={y:.2f}):\n"
                + "\n".join(
                    f"{z_top:.0f}m - {z_bottom:.0f}m: {mat_id}"
                    for z_top, z_bottom, mat_id in output_lines_sorted
                )
            )
        else:
            result = (
                f"No layer data found at x_idx={x_idx}, y_idx={y_idx} (x={x}, y={y})"
            )

        if as_string:
            return result
        else:
            print(result)

    def plot_layer_thickness(self):
        """
        Plot the thickness of each layer as a 2D map (x vs y, colored by thickness),
        with a maximum of two layers per row and a single shared colorbar outside the plot grid.
        All subplots use the same axis extent, based on the model's full x/y range.
        """
        if not self.initialized:
            print("Model not initialized. Call initialize() first.")
            return

        n_layers = len(self.volumes)
        if n_layers == 0:
            print("No volumes to plot.")
            return

        import matplotlib.pyplot as plt
        import math

        ncols = 2
        nrows = math.ceil(n_layers / ncols)

        fig, axes = plt.subplots(
            nrows, ncols, figsize=(6 * ncols, 5 * nrows), squeeze=False
        )
        axes = axes.flatten()

        # Compute the full model grid
        x_vals, y_vals = self._get_unique_xy()
        x_min, x_max = x_vals.min(), x_vals.max()
        y_min, y_max = y_vals.min(), y_vals.max()

        # Create a MultiIndex for the full grid
        full_index = pd.MultiIndex.from_product([y_vals, x_vals], names=["y", "x"])

        vmin, vmax = None, None
        pivots = []
        for i, v in enumerate(self.volumes):
            top = v["top"].reset_index()
            bottom = v["bottom"].reset_index()
            merged = top.merge(bottom, on=["x", "y"], suffixes=("_top", "_bottom"))
            merged["thickness"] = merged["z_top"] - merged["z_bottom"]
            # Set index to (y, x) for reindexing
            merged = merged.set_index(["y", "x"])
            # Reindex to the full grid, fill missing with np.nan
            merged = merged.reindex(full_index)
            # Pivot to 2D array (y vs x)
            thickness_grid = merged["thickness"].values.reshape(
                len(y_vals), len(x_vals)
            )
            pivots.append(thickness_grid)
            # Update vmin/vmax
            if vmin is None or np.nanmin(thickness_grid) < vmin:
                vmin = np.nanmin(thickness_grid)
            if vmax is None or np.nanmax(thickness_grid) > vmax:
                vmax = np.nanmax(thickness_grid)

        ims = []
        for i, thickness_grid in enumerate(pivots):
            ax = axes[i]
            im = ax.imshow(
                thickness_grid,
                extent=[x_min, x_max, y_min, y_max],
                origin="lower",
                aspect="auto",
                cmap=cm.navia_r,
                vmin=vmin,
                vmax=vmax,
            )
            ims.append(im)
            v = self.volumes[i]
            ax.set_title(
                f"Layer {i+1}: {v['material'].id if hasattr(v['material'], 'id') else str(v['material'])}"
            )
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_min, y_max)

        # Hide unused axes if n_layers is odd
        for j in range(n_layers, nrows * ncols):
            fig.delaxes(axes[j])

        # Adjust layout to make space for colorbar
        plt.tight_layout(rect=[0, 0, 0.92, 1])  # leave space on the right

        # Place colorbar outside the plot grid
        cbar = fig.colorbar(
            ims[0],
            ax=axes[:n_layers],
            orientation="vertical",
            fraction=0.03,
            pad=0.02,
            location="right",
        )
        cbar.set_label("Thickness in m")
        plt.show()
        return fig

    def __repr__(self):
        return self.print_layers_at(as_string=True)


def load_model(config_path=None):
    """
    Load the background model using the config file and layer files.
    If config_path is not given, looks for a folder called 'layers' in the current directory,
    and uses the only .ini file inside as the config. Raises an error if not found or ambiguous.
    All relative paths are resolved relative to the config file location.
    """
    import os
    import glob
    if config_path is None:
        layers_dir = os.path.join(os.getcwd(), "layers")
        if os.path.isdir(layers_dir):
            ini_files = glob.glob(os.path.join(layers_dir, "*.ini"))
            if len(ini_files) == 1:
                config_path = ini_files[0]
            elif len(ini_files) == 0:
                raise FileNotFoundError(
                    "No .ini config file found in the 'layers' directory."
                )
            else:
                raise RuntimeError(
                    f"Multiple .ini files found in 'layers' directory: {ini_files}. Please specify one explicitly."
                )
        else:
            raise ValueError(
                "No config_path provided and no 'layers' directory found in the current working directory."
            )
    config_path = os.path.abspath(config_path)
    config_dir = os.path.dirname(config_path)

    config = load_config(config_path)
    mats = materials()
    model = BackgroundModel()
    layers = json.loads(config["Model"]["layers"])
    for i, l in enumerate(layers):
        layer_path = l[0]
        # If layer_path is a string and not a number, resolve relative to config_dir
        if isinstance(layer_path, str):
            try:
                # Try to interpret as a float (constant depth)
                float(layer_path)
            except ValueError:
                # Not a number, treat as file path
                layer_path = os.path.join(config_dir, layer_path)
        mat = get_material_by_id(mats, l[1])
        model.add_layer(layer_path, mat)
    model.initialize()
    return model
