# API Reference

This document provides a reference for the public API of the `rheolopy` package.

---

## Core Rheology Computations

### `rheolopy.core`

#### `sigma_byerlee(material, z, mode)`

Compute the differential stress using Byerlee's law.

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `material` | `Material` | Material object with `fc_e`, `fc_c`, `lambda_pore`, `rho_b`. |
| `z` | `float` | Depth in meters (z >= 0, surface = 0, positive downward). |
| `mode` | `str` | `'compression'` or `'extension'`. |

**Returns:** `float` — Differential stress in Pa.

---

#### `eta_dislocation(material, temp, strain_rate)`

Compute the viscosity for dislocation creep.

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `material` | `Material` | Material with `a_disloc`, `n`, `q_disloc`. |
| `temp` | `float` | Temperature in Kelvin. |
| `strain_rate` | `float` | Strain rate in 1/s. |

**Returns:** `float` — Dislocation creep viscosity in Pa·s.

---

#### `eta_diffusion(material, temp, strain_rate)`

Compute the viscosity for diffusion creep.

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `material` | `Material` | Material with `d`, `m`, `a_diff`, `q_diff`. |
| `temp` | `float` | Temperature in Kelvin. |
| `strain_rate` | `float` | Strain rate in 1/s. |

**Returns:** `float` — Diffusion creep viscosity in Pa·s.

---

#### `eta_effective(material, temp, strain_rate)`

Compute the effective viscosity using the harmonic mean of dislocation and diffusion creep.

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `material` | `Material` | Material with properties for both creep laws. |
| `temp` | `float` | Temperature in Kelvin. |
| `strain_rate` | `float` | Strain rate in 1/s. |

**Returns:** `tuple(float, float, float)` — `(eta_eff, eta_dislocation, eta_diffusion)`, all in Pa·s.

---

#### `calc_peierls(material, temp, strain_rate)`

Compute the Peierls stress for a material at a given temperature and strain rate.

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `material` | `Material` | Material with `eps_peierls`, `sigma_peierls`, `stress_pd`. |
| `temp` | `float` | Temperature in Kelvin. |
| `strain_rate` | `float` | Strain rate in 1/s. |

**Returns:** `float` — Peierls stress in Pa. Returns `0.0` if Peierls parameters are not available.

---

#### `sigma_d(material, z, temp, strain_rate, mode, return_all, return_index, eta_min, eta_max)`

Compute the differential stress at a given depth and temperature, considering Byerlee's law, dislocation creep, diffusion creep, and optionally Peierls creep. Returns the minimum stress required for deformation (the controlling mechanism).

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `material` | `Material` | — | Material object. |
| `z` | `float` | — | Depth in meters (z >= 0, surface = 0, positive downward). |
| `temp` | `float` | — | Temperature in Kelvin. |
| `strain_rate` | `float` | `1e-17` | Strain rate in 1/s. |
| `mode` | `str` | — | `'compression'` or `'extension'`. |
| `return_all` | `bool` | `False` | If `True`, return all stress components. |
| `return_index` | `bool` | `False` | If `True`, return the index of the controlling mechanism. |
| `eta_min` | `float` | `None` | Minimum allowable effective viscosity (Pa·s). |
| `eta_max` | `float` | `None` | Maximum allowable effective viscosity (Pa·s). |

**Returns:** `float` or `tuple` — Differential stress in Pa, or tuple of all stresses if `return_all=True`.

---

#### `compute_dsigma(background, z, T, strain_rate, x_idx, y_idx, return_all, return_index)`

Compute differential stress for a `Material` or `BackgroundModel` at given depths and temperatures. Handles both scalar and array inputs.

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `background` | `Material` or `BackgroundModel` | The model or material to compute stresses for. |
| `z` | `float` or `array-like` | Depth(s) in meters (z >= 0, surface = 0, positive downward). |
| `T` | `float` or `array-like` | Temperature(s) in Kelvin. |
| `strain_rate` | `float` | Strain rate in 1/s. |
| `x_idx` | `int`, optional | x grid index (for `BackgroundModel` only). |
| `y_idx` | `int`, optional | y grid index (for `BackgroundModel` only). |
| `return_all` | `bool` | If `True`, return all stress components. |
| `return_index` | `bool` | If `True`, return the controlling mechanism index. |

**Returns:**
- Scalar input: `(compression_stress, extension_stress)` in Pa.
- Array input: `(dsigma, depths)` where `dsigma` is a concatenated array of compression (negative) and extension (positive) stresses, and `depths` is the corresponding depth array. This is ready for direct plotting as a yield strength envelope.

---

## Materials

### `rheolopy.materials`

#### `class Material`

Represents a rock material with experimentally derived rheological parameters.

**Attributes:**

| Attribute | Type | Units | Description |
|-----------|------|-------|-------------|
| `id` | `str` | — | Unique identifier. |
| `source` | `str` | — | Literature source for the parameters. |
| `type` | `str` | — | Rock type (e.g., "quartzite wet"). |
| `fc_e` | `float` | — | Friction coefficient for extension. |
| `fc_c` | `float` | — | Friction coefficient for compression. |
| `lambda_pore` | `float` | — | Pore fluid factor. |
| `rho_b` | `float` | kg/m³ | Bulk density. |
| `a_disloc` | `float` | Pa⁻ⁿ/s | Dislocation creep pre-exponential factor. |
| `n` | `float` | — | Power law exponent. |
| `q_disloc` | `float` | J/mol | Dislocation creep activation energy. |
| `a_diff` | `float` | 1/Pa/s | Diffusion creep pre-exponential factor. |
| `q_diff` | `float` | J/mol | Diffusion creep activation energy. |
| `d` | `float` | m | Grain size. |
| `m` | `float` | — | Grain size exponent. |
| `eps_peierls` | `float` | 1/s | Peierls reference strain rate. |
| `sigma_peierls` | `float` | Pa | Peierls critical stress. |
| `stress_pd` | `float` | Pa | Peierls stress constant. |

---

#### `materials(database_path='database.json')`

Load all materials from a JSON database file.

**Parameters:**

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `database_path` | `str` | `'database.json'` | Path to the JSON database. Falls back to the bundled database. |

**Returns:** `list[Material]` — Sorted list of `Material` objects.

---

#### `get_material_by_id(mats, material_id)`

Look up a material by its `id` from a list of materials.

**Parameters:**

| Name | Type | Description |
|------|------|-------------|
| `mats` | `list[Material]` | List of materials (from `materials()`). |
| `material_id` | `str` | The `id` to search for. |

**Returns:** `Material` or `None`.

---

## Geotherm

### `rheolopy.geotherm`

#### `class Geotherm`

Represents a temperature–depth profile. Loads from CSV (expects columns: depth in km, temperature in °C). Internally converts to meters and Kelvin.

##### `__init__(path=None, config=None)`

Load a geotherm from a CSV file. If `path` is `None`, resolves from config. If `config` is also `None`, uses the bundled default config and geotherm.

##### `interpolate(z)`

Interpolate temperature at given depth(s).

| Name | Type | Description |
|------|------|-------------|
| `z` | `float` or `ndarray` | Depth(s) in meters (positive downward). |

**Returns:** Temperature(s) in Kelvin.

##### `as_array()`

**Returns:** `ndarray` — Raw geotherm data as a 2-column array `(depth [m], temperature [K])`.

---

## Background Model

### `rheolopy.background`

#### `class Layer3D`

Represents a single 3D geological layer surface.

##### `__init__(csv_path, material_id)`

| Name | Type | Description |
|------|------|-------------|
| `csv_path` | `str` | Path to a space-delimited CSV with columns `x`, `y`, `z` (all in meters). |
| `material_id` | `str` or `Material` | Material for this layer. |

---

#### `class BackgroundModel`

A layered geological model composed of multiple `Layer3D` surfaces. Provides methods for building the model, querying materials at any point, and visualization.

##### `add_layer(csv_or_depth, material_id)`

Add a layer to the model.

| Name | Type | Description |
|------|------|-------------|
| `csv_or_depth` | `str`, `float`, or `int` | Path to a CSV file, or a constant depth value (m). |
| `material_id` | `str` or `Material` | The material for this layer. |

##### `initialize()`

Finalize the model after all layers have been added. Sorts layers, validates no crossings, and builds internal volumes.

##### `get_material_at(x_idx=None, y_idx=None, z=None, all_depths=False, n_z=200)`

Query the material at a specific (x, y, z) point.

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `x_idx` | `int` | midpoint | x grid index. |
| `y_idx` | `int` | midpoint | y grid index. |
| `z` | `float` | midpoint | Depth in meters (positive = below surface). |
| `all_depths` | `bool` | `False` | If `True`, return materials at all sampled depths. |
| `n_z` | `int` | `200` | Number of depth samples when `all_depths=True`. |

##### `print_layers_at(x_idx=None, y_idx=None, as_string=False, tag=False)`

Print or return the layer stack at a given grid point.

##### `plot_slice(y_index=None, x_index=None)`

Plot a vertical cross-section through the model. Returns a `matplotlib.Figure`.

##### `plot_layer_thickness()`

Plot 2D thickness maps for all layers. Returns a `matplotlib.Figure`.

---

#### `load_model(config_path=None)`

Load a `BackgroundModel` from an INI config file that specifies layer CSV paths and material IDs.

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `config_path` | `str` | auto-detect | Path to the `.ini` config. If `None`, looks for a `layers/` directory with a single `.ini` file. |

**Returns:** `BackgroundModel` — An initialized model ready for queries.

---

## Plotting

### `rheolopy.model_plot`

#### `plot_yse(model, x=None, y=None, strain_rate=None, ax=None, geotherm=None)`

Plot the yield strength envelope at a given (x, y) grid point, including layer backgrounds and a geotherm overlay.

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `model` | `BackgroundModel` | — | The initialized model. |
| `x` | `int` | midpoint | x grid index. |
| `y` | `int` | midpoint | y grid index. |
| `strain_rate` | `float` | from config | Strain rate in 1/s. |
| `ax` | `matplotlib.Axes` | new figure | Axes to plot on. |
| `geotherm` | `Geotherm` | from config | Custom geotherm to use. |

**Returns:** `matplotlib.Figure` if a new figure was created, otherwise `None`.

---

## Configuration

### `rheolopy.io_util`

#### `load_config(path=None)`

Load an INI configuration file. Falls back to the bundled `config.ini`.

**Returns:** `configparser.ConfigParser`.

#### `resolve_path(config, option, section='General')`

Resolve a relative path from a config option, relative to the config file's directory.

**Returns:** `str` — Absolute path.

---

## Bundled Data

The package ships with:

| File | Description |
|------|-------------|
| `config.ini` | Default configuration (strain rate, paths to geotherm and database). |
| `database.json` | Material database with 29 entries covering various rock types and their rheological parameters. |
| `geotherm.csv` | Default continental geotherm (McKenzie et al. 2005), columns: depth (km), temperature (°C). |
