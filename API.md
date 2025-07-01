# Rheolog API Documentation

This document describes the main classes and functions provided by the `rheology` package.

---

## Model and Layers

### BackgroundModel
Main class for building, initializing, and querying a layered 3D geological model.
Since every layer in the model corresponds to the top boundray of the geological unit, the last layer (the deepest one) has a default thickness of 10km.
**Methods:**
- `initialize()`
  - Initialize the model after adding layers.
- `add_layer(csv_or_depth, material)`
  - Add a layer from a CSV file or constant depth.
- `get_material_at(x_idx=None, y_idx=None, z=None)`
  - Get the material at a given grid index and depth.
- `plot_layer_thickness()`
  - Plot thickness of each layer as a 2D map. Returns a matplotlib Figure.
- `plot_slice(y_index=None, x_index=None)`
  - Plot a cross-section at a given x- or y-index. Returns a matplotlib Figure.

### Layer3D
Represents a single 3D layer with associated material and grid data.

**Attributes:**
- `material`: The material assigned to this layer. Has to be of `Material` type
- `data`: DataFrame with columns `x`, `y`, `z` for the layer surface.

---

## Configuration and Utilities

### load_model(config_path=None) -> BackgroundModel
Load a background model from a configuration file and associated layer files.
- `config_path` (str, optional): Path to the config file. If not provided, tries to auto-detect (it will search for a `layers` folder in your current dict and will load a .ini file, if it exists).
- **Returns:** `BackgroundModel` instance.

### load_config(path=None) -> configparser.ConfigParser
Load a configuration (INI) file.
- `path` (str, optional): Path to the config file. If not provided, loads the built-in config.
- **Returns:** ConfigParser object.

### resolve_path(config, option, section="General") -> str
Resolve a file path from the config, relative to the config file location.
- `config`: ConfigParser object.
- `option` (str): The option name in the config.
- `section` (str): The section name in the config.
- **Returns:** Absolute path as string.

---

## Materials

### Material
Class representing a rock/material with rheological properties.

**Attributes:**
- `id`: Unique identifier, used for quering the material
- `source`: What are these values based on / where are they from
- `type`: What does your material pretend to be ? (Think of it as its name)

**Byerlee:**
- `fc_e`: Friction coefficient for extension
- `fc_c`: Friction coefficient for compression
- `lambda_pore`: Pore fluid factor (aka. lambda)
- `rho_b`: Bulk density of the rock in kg/m3

**Dislocation:**
- `a_disloc`: Power law strain rate in $Pa^{-n}s^{-1}$ for dislocation creep

    Could also be given as A in $MPa$ or $GPa$. This is usually the value one get's from publications as it is the preexponential scaling factor in the power law equations.
    Example: You find a value (i.e. 1e5 MPa) you can set a_disloc = 1e5, then you also make sure your ``convert`` option is set to ``MPa``. Same for values in $GPa$
- `q_disloc`: Activation energy for dislocation in J/mol

  Usually this is given in kJ/mol in literature so make sure to convert accordingly
- `n`: Power law exponent


**Diffusion:** 
- `a_diff`: Power law strain rate in $Pa^{-n}s^{-1}$ for diffusion creep

    Could also be given as A in $MPa$ or $GPa$. This is usually the value one get's from publications as it is the preexponential scaling factor in the power law equations.
    Example: You find a value (i.e. 1e5 MPa) you can set a_diff = 1e5, then you also make sure your ``convert`` option is set to ``MPa``. Same for values in $GPa$
- `q_diff`: Activation energy for diffusion in J/mol

    Usually this is given in kJ/mol in literature so make sure to convert accordingly 
- `d`: grain size in m
- `m`: grain size exponent

    Caution: Some define this exponent as negative. Here, it is always positive defined

**Other:**
- `convert`: Either `false`, `"MPa"` or `"GPa"` depending on what values you have for `a_disloc` and `a_diff`

### materials(database_path="database.json") -> list[Material]
Load all materials from the database file.
- `database_path` (str): Path to the material database JSON.
- **Returns:** List of `Material` objects.

### get_material_by_id(mats, material_id) -> Material
Retrieve a material by its ID from a list of materials.
- `mats` (list): List of `Material` objects.
- `material_id` (str): The material ID.
- **Returns:** `Material` object.

---

### Geotherm
Represents a geotherm (temperature-depth profile) and provides interpolation and data access.

**Initialization:**
- `Geotherm(path=None, config=None)`
  - Loads a geotherm from a CSV file. If `path` is not provided, tries to resolve from the config file.

**Attributes:**
- `path`: Path to the geotherm file used.
- `data`: Numpy array with columns `[depth (m), temperature (K)]`.

**Methods:**
- `interpolate(depths_m)`
  - Interpolate temperature at given depth(s) in meters. Returns temperature(s) in Kelvin.
- `as_array()`
  - Return the raw geotherm data as a numpy array (depth [m], temp [K]).

---

## Plotting

### plot_layer_thickness(self) -> matplotlib.figure.Figure
Plot the thickness of each model layer as a 2D map.

### plot_slice(self, y_index=None, x-index=None) -> matplotlib.figure.Figure
Plot a cross-section at a given x- or y-index.

### plot_yse(model, x=None, y=None, strain_rate=None, geotherm=None) -> matplotlib.figure.Figure
Plot the yield strength envelope (YSE) at a given (x, y) location.

---

## Rheology Calculations

### sigma_byerlee(material, z, mode) -> float
Compute maximum differential stress using Byerlee’s law.

### eta_effective(material, temp, strain_rate) -> tuple
Compute the effective viscosity $\eta_{eff}$ for a material at given T and strain rate. Also returns viscosities $\eta_{disloc}$ and $\eta_{diff}$

### calc_peierls(material, temp, strain_rate) -> float
Calculate Peierls stress for olivine. This form of deformation is not well investigated for most materials. For olivine however it is. The relevant parameters for this are hardcoded (for now) as :
```python
    epsPeierls = 5.7e11 # Dorn critical strain rate
    sigmaPeierls = 8.5e9 #Peierls critical stress
    stressPD = 200e6 
```
it is automatically calculated for olivine as long as it is clear that the material is supposed to be olivine:

```python
    if "olivine" in material.type.lower():
        #calculate 
```
### compute_dsigma(background, z, T, strain_rate, x_idx=None, y_idx=None, return_all=False, return_index=False,) -> tuple
Compute the minimum differential stress according to all the provided laws.
#### compute_dsigma(background, z, T, strain_rate, x_idx=None, y_idx=None, return_all=False, return_index=False) -> tuple

Compute the minimum differential stress according to all the provided laws.

**Parameters:**
- `background` (`BackgroundModel` or `Material`): Either a `BackgroundModel` instance (for spatially varying materials) or a single `Material` object.
- `z` (float or array-like): Depth(s) at which to compute the stress. Can be a scalar or a 1D array.
- `T` (float or array-like): Temperature(s) at which to compute the stress. Must match the shape of `z`.
- `strain_rate` (float): Strain rate to use in the calculation.
- `x_idx` (int, optional): X grid index for spatial queries (used if `background` is a `BackgroundModel`).
- `y_idx` (int, optional): Y grid index for spatial queries (used if `background` is a `BackgroundModel`).
- `return_all` (bool, optional): If `True`, return all stress components from the underlying law(s). Default is `False`.
- `return_index` (bool, optional): If `True`, also return the index of the controlling mechanism. Default is `False`.

**Returns:**
- If `z` and `T` are scalars: Returns a tuple `(s_d_c, s_d_e)` where `s_d_c` is the differential stress in compression and `s_d_e` in extension.
- If `z` and `T` are arrays: Returns a tuple `(dsigma, depths)` where `dsigma` is a concatenated array of compression and extension stresses, and `depths` is the corresponding array of depths.


---

*For more details and examples, see the docstrings in the source code or the README usage section.*
