# Rheolopy

[![PyPI version](https://badge.fury.io/py/rheolopy.svg)](https://badge.fury.io/py/rheolopy)
[![CI Tests](https://github.com/tmay/rheolopy/actions/workflows/tests.yml/badge.svg)](https://github.com/tmay/rheolopy/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

**Rheolopy** is a robust Python package for geological rheology modeling and visualization. It provides tools to construct layered lithospheric models, calculate depth-dependent yield strengths, and compute effective viscosities.

Rheology describes how rocks deform under stress. In the lithosphere, this behavior transitions from brittle failure at shallow depths to ductile creep at greater depths, heavily influenced by temperature, pressure, and strain rate. `rheolopy` programmatically models these transitions by combining experimentally derived flow laws to compute **Yield Strength Envelopes (YSE)**, giving geophysicists a direct window into lithospheric strength profiles.

---

## Key Features

* **3D Geological Modeling:** Construct complex layered models from depth arrays or 3D structural CSV grids.
* **Extensible Material Database:** Includes standard rock types (e.g., wet/dry Olivine, Quartzite, Diabase) loaded seamlessly from JSON.
* **Multi-Mechanism Physics:** Models Brittle Failure (Byerlee's Law) and Ductile Flow (Dislocation, Diffusion, and Peierls creep).
* **Yield Strength Envelopes:** Built-in visualization tools to generate publication-ready YSE plots.
* **Batch Processing:** Scalable pipeline for computing rheology across massive geological datasets (see `process_rheology/`).
* **Reproducible Pipelines:** Fully configurable via `.ini` files for version-controlled experimental setups.

---

## Installation

Install the latest release directly from PyPI:

```sh
pip install rheolopy
```

For development, including testing and documentation dependencies:

```sh
git clone https://github.com/tmay/rheolopy.git
cd rheolopy
pip install -e .[dev]
```

---

## Quick Start

### Python API

You can easily compute rheological properties using the programmatic API. Here is a minimal example calculating the differential yield strength of dry Olivine at a depth of 20 km:

```python
import rheolopy
from rheolopy.materials import get_material_by_id, materials

# 1. Load the standard material database
mats = materials()
olivine = get_material_by_id(mats, "olivine_hirth_dry")

# 2. Define geological conditions
depth = 20000.0       # meters (20 km)
temp = 800.0          # Kelvin
strain_rate = 1e-15   # 1/s

# 3. Compute differential stress (Compression and Extension)
dsigma_c, dsigma_e = rheolopy.compute_dsigma(olivine, depth, temp, strain_rate)

print(f"Compression Yield Strength: {dsigma_c / 1e6:.1f} MPa")
print(f"Extension Yield Strength: {dsigma_e / 1e6:.1f} MPa")
```

### Configuration-Driven Models

For reproducible 3D models, `rheolopy` can read an entire lithospheric setup from a configuration file. A standard `config.ini` might look like:

```ini
[General]
strain_rate = 1e-16
geotherm = "geotherm.csv"
database = "database.json"

[Model]
layers = [
    [15000, "quartzite_hansen_wet"],
    [35000, "diabase_maryland_strong"],
    [100000, "olivine_hirth_dry"]
  ]
```

To load and visualize this model:

```python
from rheolopy.background import load_model

model = load_model("config.ini")
model.plot_slice()
```

---

## Physical Models & References

This package computes the minimum stress required for deformation at a given depth by comparing brittle and ductile failure criteria. We implement the following established geological flow laws:

* **Brittle Failure:** Byerlee's Law (Byerlee, 1978) combined with Anderson's theory of faulting (Sibson, 1974).
* **Dislocation & Diffusion Creep:** Modeled using standard Arrhenius relationships (e.g., Hirth & Kohlstedt, 2003).
* **Peierls Creep:** Low-temperature plasticity activated at high stresses (> 200 MPa), governed by the exponential flow law (Goetze & Evans, 1979).

> **Note:** For a complete mathematical derivation of the implemented flow laws and the exact physical formulas used in this package, please refer to our [Theoretical Background Documentation](docs/theoretical_background.pdf).

---

## Documentation

For a complete list of classes and functions, see the [API Reference](API.md).

Primary modules include:
- `rheolopy.rheolopy`: Core rheology physics and stress calculations.
- `rheolopy.background`: 3D layered background models (`BackgroundModel`).
- `rheolopy.materials`: Material definitions and JSON database loading.
- `rheolopy.geotherm`: Interpolation and handling of thermal profiles.
- `rheolopy.model_plot`: Visualization tools for models and Yield Strength Envelopes.

### Related Tools
For an interactive usage overview, check out **[Rheology Explorer](https://git.gfz-potsdam.de/tmay/rheology_explorer)**. It is a visualization tool built to compare different laboratory-derived material parameters and highlight this package's usability for 3D rheologic investigations.

---

## Requirements

- Python 3.10+
- `numpy`
- `pandas`
- `matplotlib`
- `cmcrameri`

---

## License & Citation

**License:** MIT License. See [LICENSE](LICENSE) for details.

**Citing:** If you use this package in your research, please cite it:
> *Citation will be provided upon publication.*

---

## Contributing & Support

Contributions, bug reports, and feature requests are welcome!
Please open an issue or submit a pull request on the repository.

For questions or direct support, contact [Tilman May](mailto:tmay@gfz.de) at the GFZ Helmholtz Centre for Geosciences.

<br>
<p align="center">
  <img src="docs/images/cat.png" alt="cat" width="150">
</p>