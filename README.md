# Rheology

**Rheolopy** is a Python package for geological rheology modeling and visualization. It provides tools to construct, analyze, and visualize layered lithospheric models, compute yield strength envelopes, and more.

Rheology describes how rocks deform under stress, accounting for both brittle and ductile behaviors depending on temperature, pressure, and strain rate. In this code, rheological properties are calculated using experimentally derived flow laws for different deformation mechanisms (such as brittle failure, diffusion creep, and dislocation creep). The code combines these laws to compute yield strength profiles and effective viscosities, enabling direct comparison of rock strength and deformation across varying geological conditions.

---

## Features

- Layered 3D models with flexible material assignment
- Material database with common rock types and rheological parameters
- Configurable via INI files for reproducible model setups
- Yield Strength Envelope (YSE) plotting
- Support for custom geotherms and material databases
- Easy integration with Jupyter and scripts

---

## Installation

Install from PyPI:

```sh
pip install rheolopy
```


---

## Quick Start

### 1. Module configuration

The provided `config.ini` comes with values for strain rate as well as the material database and a geotherm (McKenzie et al. 2005):

```ini
[General]
strain_rate = 1e-17
geotherm = "geotherm.csv"
database = "database.json"
```

## API Reference

For a complete list of classes and functions, see the [API documentation](API.md).

Here are the main entry points:

- `BackgroundModel`: Build and query layered models.
- `load_model(config_path)`: Load a model from a config file.
- `plot_yse(model, ...)`: Plot yield strength envelopes.
- `plot_layer_thickness()`, `plot_slice()`: Visualize model structure.
- `compute_dsigma(...)`, `sigma_byerlee(...)`, etc.: Rheology calculations.

---

### 2. Usage

For some usage overview go take a look here `Rheology Explorer` (https://git.gfz-potsdam.de/tmay/rheology_explorer). It is a tool to compare different material parameters which are derived under laboratory conditions. Also it highlights this package useability for simple 3D rheologic investigations. This will most likely help you :)

## Requirements

- Python 3.8+
- numpy
- pandas
- matplotlib
- cmcrameri

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Citing

If you use this package in your research, consider citing it :)

> Tilman May (2025). Don't worry, something is coming soon.

---

## Contributing

Contributions, bug reports, and feature requests are welcome!
Please open an issue or submit a pull request.

---

## Contact

For questions or support, contact [tmay@gfz.de].

<img src="cat.png" alt=":/" >