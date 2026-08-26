"""
rheolopy — A Python package for geological rheology modeling and visualization.

Provides tools to construct, analyze, and visualize layered lithospheric models,
compute yield strength envelopes (YSE), and calculate effective viscosities using
experimentally derived flow laws for brittle failure, dislocation creep, diffusion
creep, and Peierls creep.
"""

from .background import (
    BackgroundModel,
    Layer3D,
    load_model,
)
from .io_util import (
    load_config,
    resolve_path,
)
from .materials import (
    Material,
    materials,
    get_material_by_id,
)
from .model_plot import (
    plot_yse,
)
from .core import (
    sigma_byerlee,
    eta_dislocation,
    eta_diffusion,
    eta_effective,
    calc_peierls,
    sigma_d,
    compute_dsigma,
)

from .geotherm import Geotherm

__version__ = "0.1.0"

__all__ = [
    "BackgroundModel",
    "Layer3D",
    "load_model",
    "load_config",
    "resolve_path",
    "Material",
    "materials",
    "get_material_by_id",
    "plot_yse",
    "sigma_byerlee",
    "eta_dislocation",
    "eta_diffusion",
    "eta_effective",
    "calc_peierls",
    "sigma_d",
    "compute_dsigma",
    "Geotherm",
]
