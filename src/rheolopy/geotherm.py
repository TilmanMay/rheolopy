import numpy as np
import os
import importlib.resources
from .io_util import load_config
from typing import Optional, Union
import configparser


class Geotherm:
    """
    Represents a geotherm (temperature-depth profile).
    Can be loaded from CSV (filesystem or package data).
    Provides interpolation and unit conversion.
    """

    def __init__(self, path: Optional[str] = None, config: Optional[configparser.ConfigParser] = None):
        """
        Load geotherm from a CSV file.
        If path is None, tries to resolve from config.
        If config is None, uses the default config.
        """
        if config is None:
            config = load_config()
        if path is None:
            path = config.get("General", "geotherm")
            config_dir = os.path.dirname(config.config_path)
            path = os.path.abspath(os.path.join(config_dir, path))

        self.path = path
        self.data = self._load_geotherm(path)

    def _load_geotherm(self, path: str) -> np.ndarray:
        # Try filesystem first
        if os.path.exists(path):
            data = np.loadtxt(path, skiprows=1, delimiter=",")
        else:
            # Try package data
            try:
                with importlib.resources.files("rheolopy").joinpath(path).open(
                    "r"
                ) as f:
                    data = np.loadtxt(f, skiprows=1, delimiter=",")
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"Geotherm file '{path}' not found in filesystem or package data."
                )
        # Convert units: km to m, C to K
        data[:, 0] *= 1000
        data[:, 1] += 273.15
        return data

    def interpolate(self, depths_m: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Interpolate temperature at given depth(s) in meters.
        Returns temperature(s) in Kelvin.
        """
        return np.interp(depths_m, self.data[:, 0], self.data[:, 1])

    def as_array(self) -> np.ndarray:
        """Return the raw geotherm data as a numpy array (depth [m], temp [K])."""
        return self.data

    def __repr__(self) -> str:
        return f"<Geotherm: {len(self.data)} depth points from {self.data[0, 0]:.0f} m to {self.data[-1, 0]:.0f} m>"
