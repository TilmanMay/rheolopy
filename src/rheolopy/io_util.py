import configparser
import os
import importlib.resources
import numpy as np


class RheolopyConfig(configparser.ConfigParser):
    """ConfigParser subclass that cleanly tracks the source file path."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.config_path: str = ""


def load_config(path=None) -> RheolopyConfig:
    """
    Load an INI configuration file.
    
    If no path is provided, the bundled default config.ini is loaded.
    
    Parameters
    ----------
    path : str, optional
        Path to the INI config file.
        
    Returns
    -------
    config : RheolopyConfig
        The loaded configuration object with the `config_path` attribute set.
    """
    config = RheolopyConfig()
    if path is None:
        # Load config.ini from the package data
        with importlib.resources.files("rheolopy").joinpath("config.ini").open(
            "r"
        ) as f:
            config.read_file(f)
        config.config_path = os.path.abspath(
            importlib.resources.files("rheolopy").joinpath("config.ini")
        )
    else:
        if not os.path.isabs(path) and not os.path.exists(path):
            # Try to load from package data if not found in filesystem
            try:
                with importlib.resources.files("rheolopy").joinpath(path).open(
                    "r"
                ) as f:
                    config.read_file(f)
                config.config_path = os.path.abspath(
                    importlib.resources.files("rheolopy").joinpath(path)
                )
                return config
            except FileNotFoundError:
                pass
        config.read(path)
        config.config_path = os.path.abspath(path)
    return config


def resolve_path(config, option, section="General"):
    """
    Resolve a relative path option from the configuration file.
    
    Paths are resolved relative to the directory containing the config file.
    
    Parameters
    ----------
    config : RheolopyConfig
        The configuration object (must have `config_path` set by `load_config`).
    option : str
        The option key to resolve.
    section : str, optional
        The section containing the option. Default is "General".
        
    Returns
    -------
    abs_path : str
        The absolute path.
    """
    rel_path = config.get(section, option)
    # Strip potential quotes that might cause path issues
    rel_path = rel_path.strip("\"'")
    config_dir = os.path.dirname(config.config_path)
    abs_path = os.path.abspath(os.path.join(config_dir, rel_path))
    return abs_path
