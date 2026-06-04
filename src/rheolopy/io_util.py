import configparser
import os
import importlib.resources
import numpy as np


def load_config(path=None):
    config = configparser.ConfigParser()
    if path is None:
        # Load config.ini from the package data
        with importlib.resources.files("rheolopy").joinpath("config.ini").open(
            "r"
        ) as f:
            config.read_file(f)
        config._path = os.path.abspath(
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
                config._path = os.path.abspath(
                    importlib.resources.files("rheolopy").joinpath(path)
                )
                return config
            except FileNotFoundError:
                pass
        config.read(path)
        config._path = os.path.abspath(path)
    return config


def resolve_path(config, option, section="General"):
    rel_path = config.get(section, option)
    config_dir = os.path.dirname(config._path)
    abs_path = os.path.abspath(os.path.join(config_dir, rel_path))
    return abs_path
