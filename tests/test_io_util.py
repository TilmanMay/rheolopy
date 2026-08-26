import pytest
import os
from rheolopy.io_util import load_config, resolve_path, RheolopyConfig

def test_load_config_default():
    config = load_config()
    assert isinstance(config, RheolopyConfig)
    assert hasattr(config, "config_path")
    assert "config.ini" in config.config_path
    assert config.has_section("General")

def test_load_config_custom(tmp_path):
    config_file = tmp_path / "custom.ini"
    config_file.write_text("[General]\ngeotherm = custom.csv\n")
    
    config = load_config(str(config_file))
    assert isinstance(config, RheolopyConfig)
    assert config.get("General", "geotherm") == "custom.csv"
    assert config.config_path == str(config_file.absolute())

def test_resolve_path(tmp_path):
    config_file = tmp_path / "custom.ini"
    config_file.write_text("[General]\nfile = data.csv\n")
    
    config = load_config(str(config_file))
    resolved = resolve_path(config, "file")
    
    expected = str(tmp_path / "data.csv")
    assert resolved == expected
