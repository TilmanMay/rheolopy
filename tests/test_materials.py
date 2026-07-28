import pytest
import os
from rheolopy.materials import Material, materials, get_material_by_id

def test_material_initialization():
    mat = Material(id="test_id", type="TestRock", rho_b=2500)
    assert mat.id == "test_id"
    assert mat.type == "TestRock"
    assert mat.rho_b == 2500
    assert mat.fc_e == 0.75 # default value

def test_get_attributes():
    mat = Material(id="test_id", rho_b=2500, a_disloc=1e-15)
    attrs = mat.get_attributes()
    assert "rho_b" in attrs
    assert "a_disloc" in attrs
    assert "q_diff" not in attrs

def test_materials_load():
    # Attempt to load the default database
    mats = materials("database.json")
    assert isinstance(mats, list)
    assert len(mats) > 0
    assert isinstance(mats[0], Material)

def test_get_material_by_id():
    mats = [
        Material(id="rock1", type="Rock One"),
        Material(id="rock2", type="Rock Two")
    ]
    mat = get_material_by_id(mats, "rock2")
    assert mat is not None
    assert mat.id == "rock2"
    assert mat.type == "Rock Two"

    mat_none = get_material_by_id(mats, "nonexistent")
    assert mat_none is None
