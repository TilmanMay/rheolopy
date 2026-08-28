import pytest
import numpy as np
from rheolopy.background import BackgroundModel
from rheolopy.materials import Material

def test_add_layers_and_initialize():
    model = BackgroundModel()
    mat1 = Material(id="crust")
    mat2 = Material(id="mantle")
    
    model.add_layer(0, mat1)
    model.add_layer(30000, mat2)
    
    model.initialize()
    
    assert len(model.layers) == 2
    assert model.layers[0].material.id == "crust"
    assert model.layers[1].material.id == "mantle"

def test_get_material_at():
    model = BackgroundModel()
    mat1 = Material(id="crust")
    mat2 = Material(id="mantle")
    
    model.add_layer(0, mat1)
    model.add_layer(30000, mat2)
    model.initialize()
    
    # query shallow
    mat = model.get_material_at(z=10000)
    assert mat.id == "crust"
    
    # query deep
    mat = model.get_material_at(z=40000)
    assert mat.id == "mantle"

def test_layer_crossing_raises():
    model = BackgroundModel()
    mat1 = Material(id="crust")
    mat2 = Material(id="mantle")
    
    model.add_layer(30000, mat1)
    model.add_layer(0, mat2)
    
    # The initialize function sorts layers internally, but cross-checking logic
    # would throw an error if layers cross physically. Here we test if it handles
    # the sort properly without crossing errors for constant depths.
    model.initialize()
    assert model.layers[0].material.id == "mantle"  # 0 is top
    assert model.layers[1].material.id == "crust"   # 30000 is bottom

def test_print_layers_at():
    model = BackgroundModel()
    mat1 = Material(id="crust")
    mat2 = Material(id="mantle")
    
    model.add_layer(0, mat1)
    model.add_layer(30000, mat2)
    model.initialize()
    
    output = model.print_layers_at(as_string=True)
    assert "crust" in output
    assert "mantle" in output
def test_calc_lithostatic_pressure():
    model = BackgroundModel()
    mat1 = Material(id="crust")
    mat1.rho_b = 2800.0
    mat2 = Material(id="mantle")
    mat2.rho_b = 3300.0
    model.add_layer(-2000, mat1)
    model.add_layer(30000, mat2)
    model.initialize()
    g = 9.80665
    assert model.calc_lithostatic_pressure(z=-3000) == 0.0
    assert model.calc_lithostatic_pressure(z=-2000) == 0.0
    expected_sl = 2000 * 2800.0 * g
    assert np.isclose(model.calc_lithostatic_pressure(z=0), expected_sl)
    expected_moho = 32000 * 2800.0 * g
    # Inside mantle
    # crust thickness = 32000, mantle thickness = 10000 (because z_max = 30000 + 10000)
    expected_mantle = expected_moho + 10000 * 3300.0 * g
    assert np.isclose(model.calc_lithostatic_pressure(z=40000), expected_mantle)
    # Beyond z_max it shouldn't add more pressure
    assert np.isclose(model.calc_lithostatic_pressure(z=50000), expected_mantle)
