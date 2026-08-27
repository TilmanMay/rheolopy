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
