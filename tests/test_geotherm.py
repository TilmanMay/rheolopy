import pytest
import numpy as np
from rheolopy.geotherm import Geotherm

def test_geotherm_initialization_and_interpolation(tmp_path):
    # Create a temporary CSV file for testing
    csv_file = tmp_path / "test_geotherm.csv"
    csv_content = "Depth,Temperature\n0,0\n10,500.0\n50,1000.0"
    csv_file.write_text(csv_content)

    geo = Geotherm(path=str(csv_file))
    
    # Test interpolation at exact points
    temp = geo.interpolate(0)
    assert np.isclose(temp, 273.15)
    
    temp2 = geo.interpolate(-10000)
    assert np.isclose(temp2, 773.15)

    # Test interpolation at intermediate point
    temp_mid = geo.interpolate(-5000)
    assert np.isclose(temp_mid, (273.15 + 773.15) / 2)

def test_geotherm_as_array(tmp_path):
    csv_file = tmp_path / "test_geotherm.csv"
    csv_content = "Depth,Temperature\n0,0\n10,500.0"
    csv_file.write_text(csv_content)

    geo = Geotherm(path=str(csv_file))
    arr = geo.as_array()
    
    assert arr.shape == (2, 2)
    assert arr[0, 0] == 0
    assert arr[0, 1] == 273.15
    assert arr[1, 0] == 10000
    assert arr[1, 1] == 773.15
