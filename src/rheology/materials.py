import numpy as np
import pandas as pd
import os


class Material:
    """
    Contains a list of materials used for strength computation with exodus
    module. Available properties are

    Meta properties
    ---------------
    id : str
        unique id of the material
    source : str
        Data source
    type : str
        What does your material pretend to be ??

    Byerlee's law
    -------------
    fc_e : float
        Friction coefficient for extension
    fc_c : float
        Friction coefficient for compression
    lambda_pore : float
        Pore fluid factor
    rho_b : float
        Bulk density of the rock / kg/m3

    Dislocation creep
    -----------------
    a_disloc : float
        Preexponential scaling factor / Pa^(-n)/s
    n : float
        Power law exponent
    q_disloc : float
        Activation energy / J/mol

    Diffusion creep
    ---------------
    a_diff : float
        Preexponential scaling factor / 1/Pa/s
    q_diff : float
        Activation energy / J/mol
    a : float
        Grain size / m
    m : float
        Grain size exponent

    Dorn's law creep
    ----------------
    sigma_d : float
        Dorn's law stress / Pa
    q_d : float
        Dorn's law activation energy / J/mol
    a_d : float
        Dorn's law strain rate
    """

    def __init__(
        self,
        id,
        source="",
        type="",
        fc_e=0.75,
        fc_c=2.0,
        lambda_pore=0.36,
        rho_b=None,
        a_disloc=None,
        n=None,
        q_disloc=None,
        a_diff=None,
        q_diff=None,
        d=None,
        m=None,
        convert=None,  # now a list of tags as strings
    ):
        self.id = id
        self.source = source
        self.type = type
        self.fc_e = fc_e
        self.fc_c = fc_c
        self.lambda_pore = lambda_pore
        self.rho_b = rho_b
        # Preexponential scaling factor in Pa^(-n)/s
        self.n = n
        self.q_disloc = q_disloc
        self.q_diff = q_diff
        self.d = d
        self.m = m
        # Handle convert as a list of tags
        convert = convert or []
        if not isinstance(convert, list):
            convert = [convert]
        # Dislocation creep
        if "MPa" in convert:
            self.a_disloc = self.aToPa(a_disloc, n, 6)
        elif "GPa" in convert:
            self.a_disloc = self.aToPa(a_disloc, n, 9)
        else:
            self.a_disloc = a_disloc
        # Diffusion creep
        if "MPa" in convert:
            self.a_diff = self.aToPa(a_diff, 1, 6, m, convert)
        elif "GPa" in convert:
            self.a_diff = self.aToPa(a_diff, 1, 9, m, convert)
        else:
            self.a_diff = self.aToPa(a_diff, 1, 0, m, convert)

    def aToPa(self, A, n, u, m=0, convert=None):
        """
        Convert A from literature units to SI (Pa, m) based on stress and grain size units.
        - A: original A value
        - n: stress exponent
        - u: stress unit factor (6 for MPa, 9 for GPa)
        - m: grain size exponent (default 0 if not applicable)
        - convert: list of tags (e.g., ["MPa", "um"])
        """
        if A is None:
            return None
        A = float(A)
        n = float(n)
        m = float(m)
        # Stress unit conversion
        stress_conversion = 10.0 ** (-n * u)
        # Grain size unit conversion
        grain_conversion = 1.0
        if convert and "um" in convert:
            grain_conversion = (10.0**-6) ** m
        elif convert and "mm" in convert:
            grain_conversion = (10.0**-3) ** m
        elif convert and "nm" in convert:
            grain_conversion = (10.0**-9) ** m
        # Default is meters (no conversion)
        return A * stress_conversion * grain_conversion

    def get_attributes(self):
        """Return a list of attribute names that are not NaN."""
        non_nan_attributes = []
        for attr, value in vars(self).items():
            if isinstance(value, (int, float)) and not np.isnan(value):
                non_nan_attributes.append(attr)
        return non_nan_attributes

    def __repr__(self):
        return f"\n{self.type}, {self.source}"


def materials(database_path="database.json"):
    """
    Load materials from a JSON database file.
    If database_path is not absolute, first try to resolve relative to the current working directory.
    If not found, try to resolve relative to the project/module directory.
    """
    import os
    import pandas as pd

    # Try current working directory first
    if not os.path.isabs(database_path):
        if os.path.exists(database_path):
            resolved_path = database_path
        else:
            # Try to resolve relative to the module/project directory
            resolved_path = os.path.join(os.path.dirname(__file__), database_path)
    else:
        resolved_path = database_path

    if not os.path.exists(resolved_path):
        raise FileNotFoundError(f"Material database not found: {database_path}")

    data = pd.read_json(resolved_path)
    r = list()
    for _, row in data.iterrows():
        obj = Material(**row)
        r.append(obj)
    r.sort(key=lambda mat: getattr(mat, "type", ""))
    return r


def get_material_by_id(mats, material_id):
    material_dict = {mat.id: mat for mat in mats}
    if material_id in material_dict:
        found_material = material_dict[material_id]
        return found_material
    else:
        return None
