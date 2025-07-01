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
        convert=False,
    ):
        self.id = id
        self.source = source
        self.type = type
        self.fc_e = fc_e
        self.fc_c = fc_c
        self.lambda_pore = lambda_pore
        self.rho_b = rho_b
        # Preexponential scaling factor in Pa^(-n)/s
        if convert == "MPa":
            self.a_disloc = self.aToPa(a_disloc, n, 6)
        elif convert == "GPa":
            self.a_disloc = self.aToPa(a_disloc, n, 9)
        else:
            self.a_disloc = a_disloc
        self.n = n
        self.q_disloc = q_disloc
        if convert == "MPa":
            self.a_diff = self.aToPa(a_diff, n, 6)
        elif convert == "GPa":
            self.a_diff = self.aToPa(a_diff, n, 9)
        else:
            self.a_diff = a_diff
        self.q_diff = q_diff
        self.d = d
        self.m = m

    def aToPa(self, A, n, u):
        A = float(A)
        n = float(n)
        return A * 10.0 ** (-1.0 * n * u)

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

    if not os.path.isabs(database_path):
        # Try to resolve relative to project root
        database_path = os.path.join(os.path.dirname(__file__), database_path)
    data = pd.read_json(database_path)

    r = list()

    for _, row in data.iterrows():
        # Pass all values in the row as arguments to Material
        obj = Material(**row)
        r.append(obj)

    # Sort materials by their name (material.type)
    r.sort(key=lambda mat: getattr(mat, "type", ""))
    return r


def get_material_by_id(mats, material_id):
    material_dict = {mat.id: mat for mat in mats}
    if material_id in material_dict:
        found_material = material_dict[material_id]
        return found_material
    else:
        return None
