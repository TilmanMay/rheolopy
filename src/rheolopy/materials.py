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

    Peierls creep
    -------------
    eps_peierls : float
        Reference strain rate / 1/s
    sigma_peierls : float
        Critical stress / Pa
    stress_pd : float
        Peierls stress constant / Pa
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
        eps_peierls=None,
        sigma_peierls=None,
        stress_pd=None,
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
        self.a_disloc = a_disloc
        self.a_diff = a_diff
        self.eps_peierls = eps_peierls
        self.sigma_peierls = sigma_peierls
        self.stress_pd = stress_pd

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
