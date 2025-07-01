import numpy as np

from .materials import Material
from .background import BackgroundModel


R = 8.314472  # m2kg/s2/K/mol
g = 9.81  # m/s2


def sigma_byerlee(material, z, mode):
    """
    Compute the differential stress using Byerlee's law for a given material, depth, and mode.

    Parameters
    ----------
    material : Material
        Material object with properties 'fc_e', 'fc_c', 'lambda_p', and 'rho_b'.
    z : float
        Depth below surface in meters (z > 0).
    mode : str
        'compression' or 'extension'.

    Returns
    -------
    sigma_d : float
        Differential stress in Pascals (Pa).
    """
    if mode == "compression":
        f_f = material.fc_c
    elif mode == "extension":
        f_f = material.fc_e
    else:
        raise ValueError("Invalid parameter for mode:", mode)
    lambda_pore = material.lambda_pore
    rho_b = material.rho_b

    return f_f * rho_b * g * z * (1.0 - lambda_pore)


def eta_dislocation(material, temp, strain_rate):
    """
    Compute the viscosity for dislocation creep at specified temperature and strain rate.

    Parameters
    ----------
    material : Material
        Material object with properties 'a_disloc', 'n', 'q_disloc', and optionally 'd' (grain size).
    temp : float
        Temperature in Kelvin.
    strain_rate : float
        Reference strain rate in 1/s.

    Returns
    -------
    eta_dislocation : float
        Dislocation creep viscosity in Pa·s.
    """
    a_disloc = material.a_disloc
    n = material.n
    d = material.d  # grain size
    m = 0  # dislocation creep does not depend on grain size
    q_disloc = material.q_disloc
    F_2 = 1 / (2 ** ((n - 1) / n) * 3 ** ((n + 1) / (2 * n)))

    if a_disloc is None:
        return np.nan

    return (
        F_2
        * 1
        / ((a_disloc) ** (1 / n) * d ** (-(m) / (n)) * strain_rate ** ((n - 1) / (n)))
        * np.exp(q_disloc / n / R / temp)
    )


def eta_diffusion(material, temp, strain_rate):
    """
    Compute the viscosity for diffusion creep at specified temperature and strain rate.

    Parameters
    ----------
    material : Material
        Material object with properties 'd', 'm', 'a_diff', and 'q_diff'.
    temp : float
        Temperature in Kelvin.
    strain_rate : float
        Reference strain rate in 1/s.

    Returns
    -------
    eta_diffusion : float
        Diffusion creep viscosity in Pa·s.
    """
    d = material.d
    m = material.m
    a_diff = material.a_diff
    q_diff = material.q_diff
    n = 1

    F_2 = 3 / (2 ** ((n - 1) / n) * 3 ** ((n + 1) / (2 * n)))

    if a_diff is None:
        return np.nan
    else:
        return (
            F_2
            * 1
            / ((a_diff) ** (1 / n) * d ** (-(m) / (n)) * strain_rate ** ((n - 1) / (n)))
            * np.exp(q_diff / n / R / temp)
        )


def eta_effective(material, temp, strain_rate):
    """
    Compute the effective viscosity for a material at specified temperature and strain rate.
    Uses the harmonic mean of dislocation and diffusion creep viscosities.

    Parameters
    ----------
    material : Material
        Material object with required properties for creep laws.
    temp : float
        Temperature in Kelvin.
    strain_rate : float
        Reference strain rate in 1/s.

    Returns
    -------
    eta_eff : float
        Effective viscosity in Pa·s.
    eta_dis : float
        Dislocation creep viscosity in Pa·s.
    eta_diff : float
        Diffusion creep viscosity in Pa·s.
    """
    eta_dis = eta_dislocation(material, temp, strain_rate)
    eta_diff = eta_diffusion(material, temp, strain_rate)

    if not np.isnan(eta_dis) and not np.isnan(eta_diff):
        eta_eff = 1 / (1 / eta_dis + 1 / eta_diff)
    elif not np.isnan(eta_dis):
        eta_eff = eta_dis
    elif not np.isnan(eta_diff):
        eta_eff = eta_diff
    else:
        eta_eff = np.nan
    return eta_eff, eta_dis, eta_diff


def calc_peierls(material, temp, strain_rate):
    """
    Compute the Peierls stress for a material at given temperature and strain rate.

    Parameters
    ----------
    material : Material
        Material object with properties 'a_disloc', 'n', 'q_disloc'.
    temp : float
        Temperature in Kelvin.
    strain_rate : float
        Reference strain rate in 1/s.

    Returns
    -------
    peierls : float
        Peierls stress in Pascals (Pa).
    """
    epsPeierls = 5.7e11
    sigmaPeierls = 8.5e9
    stressPD = 200e6

    a_disloc = material.a_disloc
    n = material.n
    q_disloc = material.q_disloc

    termln = np.log(stressPD * (a_disloc / strain_rate) ** (1.0 / n))
    TPDL = q_disloc / (n * R * termln)
    QPeierls = (
        R
        * TPDL
        * np.log(epsPeierls / strain_rate)
        * ((1.0 - stressPD / sigmaPeierls) ** -2.0)
    )
    peierls = sigmaPeierls * (
        1.0 - np.sqrt((R * temp / QPeierls) * np.log(epsPeierls / strain_rate))
    )

    return peierls if peierls > 0 else 0.0


def sigma_d(
    material,
    z,
    temp,
    strain_rate=None,
    compute=None,
    mode=None,
    return_all=False,
    return_index=False,
):
    """
    Compute the differential stress for a material at a given depth, temperature,
    and strain rate. Returns the minimum stress required for deformation, considering
    Byerlee's law, dislocation creep, diffusion creep, and (optionally) Peierls creep.

    Parameters
    ----------
    material : Material
        Material object containing properties required for calculations.
    z : float
        Positive depth in meters below the surface (z > 0).
    temp : float
        Temperature in Kelvin.
    strain_rate : float, optional
        Reference strain rate in 1/s. If None, defaults to 1e-17 1/s.
    compute : list of str, optional
        List of processes to compute: 'dislocation', 'diffusion', 'peierls'.
        If None, uses ['diffusion', 'dislocation', 'peierls'] for olivine,
        otherwise ['diffusion', 'dislocation'].
    mode : str, optional
        'compression' or 'extension'.
    return_all : bool, optional
        If True, return all computed stresses (Byerlee, creep, diffusion, dislocation).
    return_index : bool, optional
        If True, return the index of the controlling mechanism.

    Returns
    -------
    Sigma : float or tuple
        Differential stress in Pa, or tuple of stresses if return_all is True.
    """
    if z < 0:
        raise ValueError("Depth must be positive. Got z =", z)
    if strain_rate is None:
        e_prime = 1e-17
    else:
        e_prime = strain_rate
    if "olivine" in material.type.lower():
        compute_default = ["diffusion", "dislocation", "peierls"]
    else:
        compute_default = ["diffusion", "dislocation"]
    if compute is None:
        compute = compute_default
    else:
        # Check the keywords
        for kwd in compute:
            if kwd not in compute_default:
                raise ValueError("Unknown compute keyword", kwd)

    s_byerlee = sigma_byerlee(material, z, mode)

    proplist = material.get_attributes()
    eta_eff, eta_dis, eta_diff = eta_effective(material, temp, e_prime)
    creep = 2 * eta_eff * e_prime
    s_disloc = 2 * eta_dis * e_prime
    s_diff = 2 * eta_diff * e_prime

    if "peierls" in compute and "a_disloc" in proplist:
        s_peierls = calc_peierls(material, temp, e_prime)
    else:
        s_peierls = 0

    if (creep > 200e6) and (s_peierls > 0) and (s_peierls < creep):
        s_creep = s_peierls
    else:
        s_creep = creep
    if return_all:
        return s_byerlee, s_creep, s_diff, s_disloc
    if return_index:
        if (creep > 200e6) and (s_peierls > 0) and (s_peierls < creep):
            return np.nanargmin([s_byerlee, s_peierls, s_diff, s_disloc])
        else:
            return np.nanargmin([s_byerlee, np.nan, s_diff, s_disloc])
    else:
        # print("s_byerlee, s_creep, s_diff", s_byerlee, s_creep, s_diff)
        return min(s_byerlee, s_creep, s_diff, s_disloc)


def compute_dsigma(
    background,
    z,
    T,
    strain_rate,
    x_idx=None,
    y_idx=None,
    return_all=False,
    return_index=False,
):
    """
    Compute differential stress for a given material or background model at depths z and
    temperatures T. Computes for both compression and extension and returns arrays.
    Handles both scalar and array input for z and T (no vectorization).

    Parameters
    ----------
    background : Material or BackgroundModel
        Material object or BackgroundModel instance.
    z : float or array-like
        Depth(s) in meters below the surface (z > 0).
    T : float or array-like
        Temperature(s) in Kelvin.
    strain_rate : float
        Reference strain rate in 1/s.
    x_idx : int, optional
        x grid index for BackgroundModel (ignored for Material).
    y_idx : int, optional
        y grid index for BackgroundModel (ignored for Material).
    return_all : bool, optional
        If True, return all computed stresses (Byerlee, creep, diffusion, dislocation).
    return_index : bool, optional
        If True, return the index of the controlling mechanism.

    Returns
    -------
    dsigma : tuple or ndarray
        For scalar input: (compression, extension) stress in Pa.
        For array input: (dsigma, depths) where dsigma is a concatenated array of
        compression and extension stresses, and depths is the corresponding array of depths.
    """
    if isinstance(background, BackgroundModel):
        # Use the correct (x_idx, y_idx) for the top boundary
        x, y = background._xy_index_to_value(x_idx, y_idx)
        top_layer = background.layers[0]
        top_z = top_layer.data[(top_layer.data["x"] == x) & (top_layer.data["y"] == y)][
            "z"
        ].values[0]
        z_for_model = -(top_z - z)
    else:
        z_for_model = z

    if np.isscalar(z):
        if isinstance(background, Material):
            mat = background
        elif isinstance(background, BackgroundModel):
            mat = background.get_material_at(x_idx=x_idx, y_idx=y_idx, z=-z_for_model)
        else:
            raise ValueError("background must be Material or BackgroundModel")
        if mat is None:
            return np.nan, np.nan
        s_d_c = -1 * sigma_d(
            mat,
            z,
            T,
            strain_rate=strain_rate,
            mode="compression",
            return_all=return_all,
            return_index=return_index,
        )
        s_d_e = sigma_d(
            mat,
            z,
            T,
            strain_rate=strain_rate,
            mode="extension",
            return_all=return_all,
            return_index=return_index,
        )
        return s_d_c, s_d_e
    else:
        s_d_c = np.empty_like(z)
        s_d_e = np.empty_like(z)
        for i in range(len(z)):
            if isinstance(background, Material):
                mat = background
                z_model = z[i]
            elif isinstance(background, BackgroundModel):
                z_model = z_for_model[i]
                mat = background.get_material_at(x_idx=x_idx, y_idx=y_idx, z=-z_model)
            else:
                raise ValueError("background must be Material or BackgroundModel")
            if mat is None:
                s_d_c[i] = np.nan
                s_d_e[i] = np.nan
            else:
                s_d_c[i] = -1 * sigma_d(
                    mat,
                    z[i],
                    T[i],
                    strain_rate=strain_rate,
                    mode="compression",
                    return_all=return_all,
                    return_index=return_index,
                )
                s_d_e[i] = sigma_d(
                    mat,
                    z[i],
                    T[i],
                    strain_rate=strain_rate,
                    mode="extension",
                    return_all=return_all,
                    return_index=return_index,
                )
        dsigma = np.concatenate((s_d_c, s_d_e[::-1]))
        depths = np.concatenate((z, z[::-1]))
        return dsigma, depths
