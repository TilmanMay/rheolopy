"""
Process an Exodus thermal model output for rheological calculations.

Supports complex 3D meshes (subduction zones, non-layer-cake geometries)
by using true element-block connectivity to assign materials at every node.

Workflow:
  1. Read configuration (keyword-based INI)
  2. Read Exodus file -> coordinates, temperature, per-element density,
     and element-block connectivity
  3. Build node -> block_id map from connectivity arrays
  4. Build per-node density from element densities
  5. Read rheology property file -> {block_name: material + strain_rate}
  6. For each 2D column (x, y):
       - Extract vertical profile z(iz), T(iz), block_ids(iz), density(iz)
       - Compute YSE on constant-resolution grid AND/OR original resolution
       - Write per-column integrated values + 3D field data + BDT
  7. Write output CSVs and augmented Exodus file

Uses the rheolopy package for physics (Byerlee's law, dislocation/diffusion
creep viscosity, Peierls creep).

Usage:
    python process_exodus.py [config_file]
    Default config: process_exodus.ini (in same directory as script)
"""

import os
import sys
import numpy as np
import netCDF4 as nc
from datetime import datetime
from multiprocessing import Pool, cpu_count

# ---------------------------------------------------------------------------
# Add rheolopy source to path
# ---------------------------------------------------------------------------
_script_dir = os.path.dirname(os.path.abspath(__file__))
_src_dir = os.path.normpath(os.path.join(_script_dir, "..", "src"))
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)

from rheolopy import materials as load_materials, get_material_by_id
from rheolopy.core import sigma_d

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
GRAV = 9.81        # m/s²
KELVIN = 273.15    # offset °C -> K
STRESS_LIM = 20e6  # Pa  – threshold for mechanical thickness


# ============================================================================
#  Data structures
# ============================================================================
class LayerProp:
    """Rheological properties for one geological block.

    Parameters
    ----------
    line : str
        Space-separated line: ``Block_Name Material_ID Strain_Rate [Common_Layers]``
    db_mats : list
        List of Material objects from the rheolopy database.
    """

    def __init__(self, line: str, db_mats: list) -> None:
        # Strip inline comments before parsing
        if "#" in line:
            line = line[:line.index("#")].strip()
        f = line.split()
        self.name: str = f[0]
        self.material_id: str = f[1]
        self.strain_rate: float = float(f[2])
        # Common_Layers is parsed for backward compat but ignored in block-aware mode
        self.common_layers: int = int(f[3]) if len(f) > 3 else 1

        self.material = get_material_by_id(db_mats, self.material_id)
        if self.material is None:
            raise ValueError(
                f"Material ID '{self.material_id}' not found in rheolopy database."
            )

        # Fallback density from the material database
        self.density = self.material.rho_b
        if self.density is None or np.isnan(self.density):
            self.density = 2700.0  # Safe fallback


class Config:
    """Parsed configuration from the keyword-based INI file."""

    def __init__(self) -> None:
        self.input_file: str = ""
        self.nx: int = 0          # 0 = auto-detect
        self.ny: int = 0          # 0 = auto-detect
        self.rheo_file: str = ""
        self.eta_low: float = 1e18
        self.eta_up: float = 1e25
        self.out_dir: str = "output"
        self.out_name: str = ""   # empty = derived from input filename
        self.resolution: float = 500.0  # metres
        self.n_workers: int = 0   # 0 = auto (cpu_count)


# ============================================================================
#  Readers
# ============================================================================
def read_config(filename: str) -> Config:
    """Read keyword-based configuration file.

    Parameters
    ----------
    filename : str
        Path to the configuration .ini file.

    Returns
    -------
    Config
        Parsed configuration object.
    """
    cfg = Config()
    with open(filename, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith("#INPUT_FILE:"):
                parts = line.split()
                cfg.input_file = parts[2]
                # nx/ny are now optional
                if len(parts) >= 5:
                    cfg.nx = int(parts[3])
                    cfg.ny = int(parts[4])
            elif line.startswith("#RHEO_FILE:"):
                cfg.rheo_file = line.split(None, 1)[1]
            elif line.startswith("#ETA_BOUNDS:"):
                parts = line.split()
                cfg.eta_low = float(parts[1])
                cfg.eta_up = float(parts[2])
            elif line.startswith("#OUT_DIR:"):
                parts = line.split()
                cfg.out_dir = parts[1]
                if len(parts) >= 3:
                    cfg.out_name = parts[2]
            elif line.startswith("#RESOLUTION:"):
                cfg.resolution = float(line.split()[1])
            elif line.startswith("#WORKERS:"):
                cfg.n_workers = int(line.split()[1])
    return cfg


# ============================================================================
#  Exodus reading and block mapping
# ============================================================================
def read_exodus(filename: str, nx: int = 0, ny: int = 0):
    """Read an Exodus (.e / netCDF) mesh file.

    Parameters
    ----------
    filename : str
        Path to the Exodus file.
    nx, ny : int, optional
        Grid dimensions. If 0, auto-detected from unique (x,y) pairs.

    Returns
    -------
    dict
        Dictionary with keys: x, y, z_all, T_all, num_nodes, nodes2D,
        nsurf, columns, num_el_blk, block_names, ds (open Dataset handle).
    """
    ds = nc.Dataset(filename, "r")

    num_nodes = ds.dimensions["num_nodes"].size
    num_el_blk = ds.dimensions["num_el_blk"].size

    x_all = np.array(ds.variables["coordx"][:], dtype=np.float64)
    y_all = np.array(ds.variables["coordy"][:], dtype=np.float64)
    z_all = np.array(ds.variables["coordz"][:], dtype=np.float64)

    # Temperature: nodal variable 1, last time step
    T_all = np.array(ds.variables["vals_nod_var1"][-1, :], dtype=np.float64)

    # --- Auto-detect nodes2D ------------------------------------------------
    if nx > 0 and ny > 0:
        nodes2D = nx * ny
    else:
        xy = np.column_stack((x_all, y_all))
        nodes2D = len(np.unique(xy, axis=0))
        print(f"    Auto-detected nodes2D = {nodes2D}")

    nsurf = num_nodes // nodes2D
    assert num_nodes == nodes2D * nsurf, (
        f"Node count mismatch: {num_nodes} != {nodes2D} * {nsurf}"
    )

    # --- Read block names ---------------------------------------------------
    block_names: list[str] = []
    if "eb_names" in ds.variables:
        raw_names = ds.variables["eb_names"][:]
    elif "name_elem_blk" in ds.variables:
        raw_names = ds.variables["name_elem_blk"][:]
    else:
        raw_names = None

    for blk in range(num_el_blk):
        if raw_names is not None:
            name = (
                raw_names[blk]
                .tobytes()
                .decode("ascii", errors="ignore")
                .strip()
                .rstrip("\x00")
            )
            if not name:
                name = f"Block_{blk + 1}"
        else:
            name = f"Block_{blk + 1}"
        block_names.append(name)

    # --- Group nodes by (x, y) columns, sorted top-to-bottom ----------------
    xy = np.column_stack((x_all, y_all))
    unique_xy, unique_indices, inverse_indices = np.unique(
        xy, axis=0, return_index=True, return_inverse=True
    )

    # Sort unique points by first appearance to preserve ordering
    order = np.argsort(unique_indices)
    map_back = np.empty_like(order)
    map_back[order] = np.arange(len(order))
    column_mapping = map_back[inverse_indices]

    columns = [[] for _ in range(nodes2D)]
    for i in range(num_nodes):
        columns[column_mapping[i]].append(i)

    x = np.zeros(nodes2D, dtype=np.float64)
    y = np.zeros(nodes2D, dtype=np.float64)

    for i in range(nodes2D):
        col_nodes = np.array(columns[i])
        sorted_idx = np.argsort(-z_all[col_nodes])  # top to bottom
        col_nodes_sorted = col_nodes[sorted_idx]
        columns[i] = col_nodes_sorted
        x[i] = x_all[col_nodes_sorted[0]]
        y[i] = y_all[col_nodes_sorted[0]]

    print(f"  Exodus file read successfully:")
    print(f"    nodes2D = {nodes2D}, nsurf = {nsurf}")
    print(f"    x range: [{x.min():.0f}, {x.max():.0f}] m")
    print(f"    y range: [{y.min():.0f}, {y.max():.0f}] m")
    print(f"    z range: [{z_all.min():.0f}, {z_all.max():.0f}] m")
    print(f"    T range: [{T_all.min():.2f}, {T_all.max():.2f}] °C")
    print(f"    Blocks: {block_names}")

    return {
        "x": x,
        "y": y,
        "z_all": z_all,
        "T_all": T_all,
        "num_nodes": num_nodes,
        "nodes2D": nodes2D,
        "nsurf": nsurf,
        "columns": columns,
        "num_el_blk": num_el_blk,
        "block_names": block_names,
        "ds": ds,  # Keep open for connectivity + density reading
    }


def build_node_block_map(
    ds: nc.Dataset,
    num_nodes: int,
    num_el_blk: int,
    z_all: np.ndarray,
) -> np.ndarray:
    """Build a node -> block_id map using element connectivity.

    For nodes shared between two vertically adjacent blocks, assigns the
    block of the element whose centroid is deeper (lower Z), which is
    physically correct for downward integration of lithostatic pressure.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Open Exodus file handle.
    num_nodes : int
        Total number of nodes.
    num_el_blk : int
        Number of element blocks.
    z_all : np.ndarray
        Z-coordinates of all nodes (metres).

    Returns
    -------
    np.ndarray
        Array of shape ``(num_nodes,)`` with block IDs (1-based).
        0 means unassigned.
    """
    node_block = np.zeros(num_nodes, dtype=np.int32)
    node_block_z = np.full(num_nodes, np.inf, dtype=np.float64)

    for blk in range(1, num_el_blk + 1):
        conn_name = f"connect{blk}"
        if conn_name not in ds.variables:
            continue

        print(f"      Block {blk}: reading connectivity...", end="\r", flush=True)
        conn = ds.variables[conn_name][:]  # shape: (n_elem, nodes_per_elem)
        conn_0 = conn - 1  # Exodus uses 1-based indexing

        # Compute element centroid Z (mean of node Zs) — fully vectorized
        elem_z = np.mean(z_all[conn_0], axis=1)  # shape: (n_elem,)

        # Broadcast centroid Z to all nodes of each element
        # conn_0 is (n_elem, 8), elem_z is (n_elem,)
        flat_nodes = conn_0.ravel()
        flat_centroids = np.repeat(elem_z, conn_0.shape[1])

        # Only update where this element's centroid is deeper (lower Z)
        mask = flat_centroids < node_block_z[flat_nodes]
        update_nodes = flat_nodes[mask]
        update_z = flat_centroids[mask]

        # For nodes updated multiple times in this block, the last write wins.
        # Since we process blocks sequentially and use "<", this is correct.
        node_block[update_nodes] = blk
        node_block_z[update_nodes] = update_z

    unassigned = np.sum(node_block == 0)
    if unassigned > 0:
        print(f"    WARNING: {unassigned} nodes unassigned to any block.")

    return node_block


def build_node_density(
    ds: nc.Dataset,
    num_nodes: int,
    num_el_blk: int,
) -> np.ndarray:
    """Build per-node density by averaging element densities to nodes.

    Parameters
    ----------
    ds : netCDF4.Dataset
        Open Exodus file handle.
    num_nodes : int
        Total number of nodes.
    num_el_blk : int
        Number of element blocks.

    Returns
    -------
    np.ndarray
        Per-node density in kg/m³. NaN where no density data is available.
    """
    # Find which element variable is "density"
    density_var_idx = -1
    if "name_elem_var" in ds.variables:
        names = ds.variables["name_elem_var"][:]
        for i in range(names.shape[0]):
            name_str = (
                names[i]
                .tobytes()
                .decode("ascii", errors="ignore")
                .strip()
                .rstrip("\x00")
            )
            # Use the first variable containing "density" (prefer exact match)
            if name_str.lower() == "density":
                density_var_idx = i + 1
                break
            elif "density" in name_str.lower() and density_var_idx < 0:
                density_var_idx = i + 1

    if density_var_idx < 0:
        print("    No density element variable found in Exodus file.")
        return np.full(num_nodes, np.nan, dtype=np.float64)

    print(f"    Building per-node density from element variable {density_var_idx}")

    node_density_sum = np.zeros(num_nodes, dtype=np.float64)
    node_density_cnt = np.zeros(num_nodes, dtype=np.int32)

    for blk in range(1, num_el_blk + 1):
        var_name = f"vals_elem_var{density_var_idx}eb{blk}"
        conn_name = f"connect{blk}"
        if var_name not in ds.variables or conn_name not in ds.variables:
            continue

        print(f"      Block {blk}: reading density...", end="\r", flush=True)
        elem_rho = ds.variables[var_name][-1, :]  # last time step
        conn = ds.variables[conn_name][:] - 1  # 0-indexed

        # Vectorized: broadcast element density to all its nodes
        flat_nodes = conn.ravel()
        flat_rho = np.repeat(elem_rho, conn.shape[1])

        np.add.at(node_density_sum, flat_nodes, flat_rho)
        np.add.at(node_density_cnt, flat_nodes, 1)

    # Average
    valid = node_density_cnt > 0
    node_density = np.full(num_nodes, np.nan, dtype=np.float64)
    node_density[valid] = node_density_sum[valid] / node_density_cnt[valid]

    n_valid = np.sum(valid)
    print(f"    Per-node density: {n_valid}/{num_nodes} nodes have density data")
    if n_valid > 0:
        print(
            f"    Density range: [{np.nanmin(node_density):.1f}, "
            f"{np.nanmax(node_density):.1f}] kg/m³"
        )

    return node_density


def read_rheology(filename: str, db_mats: list) -> dict[str, LayerProp]:
    """Read rheology file and return {block_name: LayerProp} mapping.

    Parameters
    ----------
    filename : str
        Path to the rheology .ini file.
    db_mats : list
        List of Material objects from the rheolopy database.

    Returns
    -------
    dict[str, LayerProp]
        Dictionary mapping block names to their rheological properties.
    """
    props: dict[str, LayerProp] = {}
    with open(filename, "r") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            lp = LayerProp(line, db_mats)
            props[lp.name] = lp

    print(f"  Rheology file: {len(props)} block-material mappings loaded")
    return props


# ============================================================================
#  YSE computation
# ============================================================================
def interpolate_value(v0: float, v1: float, grad: float) -> float:
    """Linear interpolation between two values."""
    return v0 + (v1 - v0) * grad


def evaluate_point(
    depth_below_surface: float,
    temp_K: float,
    prop: LayerProp,
    density: float,
    eta_low: float,
    eta_up: float,
    p_litho: float,
) -> dict:
    """Evaluate rheological properties at a single depth point.

    Parameters
    ----------
    depth_below_surface : float
        Absolute depth below the surface in metres (positive downward).
    temp_K : float
        Temperature in Kelvin.
    prop : LayerProp
        Material properties for this point.
    density : float
        Density at this point in kg/m³.
    eta_low, eta_up : float
        Viscosity clipping bounds in Pa·s.
    p_litho : float
        Accumulated lithostatic pressure in Pa.

    Returns
    -------
    dict
        Dictionary with keys: rho, brittle, ductile, peierls, yield,
        eta_diff, eta_disl, eta_eff, bdt.
    """
    mat = prop.material
    strain_rate = prop.strain_rate
    depth_z = abs(depth_below_surface)

    brittle, creep, s_diff, s_disl, peierls_val = sigma_d(
        mat, depth_z, temp_K, strain_rate, mode="compression",
        return_all=True, eta_min=eta_low, eta_max=eta_up,
        P_litho=p_litho,
    )

    dsigma = float(np.nanmin([brittle, creep]))
    
    if np.isnan(creep):
        bdt = 0
    elif np.isnan(brittle):
        bdt = 1
    else:
        bdt = 0 if brittle < creep else 1  # 0 = brittle, 1 = ductile

    eta_diff_val = (s_diff / (2.0 * strain_rate)) if not np.isnan(s_diff) else np.nan
    eta_disl_val = (s_disl / (2.0 * strain_rate)) if not np.isnan(s_disl) else np.nan
    eta_eff_val = (creep / (2.0 * strain_rate)) if not np.isnan(creep) else np.nan

    return {
        "rho": density,
        "brittle": brittle,
        "ductile": creep,
        "peierls": peierls_val if peierls_val is not None else np.nan,
        "yield": dsigma,
        "eta_diff": eta_diff_val,
        "eta_disl": eta_disl_val,
        "eta_eff": eta_eff_val,
        "bdt": bdt,
    }


def compute_yse(
    zz: np.ndarray,
    temp_profile: np.ndarray,
    block_ids: np.ndarray,
    densities: np.ndarray,
    block_props: dict[int, LayerProp],
    block_names: list[str],
    resolution: float,
    eta_bounds: tuple[float, float],
    mode: str = "CONSTANT",
) -> dict:
    """Compute the Yield Strength Envelope for a single vertical column.

    Parameters
    ----------
    zz : np.ndarray
        Z-coordinates for this column, sorted top-to-bottom.
    temp_profile : np.ndarray
        Temperature at each node (°C).
    block_ids : np.ndarray
        Block ID (1-based) for each node in this column.
    densities : np.ndarray
        Per-node density (kg/m³) for this column.
    block_props : dict[int, LayerProp]
        Mapping from block_id (1-based) to LayerProp.
    block_names : list[str]
        Block names indexed by (block_id - 1).
    resolution : float
        Vertical step size in metres (for CONSTANT mode).
    eta_bounds : tuple[float, float]
        (eta_low, eta_up) viscosity clipping bounds.
    mode : str
        ``'CONSTANT'`` or ``'ORIGINAL'``.

    Returns
    -------
    dict
        YSE results including field arrays and integrated quantities.
    """
    eta_low, eta_up = eta_bounds
    nsurf = len(zz)
    top = zz[0]
    bottom = zz[-1]

    temp_K = temp_profile + KELVIN

    if mode == "CONSTANT":
        dz = resolution
        nz = int((top - bottom) / dz)
        if nz < 2:
            nz = 2
        out_z = np.array([-iz * dz + top for iz in range(nz)])
    else:
        out_z = zz
        nz = len(out_z)

    # Output arrays
    out_T = np.zeros(nz)
    out_rho = np.zeros(nz)
    out_brittle = np.zeros(nz)
    out_ductile = np.zeros(nz)
    out_peierls = np.zeros(nz)
    out_yield = np.zeros(nz)
    out_eta_diff = np.zeros(nz)
    out_eta_disl = np.zeros(nz)
    out_eta_eff = np.zeros(nz)
    out_bdt = np.zeros(nz, dtype=int)
    out_layer_id = np.zeros(nz, dtype=int)
    out_layer_name = [""] * nz

    # Integrated quantities
    strength = 0.0
    crustal_strength = 0.0
    h_mech = 0.0
    current_h = 0.0
    te_cubed_sum = 0.0
    dsigma_prev = 0.0
    thickness_crust = 0.0
    p_litho = 0.0

    for iz in range(nz):
        zp = out_z[iz]

        # --- Locate position and interpolate temperature/block/density ------
        if mode == "ORIGINAL":
            temp = temp_K[iz]
            blk_id = block_ids[iz]
            rho = densities[iz]
        else:
            # Interpolate between original nodes
            if iz == 0:
                temp = temp_K[0]
                blk_id = block_ids[0]
                rho = densities[0]
            elif iz == nz - 1:
                temp = temp_K[-1]
                blk_id = block_ids[-1]
                rho = densities[-1]
            else:
                found = False
                for isrf in range(nsurf - 1):
                    z0, z1 = zz[isrf], zz[isrf + 1]
                    if (zp - z0) * (zp - z1) <= 0.0:
                        gradient = abs((zp - z0) / (z1 - z0)) if z1 != z0 else 0.0
                        temp = interpolate_value(
                            temp_K[isrf], temp_K[isrf + 1], gradient
                        )
                        rho = interpolate_value(
                            densities[isrf], densities[isrf + 1], gradient
                        )
                        # Use block of the node below (deeper)
                        blk_id = block_ids[isrf + 1]
                        found = True
                        break
                if not found:
                    if zp > zz[0]:
                        temp, blk_id, rho = temp_K[0], block_ids[0], densities[0]
                    else:
                        temp, blk_id, rho = temp_K[-1], block_ids[-1], densities[-1]

        # Look up the material for this block
        prop = block_props.get(blk_id)
        if prop is None:
            # Fallback: skip this node (shouldn't happen with correct .ini)
            continue

        blk_name = block_names[blk_id - 1] if 0 < blk_id <= len(block_names) else "?"
        is_crust = "mantle" not in blk_name.lower()

        # Handle NaN density (fallback to material default)
        if np.isnan(rho):
            rho = prop.density

        # --- Integrated quantities ------------------------------------------
        thickness = 0.0 if iz == 0 else abs(zp - out_z[iz - 1])

        # Lithostatic pressure integral
        p_litho += thickness * rho * GRAV

        # Evaluate rheology at this point
        pt = evaluate_point(
            zp - top, temp, prop, rho, eta_low, eta_up, p_litho
        )

        if pt["yield"] > STRESS_LIM:
            h_mech += thickness
            current_h += thickness
        else:
            if current_h > 0:
                te_cubed_sum += current_h ** 3
                current_h = 0.0

        if iz > 0:
            dsigma_median = 0.5 * (pt["yield"] + dsigma_prev)
            strength += dsigma_median * thickness
            if is_crust:
                crustal_strength += dsigma_median * thickness
                thickness_crust += thickness
        dsigma_prev = pt["yield"]

        # --- Store ----------------------------------------------------------
        out_T[iz] = temp - KELVIN
        out_rho[iz] = pt["rho"]
        out_brittle[iz] = pt["brittle"]
        out_ductile[iz] = pt["ductile"]
        out_peierls[iz] = pt["peierls"]
        out_yield[iz] = pt["yield"]
        out_eta_diff[iz] = pt["eta_diff"]
        out_eta_disl[iz] = pt["eta_disl"]
        out_eta_eff[iz] = pt["eta_eff"]
        out_bdt[iz] = pt["bdt"]
        out_layer_id[iz] = blk_id
        out_layer_name[iz] = blk_name

    # --- Compute integrated viscosities ------------------------------------
    layer_thickness = top - bottom
    nominal_sr = 1e-15
    avg_eta = (
        strength / (2.0 * nominal_sr * layer_thickness)
        if layer_thickness > 0 else eta_low
    )
    avg_eta = np.clip(avg_eta, eta_low, eta_up)

    eta_crust_avg = (
        crustal_strength / (2.0 * nominal_sr * thickness_crust)
        if thickness_crust > 0 else eta_low
    )
    eta_crust_avg = np.clip(eta_crust_avg, eta_low, eta_up)

    mantle_thickness = layer_thickness - thickness_crust
    eta_mantle_avg = (
        (strength - crustal_strength) / (2.0 * nominal_sr * mantle_thickness)
        if mantle_thickness > 0 else eta_low
    )
    eta_mantle_avg = np.clip(eta_mantle_avg, eta_low, eta_up)

    # --- Find BDT depth ----------------------------------------------------
    bdt_z = 0.0
    found_bdt = False
    for iz in range(1, nz):
        if out_bdt[iz - 1] <= 0 and out_bdt[iz] != out_bdt[iz - 1]:
            bdt_z = out_z[iz] * 1e-3
            found_bdt = True
            break

    if current_h > 0:
        te_cubed_sum += current_h ** 3
    te_decoupled = te_cubed_sum ** (1.0 / 3.0)

    return {
        "nz": nz,
        "z": out_z,
        "T": out_T,
        "rho": out_rho,
        "brittle": out_brittle,
        "ductile": out_ductile,
        "peierls": out_peierls,
        "yield": out_yield,
        "eta_diff": out_eta_diff,
        "eta_disl": out_eta_disl,
        "eta_eff": out_eta_eff,
        "bdt": out_bdt,
        "layer_id": out_layer_id,
        "layer_name": out_layer_name,
        "total_strength": strength,
        "crustal_strength": crustal_strength,
        "mechanical_thickness": h_mech,
        "te_decoupled": te_decoupled,
        "average_eta": avg_eta,
        "eta_crust_average": eta_crust_avg,
        "eta_mantle_average": eta_mantle_avg,
        "bdt_z": bdt_z,
        "found_bdt": found_bdt,
        "topo_z": top,
    }


# ============================================================================
#  Exodus writer
# ============================================================================
def write_exodus_rheology(src_path: str, dst_path: str, data_dict: dict) -> None:
    """Copy an Exodus file and append rheological nodal variables.

    Parameters
    ----------
    src_path : str
        Path to the original Exodus file.
    dst_path : str
        Path for the output Exodus file.
    data_dict : dict[str, np.ndarray]
        Mapping of variable names to full nodal arrays.
    """
    print(f"\n[+] Writing rheology back to Exodus: {dst_path}")
    new_vars = list(data_dict.keys())
    n_new = len(new_vars)

    with nc.Dataset(src_path, "r") as src, nc.Dataset(dst_path, "w") as dst:
        # Copy dimensions
        for name, dim in src.dimensions.items():
            size = len(dim) if not dim.isunlimited() else None
            if name == "num_nod_var":
                size += n_new
            dst.createDimension(name, size)

        # Copy global attributes
        dst.setncatts({k: src.getncattr(k) for k in src.ncattrs()})

        # Copy variables
        for name, var in src.variables.items():
            out_var = dst.createVariable(name, var.datatype, var.dimensions)
            out_var.setncatts({k: var.getncattr(k) for k in var.ncattrs()})

            if name == "name_nod_var":
                old_names = var[:]
                new_names = np.zeros((n_new, old_names.shape[1]), dtype="S1")
                for i, vname in enumerate(new_vars):
                    v_arr = np.array(
                        list(vname.ljust(old_names.shape[1], "\0")), dtype="S1"
                    )
                    new_names[i, :] = v_arr
                out_var[:] = np.vstack([old_names, new_names])
            else:
                out_var[:] = var[:]

        # Create new variables
        old_n_nod_var = src.dimensions["num_nod_var"].size
        n_time = (
            src.dimensions["time_step"].size
            if "time_step" in src.dimensions
            else 1
        )
        for i, vname in enumerate(new_vars):
            v_idx = old_n_nod_var + i + 1
            out_var = dst.createVariable(
                f"vals_nod_var{v_idx}", "f8", ("time_step", "num_nodes")
            )
            out_var[:] = np.tile(data_dict[vname], (n_time, 1))


# ============================================================================
#  Multiprocessing worker
# ============================================================================

# Module-level globals populated by pool initializer — avoids pickling
# large arrays for every task.
_shared: dict = {}


def _pool_init(
    z_all: np.ndarray,
    T_all: np.ndarray,
    node_block_id: np.ndarray,
    node_density: np.ndarray,
    columns: list,
    x: np.ndarray,
    y: np.ndarray,
    block_props: dict,
    block_names: list,
    resolution: float,
    eta_bounds: tuple,
) -> None:
    """Initializer for each worker process — stores shared read-only data."""
    _shared["z_all"] = z_all
    _shared["T_all"] = T_all
    _shared["node_block_id"] = node_block_id
    _shared["node_density"] = node_density
    _shared["columns"] = columns
    _shared["x"] = x
    _shared["y"] = y
    _shared["block_props"] = block_props
    _shared["block_names"] = block_names
    _shared["resolution"] = resolution
    _shared["eta_bounds"] = eta_bounds


def _process_column(pos2D: int) -> dict:
    """Process a single vertical column — called by pool workers.

    Returns a dict with all results needed for CSV writing and Exodus arrays.
    """
    z_all = _shared["z_all"]
    T_all = _shared["T_all"]
    node_block_id = _shared["node_block_id"]
    node_density_all = _shared["node_density"]
    columns = _shared["columns"]
    block_props = _shared["block_props"]
    block_names = _shared["block_names"]
    resolution = _shared["resolution"]
    eta_bounds = _shared["eta_bounds"]

    col_nodes = columns[pos2D]
    zz = z_all[col_nodes]
    T_profile = T_all[col_nodes]
    col_block_ids = node_block_id[col_nodes]
    col_densities = node_density_all[col_nodes].copy()

    # Fallback density for NaN values
    for i, (bid, rho) in enumerate(zip(col_block_ids, col_densities)):
        if np.isnan(rho) and bid in block_props:
            col_densities[i] = block_props[bid].density

    xkm = _shared["x"][pos2D] * 1e-3
    ykm = _shared["y"][pos2D] * 1e-3

    results = {"pos2D": pos2D, "xkm": xkm, "ykm": ykm, "col_nodes": col_nodes}

    for mode in ["CONSTANT", "ORIGINAL"]:
        yse = compute_yse(
            zz, T_profile, col_block_ids, col_densities,
            block_props, block_names, resolution, eta_bounds, mode=mode,
        )
        results[mode] = yse

    return results


# ============================================================================
#  Main processing loop
# ============================================================================
def main() -> None:
    t_start = datetime.now()
    print("=" * 60)
    print("  process_exodus.py  (block-aware, parallel)")
    print("=" * 60)

    import argparse
    parser = argparse.ArgumentParser(
        description="Process an Exodus thermal model for rheological calculations."
    )
    parser.add_argument(
        "config", nargs="?",
        default=os.path.join(_script_dir, "process_exodus.ini"),
        help="Path to the configuration .ini file (default: process_exodus.ini)",
    )
    parser.add_argument(
        "--workers", "-w", type=int, default=None,
        help="Number of parallel worker processes (overrides config file). "
             "Default: value from config, or all available CPU cores.",
    )
    args = parser.parse_args()
    cfg_path = args.config

    workspace = os.path.normpath(os.path.join(_script_dir, ".."))

    # --- [1] Read config ----------------------------------------------------
    print(f"\n[1] Reading config: {cfg_path}")
    cfg = read_config(cfg_path)

    # Derive output name from input filename if not set
    if not cfg.out_name:
        cfg.out_name = os.path.splitext(os.path.basename(cfg.input_file))[0]

    # CLI --workers overrides config, config overrides auto-detect
    if args.workers is not None:
        n_workers = args.workers
    elif cfg.n_workers > 0:
        n_workers = cfg.n_workers
    else:
        n_workers = cpu_count()
    print(f"    Workers: {n_workers}")

    # Load rheolopy materials database
    db_mats = load_materials()
    print(f"    Loaded {len(db_mats)} materials from rheolopy.")

    # --- [2] Read Exodus file -----------------------------------------------
    print(f"\n[2] Reading Exodus: {cfg.input_file}")
    exo_path = os.path.join(workspace, cfg.input_file)
    exo = read_exodus(exo_path, cfg.nx, cfg.ny)

    z_all = exo["z_all"]
    T_all = exo["T_all"]
    num_nodes = exo["num_nodes"]
    nodes2D = exo["nodes2D"]
    block_names = exo["block_names"]
    num_el_blk = exo["num_el_blk"]
    ds = exo["ds"]
    columns = exo["columns"]
    x = exo["x"]
    y = exo["y"]

    # --- [3] Build node -> block map ----------------------------------------
    print(f"\n[3] Building node -> block map from connectivity...")
    node_block_id = build_node_block_map(ds, num_nodes, num_el_blk, z_all)
    print(f"    Done. Block distribution:")
    for blk in range(1, num_el_blk + 1):
        n = np.sum(node_block_id == blk)
        print(f"      {block_names[blk-1]}: {n} nodes")

    # --- [4] Build per-node density -----------------------------------------
    print(f"\n[4] Building per-node density...")
    node_density = build_node_density(ds, num_nodes, num_el_blk)

    # Close the Exodus dataset (all data extracted)
    ds.close()

    # --- [5] Read rheology mapping ------------------------------------------
    print(f"\n[5] Reading rheology: {cfg.rheo_file}")
    rheo_path = os.path.join(workspace, cfg.rheo_file)
    rheo_dict = read_rheology(rheo_path, db_mats)

    # Build block_id -> LayerProp mapping
    block_props: dict[int, LayerProp] = {}
    for blk_idx, blk_name in enumerate(block_names):
        blk_id = blk_idx + 1
        if blk_name in rheo_dict:
            block_props[blk_id] = rheo_dict[blk_name]
            print(f"    Block {blk_id} ({blk_name}) -> {rheo_dict[blk_name].material_id}")
        else:
            print(f"    WARNING: Block {blk_id} ({blk_name}) has no entry in rheology file!")

    # --- [6] Set up outputs -------------------------------------------------
    out_dir = os.path.join(workspace, cfg.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n[6] Output directory: {out_dir}")

    def init_outputs(mode):
        mode_str = mode.lower()
        fn_str = os.path.join(out_dir, f"{cfg.out_name}_strength_{mode_str}.csv")
        fn_3d = os.path.join(out_dir, f"{cfg.out_name}_info3D_{mode_str}.csv")
        fn_bdt = os.path.join(out_dir, f"{cfg.out_name}_bdt_{mode_str}.csv")
        fn_thick = os.path.join(out_dir, f"{cfg.out_name}_thickness_{mode_str}.csv")

        f_s = open(fn_str, "w")
        f_3 = open(fn_3d, "w")
        f_b = open(fn_bdt, "w")
        f_t = open(fn_thick, "w")

        f_s.write(
            "x[km],y[km],log_total_strength,log_crustal_strength,strength_ratio[],"
            "mechanical_thickness[km],average_eta[log Pa_s],average_eta_crust[log Pa_s],"
            "average_eta_mantle[log Pa_s]\n"
        )
        f_3.write(
            "x[km],y[km],z[km],block_id[],block_name[],T[degC],rho[kg/m3],dsigma[MPa],"
            "log_eta_diffusion[Pa*s],log_eta_dislocation[Pa*s],log_eta_effective[Pa*s],bdt[]\n"
        )
        f_b.write("x[km],y[km],bdt[km],topography[km]\n")
        f_t.write("x[km],y[km],h_mech[km],te_decoupled[km]\n")

        return (f_s, f_3, f_b, f_t), (fn_str, fn_3d, fn_bdt, fn_thick)

    files_const, names_const = init_outputs("CONSTANT")
    files_orig, names_orig = init_outputs("ORIGINAL")

    # Arrays to hold original resolution data for Exodus nodal variables
    exo_dsigma = np.zeros(num_nodes)
    exo_eta_eff = np.zeros(num_nodes)
    exo_eta_diff = np.zeros(num_nodes)
    exo_eta_disl = np.zeros(num_nodes)
    exo_bdt = np.zeros(num_nodes)
    exo_h_mech = np.zeros(num_nodes)
    exo_te_decoupled = np.zeros(num_nodes)

    # --- [7] Process all columns in parallel --------------------------------
    eta_bounds = (cfg.eta_low, cfg.eta_up)
    print(f"\n[7] Processing {nodes2D} columns with {n_workers} workers ...")

    pool = Pool(
        processes=n_workers,
        initializer=_pool_init,
        initargs=(
            z_all, T_all, node_block_id, node_density,
            columns, x, y, block_props, block_names,
            cfg.resolution, eta_bounds,
        ),
    )

    done = 0
    for result in pool.imap_unordered(_process_column, range(nodes2D), chunksize=256):
        done += 1
        if done % 1000 == 0 or done == nodes2D:
            elapsed_so_far = (datetime.now() - t_start).total_seconds()
            rate = done / elapsed_so_far if elapsed_so_far > 0 else 0
            eta_sec = (nodes2D - done) / rate if rate > 0 else 0
            print(
                f"    {done:>8} / {nodes2D}  "
                f"({100*done/nodes2D:5.1f}%)  "
                f"{rate:.0f} col/s  "
                f"ETA {eta_sec/60:.1f} min",
                end="\r", flush=True,
            )

        xkm = result["xkm"]
        ykm = result["ykm"]
        col_nodes = result["col_nodes"]

        for mode, (f_str, f_3d, f_bdt, f_thick) in zip(
            ["CONSTANT", "ORIGINAL"], [files_const, files_orig]
        ):
            yse = result[mode]

            ts = yse["total_strength"]
            cs = yse["crustal_strength"]
            log_ts = np.log10(ts) if ts > 0 else -99
            log_cs = np.log10(cs) if cs > 0 else -99
            ratio = (cs * 100.0 / ts) if ts > 0 else 0.0

            f_str.write(
                f"{xkm:.4f},{ykm:.4f},{log_ts:.6f},{log_cs:.6f},{ratio:.4f},"
                f"{yse['mechanical_thickness'] * 1e-3:.4f},{np.log10(yse['average_eta']):.6f},"
                f"{np.log10(yse['eta_crust_average']):.6f},{np.log10(yse['eta_mantle_average']):.6f}\n"
            )

            f_thick.write(
                f"{xkm:.4f},{ykm:.4f},"
                f"{yse['mechanical_thickness'] * 1e-3:.4f},"
                f"{yse['te_decoupled'] * 1e-3:.4f}\n"
            )

            for iz in range(yse["nz"]):
                e_diff = yse["eta_diff"][iz]
                e_disl = yse["eta_disl"][iz]
                e_eff = yse["eta_eff"][iz]
                log_diff = np.log10(e_diff) if (not np.isnan(e_diff) and e_diff > 0) else np.nan
                log_disl = np.log10(e_disl) if (not np.isnan(e_disl) and e_disl > 0) else np.nan
                log_eff = np.log10(e_eff) if (not np.isnan(e_eff) and e_eff > 0) else np.nan

                f_3d.write(
                    f"{xkm:.4f},{ykm:.4f},{yse['z'][iz] * 1e-3:.4f},"
                    f"{yse['layer_id'][iz]},{yse['layer_name'][iz]},"
                    f"{yse['T'][iz]:.4f},{yse['rho'][iz]:.1f},"
                    f"{yse['yield'][iz] * 1e-6:.6f},"
                    f"{log_diff:.6f},{log_disl:.6f},{log_eff:.6f},"
                    f"{yse['bdt'][iz]}\n"
                )

                # Save into full arrays if ORIGINAL mode for Exodus output
                if mode == "ORIGINAL":
                    idx = col_nodes[iz]
                    exo_dsigma[idx] = yse["yield"][iz] * 1e-6  # MPa
                    exo_eta_eff[idx] = log_eff
                    exo_eta_diff[idx] = log_diff
                    exo_eta_disl[idx] = log_disl
                    exo_bdt[idx] = yse["bdt"][iz]
                    exo_h_mech[idx] = yse["mechanical_thickness"] * 1e-3
                    exo_te_decoupled[idx] = yse["te_decoupled"] * 1e-3

            if yse["found_bdt"]:
                f_bdt.write(
                    f"{xkm:.4f},{ykm:.4f},{yse['bdt_z']:.4f},"
                    f"{yse['topo_z'] * 1e-3:.4f}\n"
                )

    pool.close()
    pool.join()
    print()

    for fs in files_const + files_orig:
        fs.close()

    # --- [8] Write augmented Exodus file ------------------------------------
    out_exo_name = f"{cfg.out_name}_rheology.e"
    dst_exo = os.path.join(out_dir, out_exo_name)
    data_dict = {
        "dsigma_MPa": exo_dsigma,
        "log_eta_eff": exo_eta_eff,
        "log_eta_diff": exo_eta_diff,
        "log_eta_disl": exo_eta_disl,
        "bdt": exo_bdt,
        "h_mech": exo_h_mech,
        "te_decoupled": exo_te_decoupled,
    }
    write_exodus_rheology(exo_path, dst_exo, data_dict)

    elapsed = (datetime.now() - t_start).total_seconds()
    print(f"\n[Done] Elapsed: {elapsed:.1f} s ({elapsed/60:.1f} min)")


if __name__ == "__main__":
    main()
