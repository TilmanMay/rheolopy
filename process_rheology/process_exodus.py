"""
Process an Exodus thermal model output for rheological calculations.

Workflow:
  1. Read configuration (keyword-based INI)
  2. Read Exodus file -> node coordinates (x, y, z), temperature (T), and block densities
  3. Read rheology property file -> maps layer to rheolopy material ID & strain rate
  4. For each 2D column (x, y):
       - Extract vertical profile z(iz), T(iz)
       - Compute YSE on constant-resolution grid AND/OR original resolution
       - Write per-column integrated values + 3D field data + BDT
  5. Write output CSVs

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
from scipy.io import netcdf_file
from datetime import datetime

# ---------------------------------------------------------------------------
# Add rheolopy source to path
# ---------------------------------------------------------------------------
_script_dir = os.path.dirname(os.path.abspath(__file__))
_src_dir = os.path.normpath(os.path.join(_script_dir, "..", "src"))
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)

from rheolopy import materials as load_materials, get_material_by_id
from rheolopy.rheolopy import (
    sigma_byerlee,
    eta_effective as rheo_eta_effective,
    calc_peierls,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
GRAV = 9.81        # m/s²
KELVIN = 273.15    # offset °C -> K
STRESS_LIM = 20e6  # Pa  – threshold for mechanical thickness
STRESS_PD = 200e6  # Pa  – threshold for Peierls dislocation transition


# ============================================================================
#  Data structures
# ============================================================================
class LayerProp:
    """Rheological properties for one geological layer."""

    def __init__(self, line: str, db_mats: list) -> None:
        f = line.split()
        self.name: str = f[0]
        self.material_id: str = f[1]
        self.strain_rate: float = float(f[2])
        self.common_layers: int = int(f[3]) if len(f) > 3 else 1

        self.material = get_material_by_id(db_mats, self.material_id)
        if self.material is None:
            raise ValueError(f"Material ID '{self.material_id}' not found in rheolopy database.")
            
        # The density will be updated later if Exodus provides it
        self.density = self.material.rho_b
        if self.density is None or np.isnan(self.density):
            self.density = 2700.0  # Safe fallback


class Config:
    """Parsed configuration from the keyword-based INI file."""

    def __init__(self) -> None:
        self.input_file: str = ""
        self.nx: int = 0
        self.ny: int = 0
        self.rheo_file: str = ""
        self.eta_low: float = 1e18
        self.eta_up: float = 1e25
        self.out_dir: str = "output"
        self.out_name: str = "model"
        self.resolution: float = 500.0  # metres


# ============================================================================
#  Readers
# ============================================================================
def read_config(filename: str) -> Config:
    cfg = Config()
    with open(filename, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith("#INPUT_FILE:"):
                parts = line.split()
                cfg.input_file = parts[2]
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
                cfg.out_name = parts[2]
            elif line.startswith("#RESOLUTION:"):
                cfg.resolution = float(line.split()[1])
    return cfg


def read_exodus(filename: str, nx: int, ny: int):
    f = netcdf_file(filename, "r", mmap=False)

    num_nodes = f.dimensions["num_nodes"]
    num_el_blk = f.dimensions["num_el_blk"]
    nodes2D = nx * ny
    nsurf = num_nodes // nodes2D

    x_all = np.array(f.variables["coordx"].data, dtype=np.float64)
    y_all = np.array(f.variables["coordy"].data, dtype=np.float64)
    z_all = np.array(f.variables["coordz"].data, dtype=np.float64)

    # Temperature: nodal variable 1, last time step
    T_all = np.array(f.variables["vals_nod_var1"].data[-1, :], dtype=np.float64)

    # Extract block densities
    block_densities = [None] * num_el_blk
    
    # Find which element variable is density
    density_var_idx = -1
    if "name_elem_var" in f.variables:
        names = f.variables["name_elem_var"].data
        for i in range(names.shape[0]):
            name_str = names[i].tobytes().decode("ascii", errors="ignore").strip().rstrip("\x00")
            if "density" in name_str.lower():
                density_var_idx = i + 1
                break
                
    if density_var_idx > 0:
        for blk in range(1, num_el_blk + 1):
            var_name = f"vals_elem_var{density_var_idx}eb{blk}"
            if var_name in f.variables:
                data = f.variables[var_name].data[-1, :]
                if data.size > 0:
                    block_densities[blk - 1] = float(data[0])

    f.close()

    assert num_nodes == nodes2D * nsurf, f"Node count mismatch: {num_nodes} != {nodes2D} * {nsurf}"

    # Group nodes by unique (X, Y) pairs to form columns
    xy = np.column_stack((x_all, y_all))
    unique_xy, unique_indices, inverse_indices = np.unique(xy, axis=0, return_index=True, return_inverse=True)
    
    # Sort the unique points by their first appearance to preserve original column ordering
    order = np.argsort(unique_indices)
    unique_xy = unique_xy[order]
    
    # Re-map inverse_indices so they correspond to the appearance-ordered columns
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
        sorted_idx = np.argsort(-z_all[col_nodes])
        col_nodes_sorted = col_nodes[sorted_idx]
        columns[i] = col_nodes_sorted
        x[i] = x_all[col_nodes_sorted[0]]
        y[i] = y_all[col_nodes_sorted[0]]

    print(f"  Exodus file read successfully:")
    print(f"    nodes2D = {nodes2D} ({nx} x {ny}), nsurf = {nsurf}")
    print(f"    x range: [{x.min():.0f}, {x.max():.0f}] m")
    print(f"    y range: [{y.min():.0f}, {y.max():.0f}] m")
    print(f"    z range: [{z_all.min():.0f}, {z_all.max():.0f}] m")
    print(f"    T range: [{T_all.min():.2f}, {T_all.max():.2f}] °C")
    if density_var_idx > 0:
        print(f"    Found density in element variable {density_var_idx}")

    return x, y, z_all, T_all, nsurf, block_densities, num_nodes, columns


def read_rheology(filename: str, db_mats: list, block_densities: list):
    raw: list[LayerProp] = []
    with open(filename, "r") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            raw.append(LayerProp(line, db_mats))

    # Override density from Exodus if available
    for i, prop in enumerate(raw):
        if i < len(block_densities) and block_densities[i] is not None:
            prop.density = block_densities[i]

    # Expand common_layers
    props: list[LayerProp] = []
    for lp in raw:
        for _ in range(max(lp.common_layers, 1)):
            props.append(lp)

    print(f"  Rheology file: {len(raw)} raw layers -> {len(props)} expanded layers")
    return props


# ============================================================================
#  YSE computation
# ============================================================================
def interpolate_value(v0: float, v1: float, grad: float) -> float:
    return v0 + (v1 - v0) * grad

def evaluate_point(zp, temp_K, props, is_crust, eta_low, eta_up, effective_rho):
    """Evaluate rheological properties at a single depth point."""
    mat = props.material
    strain_rate = props.strain_rate
    
    # Workaround: we temporarily override the material's bulk density to use 
    # the effective integrated density for accurate lithostatic pressure
    orig_rho = mat.rho_b
    mat.rho_b = effective_rho

    depth = abs(zp)  # Using absolute Z as depth below surface
    brittle = sigma_byerlee(mat, depth, "compression")

    eta_eff_val, eta_dis_val, eta_diff_val = rheo_eta_effective(
        mat, temp_K, strain_rate
    )

    if not np.isnan(eta_eff_val):
        eta_eff_val = np.clip(eta_eff_val, eta_low, eta_up)

    creep = 2.0 * eta_eff_val * strain_rate if not np.isnan(eta_eff_val) else 1e99

    # Peierls stress (activated if the material has peierls parameters)
    peierls_val = 0.0
    if hasattr(mat, 'eps_peierls') and mat.eps_peierls is not None:
        peierls_val = calc_peierls(mat, temp_K, strain_rate)
        if peierls_val < 0.0:
            peierls_val = 0.0
        if creep > STRESS_PD and peierls_val > 0.0 and peierls_val < creep:
            creep = peierls_val

    dsigma = min(brittle, creep)
    bdt = 0 if brittle < creep else 1  # 0 = brittle, 1 = ductile
    
    # Restore original density
    mat.rho_b = orig_rho

    return {
        "rho": props.density,
        "brittle": brittle,
        "ductile": creep,
        "peierls": peierls_val,
        "yield": dsigma,
        "eta_diff": eta_diff_val if not np.isnan(eta_diff_val) else 0.0,
        "eta_disl": eta_dis_val if not np.isnan(eta_dis_val) else 0.0,
        "eta_eff": eta_eff_val if not np.isnan(eta_eff_val) else 0.0,
        "bdt": bdt,
    }


def compute_yse(
    pos2D: int,
    zz: np.ndarray,
    temp_profile: np.ndarray,
    props: list,
    resolution: float,
    eta_bounds: tuple[float, float],
    nsurf: int,
    mode: str = "CONSTANT"  # "CONSTANT" or "ORIGINAL"
):
    eta_low, eta_up = eta_bounds
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
        # ORIGINAL resolution
        out_z = zz
        nz = len(out_z)

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

        # --- Locate layer and interpolate temperature ----------------------
        if mode == "ORIGINAL":
            # Direct mapping
            layer_idx = min(iz, len(props) - 1)
            temp = temp_K[iz]
            if iz == nz - 1:
                layer_idx = len(props) - 1
        else:
            if iz == 0:
                temp = temp_K[0]
                layer_idx = 0
            elif iz == nz - 1:
                temp = temp_K[-1]
                layer_idx = len(props) - 1
            else:
                found = False
                for isrf in range(nsurf - 1):
                    z0, z1 = zz[isrf], zz[isrf + 1]
                    if (zp - z0) * (zp - z1) <= 0.0:
                        gradient = abs((zp - z0) / (z1 - z0)) if z1 != z0 else 0.0
                        temp = interpolate_value(temp_K[isrf], temp_K[isrf + 1], gradient)
                        layer_idx = isrf
                        found = True
                        break
                if not found:
                    if zp > zz[0]:
                        temp, layer_idx = temp_K[0], 0
                    else:
                        temp, layer_idx = temp_K[-1], len(props) - 1

        prop = props[layer_idx]
        is_crust = ("mantle" not in prop.name.lower())
        
        # --- Integrated quantities -----------------------------------------
        thickness = 0.0 if iz == 0 else abs(zp - out_z[iz - 1])
        
        # Compute exact lithostatic pressure (integral of rho * g * dz)
        p_litho += thickness * prop.density * GRAV
        
        # Compute effective density for sigma_byerlee
        depth = abs(zp - top)
        effective_rho = p_litho / (GRAV * depth) if depth > 0 else prop.density
        
        # Evaluate properties at this point
        pt = evaluate_point(zp - top, temp, prop, is_crust, eta_low, eta_up, effective_rho)
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

        # --- Store ---------------------------------------------------------
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
        out_layer_id[iz] = layer_idx
        out_layer_name[iz] = prop.name

    # --- Compute integrated viscosities ------------------------------------
    # (Simple weighted averages based on thickness for this column)
    layer_thickness = top - bottom
    
    # We use a nominal strain rate of 1e-15 to calculate equivalent average viscosity
    nominal_sr = 1e-15
    avg_eta = strength / (2.0 * nominal_sr * layer_thickness) if layer_thickness > 0 else eta_low
    avg_eta = np.clip(avg_eta, eta_low, eta_up)

    eta_crust_avg = crustal_strength / (2.0 * nominal_sr * thickness_crust) if thickness_crust > 0 else eta_low
    eta_crust_avg = np.clip(eta_crust_avg, eta_low, eta_up)

    mantle_thickness = layer_thickness - thickness_crust
    eta_mantle_avg = (strength - crustal_strength) / (2.0 * nominal_sr * mantle_thickness) if mantle_thickness > 0 else eta_low
    eta_mantle_avg = np.clip(eta_mantle_avg, eta_low, eta_up)

    # --- Find BDT depth ---------------------------------------------------
    bdt_z = 0.0
    found_bdt = False
    for iz in range(1, nz):
        if (out_bdt[iz - 1] <= 0
                and out_bdt[iz] != out_bdt[iz - 1]):
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


def write_exodus_rheology(src_path: str, dst_path: str, data_dict: dict):
    print(f"\n[+] Writing rheology back to Exodus: {dst_path}")
    new_vars = list(data_dict.keys())
    n_new = len(new_vars)

    with nc.Dataset(src_path, 'r') as src, nc.Dataset(dst_path, 'w') as dst:
        # Copy dimensions
        for name, dim in src.dimensions.items():
            size = len(dim) if not dim.isunlimited() else None
            if name == 'num_nod_var':
                size += n_new
            dst.createDimension(name, size)
        
        # Copy global attributes
        dst.setncatts({k: src.getncattr(k) for k in src.ncattrs()})
        
        # Copy variables
        for name, var in src.variables.items():
            out_var = dst.createVariable(name, var.datatype, var.dimensions)
            out_var.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
            
            if name == 'name_nod_var':
                old_names = var[:]
                new_names = np.zeros((n_new, old_names.shape[1]), dtype='S1')
                for i, vname in enumerate(new_vars):
                    v_arr = np.array(list(vname.ljust(old_names.shape[1], '\0')), dtype='S1')
                    new_names[i, :] = v_arr
                out_var[:] = np.vstack([old_names, new_names])
            else:
                out_var[:] = var[:]
                
        # Create new variables
        old_n_nod_var = src.dimensions['num_nod_var'].size
        n_time = src.dimensions['time_step'].size if 'time_step' in src.dimensions else 1
        for i, vname in enumerate(new_vars):
            v_idx = old_n_nod_var + i + 1
            out_var = dst.createVariable(f'vals_nod_var{v_idx}', 'f8', ('time_step', 'num_nodes'))
            out_var[:] = np.tile(data_dict[vname], (n_time, 1))


# ============================================================================
#  Main processing loop
# ============================================================================
def main() -> None:
    t_start = datetime.now()
    print("=" * 60)
    print("  process_exodus.py")
    print("=" * 60)

    if len(sys.argv) > 1:
        cfg_path = sys.argv[1]
    else:
        cfg_path = os.path.join(_script_dir, "process_exodus.ini")

    workspace = os.path.normpath(os.path.join(_script_dir, ".."))

    print(f"\n[1] Reading config: {cfg_path}")
    cfg = read_config(cfg_path)
    
    # Load rheolopy materials database
    db_mats = load_materials()
    print(f"    Loaded {len(db_mats)} materials from rheolopy.")

    print(f"\n[2] Reading Exodus: {cfg.input_file}")
    exo_path = os.path.join(workspace, cfg.input_file)
    x, y, z_all, T_all, nsurf, block_densities, num_nodes, columns = read_exodus(exo_path, cfg.nx, cfg.ny)
    nodes2D = cfg.nx * cfg.ny

    print(f"\n[3] Reading rheology: {cfg.rheo_file}")
    rheo_path = os.path.join(workspace, cfg.rheo_file)
    props = read_rheology(rheo_path, db_mats, block_densities)

    n_layers = nsurf - 1
    if len(props) != n_layers:
        print(f"  WARNING: expanded layers ({len(props)}) != nsurf-1 ({n_layers}). Adjusting.")
        while len(props) < n_layers:
            props.append(props[-1])
        props = props[:n_layers]

    out_dir = os.path.join(workspace, cfg.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n[4] Output directory: {out_dir}")

    # Helper function to open and write headers
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
        
        f_s.write("x[km],y[km],log_total_strength,log_crustal_strength,strength_ratio[],"
                  "mechanical_thickness[km],average_eta[log Pa_s],average_eta_crust[log Pa_s],"
                  "average_eta_mantle[log Pa_s]\n")
        f_3.write("x[km],y[km],z[km],layer_id[],layer_name[],T[degC],rho[kg/m3],dsigma[MPa],"
                  "log_eta_diffusion[Pa*s],log_eta_dislocation[Pa*s],log_eta_effective[Pa*s],bdt[]\n")
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

    print(f"\n[5] Processing {nodes2D} columns (constant and original resolution) ...")
    eta_bounds = (cfg.eta_low, cfg.eta_up)

    for pos2D in range(nodes2D):
        if (pos2D + 1) % 100 == 0 or pos2D == 0:
            print(f"    column {pos2D + 1} / {nodes2D}", end="\r", flush=True)

        col_nodes = columns[pos2D]
        zz = z_all[col_nodes]
        T_profile = T_all[col_nodes]

        xkm = x[pos2D] * 1e-3
        ykm = y[pos2D] * 1e-3

        for mode, (f_str, f_3d, f_bdt, f_thick) in zip(["CONSTANT", "ORIGINAL"], [files_const, files_orig]):
            yse = compute_yse(
                pos2D, zz, T_profile, props, cfg.resolution, eta_bounds, nsurf, mode=mode
            )

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
            
            f_thick.write(f"{xkm:.4f},{ykm:.4f},{yse['mechanical_thickness'] * 1e-3:.4f},{yse['te_decoupled'] * 1e-3:.4f}\n")

            for iz in range(yse["nz"]):
                e_diff, e_disl, e_eff = yse["eta_diff"][iz], yse["eta_disl"][iz], yse["eta_eff"][iz]
                log_diff = np.log10(e_diff) if e_diff > 0 else -99
                log_disl = np.log10(e_disl) if e_disl > 0 else -99
                log_eff = np.log10(e_eff) if e_eff > 0 else -99
                
                f_3d.write(
                    f"{xkm:.4f},{ykm:.4f},{yse['z'][iz] * 1e-3:.4f},{yse['layer_id'][iz]},"
                    f"{yse['layer_name'][iz]},{yse['T'][iz]:.4f},{yse['rho'][iz]:.1f},"
                    f"{yse['yield'][iz] * 1e-6:.6f},{log_diff:.6f},{log_disl:.6f},"
                    f"{log_eff:.6f},{yse['bdt'][iz]}\n"
                )
                
                # Save into full arrays if ORIGINAL mode for Exodus output
                if mode == "ORIGINAL":
                    idx = col_nodes[iz]
                    exo_dsigma[idx] = np.log10(yse['yield'][iz]) if yse['yield'][iz] > 0 else -99
                    exo_eta_eff[idx] = log_eff
                    exo_eta_diff[idx] = log_diff
                    exo_eta_disl[idx] = log_disl
                    exo_bdt[idx] = yse['bdt'][iz]
                    exo_h_mech[idx] = yse['mechanical_thickness'] * 1e-3  # Store in km for Paraview
                    exo_te_decoupled[idx] = yse['te_decoupled'] * 1e-3

            if yse["found_bdt"]:
                f_bdt.write(f"{xkm:.4f},{ykm:.4f},{yse['bdt_z']:.4f},{yse['topo_z'] * 1e-3:.4f}\n")

    for fs in files_const + files_orig:
        fs.close()
        
    # Write to duplicate Exodus file
    dst_exo = os.path.join(workspace, cfg.out_dir, "AlpsModel_rheology.e")
    data_dict = {
        'log_dsigma': exo_dsigma,
        'log_eta_eff': exo_eta_eff,
        'log_eta_diff': exo_eta_diff,
        'log_eta_disl': exo_eta_disl,
        'bdt': exo_bdt,
        'h_mech': exo_h_mech,
        'te_decoupled': exo_te_decoupled
    }
    write_exodus_rheology(exo_path, dst_exo, data_dict)

    elapsed = (datetime.now() - t_start).total_seconds()
    print(f"\n[Done] Elapsed: {elapsed:.1f} s")


if __name__ == "__main__":
    main()
