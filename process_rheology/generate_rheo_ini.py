"""
Auto-generate a rheology .ini file from an Exodus (.e) mesh.

Reads the Exodus file, extracts all element block names and their vertical
layer counts (for structured meshes), and writes a template .ini file that
can be edited before running process_exodus.py.

Usage:
    python generate_rheo_ini.py <exodus_file> [output_ini]

    If output_ini is omitted, writes to <exodus_basename>_rheology.ini in the
    same directory as the Exodus file.

Example:
    python generate_rheo_ini.py AndesModel_SACI_deep.e
    -> produces AndesModel_SACI_deep_rheology.ini
"""

import os
import sys
import numpy as np
import netCDF4 as nc

def generate_rheo_ini(exodus_path: str, output_path: str | None = None) -> str:
    """Generate a rheology .ini template from an Exodus mesh file.

    Parameters
    ----------
    exodus_path : str
        Path to the Exodus (.e) file.
    output_path : str, optional
        Path for the output .ini file.  Derived from the input name if omitted.

    Returns
    -------
    str
        The absolute path of the generated .ini file.
    """
    if output_path is None:
        base = os.path.splitext(os.path.basename(exodus_path))[0]
        output_path = os.path.join(
            os.path.dirname(os.path.abspath(exodus_path)),
            f"{base}_rheology.ini",
        )

    print(f"[1] Reading Exodus file: {exodus_path}")
    ds = nc.Dataset(exodus_path, "r")

    num_nodes = ds.dimensions["num_nodes"].size
    num_el_blk = ds.dimensions["num_el_blk"].size

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

    # --- Auto-detect nodes2D ------------------------------------------------
    x = ds.variables["coordx"][:]
    y = ds.variables["coordy"][:]
    xy = np.column_stack((x, y))
    nodes2D = len(np.unique(xy, axis=0))
    nsurf = num_nodes // nodes2D

    print(f"    Auto-detected: nodes2D = {nodes2D}, nsurf = {nsurf}")

    # --- Calculate elements per block and layers ----------------------------
    lines: list[str] = []
    total_layers = 0
    for blk in range(1, num_el_blk + 1):
        dim_name = f"num_el_in_blk{blk}"
        if dim_name not in ds.dimensions:
            continue
        n_elem = ds.dimensions[dim_name].size
        name = block_names[blk - 1]

        # For a structured mesh: layers = n_elem / nodes2D
        # For unstructured/non-layered blocks this won't be integer
        nodes2D_elem = (nodes2D - 1) if nodes2D > 1 else nodes2D  # elements vs nodes
        # Actually for hex meshes: elem_2D ≈ (nx-1)*(ny-1), but we don't know
        # nx/ny separately. Let's just report elements and let user verify.
        layers_exact = n_elem / nodes2D if nodes2D > 0 else 0

        if layers_exact > 0 and abs(layers_exact - round(layers_exact)) < 0.01:
            layers_int = int(round(layers_exact))
            layers_note = ""
        else:
            layers_int = n_elem  # Can't cleanly divide
            layers_note = "  # WARNING: non-integer layer count, check mesh structure"

        total_layers += layers_int
        lines.append((name, layers_int, layers_note, n_elem))

    ds.close()

    # --- Write the .ini file ------------------------------------------------
    print(f"[2] Writing rheology template: {output_path}")

    max_name_len = max(len(name) for name, _, _, _ in lines) if lines else 20

    with open(output_path, "w") as fh:
        fh.write(f"# Auto-generated rheology template for: {os.path.basename(exodus_path)}\n")
        fh.write(f"# Detected {num_el_blk} element blocks, {nodes2D} horizontal nodes, {nsurf} vertical nodes\n")
        fh.write("#\n")
        fh.write("# Format: Block_Name  Material_ID  Strain_Rate\n")
        fh.write("# \n")
        fh.write("# INSTRUCTIONS:\n")
        fh.write("#   1. Replace REPLACE_ME with a valid material ID from the rheolopy database.\n")
        fh.write("#      Available materials can be listed with:\n")
        fh.write("#        python -c \"from rheolopy.materials import materials; [print(m.id) for m in materials()]\"\n")
        fh.write("#   2. Adjust strain rates as needed for each block.\n")
        fh.write("#   3. Block names MUST match the Exodus file block names exactly.\n")
        fh.write("#\n")

        for name, layers, note, n_elem in lines:
            padded = name.ljust(max_name_len + 2)
            fh.write(f"{padded} REPLACE_ME    1e-16{note}\n")

    print(f"[Done] Template written with {len(lines)} blocks.")
    print(f"       Edit the file to assign materials, then run process_exodus.py.")

    return os.path.abspath(output_path)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: python {sys.argv[0]} <exodus_file> [output_ini]")
        sys.exit(1)

    exo = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else None
    generate_rheo_ini(exo, out)
