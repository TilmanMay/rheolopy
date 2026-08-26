"""
Rheology processing script
Loads data from CSV, computes rheological properties, and saves results.
"""

import os
import sys
import configparser
import numpy as np
from datetime import date

# Import rheology package
from rheolopy import materials, get_material_by_id, compute_dsigma, eta_effective


def main():
    # Get current working directory
    path = os.getcwd()

    # Create an instance of ConfigParser
    config = configparser.ConfigParser()

    # Read the config file (try script directory first, fallback to CWD)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(script_dir, "config.ini")
    if not os.path.exists(config_path):
        config_path = "process_rheology/config.ini"
    config.read(config_path)

    # Accessing parameters from the General section
    strain_rate = config.getfloat("General", "strain_rate")

    # Accessing parameters from the Settings section
    inputfile = config.get("Settings", "inputfile")
    outputfile = config.get("Settings", "outputfile")
    rheology_law = config.get("Settings", "rheology_law")
    
    # Strip "process_rheology/" from paths if we are already inside that directory
    if os.path.basename(path) == "process_rheology":
        if inputfile.startswith("process_rheology/"):
            inputfile = inputfile.replace("process_rheology/", "")
        if outputfile.startswith("process_rheology/"):
            outputfile = outputfile.replace("process_rheology/", "")

    # Load input data
    data = np.loadtxt(os.path.join(path, inputfile), delimiter=",", comments="#")

    ################################
    ### Calculating viscosity and strength

    # Load material database
    mat_dbase = materials()

    # Find the requested material by id or type
    mat = None
    for m in mat_dbase:
        if m.id == rheology_law or rheology_law in m.type:
            mat = m
            break

    if mat is None:
        print(f"ERROR: Material '{rheology_law}' not found in database!")
        print("Available materials:")
        for m in mat_dbase:
            print(f"  - id: {m.id}, type: {m.type}")
        sys.exit(1)

    print("Opted material:", rheology_law)
    print("Chosen material:", mat)

    # Initialize output arrays
    dsigma_c = np.empty(len(data))
    dsigma_e = np.empty(len(data))
    eff_vis = np.empty(len(data))

    # Loop through all points
    for i in range(len(data)):
        # Extract data: depth in km (column 2), temperature in °C (column 4)
        depth_m = data[i, 2] * 1e3  # Convert km to m
        depth_z = -abs(depth_m)     # Enforce negative depth convention
        temp_k = data[i, 4] + 273.15  # Convert °C to K

        # Compute differential stress for compression and extension
        s_d_c, s_d_e = compute_dsigma(mat, depth_z, temp_k, strain_rate)

        dsigma_c[i] = s_d_c
        dsigma_e[i] = s_d_e

        # Compute effective viscosity
        eta_eff, _, _ = eta_effective(mat, temp_k, strain_rate)
        eff_vis[i] = eta_eff

    # Add computed columns to data
    data = np.column_stack((data, dsigma_c))
    data = np.column_stack((data, dsigma_e))
    data = np.column_stack((data, np.log10(eff_vis)))

    # Create metadata for output file
    meta_data = (
        "#Created on: " + str(date.today()) + "\n"
        "#Input file is: " + str(inputfile) + "\n"
        "#Output file is: " + str(outputfile) + "\n"
        "#Material is: " + str(rheology_law) + "\n"
        "#Strain rate is: " + str(strain_rate) + "\n"
        "#\n"
    )

    # Save the output
    np.savetxt(
        str(outputfile),
        data,
        delimiter=",\t",
        header="#x(km), y(km), depth(km), Pressure(bar), Temperature(oC), Density(kg/m3), "
        "Vp(km/s), Vs(km/s), Vs_diff(%), Pseudo-melts(%), dsigma_c(Pascal), dsigma_e(Pascal), Viscosity(log10Pas)",
        comments=meta_data,
        fmt="%10.3f",
    )

    print(f"\nProcessing complete!")
    print(f"Results saved to: {outputfile}")


if __name__ == "__main__":
    main()
