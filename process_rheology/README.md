# Rheology Data Processing

This folder contains a script to process geological data using the rheology package.

## Files

- `config.ini` - Configuration file with processing parameters
- `process_data.py` - Main processing script  
- `input.csv` - Sample input data file
- `output.csv` - Generated output file (after running the script)

## Configuration

Edit `config.ini` to set your parameters:

```ini
[General]
strain_rate = 1e-17

[Settings]
inputfile = input.csv
outputfile = output.csv
rheology_law = peridotite_dry
```

### Available Materials

The `rheology_law` parameter should match a material ID or type from the database. Common options include:
- `olivine_karato_dry`
- `olivine_karato_wet`
- `olivine_hirth_dry`
- `olivine_hirth_wet`

**Note:** "peridotite_dry" is not in the default database. You may want to use `olivine_hirth_dry` or `olivine_karato_dry` instead, as peridotite is primarily composed of olivine.

## Input Data Format

The input CSV file should have the following columns (comma-separated):
1. x (km)
2. y (km)
3. depth (km)
4. Pressure (bar)
5. Temperature (°C)
6. Density (kg/m³)
7. Vp (km/s)
8. Vs (km/s)
9. Vs_diff (%)
10. Pseudo-melts (%)

## Running the Script

1. Make sure the rheology package is installed:
   ```
   pip install -e ../
   ```

2. Navigate to this directory:
   ```
   cd process_rheology
   ```

3. Run the script:
   ```
   python process_data.py
   ```

## Output

The script generates an output CSV file with the original columns plus three additional columns:
- dsigma_c (Pascal) - Differential stress for compression
- dsigma_e (Pascal) - Differential stress for extension
- Viscosity (log10 Pa·s) - Effective viscosity in log10 scale

The output file includes metadata comments with processing details.
