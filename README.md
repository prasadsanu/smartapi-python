# Gaussian Spectroscopy Suite

A comprehensive Python toolkit for automated spectroscopic calculations using Gaussian quantum chemistry software. This suite automates the generation of input files, execution of calculations, and analysis of results for various spectroscopic properties.

## Features

### Supported Spectroscopic Calculations
- **IR and Raman Spectroscopy** - Vibrational frequencies and intensities
- **UV-Visible Spectroscopy** - Electronic transitions via TD-DFT
- **NMR Spectroscopy** - Chemical shifts and coupling constants
- **Vibrational Circular Dichroism (VCD)** - Chiral vibrational spectroscopy
- **Raman Optical Activity (ROA)** - Chiral Raman spectroscopy
- **Electronic Circular Dichroism (ECD)** - Electronic chiroptical properties
- **Anharmonic Analysis** - Beyond harmonic approximation
- **Polarizability and Hyperpolarizability** - Nonlinear optical properties
- **Solvent Effects** - SCRF solvation models

### Key Capabilities
- **Automated Input Generation** - Creates properly formatted Gaussian input files
- **Job Execution** - Runs calculations automatically
- **Result Parsing** - Extracts data from Gaussian output files
- **Data Analysis** - Generates comprehensive reports and CSV files
- **Multiple Output Formats** - Text reports, JSON data, CSV files for Excel
- **Error Handling** - Robust error checking and reporting

## Installation

### Prerequisites
1. **Gaussian Software** - Gaussian 16 or Gaussian 09 must be installed and accessible
2. **Python 3.7+** - Required for the Python scripts

### Python Dependencies
Install the required Python packages:

```bash
pip install -r requirements.txt
```

Or install individually:
```bash
pip install numpy pandas matplotlib scipy
```

### Setup
1. Clone or download the files to your working directory
2. Ensure Gaussian is properly installed and the `g16` (or `g09`) command is accessible
3. Modify the Gaussian command in the script if needed (line ~130 in `gaussian_spectroscopy_suite.py`)

## Quick Start

### Basic Usage

```python
from gaussian_spectroscopy_suite import GaussianSpectroscopyCalculator

# Initialize calculator
calc = GaussianSpectroscopyCalculator("my_molecule", "analysis_directory")

# Set molecular geometry (Cartesian coordinates)
geometry = """C          0.000000    0.000000    0.000000
H          0.000000    0.000000    1.089000
H          1.026719    0.000000   -0.363000
H         -0.513360    0.889165   -0.363000
H         -0.513360   -0.889165   -0.363000"""

calc.set_molecule_geometry(geometry)

# Run complete spectroscopic analysis
calc.run_complete_spectroscopic_analysis()
```

### Custom Settings

```python
# Set custom computational parameters
calc.default_method = "M06-2X"           # DFT functional
calc.default_basis = "6-311G(d,p)"       # Basis set
calc.default_charge = 0                  # Molecular charge
calc.default_multiplicity = 1           # Spin multiplicity

# Run specific calculations
calc.generate_input_file('optimization', 'opt.gjf')
calc.generate_input_file('frequency', 'freq.gjf')
calc.generate_input_file('td_dft', 'td.gjf', additional_keywords="nstates=10")
```

## Available Calculation Types

| Calculation Type | Description | Gaussian Keywords |
|-----------------|-------------|-------------------|
| `optimization` | Geometry optimization | `opt` |
| `frequency` | Vibrational frequencies | `freq` |
| `raman` | Raman spectroscopy | `freq raman` |
| `td_dft` | UV-Visible via TD-DFT | `td` |
| `nmr` | NMR chemical shifts | `nmr` |
| `vcd` | Vibrational circular dichroism | `freq vcd` |
| `roa` | Raman optical activity | `freq roa` |
| `ecd` | Electronic circular dichroism | `td ecd` |
| `polarizability` | Linear polarizability | `polar` |
| `hyperpolarizability` | Nonlinear polarizability | `polar=gamma` |
| `anharmonic` | Anharmonic frequencies | `freq=anharmonic` |
| `solvent_freq` | Frequencies in solvent | `freq scrf=(solvent=water)` |

## Output Files

The suite generates several types of output files:

### Analysis Reports
- `{molecule_name}_analysis_report.txt` - Human-readable summary
- `{molecule_name}_detailed_results.json` - Complete data in JSON format

### CSV Data Files (for Excel/plotting)
- `{molecule_name}_frequencies.csv` - IR/Raman frequencies and intensities
- `{molecule_name}_uv_vis.csv` - UV-Visible transitions
- `{molecule_name}_nmr.csv` - NMR chemical shifts

### Gaussian Files
- `*.gjf` - Gaussian input files
- `*.out` - Gaussian output files

## Example Output

### Analysis Report Sample
```
GAUSSIAN SPECTROSCOPIC ANALYSIS REPORT
==================================================
Molecule: water_molecule
Analysis Date: 2025-01-15 14:30:22
Method: B3LYP/6-31G(d,p)

GEOMETRY OPTIMIZATION
--------------------
Final Energy: -76.419358 Hartree
Converged: Yes

VIBRATIONAL FREQUENCY ANALYSIS
------------------------------
Number of frequencies: 3
Frequencies (cm⁻¹):
    1:  1627.35 (IR intensity:  68.25)
    2:  3657.89 (IR intensity:  12.45)
    3:  3758.12 (IR intensity:  58.92)
Zero-point energy: 0.021468 Hartree

UV-VISIBLE SPECTROSCOPY (TD-DFT)
--------------------------------
Excited States:
  State 1: 167.2 nm (7.42 eV, f=0.0542)
  State 2: 132.8 nm (9.34 eV, f=0.1234)
```

## Examples

The `example_usage.py` file contains several examples:

1. **Water Molecule** - Basic complete analysis
2. **Benzene** - Custom method and basis set
3. **Methanol** - Selective calculations
4. **Acetone** - Solvent effects
5. **Template** - For your own molecules

Run examples:
```bash
python example_usage.py
```

## Advanced Usage

### Solvent Effects
```python
# Include solvation effects
calc.generate_input_file('td_dft', 'solvent.gjf', 
                        additional_keywords="scrf=(solvent=water)")
```

### High-Accuracy Calculations
```python
# Use higher-level methods
calc.default_method = "CCSD"
calc.default_basis = "cc-pVTZ"
```

### Large Molecules (ONIOM)
```python
# For very large systems, consider ONIOM
calc.generate_input_file('optimization', 'large.gjf',
                        additional_keywords="oniom(b3lyp/6-31g(d):uff)")
```

## Computational Methods Available

### DFT Functionals
- B3LYP (default)
- M06-2X
- ωB97X-D
- PBE0
- CAM-B3LYP

### Basis Sets
- 6-31G(d,p) (default)
- 6-311G(d,p)
- cc-pVDZ, cc-pVTZ
- def2-SVP, def2-TZVP

### Post-HF Methods
- MP2
- CCSD
- CCSD(T)

## Troubleshooting

### Common Issues

1. **Gaussian Not Found**
   - Ensure Gaussian is installed and `g16` command is accessible
   - Modify the command in the script if using different version

2. **Memory/Disk Issues**
   - Add memory specification: `additional_keywords="%mem=4GB"`
   - Add scratch directory: `additional_keywords="%chk=molecule.chk"`

3. **Convergence Problems**
   - Try different method: `calc.default_method = "B3LYP"`
   - Add convergence keywords: `additional_keywords="scf=xqc"`

### Performance Tips

1. **For Large Molecules**
   - Use smaller basis sets initially
   - Consider ONIOM for very large systems
   - Use parallel processing: `additional_keywords="%nprocshared=8"`

2. **For Many Calculations**
   - Run calculations in sequence to avoid memory conflicts
   - Use job scheduling systems for clusters

## Citation

If you use this suite in your research, please cite:
- Gaussian software: [Gaussian Citation](https://gaussian.com/citation/)
- This toolkit: Include reference to this repository

## License

This project is provided as-is for educational and research purposes. Please ensure you have proper licenses for Gaussian software.

## Support

For issues or questions:
1. Check the troubleshooting section
2. Review the example files
3. Consult Gaussian documentation for method-specific issues

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with detailed description

---

**Note**: This toolkit requires a valid Gaussian license. Gaussian is commercial software and must be purchased separately from Gaussian, Inc.