#!/usr/bin/env python3
"""
Example Usage of Gaussian Spectroscopy Suite
============================================
This script demonstrates various ways to use the GaussianSpectroscopyCalculator
for different types of molecules and calculations.
"""

from gaussian_spectroscopy_suite import GaussianSpectroscopyCalculator

def example_water_molecule():
    """
    Example 1: Complete spectroscopic analysis of water molecule
    This demonstrates the basic usage with a simple molecule
    """
    print("="*60)
    print("EXAMPLE 1: Water Molecule Analysis")
    print("="*60)
    
    # Initialize calculator for water
    calc = GaussianSpectroscopyCalculator("water", "water_analysis")
    
    # Water molecule geometry (optimized structure)
    water_geometry = """O          0.000000    0.000000    0.117176
H          0.000000    0.757480   -0.468706
H          0.000000   -0.757480   -0.468706"""
    
    calc.set_molecule_geometry(water_geometry)
    
    # Run complete analysis
    calc.run_complete_spectroscopic_analysis()

def example_benzene_molecule():
    """
    Example 2: Benzene molecule with custom method
    This shows how to use different computational methods
    """
    print("="*60)
    print("EXAMPLE 2: Benzene Molecule with Custom Settings")
    print("="*60)
    
    # Initialize with custom settings
    calc = GaussianSpectroscopyCalculator("benzene", "benzene_analysis")
    
    # Set custom method and basis set
    calc.default_method = "M06-2X"  # Different DFT functional
    calc.default_basis = "6-311G(d,p)"  # Larger basis set
    
    # Benzene molecule geometry
    benzene_geometry = """C          0.000000    1.396000    0.000000
C          1.209000    0.698000    0.000000
C          1.209000   -0.698000    0.000000
C          0.000000   -1.396000    0.000000
C         -1.209000   -0.698000    0.000000
C         -1.209000    0.698000    0.000000
H          0.000000    2.480000    0.000000
H          2.148000    1.240000    0.000000
H          2.148000   -1.240000    0.000000
H          0.000000   -2.480000    0.000000
H         -2.148000   -1.240000    0.000000
H         -2.148000    1.240000    0.000000"""
    
    calc.set_molecule_geometry(benzene_geometry)
    
    # Run complete analysis
    calc.run_complete_spectroscopic_analysis()

def example_selective_calculations():
    """
    Example 3: Running only specific calculations
    This shows how to run individual calculation types
    """
    print("="*60)
    print("EXAMPLE 3: Selective Calculations for Methanol")
    print("="*60)
    
    calc = GaussianSpectroscopyCalculator("methanol", "methanol_selective")
    
    # Methanol geometry
    methanol_geometry = """C          0.659000    0.000000    0.000000
O         -0.759000    0.000000    0.000000
H          1.026000    0.000000    1.026000
H          1.026000    0.889000   -0.513000
H          1.026000   -0.889000   -0.513000
H         -1.124000    0.000000    0.889000"""
    
    calc.set_molecule_geometry(methanol_geometry)
    
    # Run only specific calculations
    print("\nRunning geometry optimization...")
    opt_input = calc.generate_input_file('optimization', 'methanol_opt.gjf')
    calc.run_gaussian_calculation(opt_input, 'methanol_opt.out')
    
    print("\nRunning frequency analysis...")
    freq_input = calc.generate_input_file('frequency', 'methanol_freq.gjf')
    calc.run_gaussian_calculation(freq_input, 'methanol_freq.out')
    
    print("\nRunning TD-DFT for UV-Vis...")
    td_input = calc.generate_input_file('td_dft', 'methanol_td.gjf', 
                                       additional_keywords="nstates=10")
    calc.run_gaussian_calculation(td_input, 'methanol_td.out')

def example_solvent_effects():
    """
    Example 4: Studying solvent effects
    This demonstrates how to include solvation in calculations
    """
    print("="*60)
    print("EXAMPLE 4: Solvent Effects on Acetone")
    print("="*60)
    
    calc = GaussianSpectroscopyCalculator("acetone_solvent", "acetone_solvent_analysis")
    
    # Acetone geometry
    acetone_geometry = """C          0.000000    0.000000    0.000000
C          1.507000    0.000000    0.000000
C         -1.507000    0.000000    0.000000
O          0.000000    1.210000    0.000000
H          2.141000   -0.889000    0.000000
H          2.141000    0.889000    0.000000
H         -2.141000    0.889000    0.000000
H         -2.141000   -0.889000    0.000000"""
    
    calc.set_molecule_geometry(acetone_geometry)
    
    # Generate input files with different solvents
    print("\nGenerating gas phase calculation...")
    gas_input = calc.generate_input_file('td_dft', 'acetone_gas.gjf', 
                                        additional_keywords="nstates=5")
    
    print("\nGenerating water solvation calculation...")
    water_input = calc.generate_input_file('td_dft', 'acetone_water.gjf',
                                         additional_keywords="nstates=5 scrf=(solvent=water)")
    
    print("\nGenerating acetonitrile solvation calculation...")
    acn_input = calc.generate_input_file('td_dft', 'acetone_acn.gjf',
                                       additional_keywords="nstates=5 scrf=(solvent=acetonitrile)")
    
    # Note: In practice, you would run these calculations and compare results

def example_custom_molecule():
    """
    Example 5: Template for your own molecule
    Copy and modify this function for your specific molecule
    """
    print("="*60)
    print("EXAMPLE 5: Template for Custom Molecule")
    print("="*60)
    
    # Initialize calculator with your molecule name
    calc = GaussianSpectroscopyCalculator("your_molecule_name", "your_analysis_directory")
    
    # Set your computational parameters
    calc.default_method = "B3LYP"        # Choose your method
    calc.default_basis = "6-31G(d,p)"    # Choose your basis set
    calc.default_charge = 0              # Set molecular charge
    calc.default_multiplicity = 1       # Set spin multiplicity
    
    # Your molecule geometry (replace with your coordinates)
    your_geometry = """C          0.000000    0.000000    0.000000
H          0.000000    0.000000    1.089000
H          1.026719    0.000000   -0.363000
H         -0.513360    0.889165   -0.363000
H         -0.513360   -0.889165   -0.363000"""
    
    calc.set_molecule_geometry(your_geometry)
    
    # Run the analysis you need
    print("\nTo run complete analysis, uncomment the next line:")
    # calc.run_complete_spectroscopic_analysis()
    
    print("Or run specific calculations as needed:")
    print("- calc.generate_input_file('optimization', 'opt.gjf')")
    print("- calc.generate_input_file('frequency', 'freq.gjf')")
    print("- calc.generate_input_file('td_dft', 'td.gjf')")

def main():
    """
    Run all examples
    """
    print("Gaussian Spectroscopy Suite - Examples")
    print("="*60)
    print("This script demonstrates various usage patterns.")
    print("Uncomment the examples you want to run.\n")
    
    # Uncomment the examples you want to run:
    
    # Example 1: Basic water molecule analysis
    # example_water_molecule()
    
    # Example 2: Benzene with custom settings
    # example_benzene_molecule()
    
    # Example 3: Selective calculations
    # example_selective_calculations()
    
    # Example 4: Solvent effects
    # example_solvent_effects()
    
    # Example 5: Template for your molecule
    example_custom_molecule()
    
    print("\nTo run actual calculations:")
    print("1. Make sure Gaussian is installed and accessible")
    print("2. Uncomment the desired examples above")
    print("3. Run this script: python example_usage.py")

if __name__ == "__main__":
    main()