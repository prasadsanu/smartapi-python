#!/usr/bin/env python3
"""
Gaussian Spectroscopy Suite
===========================
A comprehensive Python script for automated spectroscopic calculations using Gaussian.
This script can generate input files, run calculations, and analyze results for various
spectroscopic properties.

Author: AI Assistant
Date: 2025
"""

import os
import sys
import subprocess
import re
import numpy as np
import pandas as pd
from pathlib import Path
import json
import time
from datetime import datetime

class GaussianSpectroscopyCalculator:
    """
    Main class for handling Gaussian spectroscopic calculations.
    This class manages input generation, job execution, and result analysis.
    """
    
    def __init__(self, molecule_name="molecule", base_dir="gaussian_calculations"):
        """
        Initialize the calculator with molecule name and working directory.
        
        Args:
            molecule_name (str): Name of the molecule being studied
            base_dir (str): Base directory for all calculations
        """
        self.molecule_name = molecule_name
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(exist_ok=True)
        
        # Dictionary to store all calculation results
        self.results = {}
        
        # Standard calculation parameters
        self.default_method = "B3LYP"
        self.default_basis = "6-31G(d,p)"
        self.default_charge = 0
        self.default_multiplicity = 1
        
        print(f"Initialized Gaussian Spectroscopy Calculator for {molecule_name}")
        print(f"Working directory: {self.base_dir.absolute()}")
    
    def set_molecule_geometry(self, geometry):
        """
        Set the molecular geometry for calculations.
        
        Args:
            geometry (str): Molecular geometry in Gaussian format
                          (either Cartesian coordinates or Z-matrix)
        """
        self.geometry = geometry
        print("Molecular geometry set successfully")
    
    def generate_input_file(self, calc_type, filename, additional_keywords="", 
                          method=None, basis=None, charge=None, multiplicity=None):
        """
        Generate Gaussian input file for specific calculation type.
        
        Args:
            calc_type (str): Type of calculation (opt, freq, td, etc.)
            filename (str): Output filename for the input file
            additional_keywords (str): Additional Gaussian keywords
            method (str): Quantum mechanical method (default: B3LYP)
            basis (str): Basis set (default: 6-31G(d,p))
            charge (int): Molecular charge (default: 0)
            multiplicity (int): Spin multiplicity (default: 1)
        """
        # Use defaults if not specified
        method = method or self.default_method
        basis = basis or self.default_basis
        charge = charge if charge is not None else self.default_charge
        multiplicity = multiplicity or self.default_multiplicity
        
        # Define calculation-specific keywords
        calc_keywords = {
            'optimization': 'opt',
            'frequency': 'freq',
            'td_dft': 'td',
            'nmr': 'nmr',
            'raman': 'freq raman',
            'vcd': 'freq vcd',
            'roa': 'freq roa',
            'ecd': 'td ecd',
            'uv_vis': 'td',
            'polarizability': 'polar',
            'hyperpolarizability': 'polar=gamma',
            'anharmonic': 'freq=anharmonic',
            'solvent_freq': 'freq scrf=(solvent=water)',
            'excited_opt': 'td opt',
            'irc': 'irc'
        }
        
        keywords = calc_keywords.get(calc_type, calc_type)
        if additional_keywords:
            keywords += f" {additional_keywords}"
        
        # Create input file content
        input_content = f"""#p {method}/{basis} {keywords}
        
{self.molecule_name} - {calc_type} calculation

{charge} {multiplicity}
{self.geometry}

"""
        
        # Write input file
        input_path = self.base_dir / filename
        with open(input_path, 'w') as f:
            f.write(input_content)
        
        print(f"Generated input file: {input_path}")
        return input_path
    
    def run_gaussian_calculation(self, input_file, output_file=None):
        """
        Execute Gaussian calculation.
        
        Args:
            input_file (str/Path): Path to Gaussian input file
            output_file (str/Path): Path to output file (optional)
        
        Returns:
            bool: True if calculation completed successfully
        """
        input_path = Path(input_file)
        
        if output_file is None:
            output_file = input_path.with_suffix('.out')
        
        try:
            # Run Gaussian calculation
            # Note: Adjust the command based on your Gaussian installation
            cmd = f"g16 < {input_path} > {output_file}"
            print(f"Running: {cmd}")
            
            # Execute the command
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0:
                print(f"Calculation completed successfully: {output_file}")
                return True
            else:
                print(f"Calculation failed with error: {result.stderr}")
                return False
                
        except Exception as e:
            print(f"Error running calculation: {e}")
            return False
    
    def parse_optimization_results(self, output_file):
        """
        Parse geometry optimization results from Gaussian output.
        
        Args:
            output_file (str/Path): Path to Gaussian output file
            
        Returns:
            dict: Parsed optimization data
        """
        results = {}
        
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Extract final energy
            energy_pattern = r'SCF Done:\s+E\([^)]+\)\s+=\s+([-\d.]+)'
            energy_matches = re.findall(energy_pattern, content)
            if energy_matches:
                results['final_energy'] = float(energy_matches[-1])
            
            # Extract optimization status
            if "Optimization completed" in content:
                results['converged'] = True
            else:
                results['converged'] = False
            
            # Extract final geometry
            geometry_pattern = r'Standard orientation:.*?\n(.*?)(?=\n\s*-{5}|\n\s*Rotational constants)'
            geometry_match = re.search(geometry_pattern, content, re.DOTALL)
            if geometry_match:
                results['final_geometry'] = geometry_match.group(1)
            
            print(f"Parsed optimization results from {output_file}")
            
        except Exception as e:
            print(f"Error parsing optimization results: {e}")
            results['error'] = str(e)
        
        return results
    
    def parse_frequency_results(self, output_file):
        """
        Parse vibrational frequency analysis results.
        
        Args:
            output_file (str/Path): Path to Gaussian output file
            
        Returns:
            dict: Parsed frequency data including IR and Raman
        """
        results = {}
        
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Extract vibrational frequencies
            freq_pattern = r'Frequencies --\s+([\d\s.-]+)'
            freq_matches = re.findall(freq_pattern, content)
            frequencies = []
            for match in freq_matches:
                freqs = [float(x) for x in match.split()]
                frequencies.extend(freqs)
            results['frequencies'] = frequencies
            
            # Extract IR intensities
            ir_pattern = r'IR Inten\s+--\s+([\d\s.-]+)'
            ir_matches = re.findall(ir_pattern, content)
            ir_intensities = []
            for match in ir_matches:
                intensities = [float(x) for x in match.split()]
                ir_intensities.extend(intensities)
            results['ir_intensities'] = ir_intensities
            
            # Extract Raman activities (if present)
            raman_pattern = r'Raman Activ --\s+([\d\s.-]+)'
            raman_matches = re.findall(raman_pattern, content)
            raman_activities = []
            for match in raman_matches:
                activities = [float(x) for x in match.split()]
                raman_activities.extend(activities)
            if raman_activities:
                results['raman_activities'] = raman_activities
            
            # Extract thermochemical data
            if "Zero-point correction=" in content:
                zpe_match = re.search(r'Zero-point correction=\s+([\d.-]+)', content)
                if zpe_match:
                    results['zero_point_energy'] = float(zpe_match.group(1))
            
            if "Thermal correction to Enthalpy=" in content:
                enthalpy_match = re.search(r'Thermal correction to Enthalpy=\s+([\d.-]+)', content)
                if enthalpy_match:
                    results['enthalpy_correction'] = float(enthalpy_match.group(1))
            
            if "Thermal correction to Gibbs Free Energy=" in content:
                gibbs_match = re.search(r'Thermal correction to Gibbs Free Energy=\s+([\d.-]+)', content)
                if gibbs_match:
                    results['gibbs_correction'] = float(gibbs_match.group(1))
            
            print(f"Parsed frequency results from {output_file}")
            
        except Exception as e:
            print(f"Error parsing frequency results: {e}")
            results['error'] = str(e)
        
        return results
    
    def parse_td_dft_results(self, output_file):
        """
        Parse TD-DFT results for UV-Vis and electronic transitions.
        
        Args:
            output_file (str/Path): Path to Gaussian output file
            
        Returns:
            dict: Parsed TD-DFT data
        """
        results = {}
        
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Extract excited states
            excited_states = []
            state_pattern = r'Excited State\s+(\d+):\s+\S+\s+([\d.]+)\s+eV\s+([\d.]+)\s+nm\s+f=([\d.]+)'
            state_matches = re.findall(state_pattern, content)
            
            for match in state_matches:
                state_num, energy_ev, wavelength_nm, osc_strength = match
                excited_states.append({
                    'state': int(state_num),
                    'energy_ev': float(energy_ev),
                    'wavelength_nm': float(wavelength_nm),
                    'oscillator_strength': float(osc_strength)
                })
            
            results['excited_states'] = excited_states
            
            # Extract ground state dipole moment
            dipole_pattern = r'Dipole moment \(field-independent basis, Debye\):\s*\n\s*X=\s*([-\d.]+)\s*Y=\s*([-\d.]+)\s*Z=\s*([-\d.]+)\s*Tot=\s*([-\d.]+)'
            dipole_match = re.search(dipole_pattern, content)
            if dipole_match:
                results['ground_state_dipole'] = {
                    'x': float(dipole_match.group(1)),
                    'y': float(dipole_match.group(2)),
                    'z': float(dipole_match.group(3)),
                    'total': float(dipole_match.group(4))
                }
            
            print(f"Parsed TD-DFT results from {output_file}")
            
        except Exception as e:
            print(f"Error parsing TD-DFT results: {e}")
            results['error'] = str(e)
        
        return results
    
    def parse_nmr_results(self, output_file):
        """
        Parse NMR chemical shifts and coupling constants.
        
        Args:
            output_file (str/Path): Path to Gaussian output file
            
        Returns:
            dict: Parsed NMR data
        """
        results = {}
        
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Extract NMR shielding tensors
            shielding_pattern = r'(\d+)\s+([A-Z][a-z]?)\s+Isotropic =\s+([-\d.]+)\s+Anisotropy =\s+([-\d.]+)'
            shielding_matches = re.findall(shielding_pattern, content)
            
            chemical_shifts = []
            for match in shielding_matches:
                atom_num, element, isotropic, anisotropy = match
                chemical_shifts.append({
                    'atom_number': int(atom_num),
                    'element': element,
                    'isotropic_shielding': float(isotropic),
                    'anisotropy': float(anisotropy)
                })
            
            results['chemical_shifts'] = chemical_shifts
            
            print(f"Parsed NMR results from {output_file}")
            
        except Exception as e:
            print(f"Error parsing NMR results: {e}")
            results['error'] = str(e)
        
        return results
    
    def run_complete_spectroscopic_analysis(self):
        """
        Run a complete suite of spectroscopic calculations.
        This is the main method that coordinates all calculations.
        """
        print("\n" + "="*60)
        print("STARTING COMPLETE SPECTROSCOPIC ANALYSIS")
        print("="*60)
        
        if not hasattr(self, 'geometry'):
            print("ERROR: No molecular geometry set. Use set_molecule_geometry() first.")
            return
        
        # Dictionary to track calculation status
        calculations = {
            'optimization': {'input': 'opt.gjf', 'output': 'opt.out'},
            'frequency': {'input': 'freq.gjf', 'output': 'freq.out'},
            'raman': {'input': 'raman.gjf', 'output': 'raman.out'},
            'td_dft': {'input': 'td.gjf', 'output': 'td.out'},
            'nmr': {'input': 'nmr.gjf', 'output': 'nmr.out'},
            'vcd': {'input': 'vcd.gjf', 'output': 'vcd.out'},
            'polarizability': {'input': 'polar.gjf', 'output': 'polar.out'},
            'anharmonic': {'input': 'anharm.gjf', 'output': 'anharm.out'}
        }
        
        # Generate all input files
        print("\nGenerating input files...")
        for calc_type, files in calculations.items():
            try:
                self.generate_input_file(calc_type, files['input'])
            except Exception as e:
                print(f"Error generating {calc_type} input: {e}")
        
        # Run calculations (in sequence to avoid conflicts)
        print("\nRunning calculations...")
        for calc_type, files in calculations.items():
            print(f"\nStarting {calc_type} calculation...")
            input_path = self.base_dir / files['input']
            output_path = self.base_dir / files['output']
            
            if input_path.exists():
                success = self.run_gaussian_calculation(input_path, output_path)
                calculations[calc_type]['success'] = success
                
                # Parse results if calculation was successful
                if success and output_path.exists():
                    if calc_type == 'optimization':
                        self.results[calc_type] = self.parse_optimization_results(output_path)
                    elif calc_type in ['frequency', 'raman', 'anharmonic']:
                        self.results[calc_type] = self.parse_frequency_results(output_path)
                    elif calc_type == 'td_dft':
                        self.results[calc_type] = self.parse_td_dft_results(output_path)
                    elif calc_type == 'nmr':
                        self.results[calc_type] = self.parse_nmr_results(output_path)
                    else:
                        # For other calculations, store basic info
                        self.results[calc_type] = {'completed': True}
        
        # Generate comprehensive report
        self.generate_analysis_report()
        
        print("\n" + "="*60)
        print("SPECTROSCOPIC ANALYSIS COMPLETED")
        print("="*60)
    
    def generate_analysis_report(self):
        """
        Generate a comprehensive analysis report of all spectroscopic data.
        Creates both detailed JSON output and simplified summary files.
        """
        print("\nGenerating analysis report...")
        
        # Create detailed JSON report
        json_report_path = self.base_dir / f"{self.molecule_name}_detailed_results.json"
        with open(json_report_path, 'w') as f:
            json.dump(self.results, f, indent=2, default=str)
        
        # Create simplified text report
        report_path = self.base_dir / f"{self.molecule_name}_analysis_report.txt"
        
        with open(report_path, 'w') as f:
            f.write(f"GAUSSIAN SPECTROSCOPIC ANALYSIS REPORT\n")
            f.write(f"{'='*50}\n")
            f.write(f"Molecule: {self.molecule_name}\n")
            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Method: {self.default_method}/{self.default_basis}\n\n")
            
            # Optimization results
            if 'optimization' in self.results:
                opt_data = self.results['optimization']
                f.write("GEOMETRY OPTIMIZATION\n")
                f.write("-"*20 + "\n")
                if 'final_energy' in opt_data:
                    f.write(f"Final Energy: {opt_data['final_energy']:.6f} Hartree\n")
                if 'converged' in opt_data:
                    f.write(f"Converged: {'Yes' if opt_data['converged'] else 'No'}\n")
                f.write("\n")
            
            # Frequency analysis results
            if 'frequency' in self.results:
                freq_data = self.results['frequency']
                f.write("VIBRATIONAL FREQUENCY ANALYSIS\n")
                f.write("-"*30 + "\n")
                if 'frequencies' in freq_data:
                    f.write(f"Number of frequencies: {len(freq_data['frequencies'])}\n")
                    f.write("Frequencies (cm⁻¹):\n")
                    for i, freq in enumerate(freq_data['frequencies'][:10]):  # Show first 10
                        intensity = freq_data.get('ir_intensities', [0]*(i+1))[i] if i < len(freq_data.get('ir_intensities', [])) else 0
                        f.write(f"  {i+1:3d}: {freq:8.2f} (IR intensity: {intensity:6.2f})\n")
                    if len(freq_data['frequencies']) > 10:
                        f.write(f"  ... and {len(freq_data['frequencies'])-10} more frequencies\n")
                
                if 'zero_point_energy' in freq_data:
                    f.write(f"Zero-point energy: {freq_data['zero_point_energy']:.6f} Hartree\n")
                f.write("\n")
            
            # TD-DFT results
            if 'td_dft' in self.results:
                td_data = self.results['td_dft']
                f.write("UV-VISIBLE SPECTROSCOPY (TD-DFT)\n")
                f.write("-"*30 + "\n")
                if 'excited_states' in td_data:
                    f.write("Excited States:\n")
                    for state in td_data['excited_states'][:5]:  # Show first 5 states
                        f.write(f"  State {state['state']}: {state['wavelength_nm']:.1f} nm "
                               f"({state['energy_ev']:.2f} eV, f={state['oscillator_strength']:.4f})\n")
                f.write("\n")
            
            # NMR results
            if 'nmr' in self.results:
                nmr_data = self.results['nmr']
                f.write("NMR CHEMICAL SHIFTS\n")
                f.write("-"*20 + "\n")
                if 'chemical_shifts' in nmr_data:
                    f.write("Chemical Shifts (ppm from TMS):\n")
                    for shift in nmr_data['chemical_shifts'][:10]:  # Show first 10
                        # Convert shielding to chemical shift (approximate)
                        chem_shift = 200 - shift['isotropic_shielding']  # Rough conversion
                        f.write(f"  Atom {shift['atom_number']} ({shift['element']}): {chem_shift:.2f} ppm\n")
                f.write("\n")
            
            # Summary
            f.write("CALCULATION SUMMARY\n")
            f.write("-"*20 + "\n")
            for calc_type in self.results:
                status = "✓ Completed" if calc_type in self.results else "✗ Failed"
                f.write(f"{calc_type.capitalize()}: {status}\n")
        
        # Create CSV files for easy data analysis
        self.create_csv_outputs()
        
        print(f"Analysis report saved to: {report_path}")
        print(f"Detailed JSON data saved to: {json_report_path}")
    
    def create_csv_outputs(self):
        """
        Create CSV files for different types of spectroscopic data.
        This makes the data easy to import into Excel or other analysis tools.
        """
        # IR/Raman frequencies CSV
        if 'frequency' in self.results:
            freq_data = self.results['frequency']
            if 'frequencies' in freq_data:
                df_freq = pd.DataFrame({
                    'Frequency_cm-1': freq_data['frequencies'],
                    'IR_Intensity': freq_data.get('ir_intensities', [0]*len(freq_data['frequencies'])),
                    'Raman_Activity': freq_data.get('raman_activities', [0]*len(freq_data['frequencies']))
                })
                df_freq.to_csv(self.base_dir / f"{self.molecule_name}_frequencies.csv", index=False)
        
        # UV-Vis transitions CSV
        if 'td_dft' in self.results and 'excited_states' in self.results['td_dft']:
            df_uv = pd.DataFrame(self.results['td_dft']['excited_states'])
            df_uv.to_csv(self.base_dir / f"{self.molecule_name}_uv_vis.csv", index=False)
        
        # NMR chemical shifts CSV
        if 'nmr' in self.results and 'chemical_shifts' in self.results['nmr']:
            df_nmr = pd.DataFrame(self.results['nmr']['chemical_shifts'])
            df_nmr.to_csv(self.base_dir / f"{self.molecule_name}_nmr.csv", index=False)
        
        print("CSV data files created for easy analysis")

def main():
    """
    Main function demonstrating how to use the Gaussian Spectroscopy Calculator.
    This example uses water molecule for demonstration.
    """
    print("Gaussian Spectroscopy Suite - Demo")
    print("="*40)
    
    # Initialize calculator
    calc = GaussianSpectroscopyCalculator("water_molecule")
    
    # Set molecular geometry (water molecule example)
    water_geometry = """O          0.000000    0.000000    0.117176
H          0.000000    0.757480   -0.468706
H          0.000000   -0.757480   -0.468706"""
    
    calc.set_molecule_geometry(water_geometry)
    
    # Run complete spectroscopic analysis
    calc.run_complete_spectroscopic_analysis()
    
    print("\nDemo completed! Check the 'gaussian_calculations' directory for results.")
    print("\nTo use with your own molecule:")
    print("1. Replace the geometry with your molecule's coordinates")
    print("2. Adjust the method and basis set if needed")
    print("3. Make sure Gaussian is properly installed and accessible")

if __name__ == "__main__":
    # Check if pandas is available
    try:
        import pandas as pd
    except ImportError:
        print("Warning: pandas not found. CSV output will be limited.")
        print("Install with: pip install pandas")
    
    main()