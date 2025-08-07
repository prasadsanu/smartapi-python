#!/usr/bin/env python3
"""
Gaussian Spectroscopy Suite - Fixed Version
==========================================
A comprehensive Python script for automated spectroscopic calculations using Gaussian.
This script can generate input files, run calculations, and analyze results for various
spectroscopic properties.

Fixed issues:
- Improved Raman and VCD calculation keywords
- Better error handling
- More robust output parsing
- Optional calculation execution

Author: AI Assistant
Date: 2025
"""

import os
import sys
import subprocess
import re
import numpy as np
try:
    import pandas as pd
except ImportError:
    pd = None  # pandas is optional, CSV export disabled if missing

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
        
        # Gaussian executable command (modify if needed)
        self.gaussian_command = "g16"  # Change to "g09" or full path if needed
        
        print(f"Initialized Gaussian Spectroscopy Calculator for {molecule_name}")
        print(f"Working directory: {self.base_dir.absolute()}")
    
    def set_molecule_geometry(self, geometry):
        """
        Set the molecular geometry for calculations.
        
        Args:
            geometry (str): Molecular geometry in Gaussian format
                          (either Cartesian coordinates or Z-matrix)
        """
        self.geometry = geometry.strip()
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
        
        # Define calculation-specific keywords with improved formatting
        calc_keywords = {
            'optimization': 'opt',
            'frequency': 'freq',
            'td_dft': 'td=(nstates=10)',
            'nmr': 'nmr',
            'raman': 'freq=(raman)',  # Fixed: proper parentheses
            'vcd': 'freq=(vcd)',      # Fixed: proper parentheses
            'roa': 'freq=(roa)',      # Fixed: proper parentheses
            'ecd': 'td=(ecd,nstates=10)',  # Fixed: proper format
            'uv_vis': 'td=(nstates=10)',
            'polarizability': 'polar',
            'hyperpolarizability': 'polar=(gamma)',
            'anharmonic': 'freq=(anharmonic)',  # Fixed: proper parentheses
            'solvent_freq': 'freq scrf=(solvent=water)',
            'excited_opt': 'td=(nstates=5) opt',
            'irc': 'irc=(calcfc,maxpoints=20)',
            'opt_freq': 'opt freq'  # Combined optimization and frequency
        }
        
        keywords = calc_keywords.get(calc_type, calc_type)
        if additional_keywords:
            keywords += f" {additional_keywords}"
        
        # Add memory and processor specifications for stability
        memory_line = "%mem=2GB"
        proc_line = "%nprocshared=2"
        
        # Create input file content with proper formatting
        input_content = f"""{memory_line}
{proc_line}
#p {method}/{basis} {keywords}

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
    
    def run_gaussian_calculation(self, input_file, output_file=None, timeout=3600):
        """
        Execute Gaussian calculation with timeout and better error handling.
        
        Args:
            input_file (str/Path): Path to Gaussian input file
            output_file (str/Path): Path to output file (optional)
            timeout (int): Maximum time to wait for calculation (seconds)
        
        Returns:
            bool: True if calculation completed successfully
        """
        input_path = Path(input_file)
        
        if output_file is None:
            output_file = input_path.with_suffix('.out')
        
        try:
            # Check if Gaussian is available
            check_cmd = f"which {self.gaussian_command}"
            check_result = subprocess.run(check_cmd, shell=True, capture_output=True, text=True)
            
            if check_result.returncode != 0:
                print(f"Warning: {self.gaussian_command} not found in PATH. Skipping calculation.")
                print("To run calculations, ensure Gaussian is properly installed and accessible.")
                return False
            
            # Run Gaussian calculation
            cmd = f"{self.gaussian_command} < {input_path} > {output_file}"
            print(f"Running: {cmd}")
            
            # Execute the command with timeout
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
            
            if result.returncode == 0:
                print(f"Calculation completed successfully: {output_file}")
                
                # Check if output file contains "Normal termination"
                if Path(output_file).exists():
                    with open(output_file, 'r') as f:
                        content = f.read()
                    if "Normal termination" in content:
                        return True
                    else:
                        print(f"Warning: Calculation may not have completed normally")
                        return False
                return True
            else:
                print(f"Calculation failed with return code {result.returncode}")
                if result.stderr:
                    print(f"Error output: {result.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            print(f"Calculation timed out after {timeout} seconds")
            return False
        except Exception as e:
            print(f"Error running calculation: {e}")
            return False
    
    def parse_optimization_results(self, output_file):
        """
        Parse geometry optimization results from Gaussian output.
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
            results['converged'] = "Optimization completed" in content
            
            # Extract number of optimization steps
            step_pattern = r'Step number\s+(\d+)'
            step_matches = re.findall(step_pattern, content)
            if step_matches:
                results['optimization_steps'] = int(step_matches[-1])
            
            # Extract final geometry coordinates
            geom_pattern = r'Standard orientation:.*?\n.*?\n.*?\n.*?\n.*?\n(.*?)(?=\n\s*-{5})'
            geom_match = re.search(geom_pattern, content, re.DOTALL)
            if geom_match:
                geom_lines = geom_match.group(1).strip().split('\n')
                coordinates = []
                for line in geom_lines:
                    if line.strip() and not line.startswith('-'):
                        parts = line.split()
                        if len(parts) >= 6:
                            coordinates.append({
                                'atom_number': int(parts[0]),
                                'atomic_number': int(parts[1]),
                                'x': float(parts[3]),
                                'y': float(parts[4]),
                                'z': float(parts[5])
                            })
                results['final_coordinates'] = coordinates
            
            print(f"Parsed optimization results from {output_file}")
            
        except Exception as e:
            print(f"Error parsing optimization results: {e}")
            results['error'] = str(e)
        
        return results
    
    def parse_frequency_results(self, output_file):
        """
        Parse vibrational frequency analysis results with improved patterns.
        """
        results = {}
        
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Extract vibrational frequencies
            freq_pattern = r'Frequencies\s+--\s+([\d\s.-]+)'
            freq_matches = re.findall(freq_pattern, content)
            frequencies = []
            for match in freq_matches:
                freqs = [float(x) for x in match.split() if x.replace('.', '').replace('-', '').isdigit()]
                frequencies.extend(freqs)
            results['frequencies'] = frequencies
            
            # Extract IR intensities
            ir_pattern = r'IR Inten\s+--\s+([\d\s.-]+)'
            ir_matches = re.findall(ir_pattern, content)
            ir_intensities = []
            for match in ir_matches:
                intensities = [float(x) for x in match.split() if x.replace('.', '').replace('-', '').isdigit()]
                ir_intensities.extend(intensities)
            results['ir_intensities'] = ir_intensities
            
            # Extract Raman activities with improved pattern
            raman_pattern = r'Raman Activ\s+--\s+([\d\s.-]+)'
            raman_matches = re.findall(raman_pattern, content)
            raman_activities = []
            for match in raman_matches:
                activities = [float(x) for x in match.split() if x.replace('.', '').replace('-', '').isdigit()]
                raman_activities.extend(activities)
            if raman_activities:
                results['raman_activities'] = raman_activities
            
            # Extract VCD intensities (rotational strengths)
            vcd_pattern = r'Rot\. str\.\s+--\s+([\d\s.-]+)'
            vcd_matches = re.findall(vcd_pattern, content)
            vcd_intensities = []
            for match in vcd_matches:
                intensities = [float(x) for x in match.split() if x.replace('.', '').replace('-', '').isdigit()]
                vcd_intensities.extend(intensities)
            if vcd_intensities:
                results['vcd_intensities'] = vcd_intensities
            
            # Extract ROA intensities
            roa_pattern = r'ROA Inten\s+--\s+([\d\s.-]+)'
            roa_matches = re.findall(roa_pattern, content)
            roa_intensities = []
            for match in roa_matches:
                intensities = [float(x) for x in match.split() if x.replace('.', '').replace('-', '').isdigit()]
                roa_intensities.extend(intensities)
            if roa_intensities:
                results['roa_intensities'] = roa_intensities
            
            # Extract thermochemical data
            thermo_patterns = {
                'zero_point_energy': r'Zero-point correction=\s+([\d.-]+)',
                'enthalpy_correction': r'Thermal correction to Enthalpy=\s+([\d.-]+)',
                'gibbs_correction': r'Thermal correction to Gibbs Free Energy=\s+([\d.-]+)',
                'electronic_energy': r'Sum of electronic and zero-point Energies=\s+([-\d.]+)',
                'enthalpy': r'Sum of electronic and thermal Enthalpies=\s+([-\d.]+)',
                'gibbs_energy': r'Sum of electronic and thermal Free Energies=\s+([-\d.]+)'
            }
            
            for key, pattern in thermo_patterns.items():
                match = re.search(pattern, content)
                if match:
                    results[key] = float(match.group(1))
            
            print(f"Parsed frequency results from {output_file}")
            
        except Exception as e:
            print(f"Error parsing frequency results: {e}")
            results['error'] = str(e)
        
        return results
    
    def parse_td_dft_results(self, output_file):
        """
        Parse TD-DFT results for UV-Vis and electronic transitions.
        """
        results = {}
        
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Extract excited states with improved pattern
            excited_states = []
            state_pattern = r'Excited State\s+(\d+):\s+(\S+)\s+([\d.]+)\s+eV\s+([\d.]+)\s+nm\s+f=([\d.]+)'
            state_matches = re.findall(state_pattern, content)
            
            for match in state_matches:
                state_num, symmetry, energy_ev, wavelength_nm, osc_strength = match
                excited_states.append({
                    'state': int(state_num),
                    'symmetry': symmetry,
                    'energy_ev': float(energy_ev),
                    'wavelength_nm': float(wavelength_nm),
                    'oscillator_strength': float(osc_strength)
                })
            
            results['excited_states'] = excited_states
            
            # Extract ground state dipole moment
            dipole_patterns = [
                r'Dipole moment \(field-independent basis, Debye\):\s*\n\s*X=\s*([-\d.]+)\s*Y=\s*([-\d.]+)\s*Z=\s*([-\d.]+)\s*Tot=\s*([-\d.]+)',
                r'X=\s*([-\d.]+)\s*Y=\s*([-\d.]+)\s*Z=\s*([-\d.]+)\s*Tot=\s*([-\d.]+)'
            ]
            
            for pattern in dipole_patterns:
                dipole_match = re.search(pattern, content)
                if dipole_match:
                    results['ground_state_dipole'] = {
                        'x': float(dipole_match.group(1)),
                        'y': float(dipole_match.group(2)),
                        'z': float(dipole_match.group(3)),
                        'total': float(dipole_match.group(4))
                    }
                    break
            
            # Extract ECD data if present
            ecd_pattern = r'R\(length\)\s+([-\d.]+)'
            ecd_matches = re.findall(ecd_pattern, content)
            if ecd_matches:
                results['ecd_rotational_strengths'] = [float(x) for x in ecd_matches]
            
            print(f"Parsed TD-DFT results from {output_file}")
            
        except Exception as e:
            print(f"Error parsing TD-DFT results: {e}")
            results['error'] = str(e)
        
        return results
    
    def parse_nmr_results(self, output_file):
        """
        Parse NMR chemical shifts and coupling constants.
        """
        results = {}
        
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Extract NMR shielding tensors with improved pattern
            shielding_pattern = r'(\d+)\s+([A-Z][a-z]?)\s+Isotropic\s+=\s+([-\d.]+)\s+Anisotropy\s+=\s+([-\d.]+)'
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
            
            # Extract spin-spin coupling constants if present
            coupling_pattern = r'(\d+)\s+(\d+)\s+([-\d.]+)\s+Hz'
            coupling_matches = re.findall(coupling_pattern, content)
            if coupling_matches:
                coupling_constants = []
                for match in coupling_matches:
                    atom1, atom2, coupling = match
                    coupling_constants.append({
                        'atom1': int(atom1),
                        'atom2': int(atom2),
                        'coupling_hz': float(coupling)
                    })
                results['coupling_constants'] = coupling_constants
            
            print(f"Parsed NMR results from {output_file}")
            
        except Exception as e:
            print(f"Error parsing NMR results: {e}")
            results['error'] = str(e)
        
        return results
    
    def run_complete_spectroscopic_analysis(self, run_calculations=True):
        """
        Run a complete suite of spectroscopic calculations.
        
        Args:
            run_calculations (bool): If False, only generate input files without running
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
            'polarizability': {'input': 'polar.gjf', 'output': 'polar.out'}
        }
        
        # Generate all input files
        print("\nGenerating input files...")
        for calc_type, files in calculations.items():
            try:
                self.generate_input_file(calc_type, files['input'])
            except Exception as e:
                print(f"Error generating {calc_type} input: {e}")
        
        if not run_calculations:
            print("\nInput files generated. Set run_calculations=True to execute them.")
            return
        
        # Run calculations (in sequence to avoid conflicts)
        print("\nRunning calculations...")
        for calc_type, files in calculations.items():
            print(f"\n{'='*40}")
            print(f"Starting {calc_type.upper()} calculation...")
            print('='*40)
            
            input_path = self.base_dir / files['input']
            output_path = self.base_dir / files['output']
            
            if input_path.exists():
                success = self.run_gaussian_calculation(input_path, output_path)
                calculations[calc_type]['success'] = success
                
                # Parse results if calculation was successful and output exists
                if success and output_path.exists():
                    try:
                        if calc_type == 'optimization':
                            self.results[calc_type] = self.parse_optimization_results(output_path)
                        elif calc_type in ['frequency', 'raman', 'vcd']:
                            self.results[calc_type] = self.parse_frequency_results(output_path)
                        elif calc_type == 'td_dft':
                            self.results[calc_type] = self.parse_td_dft_results(output_path)
                        elif calc_type == 'nmr':
                            self.results[calc_type] = self.parse_nmr_results(output_path)
                        else:
                            # For other calculations, store basic info
                            self.results[calc_type] = {'completed': True}
                    except Exception as e:
                        print(f"Error parsing {calc_type} results: {e}")
                        self.results[calc_type] = {'error': str(e)}
                else:
                    self.results[calc_type] = {'error': 'Calculation failed or output not found'}
            else:
                print(f"Input file not found: {input_path}")
        
        # Generate comprehensive report
        self.generate_analysis_report()
        
        print("\n" + "="*60)
        print("SPECTROSCOPIC ANALYSIS COMPLETED")
        print("="*60)
        print(f"Check the '{self.base_dir}' directory for all results.")
    
    def generate_analysis_report(self):
        """
        Generate a comprehensive analysis report of all spectroscopic data.
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
            if 'optimization' in self.results and 'error' not in self.results['optimization']:
                opt_data = self.results['optimization']
                f.write("GEOMETRY OPTIMIZATION\n")
                f.write("-"*20 + "\n")
                if 'final_energy' in opt_data:
                    f.write(f"Final Energy: {opt_data['final_energy']:.6f} Hartree\n")
                if 'converged' in opt_data:
                    f.write(f"Converged: {'Yes' if opt_data['converged'] else 'No'}\n")
                if 'optimization_steps' in opt_data:
                    f.write(f"Optimization Steps: {opt_data['optimization_steps']}\n")
                f.write("\n")
            
            # Frequency analysis results
            if 'frequency' in self.results and 'error' not in self.results['frequency']:
                freq_data = self.results['frequency']
                f.write("VIBRATIONAL FREQUENCY ANALYSIS\n")
                f.write("-"*30 + "\n")
                if 'frequencies' in freq_data and freq_data['frequencies']:
                    f.write(f"Number of frequencies: {len(freq_data['frequencies'])}\n")
                    f.write("Frequencies (cm⁻¹):\n")
                    for i, freq in enumerate(freq_data['frequencies'][:10]):  # Show first 10
                        intensity = freq_data.get('ir_intensities', [0]*(i+1))[i] if i < len(freq_data.get('ir_intensities', [])) else 0
                        f.write(f"  {i+1:3d}: {freq:8.2f} (IR intensity: {intensity:8.2f})\n")
                    if len(freq_data['frequencies']) > 10:
                        f.write(f"  ... and {len(freq_data['frequencies'])-10} more frequencies\n")
                
                if 'zero_point_energy' in freq_data:
                    f.write(f"Zero-point energy: {freq_data['zero_point_energy']:.6f} Hartree\n")
                f.write("\n")
            
            # Raman results
            if 'raman' in self.results and 'error' not in self.results['raman']:
                raman_data = self.results['raman']
                f.write("RAMAN SPECTROSCOPY\n")
                f.write("-"*18 + "\n")
                if 'raman_activities' in raman_data and raman_data['raman_activities']:
                    f.write("Raman Activities:\n")
                    frequencies = raman_data.get('frequencies', [])
                    for i, activity in enumerate(raman_data['raman_activities'][:10]):
                        freq = frequencies[i] if i < len(frequencies) else 0
                        f.write(f"  {freq:8.2f} cm⁻¹: {activity:8.2f} Å⁴/amu\n")
                f.write("\n")
            
            # VCD results
            if 'vcd' in self.results and 'error' not in self.results['vcd']:
                vcd_data = self.results['vcd']
                f.write("VIBRATIONAL CIRCULAR DICHROISM (VCD)\n")
                f.write("-"*36 + "\n")
                if 'vcd_intensities' in vcd_data and vcd_data['vcd_intensities']:
                    f.write("VCD Rotational Strengths:\n")
                    frequencies = vcd_data.get('frequencies', [])
                    for i, intensity in enumerate(vcd_data['vcd_intensities'][:10]):
                        freq = frequencies[i] if i < len(frequencies) else 0
                        f.write(f"  {freq:8.2f} cm⁻¹: {intensity:12.6f} 10⁻⁴⁴ esu²·cm²\n")
                else:
                    f.write("VCD intensities not found (may be zero for symmetric molecules)\n")
                f.write("\n")
            
            # TD-DFT results
            if 'td_dft' in self.results and 'error' not in self.results['td_dft']:
                td_data = self.results['td_dft']
                f.write("UV-VISIBLE SPECTROSCOPY (TD-DFT)\n")
                f.write("-"*30 + "\n")
                if 'excited_states' in td_data and td_data['excited_states']:
                    f.write("Excited States:\n")
                    for state in td_data['excited_states'][:5]:  # Show first 5 states
                        f.write(f"  State {state['state']}: {state['wavelength_nm']:.1f} nm "
                               f"({state['energy_ev']:.2f} eV, f={state['oscillator_strength']:.4f})\n")
                f.write("\n")
            
            # NMR results
            if 'nmr' in self.results and 'error' not in self.results['nmr']:
                nmr_data = self.results['nmr']
                f.write("NMR CHEMICAL SHIFTS\n")
                f.write("-"*20 + "\n")
                if 'chemical_shifts' in nmr_data and nmr_data['chemical_shifts']:
                    f.write("Chemical Shifts (approximate ppm from TMS):\n")
                    for shift in nmr_data['chemical_shifts'][:10]:  # Show first 10
                        # Convert shielding to chemical shift (approximate)
                        chem_shift = 200 - shift['isotropic_shielding']  # Rough conversion
                        f.write(f"  Atom {shift['atom_number']} ({shift['element']}): {chem_shift:.2f} ppm\n")
                f.write("\n")
            
            # Summary
            f.write("CALCULATION SUMMARY\n")
            f.write("-"*20 + "\n")
            for calc_type in ['optimization', 'frequency', 'raman', 'vcd', 'td_dft', 'nmr', 'polarizability']:
                if calc_type in self.results:
                    if 'error' not in self.results[calc_type]:
                        status = "✓ Completed successfully"
                    else:
                        status = f"✗ Failed: {self.results[calc_type]['error']}"
                else:
                    status = "- Not attempted"
                f.write(f"{calc_type.capitalize()}: {status}\n")
        
        # Create CSV files for easy data analysis
        self.create_csv_outputs()
        
        print(f"Analysis report saved to: {report_path}")
        print(f"Detailed JSON data saved to: {json_report_path}")
    
    def create_csv_outputs(self):
        """
        Create CSV files for different types of spectroscopic data.
        """
        if pd is None:
            print("Pandas not installed; skipping CSV export.")
            return
        
        try:
            # Combined frequencies CSV with all spectroscopic data
            freq_data_combined = {}
            
            # Get frequency data from different calculations
            for calc_type in ['frequency', 'raman', 'vcd']:
                if calc_type in self.results and 'frequencies' in self.results[calc_type]:
                    if not freq_data_combined:  # First time, set up the base
                        freq_data_combined['Frequency_cm-1'] = self.results[calc_type]['frequencies']
                    
                    # Add IR intensities
                    if 'ir_intensities' in self.results[calc_type]:
                        freq_data_combined['IR_Intensity'] = self.results[calc_type]['ir_intensities']
                    
                    # Add Raman activities
                    if 'raman_activities' in self.results[calc_type]:
                        freq_data_combined['Raman_Activity'] = self.results[calc_type]['raman_activities']
                    
                    # Add VCD intensities
                    if 'vcd_intensities' in self.results[calc_type]:
                        freq_data_combined['VCD_Intensity'] = self.results[calc_type]['vcd_intensities']
            
            # Create frequencies CSV if we have data
            if freq_data_combined:
                # Ensure all arrays have the same length
                max_len = max(len(v) for v in freq_data_combined.values())
                for key, values in freq_data_combined.items():
                    if len(values) < max_len:
                        freq_data_combined[key].extend([0.0] * (max_len - len(values)))
                
                df_freq = pd.DataFrame(freq_data_combined)
                df_freq.to_csv(self.base_dir / f"{self.molecule_name}_frequencies.csv", index=False)
                print(f"Frequencies CSV saved: {self.molecule_name}_frequencies.csv")
            
            # UV-Vis transitions CSV
            if 'td_dft' in self.results and 'excited_states' in self.results['td_dft']:
                df_uv = pd.DataFrame(self.results['td_dft']['excited_states'])
                df_uv.to_csv(self.base_dir / f"{self.molecule_name}_uv_vis.csv", index=False)
                print(f"UV-Vis CSV saved: {self.molecule_name}_uv_vis.csv")
            
            # NMR chemical shifts CSV
            if 'nmr' in self.results and 'chemical_shifts' in self.results['nmr']:
                df_nmr = pd.DataFrame(self.results['nmr']['chemical_shifts'])
                df_nmr.to_csv(self.base_dir / f"{self.molecule_name}_nmr.csv", index=False)
                print(f"NMR CSV saved: {self.molecule_name}_nmr.csv")
            
        except Exception as e:
            print(f"Error creating CSV files: {e}")

def main():
    """
    Main function demonstrating how to use the Gaussian Spectroscopy Calculator.
    """
    print("Gaussian Spectroscopy Suite - Fixed Version")
    print("="*50)
    
    # Initialize calculator
    calc = GaussianSpectroscopyCalculator("water_molecule")
    
    # Set molecular geometry (water molecule example)
    water_geometry = """O          0.000000    0.000000    0.117176
H          0.000000    0.757480   -0.468706
H          0.000000   -0.757480   -0.468706"""
    
    calc.set_molecule_geometry(water_geometry)
    
    # Option 1: Just generate input files (for testing)
    print("\nGenerating input files only (set run_calculations=True to execute)...")
    calc.run_complete_spectroscopic_analysis(run_calculations=False)
    
    # Option 2: Run complete analysis (uncomment to run actual calculations)
    # calc.run_complete_spectroscopic_analysis(run_calculations=True)
    
    print(f"\nDemo completed! Check the '{calc.base_dir}' directory for results.")
    print("\nTo use with your own molecule:")
    print("1. Replace the geometry with your molecule's coordinates")
    print("2. Adjust the method and basis set if needed")
    print("3. Make sure Gaussian is properly installed and accessible")
    print("4. Set run_calculations=True to execute the calculations")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nCalculation interrupted by user.")
    except Exception as e:
        print(f"Fatal error occurred: {e}")
        import traceback
        traceback.print_exc()