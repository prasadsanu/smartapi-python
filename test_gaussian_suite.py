#!/usr/bin/env python3
"""
Test Script for Gaussian Spectroscopy Suite
===========================================
This script tests the fixed version of the Gaussian spectroscopy calculator,
particularly focusing on Raman and VCD calculations.
"""

from gaussian_spectroscopy_suite_fixed import GaussianSpectroscopyCalculator
import os
from pathlib import Path

def test_input_generation():
    """Test that input files are generated correctly"""
    print("="*50)
    print("TESTING INPUT FILE GENERATION")
    print("="*50)
    
    # Initialize calculator
    calc = GaussianSpectroscopyCalculator("test_molecule", "test_calculations")
    
    # Set a simple molecule (methane)
    methane_geometry = """C          0.000000    0.000000    0.000000
H          0.629118    0.629118    0.629118
H         -0.629118   -0.629118    0.629118
H         -0.629118    0.629118   -0.629118
H          0.629118   -0.629118   -0.629118"""
    
    calc.set_molecule_geometry(methane_geometry)
    
    # Test individual input file generation
    test_calculations = ['optimization', 'frequency', 'raman', 'vcd', 'td_dft', 'nmr']
    
    for calc_type in test_calculations:
        try:
            input_file = calc.generate_input_file(calc_type, f'test_{calc_type}.gjf')
            print(f"✓ {calc_type}: Generated successfully")
            
            # Check if file exists and has content
            if input_file.exists():
                with open(input_file, 'r') as f:
                    content = f.read()
                    if len(content) > 100:  # Basic content check
                        print(f"  File size: {len(content)} characters")
                    else:
                        print(f"  Warning: File seems too small")
            
        except Exception as e:
            print(f"✗ {calc_type}: Failed - {e}")
    
    return calc

def test_keyword_formats():
    """Test that Gaussian keywords are formatted correctly"""
    print("\n" + "="*50)
    print("TESTING KEYWORD FORMATS")
    print("="*50)
    
    calc = GaussianSpectroscopyCalculator("keyword_test", "keyword_test")
    calc.set_molecule_geometry("C 0 0 0\nH 1 0 0")
    
    # Test specific problematic keywords
    test_cases = {
        'raman': 'freq=(raman)',
        'vcd': 'freq=(vcd)',
        'roa': 'freq=(roa)',
        'ecd': 'td=(ecd,nstates=10)',
        'anharmonic': 'freq=(anharmonic)'
    }
    
    for calc_type, expected_keyword in test_cases.items():
        input_file = calc.generate_input_file(calc_type, f'keyword_{calc_type}.gjf')
        
        with open(input_file, 'r') as f:
            content = f.read()
        
        if expected_keyword in content:
            print(f"✓ {calc_type}: Keyword '{expected_keyword}' found correctly")
        else:
            print(f"✗ {calc_type}: Expected keyword '{expected_keyword}' not found")
            print(f"  Content preview: {content[:200]}...")

def test_with_your_molecule():
    """Template for testing with your own molecule"""
    print("\n" + "="*50)
    print("TESTING WITH CUSTOM MOLECULE")
    print("="*50)
    
    # Initialize calculator with your molecule name
    calc = GaussianSpectroscopyCalculator("your_molecule", "your_molecule_test")
    
    # Replace this with your actual molecular geometry
    your_geometry = """C          0.000000    0.000000    0.000000
H          0.000000    0.000000    1.089000
H          1.026719    0.000000   -0.363000
H         -0.513360    0.889165   -0.363000
H         -0.513360   -0.889165   -0.363000"""
    
    calc.set_molecule_geometry(your_geometry)
    
    # Test generating all input files
    print("Generating all spectroscopic input files...")
    calc.run_complete_spectroscopic_analysis(run_calculations=False)
    
    print("Input files generated successfully!")
    print(f"Check directory: {calc.base_dir}")
    
    return calc

def show_generated_files(calc):
    """Show what files were generated"""
    print("\n" + "="*50)
    print("GENERATED FILES")
    print("="*50)
    
    if calc.base_dir.exists():
        files = list(calc.base_dir.glob("*.gjf"))
        if files:
            print("Input files (.gjf):")
            for file in sorted(files):
                size = file.stat().st_size
                print(f"  {file.name:20} ({size} bytes)")
        
        output_files = list(calc.base_dir.glob("*.out"))
        if output_files:
            print("\nOutput files (.out):")
            for file in sorted(output_files):
                size = file.stat().st_size
                print(f"  {file.name:20} ({size} bytes)")
        
        other_files = list(calc.base_dir.glob("*.txt")) + list(calc.base_dir.glob("*.json")) + list(calc.base_dir.glob("*.csv"))
        if other_files:
            print("\nAnalysis files:")
            for file in sorted(other_files):
                size = file.stat().st_size
                print(f"  {file.name:20} ({size} bytes)")
    else:
        print("No files found - directory doesn't exist")

def show_sample_input_file(calc):
    """Show content of a sample input file"""
    print("\n" + "="*50)
    print("SAMPLE INPUT FILE CONTENT")
    print("="*50)
    
    # Look for a Raman input file to show
    raman_file = calc.base_dir / "raman.gjf"
    if raman_file.exists():
        print(f"Content of {raman_file}:")
        print("-" * 30)
        with open(raman_file, 'r') as f:
            print(f.read())
        print("-" * 30)
    else:
        # Show any .gjf file
        gjf_files = list(calc.base_dir.glob("*.gjf"))
        if gjf_files:
            sample_file = gjf_files[0]
            print(f"Content of {sample_file}:")
            print("-" * 30)
            with open(sample_file, 'r') as f:
                print(f.read())
            print("-" * 30)

def main():
    """Run all tests"""
    print("Gaussian Spectroscopy Suite - Test Script")
    print("="*60)
    
    # Test 1: Basic input generation
    calc1 = test_input_generation()
    
    # Test 2: Keyword formats
    test_keyword_formats()
    
    # Test 3: With custom molecule
    calc2 = test_with_your_molecule()
    
    # Show what was generated
    show_generated_files(calc2)
    show_sample_input_file(calc2)
    
    print("\n" + "="*60)
    print("TESTING COMPLETED")
    print("="*60)
    print("\nNext steps:")
    print("1. Check the generated .gjf files look correct")
    print("2. If you have Gaussian installed, set run_calculations=True")
    print("3. Run the actual calculations to test parsing")
    print("4. Check the analysis reports and CSV files")

if __name__ == "__main__":
    main()