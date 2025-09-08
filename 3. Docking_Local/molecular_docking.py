#!/usr/bin/env python3
"""
Molecular Docking Workflow: 2M4J (Beta-amyloid fibril) + R-flurbiprofen
This script performs molecular docking using AutoDock Vina
"""

import os
import subprocess
import sys
from pathlib import Path

def run_command(cmd, description=""):
    """Run shell command and handle errors"""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Command: {cmd}")
    print('='*60)
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"ERROR: {description} failed")
        print(f"STDERR: {result.stderr}")
        return False
    else:
        print(f"SUCCESS: {description} completed")
        if result.stdout:
            print(f"Output: {result.stdout}")
        return True

def prepare_receptor(pdb_file, output_pdbqt):
    """Convert PDB receptor to PDBQT format using Open Babel"""
    cmd = f"obabel {pdb_file} -O {output_pdbqt} -p 7.4"
    return run_command(cmd, f"Converting receptor {pdb_file} to PDBQT")

def prepare_ligand(mol2_file, output_pdbqt):
    """Convert MOL2 ligand to PDBQT format using Open Babel"""
    cmd = f"obabel {mol2_file} -O {output_pdbqt} --gen3d"
    return run_command(cmd, f"Converting ligand {mol2_file} to PDBQT")

def calculate_binding_site_center(pdb_file):
    """Calculate the center of mass of the protein for binding site"""
    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
        
        x_coords, y_coords, z_coords = [], [], []
        
        for line in lines:
            if line.startswith('ATOM'):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                x_coords.append(x)
                y_coords.append(y)
                z_coords.append(z)
        
        center_x = sum(x_coords) / len(x_coords)
        center_y = sum(y_coords) / len(y_coords)
        center_z = sum(z_coords) / len(z_coords)
        
        return center_x, center_y, center_z
    
    except Exception as e:
        print(f"Error calculating binding site center: {e}")
        return 0, 0, 0

def create_vina_config(receptor_pdbqt, ligand_pdbqt, center_x, center_y, center_z, 
                      size_x=20, size_y=20, size_z=20, exhaustiveness=8):
    """Create AutoDock Vina configuration file"""
    
    config_content = f"""# AutoDock Vina configuration file
# Receptor and ligand files
receptor = {receptor_pdbqt}
ligand = {ligand_pdbqt}

# Output file
out = docking_result.pdbqt

# Search space (binding site)
center_x = {center_x:.3f}
center_y = {center_y:.3f}
center_z = {center_z:.3f}

# Search space size (Angstroms)
size_x = {size_x}
size_y = {size_y}
size_z = {size_z}

# Docking parameters
exhaustiveness = {exhaustiveness}
num_modes = 9
energy_range = 3

# CPU usage
cpu = 4
"""
    
    with open('vina_config.txt', 'w') as f:
        f.write(config_content)
    
    print(f"Created Vina configuration file: vina_config.txt")
    print(f"Binding site center: ({center_x:.3f}, {center_y:.3f}, {center_z:.3f})")
    return 'vina_config.txt'

def run_vina_docking(config_file):
    """Run AutoDock Vina docking"""
    cmd = f"vina --config {config_file} > docking_log.txt 2>&1"
    return run_command(cmd, "AutoDock Vina molecular docking")

def analyze_results():
    """Analyze docking results"""
    print(f"\n{'='*60}")
    print("DOCKING RESULTS ANALYSIS")
    print('='*60)
    
    # Check if result files exist
    if os.path.exists('docking_result.pdbqt'):
        print("✓ Docking result file created: docking_result.pdbqt")
        
        # Parse the results
        try:
            with open('docking_log.txt', 'r') as f:
                log_content = f.read()
                print("\nDocking Log Summary:")
                print("-" * 40)
                
                # Extract binding energies
                lines = log_content.split('\n')
                for line in lines:
                    if 'REMARK VINA RESULT:' in line:
                        parts = line.split()
                        if len(parts) >= 4:
                            energy = parts[3]
                            print(f"Binding Energy: {energy} kcal/mol")
                            break
                
                print(f"\nFull log saved to: docking_log.txt")
                
        except Exception as e:
            print(f"Error reading log file: {e}")
    else:
        print("✗ Docking result file not found")
    
    if os.path.exists('docking_log.txt'):
        print("✓ Docking log file created: docking_log.txt")
    else:
        print("✗ Docking log file not found")

def main():
    """Main docking workflow"""
    print("="*80)
    print("MOLECULAR DOCKING WORKFLOW")
    print("Receptor: 2M4J (Beta-amyloid fibril from Alzheimer's disease)")
    print("Ligand: R-flurbiprofen")
    print("="*80)
    
    # File paths
    receptor_pdb = "../1. Process/2m4j_single.pdb"
    ligand_mol2 = "../1. Process/R_flurbiprofen_NCI.mol2"
    
    receptor_pdbqt = "2m4j_receptor.pdbqt"
    ligand_pdbqt = "R_flurbiprofen.pdbqt"
    
    # Check if input files exist
    if not os.path.exists(receptor_pdb):
        print(f"ERROR: Receptor file not found: {receptor_pdb}")
        return False
    
    if not os.path.exists(ligand_mol2):
        print(f"ERROR: Ligand file not found: {ligand_mol2}")
        return False
    
    # Step 1: Install Open Babel if needed
    print("\nChecking for Open Babel...")
    result = subprocess.run("which obabel", shell=True, capture_output=True)
    if result.returncode != 0:
        print("Installing Open Babel...")
        if not run_command("brew install open-babel", "Installing Open Babel"):
            print("Failed to install Open Babel. Please install manually.")
            return False
    
    # Step 2: Prepare receptor (PDB to PDBQT)
    if not prepare_receptor(receptor_pdb, receptor_pdbqt):
        return False
    
    # Step 3: Prepare ligand (MOL2 to PDBQT)
    if not prepare_ligand(ligand_mol2, ligand_pdbqt):
        return False
    
    # Step 4: Calculate binding site center
    center_x, center_y, center_z = calculate_binding_site_center(receptor_pdb)
    
    # Step 5: Create Vina configuration
    config_file = create_vina_config(receptor_pdbqt, ligand_pdbqt, 
                                   center_x, center_y, center_z)
    
    # Step 6: Run docking
    if not run_vina_docking(config_file):
        return False
    
    # Step 7: Analyze results
    analyze_results()
    
    print(f"\n{'='*80}")
    print("DOCKING WORKFLOW COMPLETED!")
    print("="*80)
    print("Files generated:")
    print(f"  - {receptor_pdbqt} (receptor in PDBQT format)")
    print(f"  - {ligand_pdbqt} (ligand in PDBQT format)")
    print(f"  - {config_file} (Vina configuration)")
    print(f"  - docking_result.pdbqt (docking poses)")
    print(f"  - docking_log.txt (docking log)")
    print("="*80)
    
    return True

if __name__ == "__main__":
    success = main()
    if not success:
        sys.exit(1)
