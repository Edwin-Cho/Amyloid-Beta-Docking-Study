# Docking Local - Molecular Docking Execution

> **Language**: [English] | [한국어]([KOR]%20Docking_Local.md)

This folder contains scripts for executing molecular docking in a local environment and analyzing the results.

## 📁 File Structure

### Input Files (Materials/)
- `2m4j.pdbqt` - Prepared 2M4J receptor structure
- `R_flurbiprofen.pdbqt` - Prepared R-flurbiprofen ligand structure

### Execution Scripts
- `molecular_docking.py` - Complete docking workflow script
- `final_analysis.py` - Comprehensive results analysis script

### Result Files
- `final_docking_result.pdbqt` - Final docking results (9 binding poses)

## 🚀 Execution Methods

### 1. Environment Setup

```bash
# Verify AutoDock Vina installation
vina --help

# Install Python packages
pip install numpy matplotlib biopython

# Grant execution permissions
chmod +x *.py
```

### 2. Complete Docking Workflow Execution

```bash
# Execute entire docking process (receptor preparation + docking + analysis)
python3 molecular_docking.py
```

**This script performs the following:**
- Receptor file validation and preparation
- Ligand file validation and preparation
- Binding site center calculation
- AutoDock Vina docking execution
- Results analysis and visualization

### 3. Results Analysis Only

```bash
# Comprehensive analysis of existing docking results
python3 final_analysis.py
```

**Analysis Content:**
- Binding energy and RMSD analysis
- Binding affinity evaluation
- Biological significance interpretation
- Future research directions

## 📊 Key Results

### 🎯 Docking Performance
- **Best Binding Energy**: -4.814 kcal/mol
- **Binding Affinity**: Moderate (10-100 μM range)
- **Number of Binding Poses**: 9 diverse conformations
- **Target**: Alzheimer's disease beta-amyloid fibril
- **Ligand**: R-flurbiprofen (anti-inflammatory drug)

### 📈 Statistical Summary
- **Average Binding Energy**: -4.676 ± 0.099 kcal/mol
- **Energy Range**: 0.282 kcal/mol
- **Consistency**: All poses show moderate-level affinity

## 🔧 Detailed Script Descriptions

### `molecular_docking.py`
**Function**: Complete molecular docking pipeline
```bash
python3 molecular_docking.py
```
**Execution Steps**:
1. Input file validation
2. Receptor/ligand preparation
3. Binding site calculation
4. Vina configuration file generation
5. Docking execution
6. Results parsing and analysis

### `final_analysis.py`
**Function**: Comprehensive docking results analysis
```bash
python3 final_analysis.py
```
**Analysis Items**:
- Energy analysis by binding pose
- Statistical summary
- Biological interpretation
- Therapeutic implications
- Research recommendations

## ⚙️ Configuration Options

### Vina Parameters (modifiable in molecular_docking.py)
```python
# Docking parameters
exhaustiveness = 32        # Search precision
num_modes = 20            # Number of poses to generate
energy_range = 4          # Energy range (kcal/mol)

# Search space
center_x, center_y, center_z = auto_calculated
size_x = size_y = size_z = 25  # Search box size (Å)
```

## 📋 Output File Descriptions

### `final_docking_result.pdbqt`
- **Format**: AutoDock Vina PDBQT
- **Content**: 9 binding poses with energy information
- **Structure**: MODEL 1-9, each including binding energy and RMSD

### Console Output
- Real-time progress updates
- Binding energy rankings
- Statistical analysis results
- Biological interpretation

## 🔍 Results Interpretation

### Binding Affinity Classification
- **Strong Binding**: < -7.0 kcal/mol
- **Good Binding**: -5.0 ~ -7.0 kcal/mol
- **Moderate Binding**: -3.0 ~ -5.0 kcal/mol ← **Current Results**
- **Weak Binding**: > -3.0 kcal/mol

### Biological Significance
The moderate binding between R-flurbiprofen and beta-amyloid fibrils suggests meaningful interactions for Alzheimer's disease treatment, showing potential for fibril formation inhibition or stability disruption.

## ⚠️ Troubleshooting

### Vina Execution Errors
```bash
# Verify Vina installation
which vina
vina --version

# Path issues
export PATH="/usr/local/bin:$PATH"
```

### Memory Shortage
```bash
# Reduce exhaustiveness value
exhaustiveness = 16  # Decrease from default 32
```

### File Permission Errors
```bash
chmod 755 *.py
chmod 644 *.pdbqt
```

## 📈 Performance Optimization

### Quick Testing
```python
# Modify in molecular_docking.py
exhaustiveness = 8
num_modes = 5
```

### High-Precision Docking
```python
exhaustiveness = 64
num_modes = 50
energy_range = 6
```

## 🎯 Next Steps

1. **Visualization**: Examine binding poses with PyMOL/ChimeraX
2. **Validation**: Molecular dynamics simulations
3. **Optimization**: Ligand structure improvement
4. **Experiments**: Biochemical binding analysis
