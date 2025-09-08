# Preprocessing & Transform

> **Language**: [English] | [한국어]([KOR]%20Processing_Transform.md)

This folder contains PDB file preprocessing and structural transformation scripts for molecular docking.

## 📁 File Structure

### Input Files
- `PDB 1IYT.pdb` - Original 1IYT PDB file
- `PDB 2M4J.pdb` - Original 2M4J PDB file
- `Conformer3D_COMPOUND_CID_92337.sdf` - R-flurbiprofen SDF file

### Python Scripts
- `create_single_model.py` - Extract single model from multi-model PDB
- `fix_pdb_headers.py` - Fix PDB header information
- `prepare_receptor.py` - Receptor preparation and PDBQT conversion
- `preprocess_pdb.py` - PDB file preprocessing
- `transform.py` - Structure transformation utilities

### Shell Scripts
- `prepare_1IYT_receptor.sh` - Automated 1IYT receptor preparation
- `prepare_2M4J_receptor.sh` - Automated 2M4J receptor preparation

## 🚀 Execution Methods

### 1. Environment Setup

```bash
# Install required Python packages
pip install biopython numpy pandas

# Install Open Babel (macOS)
brew install open-babel

# AutoDock Tools installation required (MGLTools)
# https://ccsb.scripps.edu/mgltools/downloads/
```

### 2. Single Model PDB Generation

```bash
# Extract only the first model from multi-model PDB
python create_single_model.py PDB\ 1IYT.pdb 1iyt_single.pdb
python create_single_model.py PDB\ 2M4J.pdb 2m4j_single.pdb
```

### 3. PDB Header Modification

```bash
# Clean up PDB file header information
python fix_pdb_headers.py 1iyt_single.pdb 1iyt_fixed.pdb
python fix_pdb_headers.py 2m4j_single.pdb 2m4j_fixed.pdb
```

### 4. Receptor Preparation (Automated)

```bash
# 1IYT receptor preparation
chmod +x prepare_1IYT_receptor.sh
./prepare_1IYT_receptor.sh

# 2M4J receptor preparation
chmod +x prepare_2M4J_receptor.sh
./prepare_2M4J_receptor.sh
```

### 5. Manual Receptor Preparation

```bash
# Use Python scripts directly
python prepare_receptor.py --input 1iyt_fixed.pdb --output 1iyt_receptor.pdbqt
python prepare_receptor.py --input 2m4j_fixed.pdb --output 2m4j_receptor.pdbqt
```

### 6. Complete Preprocessing Pipeline

```bash
# Execute all preprocessing steps sequentially
python preprocess_pdb.py --target 1IYT
python preprocess_pdb.py --target 2M4J
```

## 📋 Script Descriptions

### `create_single_model.py`
- **Purpose**: Extract only the first model from NMR multi-model structures
- **Usage**: `python create_single_model.py <input.pdb> <output.pdb>`
- **Output**: Single model PDB file

### `fix_pdb_headers.py`
- **Purpose**: Clean up and standardize PDB file header information
- **Usage**: `python fix_pdb_headers.py <input.pdb> <output.pdb>`
- **Function**: Clean HEADER, TITLE, REMARK lines

### `prepare_receptor.py`
- **Purpose**: Generate receptor PDBQT files for docking
- **Usage**: `python prepare_receptor.py --input <pdb> --output <pdbqt>`
- **Function**: Add hydrogens, calculate charges, PDBQT conversion

### `preprocess_pdb.py`
- **Purpose**: Execute complete preprocessing pipeline
- **Usage**: `python preprocess_pdb.py --target <1IYT|2M4J>`
- **Function**: Automate all preprocessing steps

### Shell Scripts
- **Purpose**: Complete automation of receptor preparation process
- **Included Steps**: 
  1. Single model extraction
  2. Header modification
  3. PDBQT conversion
  4. Validation

## ⚠️ Important Notes

1. **File Paths**: Enclose filenames with spaces in quotes
2. **Permission Settings**: Grant execution permissions before running shell scripts
3. **Dependencies**: Open Babel and MGLTools must be installed on the system
4. **Memory**: Sufficient memory required for processing large PDB files

## 🔧 Troubleshooting

### Open Babel Errors
```bash
# Reinstall Open Babel on macOS
brew uninstall open-babel
brew install open-babel
```

### MGLTools Path Issues
```bash
# Add MGLTools path to environment variables
export PATH="/Applications/MGLTools-1.5.7/bin:$PATH"
```

### Permission Errors
```bash
# Grant execution permissions to scripts
chmod +x *.sh
chmod +x *.py
```

## 📤 Output Files

Upon successful execution, the following files will be generated:
- `1iyt_single.pdb` - Single model 1IYT structure
- `2m4j_single.pdb` - Single model 2M4J structure
- `1iyt_receptor.pdbqt` - 1IYT receptor for docking
- `2m4j_receptor.pdbqt` - 2M4J receptor for docking
