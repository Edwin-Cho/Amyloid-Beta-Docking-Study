# SwissDock Molecular Docking Project
## Comparative Docking Study of R-flurbiprofen with Amyloid-Beta Structures

> **Language**: [English] | [한국어]([KOR]%20Overview.md)

This project conducts molecular docking research between R-flurbiprofen, a candidate therapeutic compound for Alzheimer's disease, and two different amyloid-beta structures (1IYT, 2M4J).

---

## 🎯 Project Overview

### Research Objectives
- Compare binding characteristics between soluble amyloid-beta (1IYT) and fibril form (2M4J)
- Evaluate therapeutic potential of R-flurbiprofen
- Propose stage-specific treatment strategies for disease progression

### Key Findings
- **2M4J Best Binding Energy**: -4.814 kcal/mol (moderate affinity)
- **Structural Differences**: Clear distinction in binding characteristics between soluble vs fibril forms
- **Therapeutic Implications**: Confirmed need for differentiated approaches by disease stage

---

## 📁 Project Structure

```
SwissDock Project/
├── 0. Results/                    # 📊 Final results and analysis reports
│   ├── [EN] Comparative_Docking_Analysis.md
│   ├── [KOR] Comparative_Docking_Analysis.md
│   ├── 1IYT_Results/             # SwissDock results
│   └── 2M4J_Results/             # AutoDock Vina results
├── 1. Materials/                  # 🧬 Original structure files
│   ├── 1iyt_single.pdb
│   ├── 2m4j_single.pdb
│   ├── R_flurbiprofen_NCI.mol2
│   ├── [KOR] Materials.md
│   └── [EN] Materials.md
├── 2. Preprocessing_Transform/    # ⚙️ Structure preprocessing scripts
│   ├── *.py (5 Python scripts)
│   ├── *.sh (2 Shell scripts)
│   ├── [KOR] Processing_Transform.md
│   └── [EN] Processing_Transform.md
├── 3. Docking_Local/             # 🎯 Local docking execution
│   ├── molecular_docking.py
│   ├── final_analysis.py
│   ├── [KOR] Docking_Local.md
│   ├── [EN] Docking_Local.md
│   └── Materials/
└── 4. Archives/                  # 📚 Reference materials
```

---

## 🚀 Quick Start Guide

### 1. Environment Setup

```bash
# Install Python packages
pip install biopython numpy pandas matplotlib

# Install AutoDock Vina (macOS)
brew install autodock-vina

# Install Open Babel
brew install open-babel

# Install MGLTools (manual)
# https://ccsb.scripps.edu/mgltools/downloads/
```

### 2. Complete Workflow Execution

```bash
# Step 1: Structure preprocessing
cd "2. Preprocessing_Transform"
./prepare_2M4J_receptor.sh

# Step 2: Molecular docking execution
cd "../3. Docking_Local"
python3 molecular_docking.py

# Step 3: Results analysis
python3 final_analysis.py

# Step 4: View results
cd "../0. Results"
open "[EN] Comparative_Docking_Analysis.md"
```

---

## 📋 Detailed Guides by Folder

### 🧬 1. Materials/
**Original structure file repository**
- 1IYT: Soluble amyloid-beta (1-42)
- 2M4J: Amyloid-beta fibril (1-40)  
- R-flurbiprofen: Candidate therapeutic ligand

📖 [Korean Guide](1.%20Materials/[KOR]%20Materials.md) | [English Guide](1.%20Materials/[EN]%20Materials.md)

### ⚙️ 2. Preprocessing_Transform/
**Structure preprocessing and transformation scripts**

```bash
# Automated preprocessing (recommended)
./prepare_2M4J_receptor.sh

# Manual step-by-step execution
python create_single_model.py "PDB 2M4J.pdb" 2m4j_single.pdb
python fix_pdb_headers.py 2m4j_single.pdb 2m4j_fixed.pdb
python prepare_receptor.py --input 2m4j_fixed.pdb --output 2m4j.pdbqt
```

📖 [Korean Guide](2.%20Preprocessing_Transform/[KOR]%20Processing_Transform.md) | [English Guide](2.%20Preprocessing_Transform/[EN]%20Processing_Transform.md)

### 🎯 3. Docking_Local/
**Molecular docking execution and analysis**

```bash
# Complete docking workflow
python3 molecular_docking.py

# Results analysis only
python3 final_analysis.py
```

**Key Results:**
- Best binding energy: -4.814 kcal/mol
- 9 binding poses generated
- Moderate-level binding affinity (10-100 μM)

📖 [Korean Guide](3.%20Docking_Local/[KOR]%20Docking_Local.md) | [English Guide](3.%20Docking_Local/[EN]%20Docking_Local.md)

### 📊 0. Results/
**Final analysis results and reports**
- English/Korean comparative analysis reports
- Individual docking result files
- Statistical analysis and biological interpretation

📖 [Detailed Guide](0.%20Results/README.md)

---

## 🔬 Scientific Contributions

### Novel Discoveries
1. **Structure-specific Binding Differences**: Clear distinction between soluble vs fibril forms
2. **Treatment Strategies**: Stage-specific personalized approaches
3. **Drug Design**: R-flurbiprofen optimization directions

### Clinical Implications
- Stage-specific treatment strategies for Alzheimer's disease progression
- Guidelines for amyloid-targeting therapeutic development
- Potential for personalized medicine approaches

---

## 📈 Key Results Summary

| Structure | Docking Tool | Best Binding Energy | Binding Characteristics | Therapeutic Approach |
|-----------|--------------|---------------------|------------------------|---------------------|
| 1IYT | SwissDock | Various modes | Flexible binding | Preventive |
| 2M4J | AutoDock Vina | -4.814 kcal/mol | Restricted binding | Therapeutic |

---

## 🛠️ Troubleshooting

### Common Errors

**Vina Execution Errors**
```bash
# Verify Vina installation
which vina
vina --version
```

**Python Package Errors**
```bash
pip install --upgrade biopython numpy pandas matplotlib
```

**File Permission Errors**
```bash
chmod +x *.sh
chmod +x *.py
```

---

## 📚 References

### Key Publications
1. **1IYT**: Crescenzi et al. (2002) *Eur J Biochem* 269: 5642-5648
2. **2M4J**: Lu et al. (2013) *Cell* 154: 1257-1268
3. **R-flurbiprofen**: Wilcock et al. (2008) *Lancet Neurol* 7: 483-493

### Tools and Databases
- [AutoDock Vina](https://vina.scripps.edu/)
- [SwissDock](http://www.swissdock.ch/)
- [Protein Data Bank](https://www.rcsb.org/)

---

**Last Updated**: September 8, 2025
