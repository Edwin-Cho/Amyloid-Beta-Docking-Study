# Materials - Original Structure Files

> **Language**: [English] | [한국어]([KOR]%20Materials.md)

This folder contains the original structure files used for molecular docking research.

## 📁 File Structure

### Protein Structures (PDB Files)
- `1iyt_single.pdb` - 1IYT single model structure (soluble amyloid-beta peptide)
- `2m4j_single.pdb` - 2M4J single model structure (amyloid-beta fibril)

### Ligand Structure (MOL2 File)
- `R_flurbiprofen_NCI.mol2` - R-flurbiprofen 3D structure (NCI database)

## 📋 File Descriptions

### 1IYT - Soluble Amyloid-Beta (1-42)
- **PDB ID**: 1IYT
- **Classification**: Protein binding
- **Experimental Method**: Solution NMR
- **Structure Type**: Monomeric peptide in apolar microenvironment
- **Residues**: 42 amino acids (Aβ1-42)
- **Conformation**: Alpha-helix rich, membrane-associated form
- **Biological Context**: Early-stage soluble amyloid species
- **Deposition Date**: 2002-09-06

### 2M4J - Amyloid-Beta Fibril (1-40)
- **PDB ID**: 2M4J
- **Classification**: Protein fibril
- **Experimental Method**: Solid-state NMR
- **Structure Type**: Fibrillar aggregate from patient brain tissue
- **Residues**: 40 amino acids (Aβ1-40)
- **Conformation**: Beta-sheet rich fibril structure
- **Biological Context**: Late-stage pathological amyloid deposits
- **Deposition Date**: 2013-02-05

### R-flurbiprofen
- **Chemical Formula**: C₁₅H₁₃FO₂
- **Molecular Weight**: 244.26 g/mol
- **IUPAC Name**: (2R)-2-(3-fluoro-4-phenylphenyl)propanoic acid
- **Drug Classification**: Nonsteroidal anti-inflammatory drug (NSAID)
- **Therapeutic Interest**: Alzheimer's disease treatment candidate
- **Mode of Action**: Amyloid aggregation inhibition

## 🔍 Research Purpose

These files are used for the following research objectives:

1. **Structural Comparison**: Soluble vs fibril forms of amyloid-beta
2. **Docking Studies**: Analysis of R-flurbiprofen binding characteristics
3. **Drug Development**: Alzheimer's disease therapeutic strategy research
4. **Structure-Activity Relationships**: Understanding molecular interactions

## 🚀 Next Steps

These original files are processed in the following folders:

1. **2. Preprocessing_Transform/**: Structure preprocessing and transformation
2. **3. Docking_Local/**: Molecular docking execution
3. **0. Results/**: Results analysis and reports

## 📊 Structural Features

### Key Differences Between 1IYT vs 2M4J

| Property | 1IYT (Soluble) | 2M4J (Fibril) |
|----------|----------------|---------------|
| Structural State | Monomer | Aggregate |
| Main Secondary Structure | α-helix | β-sheet |
| Flexibility | High | Low |
| Membrane Binding | Present | Absent |
| Pathological Stage | Early | Late |
| Therapeutic Target | Preventive | Therapeutic |

## ⚠️ File Processing Considerations

1. **File Integrity**: Use copies of original files, do not modify originals
2. **Coordinate System**: Verify consistency of PDB coordinate systems
3. **Missing Atoms**: Structural completeness validation required
4. **Hydrogen Atoms**: Most experimental structures lack hydrogens (added during preprocessing)

## 📖 References

### 1IYT Related
- Crescenzi, O. et al. (2002) "Solution structure of the Alzheimer amyloid β-peptide (1-42) in an apolar microenvironment" *Eur J Biochem* 269: 5642-5648

### 2M4J Related
- Lu, J.X. et al. (2013) "Molecular structure of β-amyloid fibrils in Alzheimer's disease brain tissue" *Cell* 154: 1257-1268

### R-flurbiprofen Related
- Wilcock, G.K. et al. (2008) "Efficacy and safety of tarenflurbil in mild to moderate Alzheimer's disease" *Lancet Neurol* 7: 483-493
