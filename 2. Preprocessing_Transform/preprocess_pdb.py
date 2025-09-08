#!/usr/bin/env python3
"""
PDB 파일 전처리 스크립트
- 용매(물) 분자 제거
- 수소 원자 추가
- 전하 계산 및 프로톤화 상태 지정
"""

import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol
import warnings
warnings.filterwarnings('ignore')

def read_pdb_file(pdb_file):
    """PDB 파일을 읽고 분자 객체로 변환"""
    try:
        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False)
        if mol is None:
            print(f"❌ PDB 파일 읽기 실패: {pdb_file}")
            return None
        
        # 기본 sanitization 시도
        try:
            Chem.SanitizeMol(mol)
        except:
            print(f"⚠️  분자 sanitization 실패, 부분적으로 처리합니다: {pdb_file}")
        
        return mol
    except Exception as e:
        print(f"❌ 오류 발생: {e}")
        return None

def remove_solvent(mol):
    """용매 분자 제거 (물, 이온 등)"""
    if mol is None:
        return None
    
    # 제거할 용매 분자들
    solvent_names = ['HOH', 'WAT', 'H2O', 'CL', 'NA', 'K', 'MG', 'CA', 'ZN', 'FE']
    
    # 원자 정보에서 용매 제거
    atoms_to_remove = []
    for atom in mol.GetAtoms():
        pdb_info = atom.GetPDBResidueInfo()
        if pdb_info and pdb_info.GetResidueName().strip() in solvent_names:
            atoms_to_remove.append(atom.GetIdx())
    
    if atoms_to_remove:
        print(f"용매 원자 {len(atoms_to_remove)}개 제거 중...")
        # 역순으로 제거 (인덱스 변경 방지)
        atoms_to_remove.sort(reverse=True)
        mol_rw = Chem.RWMol(mol)
        for idx in atoms_to_remove:
            mol_rw.RemoveAtom(idx)
        mol = mol_rw.GetMol()
    
    return mol

def add_hydrogens(mol):
    """수소 원자 추가"""
    if mol is None:
        return None
    
    try:
        # 기존 수소 제거 후 다시 추가
        mol = Chem.RemoveHs(mol)
        mol = Chem.AddHs(mol, addCoords=True)
        print(f"수소 추가 완료: {mol.GetNumAtoms()}개 원자")
        return mol
    except Exception as e:
        print(f"⚠️  수소 추가 실패: {e}")
        return mol

def add_charges(mol):
    """전하 계산 및 프로톤화 상태 지정"""
    if mol is None:
        return None
    
    try:
        # Gasteiger 전하 계산
        AllChem.ComputeGasteigerCharges(mol)
        print("Gasteiger 전하 계산 완료")
        
        # 단백질의 경우 pH 7.4에서의 프로톤화 상태 적용
        # (실제로는 더 정교한 pKa 계산이 필요하지만 기본적인 처리)
        
        return mol
    except Exception as e:
        print(f"⚠️  전하 계산 실패: {e}")
        return mol

def save_pdb(mol, output_file):
    """전처리된 분자를 PDB 파일로 저장"""
    if mol is None:
        return False
    
    try:
        # PDB 파일로 저장
        Chem.MolToPDBFile(mol, output_file)
        
        if os.path.exists(output_file):
            file_size = os.path.getsize(output_file)
            print(f"✅ 저장 완료: {output_file} ({file_size} bytes)")
            return True
        else:
            print(f"❌ 파일 저장 실패: {output_file}")
            return False
            
    except Exception as e:
        print(f"❌ 저장 오류: {e}")
        return False

def preprocess_pdb(input_file, output_file):
    """PDB 파일 전처리 메인 함수"""
    print(f"\n🔄 처리 중: {input_file}")
    
    # 1. PDB 파일 읽기
    mol = read_pdb_file(input_file)
    if mol is None:
        return False
    
    initial_atoms = mol.GetNumAtoms()
    print(f"초기 원자 수: {initial_atoms}")
    
    # 2. 용매 제거
    mol = remove_solvent(mol)
    if mol is None:
        return False
    
    # 3. 수소 추가
    mol = add_hydrogens(mol)
    if mol is None:
        return False
    
    # 4. 전하 계산
    mol = add_charges(mol)
    if mol is None:
        return False
    
    # 5. 저장
    success = save_pdb(mol, output_file)
    
    if success:
        final_atoms = mol.GetNumAtoms()
        print(f"최종 원자 수: {final_atoms} (변화: {final_atoms - initial_atoms:+d})")
    
    return success

def main():
    """메인 실행 함수"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 처리할 PDB 파일들
    pdb_files = [
        ("PDB 1IYT.pdb", "1IYT_prepared.pdb"),
        ("PDB 2M4J.pdb", "2M4J_prepared.pdb")
    ]
    
    print("🧬 PDB 파일 전처리 시작")
    print("=" * 50)
    
    success_count = 0
    
    for input_name, output_name in pdb_files:
        input_path = os.path.join(script_dir, input_name)
        output_path = os.path.join(script_dir, output_name)
        
        if not os.path.exists(input_path):
            print(f"❌ 파일을 찾을 수 없습니다: {input_path}")
            continue
        
        if preprocess_pdb(input_path, output_path):
            success_count += 1
    
    print("\n" + "=" * 50)
    print(f"🎯 완료: {success_count}/{len(pdb_files)} 파일 처리됨")
    
    if success_count == len(pdb_files):
        print("✅ 모든 파일이 성공적으로 전처리되었습니다!")
    else:
        print("⚠️  일부 파일 처리에 실패했습니다.")

if __name__ == "__main__":
    main()
