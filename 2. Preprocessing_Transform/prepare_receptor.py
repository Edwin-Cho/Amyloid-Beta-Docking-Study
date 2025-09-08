#!/usr/bin/env python3
"""
AutoDock용 단백질 수용체 전처리 스크립트
- 용매 제거
- 수소 추가
- 전하 계산
- PDBQT 형식으로 저장
"""

import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import warnings
warnings.filterwarnings('ignore')

def clean_pdb_for_autodock(input_pdb, output_pdb):
    """AutoDock용 PDB 파일 정리"""
    
    print(f"🔄 정리 중: {input_pdb}")
    
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
    
    cleaned_lines = []
    atom_count = 0
    
    # 필요한 라인만 유지
    for line in lines:
        if line.startswith('ATOM'):
            # 용매 분자 제거
            residue_name = line[17:20].strip()
            if residue_name not in ['HOH', 'WAT', 'H2O', 'CL', 'NA', 'K', 'MG', 'CA', 'ZN', 'FE']:
                # 전하 정보 정리 (AutoDock 호환성)
                cleaned_line = line[:54] + '  1.00  0.00' + line[66:]
                cleaned_lines.append(cleaned_line)
                atom_count += 1
        elif line.startswith(('HEADER', 'TITLE', 'COMPND', 'MODEL', 'ENDMDL', 'END')):
            cleaned_lines.append(line)
    
    # 정리된 PDB 저장
    with open(output_pdb, 'w') as f:
        f.writelines(cleaned_lines)
    
    print(f"✅ 정리 완료: {atom_count}개 원자, {output_pdb}")
    return output_pdb

def create_autodock_script(pdb_file, output_name):
    """AutoDockTools 명령어 스크립트 생성"""
    
    script_content = f"""#!/bin/bash
# AutoDockTools를 사용한 수용체 전처리 스크립트

echo "🧬 AutoDockTools로 수용체 전처리 시작: {pdb_file}"

# 1. 수용체 전처리 (Python 2.7 환경에서 실행)
python2.7 -c "
import sys
sys.path.append('/usr/local/lib/python2.7/site-packages/AutoDockTools')
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation

# 수용체 전처리
prep = AD4ReceptorPreparation()
prep.prepare_receptor('{pdb_file}', outputfilename='{output_name}.pdbqt', 
                     repairs='bonds_hydrogens', 
                     charges_to_add='gasteiger',
                     cleanup='nphs_lps_waters_nonstdres')
print('✅ 수용체 전처리 완료: {output_name}.pdbqt')
"

# 대안: MGLTools의 prepare_receptor4.py 사용
if [ -f "/usr/local/bin/prepare_receptor4.py" ]; then
    echo "🔧 prepare_receptor4.py 사용"
    python2.7 /usr/local/bin/prepare_receptor4.py -r {pdb_file} -o {output_name}.pdbqt -A hydrogens -U nphs_lps_waters_nonstdres
elif [ -f "/opt/mgltools/bin/prepare_receptor4.py" ]; then
    echo "🔧 prepare_receptor4.py 사용 (opt 경로)"
    python2.7 /opt/mgltools/bin/prepare_receptor4.py -r {pdb_file} -o {output_name}.pdbqt -A hydrogens -U nphs_lps_waters_nonstdres
else
    echo "⚠️  AutoDockTools가 설치되지 않았습니다."
    echo "다음 명령어로 설치하세요:"
    echo "conda install -c bioconda autodock-vina"
    echo "또는 MGLTools를 다운로드하세요: http://mgltools.scripps.edu/"
fi

echo "🎯 완료!"
"""
    
    script_file = f"prepare_{output_name}.sh"
    with open(script_file, 'w') as f:
        f.write(script_content)
    
    # 실행 권한 부여
    os.chmod(script_file, 0o755)
    
    print(f"📝 AutoDock 스크립트 생성: {script_file}")
    return script_file

def create_manual_pdbqt(pdb_file, output_pdbqt):
    """수동으로 PDBQT 형식 생성 (AutoDockTools 없을 때)"""
    
    print(f"🔧 수동 PDBQT 변환: {pdb_file}")
    
    # 원자 타입 매핑 (AutoDock 4 기준)
    atom_type_map = {
        'C': 'C', 'N': 'N', 'O': 'O', 'S': 'S', 'P': 'P',
        'F': 'F', 'CL': 'Cl', 'BR': 'Br', 'I': 'I',
        'H': 'HD'  # 극성 수소
    }
    
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    pdbqt_lines = []
    
    for line in lines:
        if line.startswith('ATOM'):
            # PDB에서 PDBQT로 변환
            atom_name = line[12:16].strip()
            element = line[76:78].strip()
            
            if not element:
                # 원소 정보가 없으면 원자명에서 추출
                element = atom_name[0]
            
            # AutoDock 원자 타입 결정
            autodock_type = atom_type_map.get(element.upper(), 'C')
            
            # PDBQT 형식으로 변환
            pdbqt_line = line[:77] + f"{autodock_type:>2}" + line[79:]
            pdbqt_lines.append(pdbqt_line)
        elif line.startswith(('ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF')):
            pdbqt_lines.append(line)
        elif line.startswith(('HEADER', 'COMPND', 'MODEL', 'ENDMDL', 'END')):
            pdbqt_lines.append(line)
    
    # PDBQT 파일 저장
    with open(output_pdbqt, 'w') as f:
        f.writelines(pdbqt_lines)
    
    print(f"✅ PDBQT 변환 완료: {output_pdbqt}")
    return output_pdbqt

def main():
    """메인 실행 함수"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 처리할 파일들
    files = [
        ("1IYT_prepared.pdb", "1IYT_receptor"),
        ("2M4J_prepared.pdb", "2M4J_receptor")
    ]
    
    print("🎯 AutoDock용 수용체 전처리")
    print("=" * 50)
    
    for pdb_name, output_name in files:
        pdb_path = os.path.join(script_dir, pdb_name)
        
        if not os.path.exists(pdb_path):
            print(f"❌ 파일 없음: {pdb_path}")
            continue
        
        # 1. PDB 정리
        cleaned_pdb = os.path.join(script_dir, f"{output_name}_clean.pdb")
        clean_pdb_for_autodock(pdb_path, cleaned_pdb)
        
        # 2. AutoDock 스크립트 생성
        script_file = create_autodock_script(cleaned_pdb, output_name)
        
        # 3. 수동 PDBQT 생성 (백업용)
        pdbqt_file = os.path.join(script_dir, f"{output_name}.pdbqt")
        create_manual_pdbqt(cleaned_pdb, pdbqt_file)
        
        print(f"📋 사용법:")
        print(f"   ./{script_file}  # AutoDockTools 사용")
        print(f"   또는 {pdbqt_file} 직접 사용")
        print()

if __name__ == "__main__":
    main()
