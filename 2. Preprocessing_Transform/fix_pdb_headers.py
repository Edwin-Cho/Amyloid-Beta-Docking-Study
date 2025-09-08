#!/usr/bin/env python3
"""
PDB 파일에 표준 헤더 추가하여 siwss doc 호환성 개선
"""

import os
from datetime import datetime

def add_pdb_header(input_file, output_file, pdb_id, title):
    """표준 PDB 헤더를 추가하여 파일 수정"""
    
    # 현재 날짜
    date_str = datetime.now().strftime("%d-%b-%y").upper()
    
    # 표준 PDB 헤더 생성
    header_lines = [
        f"HEADER    PROTEIN                                 {date_str}   {pdb_id.upper()}              \n",
        f"TITLE     {title}\n",
        f"COMPND    MOL_ID: 1;\n",
        f"COMPND   2 MOLECULE: {title};\n",
        f"COMPND   3 CHAIN: A;\n",
        f"SOURCE    MOL_ID: 1;\n",
        f"SOURCE   2 ORGANISM_SCIENTIFIC: HOMO SAPIENS;\n",
        f"SOURCE   3 ORGANISM_COMMON: HUMAN;\n",
        f"KEYWDS    PROTEIN, ALZHEIMER DISEASE, AMYLOID\n",
        f"EXPDTA    SOLUTION NMR\n",
        f"REMARK   2 RESOLUTION. NOT APPLICABLE.\n"
    ]
    
    # 기존 파일 읽기
    with open(input_file, 'r') as f:
        original_lines = f.readlines()
    
    # 새 파일 작성
    with open(output_file, 'w') as f:
        # 헤더 추가
        f.writelines(header_lines)
        
        # 기존 내용 추가 (COMPND 라인 제외)
        for line in original_lines:
            if not line.startswith('COMPND'):
                f.write(line)
        
        # END 라인 추가 (없으면)
        if not any(line.startswith('END') for line in original_lines):
            f.write("END\n")
    
    print(f"✅ 헤더 추가 완료: {output_file}")

def main():
    """메인 실행 함수"""
    
    files_to_fix = [
        ("1iyt.pdb", "1iyt", "ALZHEIMER'S DISEASE AMYLOID BETA-PEPTIDE"),
        ("2m4j.pdb", "2m4j", "BETA-AMYLOID FIBRIL FROM ALZHEIMER'S DISEASE")
    ]
    
    print("🔧 PDB 헤더 수정 중...")
    
    for filename, pdb_id, title in files_to_fix:
        if os.path.exists(filename):
            # 백업 생성
            backup_name = f"{filename}.backup"
            os.rename(filename, backup_name)
            
            # 헤더 추가하여 새 파일 생성
            add_pdb_header(backup_name, filename, pdb_id, title)
            
            print(f"📁 백업: {backup_name}")
        else:
            print(f"❌ 파일 없음: {filename}")
    
    print("\n🎯 siwss doc에서 다시 시도해보세요!")
    print("PDB ID 입력: 1iyt 또는 2m4j")

if __name__ == "__main__":
    main()
