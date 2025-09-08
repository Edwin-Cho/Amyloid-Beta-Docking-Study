#!/usr/bin/env python3
"""
멀티 모델 PDB를 단일 모델로 변환
"""

import os

def extract_first_model(input_file, output_file):
    """첫 번째 모델만 추출하여 새 PDB 파일 생성"""
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    output_lines = []
    in_first_model = False
    model_count = 0
    
    for line in lines:
        # 헤더 정보는 항상 포함
        if line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS', 'EXPDTA', 'REMARK')):
            output_lines.append(line)
        
        # 첫 번째 MODEL 시작
        elif line.startswith('MODEL'):
            model_count += 1
            if model_count == 1:
                in_first_model = True
                # MODEL 라인 제거 (단일 모델이므로 불필요)
                continue
        
        # ENDMDL 만나면 첫 번째 모델 종료
        elif line.startswith('ENDMDL'):
            if model_count == 1:
                in_first_model = False
                break
        
        # 첫 번째 모델의 ATOM 라인들
        elif in_first_model and line.startswith('ATOM'):
            output_lines.append(line)
    
    # END 라인 추가
    output_lines.append('END\n')
    
    # 새 파일 저장
    with open(output_file, 'w') as f:
        f.writelines(output_lines)
    
    atom_count = sum(1 for line in output_lines if line.startswith('ATOM'))
    print(f"✅ 단일 모델 생성: {output_file} ({atom_count}개 원자)")

def main():
    """메인 실행 함수"""
    
    files_to_convert = [
        ("1iyt.pdb", "1iyt_single.pdb"),
        ("2m4j.pdb", "2m4j_single.pdb")
    ]
    
    print("🔧 단일 모델 PDB 파일 생성 중...")
    
    for input_file, output_file in files_to_convert:
        if os.path.exists(input_file):
            extract_first_model(input_file, output_file)
        else:
            print(f"❌ 파일 없음: {input_file}")
    
    print("\n🎯 단일 모델 파일들:")
    print("- 1iyt_single.pdb")
    print("- 2m4j_single.pdb")
    print("\n이 파일들을 siwss doc에 업로드하세요!")

if __name__ == "__main__":
    main()
