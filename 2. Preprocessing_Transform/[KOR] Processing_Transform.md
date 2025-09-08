# Preprocessing & Transform

> **Language**: [English]([EN]%20Processing_Transform.md) | [한국어]

이 폴더는 PDB 파일 전처리 및 분자 도킹을 위한 구조 변환 스크립트들을 포함합니다.

## 📁 파일 구조

### 입력 파일
- `PDB 1IYT.pdb` - 1IYT 원본 PDB 파일
- `PDB 2M4J.pdb` - 2M4J 원본 PDB 파일
- `Conformer3D_COMPOUND_CID_92337.sdf` - R-flurbiprofen SDF 파일

### Python 스크립트
- `create_single_model.py` - 다중 모델 PDB에서 단일 모델 추출
- `fix_pdb_headers.py` - PDB 헤더 정보 수정
- `prepare_receptor.py` - 수용체 준비 및 PDBQT 변환
- `preprocess_pdb.py` - PDB 파일 전처리
- `transform.py` - 구조 변환 유틸리티

### Shell 스크립트
- `prepare_1IYT_receptor.sh` - 1IYT 수용체 자동 준비
- `prepare_2M4J_receptor.sh` - 2M4J 수용체 자동 준비

## 🚀 실행 방법

### 1. 환경 설정

```bash
# 필요한 Python 패키지 설치
pip install biopython numpy pandas

# Open Babel 설치 (macOS)
brew install open-babel

# AutoDock Tools 설치 필요 (MGLTools)
# https://ccsb.scripps.edu/mgltools/downloads/
```

### 2. 단일 모델 PDB 생성

```bash
# 다중 모델 PDB에서 첫 번째 모델만 추출
python create_single_model.py PDB\ 1IYT.pdb 1iyt_single.pdb
python create_single_model.py PDB\ 2M4J.pdb 2m4j_single.pdb
```

### 3. PDB 헤더 수정

```bash
# PDB 파일의 헤더 정보 정리
python fix_pdb_headers.py 1iyt_single.pdb 1iyt_fixed.pdb
python fix_pdb_headers.py 2m4j_single.pdb 2m4j_fixed.pdb
```

### 4. 수용체 준비 (자동화)

```bash
# 1IYT 수용체 준비
chmod +x prepare_1IYT_receptor.sh
./prepare_1IYT_receptor.sh

# 2M4J 수용체 준비
chmod +x prepare_2M4J_receptor.sh
./prepare_2M4J_receptor.sh
```

### 5. 수동 수용체 준비

```bash
# Python 스크립트를 직접 사용
python prepare_receptor.py --input 1iyt_fixed.pdb --output 1iyt_receptor.pdbqt
python prepare_receptor.py --input 2m4j_fixed.pdb --output 2m4j_receptor.pdbqt
```

### 6. 전체 전처리 파이프라인

```bash
# 모든 전처리 단계를 순차적으로 실행
python preprocess_pdb.py --target 1IYT
python preprocess_pdb.py --target 2M4J
```

## 📋 스크립트 설명

### `create_single_model.py`
- **목적**: NMR 구조의 다중 모델에서 첫 번째 모델만 추출
- **사용법**: `python create_single_model.py <input.pdb> <output.pdb>`
- **출력**: 단일 모델 PDB 파일

### `fix_pdb_headers.py`
- **목적**: PDB 파일의 헤더 정보 정리 및 표준화
- **사용법**: `python fix_pdb_headers.py <input.pdb> <output.pdb>`
- **기능**: HEADER, TITLE, REMARK 라인 정리

### `prepare_receptor.py`
- **목적**: 도킹을 위한 수용체 PDBQT 파일 생성
- **사용법**: `python prepare_receptor.py --input <pdb> --output <pdbqt>`
- **기능**: 수소 추가, 전하 계산, PDBQT 변환

### `preprocess_pdb.py`
- **목적**: 전체 전처리 파이프라인 실행
- **사용법**: `python preprocess_pdb.py --target <1IYT|2M4J>`
- **기능**: 모든 전처리 단계 자동화

### Shell 스크립트
- **목적**: 수용체 준비 과정 완전 자동화
- **포함 단계**: 
  1. 단일 모델 추출
  2. 헤더 수정
  3. PDBQT 변환
  4. 검증

## ⚠️ 주의사항

1. **파일 경로**: 공백이 포함된 파일명은 따옴표로 감싸기
2. **권한 설정**: Shell 스크립트 실행 전 실행 권한 부여 필요
3. **의존성**: Open Babel과 MGLTools가 시스템에 설치되어 있어야 함
4. **메모리**: 대용량 PDB 파일 처리 시 충분한 메모리 필요

## 🔧 문제 해결

### Open Babel 오류
```bash
# macOS에서 Open Babel 재설치
brew uninstall open-babel
brew install open-babel
```

### MGLTools 경로 문제
```bash
# MGLTools 경로를 환경변수에 추가
export PATH="/Applications/MGLTools-1.5.7/bin:$PATH"
```

### 권한 오류
```bash
# 스크립트 실행 권한 부여
chmod +x *.sh
chmod +x *.py
```

## 📤 출력 파일

성공적으로 실행되면 다음 파일들이 생성됩니다:
- `1iyt_single.pdb` - 단일 모델 1IYT 구조
- `2m4j_single.pdb` - 단일 모델 2M4J 구조
- `1iyt_receptor.pdbqt` - 도킹용 1IYT 수용체
- `2m4j_receptor.pdbqt` - 도킹용 2M4J 수용체
