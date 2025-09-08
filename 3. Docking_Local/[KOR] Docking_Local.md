# Docking Local - Molecular Docking Execution

> **Language**: [English]([EN]%20Docking_Local.md) | [한국어]

이 폴더는 로컬 환경에서 분자 도킹을 실행하고 결과를 분석하는 스크립트들을 포함합니다.

## 📁 파일 구조

### 입력 파일 (Materials/)
- `2m4j.pdbqt` - 준비된 2M4J 수용체 구조
- `R_flurbiprofen.pdbqt` - 준비된 R-flurbiprofen 리간드 구조

### 실행 스크립트
- `molecular_docking.py` - 완전한 도킹 워크플로우 스크립트
- `final_analysis.py` - 종합적인 결과 분석 스크립트

### 결과 파일
- `final_docking_result.pdbqt` - 최종 도킹 결과 (9개 결합 포즈)

## 🚀 실행 방법

### 1. 환경 설정

```bash
# AutoDock Vina 설치 확인
vina --help

# Python 패키지 설치
pip install numpy matplotlib biopython

# 실행 권한 부여
chmod +x *.py
```

### 2. 완전한 도킹 워크플로우 실행

```bash
# 전체 도킹 프로세스 실행 (수용체 준비 + 도킹 + 분석)
python3 molecular_docking.py
```

**이 스크립트는 다음을 수행합니다:**
- 수용체 파일 검증 및 준비
- 리간드 파일 검증 및 준비
- 결합 부위 중심 계산
- AutoDock Vina 도킹 실행
- 결과 분석 및 시각화

### 3. 결과 분석만 실행

```bash
# 기존 도킹 결과에 대한 종합 분석
python3 final_analysis.py
```

**분석 내용:**
- 결합 에너지 및 RMSD 분석
- 결합 친화도 평가
- 생물학적 의미 해석
- 향후 연구 방향 제시

## 📊 주요 결과

### 🎯 도킹 성과
- **최고 결합 에너지**: -4.814 kcal/mol
- **결합 친화도**: 중간 정도 (10-100 μM 범위)
- **결합 포즈 수**: 9개의 다양한 형태
- **표적**: 알츠하이머병 베타 아밀로이드 피브릴
- **리간드**: R-flurbiprofen (항염제)

### 📈 통계 요약
- **평균 결합 에너지**: -4.676 ± 0.099 kcal/mol
- **에너지 범위**: 0.282 kcal/mol
- **일관성**: 모든 포즈가 중간 수준의 친화도

## 🔧 스크립트 상세 설명

### `molecular_docking.py`
**기능**: 완전한 분자 도킹 파이프라인
```bash
python3 molecular_docking.py
```
**수행 단계**:
1. 입력 파일 검증
2. 수용체/리간드 준비
3. 결합 부위 계산
4. Vina 설정 파일 생성
5. 도킹 실행
6. 결과 파싱 및 분석

### `final_analysis.py`
**기능**: 도킹 결과 종합 분석
```bash
python3 final_analysis.py
```
**분석 항목**:
- 결합 포즈별 에너지 분석
- 통계적 요약
- 생물학적 해석
- 치료적 함의
- 연구 권장사항

## ⚙️ 설정 옵션

### Vina 매개변수 (molecular_docking.py에서 수정 가능)
```python
# 도킹 매개변수
exhaustiveness = 32        # 탐색 정밀도
num_modes = 20            # 생성할 포즈 수
energy_range = 4          # 에너지 범위 (kcal/mol)

# 검색 공간
center_x, center_y, center_z = auto_calculated
size_x = size_y = size_z = 25  # 검색 박스 크기 (Å)
```

## 📋 출력 파일 설명

### `final_docking_result.pdbqt`
- **형식**: AutoDock Vina PDBQT
- **내용**: 9개 결합 포즈와 에너지 정보
- **구조**: MODEL 1-9, 각각 결합 에너지와 RMSD 포함

### 콘솔 출력
- 실시간 진행 상황
- 결합 에너지 순위
- 통계적 분석 결과
- 생물학적 해석

## 🔍 결과 해석

### 결합 친화도 분류
- **강한 결합**: < -7.0 kcal/mol
- **좋은 결합**: -5.0 ~ -7.0 kcal/mol
- **중간 결합**: -3.0 ~ -5.0 kcal/mol ← **현재 결과**
- **약한 결합**: > -3.0 kcal/mol

### 생물학적 의미
R-flurbiprofen과 베타 아밀로이드 피브릴 간의 중간 정도 결합은 알츠하이머병 치료에서 의미 있는 상호작용을 시사하며, 피브릴 형성 억제나 안정성 방해 가능성을 보여줍니다.

## ⚠️ 문제 해결

### Vina 실행 오류
```bash
# Vina 설치 확인
which vina
vina --version

# 경로 문제 시
export PATH="/usr/local/bin:$PATH"
```

### 메모리 부족
```bash
# exhaustiveness 값 감소
exhaustiveness = 16  # 기본값 32에서 감소
```

### 파일 권한 오류
```bash
chmod 755 *.py
chmod 644 *.pdbqt
```

## 📈 성능 최적화

### 빠른 테스트
```python
# molecular_docking.py에서 수정
exhaustiveness = 8
num_modes = 5
```

### 고정밀 도킹
```python
exhaustiveness = 64
num_modes = 50
energy_range = 6
```

## 🎯 다음 단계

1. **시각화**: PyMOL/ChimeraX로 결합 포즈 확인
2. **검증**: 분자 동역학 시뮬레이션
3. **최적화**: 리간드 구조 개선
4. **실험**: 생화학적 결합 분석
