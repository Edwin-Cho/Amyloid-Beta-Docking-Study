# SwissDock Molecular Docking Project
## R-flurbiprofen과 Amyloid-Beta 구조 비교 도킹 연구

> **Language**: [English]([EN]%20Overview.md) | [한국어]

이 프로젝트는 알츠하이머병 치료 후보 물질인 R-flurbiprofen과 두 가지 아밀로이드 베타 구조(1IYT, 2M4J) 간의 분자 도킹 연구를 수행합니다.

---

## 🎯 프로젝트 개요

### 연구 목적
- 용해성 아밀로이드 베타(1IYT)와 피브릴 형태(2M4J)의 결합 특성 비교
- R-flurbiprofen의 치료적 잠재력 평가
- 질병 진행 단계별 맞춤형 치료 전략 제시

### 주요 발견
- **2M4J 최고 결합 에너지**: -4.814 kcal/mol (중간 친화도)
- **구조적 차이**: 용해성 vs 피브릴 형태의 명확한 결합 특성 구분
- **치료적 함의**: 질병 단계별 차별화된 접근법 필요성 확인

---

## 📁 프로젝트 구조

```
SwissDock Project/
├── 0. Results/                    # 📊 최종 결과 및 분석 보고서
│   ├── [EN] Comparative_Docking_Analysis.md
│   ├── [KOR] Comparative_Docking_Analysis.md
│   ├── 1IYT_Results/             # SwissDock 결과
│   └── 2M4J_Results/             # AutoDock Vina 결과
├── 1. Materials/                  # 🧬 원본 구조 파일
│   ├── 1iyt_single.pdb
│   ├── 2m4j_single.pdb
│   ├── R_flurbiprofen_NCI.mol2
│   ├── [KOR] Materials.md
│   └── [EN] Materials.md
├── 2. Preprocessing_Transform/    # ⚙️ 구조 전처리 스크립트
│   ├── *.py (5개 Python 스크립트)
│   ├── *.sh (2개 Shell 스크립트)
│   ├── [KOR] Processing_Transform.md
│   └── [EN] Processing_Transform.md
├── 3. Docking_Local/             # 🎯 로컬 도킹 실행
│   ├── molecular_docking.py
│   ├── final_analysis.py
│   ├── [KOR] Docking_Local.md
│   ├── [EN] Docking_Local.md
│   └── Materials/
└── 4. Archives/                  # 📚 참고 자료
```

---

## 🚀 빠른 시작 가이드

### 1. 환경 설정

```bash
# Python 패키지 설치
pip install biopython numpy pandas matplotlib

# AutoDock Vina 설치 (macOS)
brew install autodock-vina

# Open Babel 설치
brew install open-babel

# MGLTools 설치 (수동)
# https://ccsb.scripps.edu/mgltools/downloads/
```

### 2. 전체 워크플로우 실행

```bash
# 1단계: 구조 전처리
cd "2. Preprocessing_Transform"
./prepare_2M4J_receptor.sh

# 2단계: 분자 도킹 실행
cd "../3. Docking_Local"
python3 molecular_docking.py

# 3단계: 결과 분석
python3 final_analysis.py

# 4단계: 결과 확인
cd "../0. Results"
open "[KOR] Comparative_Docking_Analysis.md"
```

---

## 📋 각 폴더별 상세 가이드

### 🧬 1. Materials/
**원본 구조 파일 저장소**
- 1IYT: 용해성 아밀로이드 베타 (1-42)
- 2M4J: 아밀로이드 베타 피브릴 (1-40)  
- R-flurbiprofen: 치료 후보 리간드

📖 [한국어 가이드](1.%20Materials/[KOR]%20Materials.md) | [English Guide](1.%20Materials/[EN]%20Materials.md)

### ⚙️ 2. Preprocessing_Transform/
**구조 전처리 및 변환 스크립트**

```bash
# 자동 전처리 (권장)
./prepare_2M4J_receptor.sh

# 수동 단계별 실행
python create_single_model.py "PDB 2M4J.pdb" 2m4j_single.pdb
python fix_pdb_headers.py 2m4j_single.pdb 2m4j_fixed.pdb
python prepare_receptor.py --input 2m4j_fixed.pdb --output 2m4j.pdbqt
```

📖 [한국어 가이드](2.%20Preprocessing_Transform/[KOR]%20Processing_Transform.md) | [English Guide](2.%20Preprocessing_Transform/[EN]%20Processing_Transform.md)

### 🎯 3. Docking_Local/
**분자 도킹 실행 및 분석**

```bash
# 완전한 도킹 워크플로우
python3 molecular_docking.py

# 결과 분석만 실행
python3 final_analysis.py
```

**주요 결과:**
- 최고 결합 에너지: -4.814 kcal/mol
- 9개 결합 포즈 생성
- 중간 수준 결합 친화도 (10-100 μM)

📖 [한국어 가이드](3.%20Docking_Local/[KOR]%20Docking_Local.md) | [English Guide](3.%20Docking_Local/[EN]%20Docking_Local.md)

### 📊 0. Results/
**최종 분석 결과 및 보고서**
- 영문/한국어 비교 분석 보고서
- 개별 도킹 결과 파일
- 통계 분석 및 생물학적 해석

📖 [상세 가이드](0.%20Results/README.md)

---

## 🔬 과학적 기여

### 새로운 발견
1. **구조별 결합 차이**: 용해성 vs 피브릴 형태의 명확한 구분
2. **치료 전략**: 질병 단계별 맞춤형 접근법
3. **약물 설계**: R-flurbiprofen 최적화 방향 제시

### 임상적 함의
- 알츠하이머병 진행 단계별 치료 전략
- 아밀로이드 표적 치료제 개발 지침
- 개인 맞춤형 치료 가능성

---

## 📈 주요 결과 요약

| 구조 | 도킹 도구 | 최고 결합 에너지 | 결합 특성 | 치료적 접근 |
|------|-----------|------------------|-----------|-------------|
| 1IYT | SwissDock | 다양한 모드 | 유연한 결합 | 예방적 |
| 2M4J | AutoDock Vina | -4.814 kcal/mol | 제한된 결합 | 치료적 |

---

## 🛠️ 문제 해결

### 일반적인 오류

**Vina 실행 오류**
```bash
# Vina 설치 확인
which vina
vina --version
```

**Python 패키지 오류**
```bash
pip install --upgrade biopython numpy pandas matplotlib
```

**파일 권한 오류**
```bash
chmod +x *.sh
chmod +x *.py
```

---

## 📚 참고 문헌

### 주요 논문
1. **1IYT**: Crescenzi et al. (2002) *Eur J Biochem* 269: 5642-5648
2. **2M4J**: Lu et al. (2013) *Cell* 154: 1257-1268
3. **R-flurbiprofen**: Wilcock et al. (2008) *Lancet Neurol* 7: 483-493

### 도구 및 데이터베이스
- [AutoDock Vina](https://vina.scripps.edu/)
- [SwissDock](http://www.swissdock.ch/)
- [Protein Data Bank](https://www.rcsb.org/)

---

**마지막 업데이트**: 2026년 02월 25일