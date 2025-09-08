# Materials - 원본 구조 파일

> **Language**: [English]([EN]%20Materials.md) | [한국어]

이 폴더는 분자 도킹 연구에 사용되는 원본 구조 파일들을 포함합니다.

## 📁 파일 구조

### 단백질 구조 (PDB 파일)
- `1iyt_single.pdb` - 1IYT 단일 모델 구조 (용해성 아밀로이드 베타 펩타이드)
- `2m4j_single.pdb` - 2M4J 단일 모델 구조 (아밀로이드 베타 피브릴)

### 리간드 구조 (MOL2 파일)
- `R_flurbiprofen_NCI.mol2` - R-flurbiprofen 3D 구조 (NCI 데이터베이스)

## 📋 파일 설명

### 1IYT - 용해성 아밀로이드 베타 (1-42)
- **PDB ID**: 1IYT
- **분류**: 단백질 결합
- **실험 방법**: 용액 NMR
- **구조 유형**: 단량체 펩타이드 (무극성 미세환경)
- **잔기 수**: 42개 아미노산 (Aβ1-42)
- **형태**: 알파 헬릭스가 풍부한 막 결합 형태
- **생물학적 맥락**: 초기 단계 용해성 아밀로이드 종
- **등록일**: 2002-09-06

### 2M4J - 아밀로이드 베타 피브릴 (1-40)
- **PDB ID**: 2M4J
- **분류**: 단백질 피브릴
- **실험 방법**: 고체상태 NMR
- **구조 유형**: 환자 뇌조직 유래 피브릴 응집체
- **잔기 수**: 40개 아미노산 (Aβ1-40)
- **형태**: 베타 시트가 풍부한 피브릴 구조
- **생물학적 맥락**: 후기 단계 병리학적 아밀로이드 침착
- **등록일**: 2013-02-05

### R-flurbiprofen
- **화학식**: C₁₅H₁₃FO₂
- **분자량**: 244.26 g/mol
- **IUPAC명**: (2R)-2-(3-fluoro-4-phenylphenyl)propanoic acid
- **약물 분류**: 비스테로이드성 항염제 (NSAID)
- **치료적 관심**: 알츠하이머병 치료 후보 물질
- **작용 기전**: 아밀로이드 응집 억제

## 🔍 사용 목적

이 파일들은 다음 연구 목적으로 사용됩니다:

1. **구조 비교**: 용해성 vs 피브릴 형태의 아밀로이드 베타
2. **도킹 연구**: R-flurbiprofen의 결합 특성 분석
3. **치료제 개발**: 알츠하이머병 치료 전략 연구
4. **구조-활성 관계**: 분자 상호작용 이해

## 🚀 다음 단계

이 원본 파일들은 다음 폴더에서 처리됩니다:

1. **2. Preprocessing_Transform/**: 구조 전처리 및 변환
2. **3. Docking_Local/**: 분자 도킹 실행
3. **0. Results/**: 결과 분석 및 보고서

## 📊 구조적 특징

### 1IYT vs 2M4J 주요 차이점

| 특성 | 1IYT (용해성) | 2M4J (피브릴) |
|------|---------------|---------------|
| 구조 상태 | 단량체 | 응집체 |
| 주요 2차 구조 | α-헬릭스 | β-시트 |
| 유연성 | 높음 | 낮음 |
| 막 결합 | 있음 | 없음 |
| 병리학적 단계 | 초기 | 후기 |
| 치료 표적 | 예방적 | 치료적 |

## ⚠️ 파일 처리 주의사항

1. **파일 무결성**: 원본 파일은 수정하지 않고 복사본 사용
2. **좌표계**: PDB 파일의 좌표계 일관성 확인
3. **누락 원자**: 구조 완전성 검증 필요
4. **수소 원자**: 대부분 실험 구조에서 누락됨 (전처리에서 추가)

## 📖 참고 문헌

### 1IYT 관련
- Crescenzi, O. et al. (2002) "Solution structure of the Alzheimer amyloid β-peptide (1-42) in an apolar microenvironment" *Eur J Biochem* 269: 5642-5648

### 2M4J 관련
- Lu, J.X. et al. (2013) "Molecular structure of β-amyloid fibrils in Alzheimer's disease brain tissue" *Cell* 154: 1257-1268

### R-flurbiprofen 관련
- Wilcock, G.K. et al. (2008) "Efficacy and safety of tarenflurbil in mild to moderate Alzheimer's disease" *Lancet Neurol* 7: 483-493
