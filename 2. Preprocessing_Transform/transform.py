import os
from rdkit import Chem
from rdkit.Chem import AllChem

# 현재 스크립트의 디렉토리를 기준으로 파일 경로 설정
script_dir = os.path.dirname(os.path.abspath(__file__))
in_sdf = os.path.join(script_dir, "Conformer3D_COMPOUND_CID_92337.sdf")
out_mol2 = os.path.join(script_dir, "R_flurbiprofen.mol2")

try:
    # 입력 파일 존재 확인
    if not os.path.exists(in_sdf):
        raise FileNotFoundError(f"입력 파일을 찾을 수 없습니다: {in_sdf}")
    
    print(f"SDF 파일 읽는 중: {in_sdf}")
    
    # SDF 읽기
    suppl = Chem.SDMolSupplier(in_sdf, removeHs=False)
    
    if len(suppl) == 0:
        raise ValueError("SDF 파일에 분자가 없습니다.")
    
    mol = suppl[0]
    if mol is None:
        raise ValueError("SDF 파싱 실패: 파일 형식을 확인하세요.")
    
    print(f"분자 정보: {mol.GetNumAtoms()}개 원자, {mol.GetNumBonds()}개 결합")
    
    # H 추가
    mol = Chem.AddHs(mol)
    print(f"수소 추가 후: {mol.GetNumAtoms()}개 원자")
    
    # 3D 좌표가 없으면 내장법으로 생성
    if mol.GetNumConformers() == 0:
        print("3D 좌표 생성 중...")
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)
    else:
        print("기존 3D 좌표 사용")
    
    # 부분전하 계산(도킹용 MOL2에 유용)
    print("Gasteiger 전하 계산 중...")
    AllChem.ComputeGasteigerCharges(mol)
    
    # MOL2로 저장 (OpenBabel 사용)
    print(f"MOL2 파일로 저장 중: {out_mol2}")
    
    try:
        # OpenBabel을 사용하여 MOL2 변환
        from openbabel import pybel
        
        # 임시 SDF 파일 생성
        temp_sdf = out_mol2.replace('.mol2', '_temp.sdf')
        writer = Chem.SDWriter(temp_sdf)
        writer.write(mol)
        writer.close()
        
        # OpenBabel로 SDF를 MOL2로 변환
        mol_ob = pybel.readfile("sdf", temp_sdf).__next__()
        mol_ob.write("mol2", out_mol2, overwrite=True)
        
        # 임시 파일 삭제
        os.remove(temp_sdf)
        print("✅ OpenBabel을 사용하여 MOL2 변환 완료")
        
    except (ImportError, Exception) as e:
        print(f"⚠️  OpenBabel MOL2 변환 실패: {e}")
        print("SDF 형식으로 저장합니다.")
        sdf_output = out_mol2.replace('.mol2', '_converted.sdf')
        writer = Chem.SDWriter(sdf_output)
        writer.write(mol)
        writer.close()
        print(f"대신 SDF 파일로 저장됨: {sdf_output}")
        out_mol2 = sdf_output
    
    if os.path.exists(out_mol2):
        print(f"✅ 변환 완료: {out_mol2}")
        print(f"파일 크기: {os.path.getsize(out_mol2)} bytes")
    else:
        raise RuntimeError("MOL2 파일 생성에 실패했습니다.")

except Exception as e:
    print(f"❌ 오류 발생: {e}")
    raise