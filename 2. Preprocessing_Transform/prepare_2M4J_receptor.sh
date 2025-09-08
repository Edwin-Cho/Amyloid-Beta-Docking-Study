#!/bin/bash
# AutoDockTools를 사용한 수용체 전처리 스크립트

echo "🧬 AutoDockTools로 수용체 전처리 시작: /Users/edwinrcho/Downloads/siwss doc/2M4J_receptor_clean.pdb"

# 1. 수용체 전처리 (Python 2.7 환경에서 실행)
python2.7 -c "
import sys
sys.path.append('/usr/local/lib/python2.7/site-packages/AutoDockTools')
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation

# 수용체 전처리
prep = AD4ReceptorPreparation()
prep.prepare_receptor('/Users/edwinrcho/Downloads/siwss doc/2M4J_receptor_clean.pdb', outputfilename='2M4J_receptor.pdbqt', 
                     repairs='bonds_hydrogens', 
                     charges_to_add='gasteiger',
                     cleanup='nphs_lps_waters_nonstdres')
print('✅ 수용체 전처리 완료: 2M4J_receptor.pdbqt')
"

# 대안: MGLTools의 prepare_receptor4.py 사용
if [ -f "/usr/local/bin/prepare_receptor4.py" ]; then
    echo "🔧 prepare_receptor4.py 사용"
    python2.7 /usr/local/bin/prepare_receptor4.py -r /Users/edwinrcho/Downloads/siwss doc/2M4J_receptor_clean.pdb -o 2M4J_receptor.pdbqt -A hydrogens -U nphs_lps_waters_nonstdres
elif [ -f "/opt/mgltools/bin/prepare_receptor4.py" ]; then
    echo "🔧 prepare_receptor4.py 사용 (opt 경로)"
    python2.7 /opt/mgltools/bin/prepare_receptor4.py -r /Users/edwinrcho/Downloads/siwss doc/2M4J_receptor_clean.pdb -o 2M4J_receptor.pdbqt -A hydrogens -U nphs_lps_waters_nonstdres
else
    echo "⚠️  AutoDockTools가 설치되지 않았습니다."
    echo "다음 명령어로 설치하세요:"
    echo "conda install -c bioconda autodock-vina"
    echo "또는 MGLTools를 다운로드하세요: http://mgltools.scripps.edu/"
fi

echo "🎯 완료!"
