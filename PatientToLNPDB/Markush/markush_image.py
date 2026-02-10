import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

# 1. 경로 설정
current_dir = os.path.dirname(os.path.abspath(__file__))
# 최신 통합 엑셀 파일명으로 변경
file_path = os.path.join(current_dir, "WO2021021634_Unified_MultiAP.xlsx")
output_dir = os.path.join(current_dir, "structure_svgs")

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 2. 데이터 로드
if not os.path.exists(file_path):
    print(f"오류: '{file_path}' 파일을 찾을 수 없습니다.")
else:
    df = pd.read_excel(file_path)
    print(f"이미지 저장 시작: {output_dir}")

    # 3. 데이터 순회 및 이미지 생성
    for i, row in df.iterrows():
        label = str(row['Label']) # Label 컬럼 활용 (예: R1(Candidate))
        cxsmiles = str(row['CXSMILES']) # 중괄호가 제거된 CXSMILES 컬럼 사용

        if not cxsmiles or cxsmiles == 'nan':
            continue

        # CXSMILES 로드 (RDKit은 |...| 메타데이터를 자동으로 인식함)
        mol = Chem.MolFromSmiles(cxsmiles)
        
        # 실패 시 SMILES 부분만 추출 시도
        if not mol:
            pure_smiles = cxsmiles.split('|')[0].strip()
            mol = Chem.MolFromSmiles(pure_smiles)

        if mol:
            # SVG 드로잉 설정
            d2d = rdMolDraw2D.MolDraw2DSVG(400, 400)
            d2d.DrawMolecule(mol)
            d2d.FinishDrawing()
            
            # 파일명에 Label 포함 (특수문자 제거)
            clean_label = "".join(c for c in label if c.isalnum() or c in (' ', '_')).rstrip()
            filename = f"{i:03d}_{clean_label}.svg"
            save_path = os.path.join(output_dir, filename)
            
            with open(save_path, "w") as f:
                f.write(d2d.GetDrawingText())
            
            print(f"성공 [{i}]: {filename} 저장 완료")
        else:
            print(f"실패 [{i}]: {label} - 구조를 인식할 수 없습니다.")

    print(f"\n✅ 작업 완료! SVG 파일 확인 경로:\n{output_dir}")