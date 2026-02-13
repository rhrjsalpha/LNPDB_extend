import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import itertools

# 1. 기초 데이터 설정 (이미지 및 CXSMILES 기반 하드코딩)
# Core: C1CCC(C1)c2ccccc2 (사이클로펜틸-벤젠)
core_smi = "C1CCC(C1)C1=CC=CC=C1"
core = Chem.MolFromSmiles(core_smi)

# R1 후보군 (Halogens)
r1_candidates = [
    {"name": "F", "mol": Chem.MolFromSmiles("F*")},
    {"name": "Cl", "mol": Chem.MolFromSmiles("Cl*")},
    {"name": "Br", "mol": Chem.MolFromSmiles("Br*")},
    {"name": "I", "mol": Chem.MolFromSmiles("I*")}
]

# R2 후보군 (Alkyls)
r2_candidates = [
    {"name": "Ethyl", "mol": Chem.MolFromSmiles("CC*")},
    {"name": "Propyl", "mol": Chem.MolFromSmiles("CCC*")},
    {"name": "Butyl", "mol": Chem.MolFromSmiles("CCCC*")}
]

# 위치 변동 인덱스 (CXSMILES m: 태그 기준 매핑)
# 주의: RDKit 인덱스는 0부터 시작하므로 필요시 조정 필요
r1_positions = [9, 10, 11]  # 벤젠 고리의 ortho, meta, para 위치
r2_positions = [0, 1, 2]     # 사이클로펜탄 고리의 상단/측면 위치

def assemble_fixed_position(core_mol, r1_info, r1_pos, r2_info, r2_pos):
    """
    특정 위치가 지정된 상태로 R1과 R2를 강제 결합합니다.
    """
    # 1. R1 결합
    combo1 = Chem.RWMol(Chem.CombineMols(core_mol, r1_info['mol']))
    # R1의 '*' 원자 인덱스 찾기
    r1_star = [a.GetIdx() for a in combo1.GetAtoms() if a.GetSymbol() == '*'][0]
    r1_neighbor = combo1.GetAtomWithIdx(r1_star).GetNeighbors()[0].GetIdx()
    
    # 지정된 위치(r1_pos)와 결합 후 '*' 제거
    combo1.AddBond(r1_pos, r1_neighbor, Chem.rdchem.BondType.SINGLE)
    combo1.RemoveAtom(r1_star)
    
    # 2. R2 결합 (이미 R1이 붙은 분자에서 진행)
    final_mw = Chem.RWMol(Chem.CombineMols(combo1.GetMol(), r2_info['mol']))
    r2_star = [a.GetIdx() for a in final_mw.GetAtoms() if a.GetSymbol() == '*'][0]
    r2_neighbor = final_mw.GetAtomWithIdx(r2_star).GetNeighbors()[0].GetIdx()
    
    final_mw.AddBond(r2_pos, r2_neighbor, Chem.rdchem.BondType.SINGLE)
    final_mw.RemoveAtom(r2_star)
    
    res = final_mw.GetMol()
    Chem.SanitizeMol(res)
    return res

# 2. 모든 조합 실행 (4 R1 x 3 R2 x 3 Pos1 x 3 Pos2 = 108개 가능성)
output_dir = "position_variation_results"
if not os.path.exists(output_dir): os.makedirs(output_dir)

count = 0
for r1, p1, r2, p2 in itertools.product(r1_candidates, r1_positions, r2_candidates, r2_positions):
    try:
        product = assemble_fixed_position(core, r1, p1, r2, p2)
        count += 1
        
        # 이미지 저장
        AllChem.Compute2DCoords(product)
        d2d = rdMolDraw2D.MolDraw2DSVG(400, 400)
        d2d.DrawMolecule(product)
        d2d.FinishDrawing()
        
        filename = f"res_{count:03d}_{r1['name']}_at{p1}_{r2['name']}_at{p2}.svg"
        with open(os.path.join(output_dir, filename), "w") as f:
            f.write(d2d.GetDrawingText())
            
    except Exception as e:
        print(f"Error at combination {count}: {e}")

print(f"✅ 총 {count}개의 위치 변동 화합물 생성 완료!")