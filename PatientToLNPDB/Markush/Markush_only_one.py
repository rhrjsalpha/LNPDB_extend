import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import itertools
import os

# ======================================================
# 1. 경로 설정 및 데이터 로드
# ======================================================
current_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(current_dir, "marvinfile.xlsx")
image_out_dir = os.path.join(current_dir, "marvin")
os.makedirs(image_out_dir, exist_ok=True)

df = pd.read_excel(file_path)

# ======================================================
# --- [도구 함수군] ---
# ======================================================

def get_actual_label(atom):
    """원자의 속성(atomLabel, molFileAlias 등)을 뒤져 라벨 추출"""
    props = atom.GetPropsAsDict()
    candidates = ['atomLabel', '_label', '_APLabel', 'molFileAlias', '_AtomLabel']
    for key in candidates:
        if key in props:
            return str(props[key])
    return ""

def get_mol_labels(mol):
    """현재 분자가 가진 모든 별표(*) 라벨 목록 반환"""
    return [get_actual_label(a) for a in mol.GetAtoms() if a.GetSymbol() == '*']

def join_by_label_debug(parent, child, p_label, c_label, step_name=""):
    """결합 시도 전후 상태를 체크하며 안전하게 결합을 수행하는 함수"""
    if parent is None or child is None:
        return None

    mw = Chem.RWMol(Chem.CombineMols(parent, child))
    p_idx, c_idx = -1, -1
    p_dummy, c_dummy = -1, -1

    # parent/child 더미 탐색
    for atom in mw.GetAtoms():
        if atom.GetSymbol() != '*':
            continue
        val = get_actual_label(atom)

        if val == p_label and p_idx == -1:
            if atom.GetNeighbors():
                p_idx = atom.GetNeighbors()[0].GetIdx()
                p_dummy = atom.GetIdx()

        elif val == c_label and c_idx == -1:
            if atom.GetNeighbors():
                c_idx = atom.GetNeighbors()[0].GetIdx()
                c_dummy = atom.GetIdx()

    if p_idx != -1 and c_idx != -1:
        mw.AddBond(p_idx, c_idx, Chem.rdchem.BondType.SINGLE)

        # 더미 제거 (큰 idx부터 제거)
        for idx in sorted([p_dummy, c_dummy], reverse=True):
            mw.RemoveAtom(idx)

        res = mw.GetMol()
        try:
            Chem.SanitizeMol(res)
            return res
        except Exception as e:
            print(f"  [❌ Sanitize 실패] {step_name}: {e}")
            return None

    # 결합 실패
    return None

# ======================================================
# 2. R-그룹 분류 (엑셀 기반 유지)
# ======================================================
groups = {}
for _, row in df[df['Category'] == 'R_GROUP'].iterrows():
    label_key = str(row['Label']).split('(')[0]  # 'R1', 'R2', 'R3' ...
    mol = Chem.MolFromSmiles(str(row['CXSMILES']))
    if mol:
        groups.setdefault(label_key, []).append(mol)

# ======================================================
# 3. Core 로드
# ======================================================
core_row = df[df['Category'] == 'CORE'].iloc[0]
core_mol = Chem.MolFromSmiles(str(core_row['CXSMILES']))
if core_mol is None:
    raise ValueError("❌ CORE CXSMILES를 RDKit이 파싱하지 못했습니다.")

# (선택) 라벨 확인 출력
print("\n--- [라벨 점검] ---")
print("CORE labels:", get_mol_labels(core_mol))
for k in ['R1', 'R2', 'R3']:
    if k in groups and len(groups[k]) > 0:
        print(f"{k} labels(example):", get_mol_labels(groups[k][0]))

# ======================================================
# 4. Enumeration: CORE + R1 + R2 + R3
#    (요청: R3까지만)
# ======================================================
r1_pool = groups.get('R1', [])
r2_pool = groups.get('R2', [])
r3_pool = groups.get('R3', [])

print(f"\n--- [후보군 크기] ---")
print(f"R1: {len(r1_pool)}개, R2: {len(r2_pool)}개, R3: {len(r3_pool)}개")

final_smiles = set()
count = 0

print(f"\n--- [조립 시작] CORE -> R1 -> R2 -> R3 ---")

for r1 in r1_pool:
    inter1 = join_by_label_debug(core_mol, r1, '_R1', '_R1_1', "CORE+R1")
    if not inter1:
        continue

    for r2 in r2_pool:
        inter2 = join_by_label_debug(inter1, r2, '_R2', '_R2_1', "CORE+R1+R2")
        if not inter2:
            continue

        for r3 in r3_pool:
            final = join_by_label_debug(inter2, r3, '_R3', '_R3_1', "CORE+R1+R2+R3")
            if not final:
                continue

            smi = Chem.MolToSmiles(final, isomericSmiles=True)
            if smi in final_smiles:
                continue

            final_smiles.add(smi)
            count += 1

            # SVG 저장
            d2d = rdMolDraw2D.MolDraw2DSVG(600, 600)
            d2d.DrawMolecule(final)
            d2d.FinishDrawing()

            with open(os.path.join(image_out_dir, f"result_{count:03d}.svg"), "w") as f:
                f.write(d2d.GetDrawingText())

print(f"\n✅ 완료! 총 {len(final_smiles)}개의 고유 화합물 생성")
print(f"📁 SVG 저장 폴더: {image_out_dir}")
