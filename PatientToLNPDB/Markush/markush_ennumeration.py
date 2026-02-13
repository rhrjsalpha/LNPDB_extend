import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import itertools
import os

# 1. 경로 설정 및 데이터 로드
current_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(current_dir, "WO2021021634_Unified_MultiAP.xlsx")
image_out_dir = os.path.join(current_dir, "final_enumeration_svgs")
if not os.path.exists(image_out_dir): os.makedirs(image_out_dir)

df = pd.read_excel(file_path)

# --- [도구 함수군] ---

def get_actual_label(atom):
    """원자의 속성(atomLabel, molFileAlias 등)을 뒤져 라벨 추출"""
    props = atom.GetPropsAsDict()
    candidates = ['atomLabel', '_label', '_APLabel', 'molFileAlias', '_AtomLabel']
    for key in candidates:
        if key in props: return str(props[key])
    return ""

def get_mol_labels(mol):
    """현재 분자가 가진 모든 별표(*) 라벨 목록 반환"""
    return [get_actual_label(a) for a in mol.GetAtoms() if a.GetSymbol() == '*']

def join_by_label_debug(parent, child, p_label, c_label, step_name=""):
    """결합 시도 전후 상태를 체크하며 안전하게 결합을 수행하는 함수"""
    if parent is None or child is None: return None
    
    mw = Chem.RWMol(Chem.CombineMols(parent, child))
    p_idx, c_idx = -1, -1
    p_dummy, c_dummy = -1, -1

    for atom in mw.GetAtoms():
        if atom.GetSymbol() == '*':
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
        for idx in sorted([p_dummy, c_dummy], reverse=True): mw.RemoveAtom(idx)
        res = mw.GetMol()
        try:
            Chem.SanitizeMol(res)
            return res
        except Exception as e:
            print(f"  [❌ Sanitize 실패] {step_name}: {e}")
            return None
    return None

def decompose_and_insert_bridge(scaffold, bridge_frag, target_label):
    """
    R5 원자를 제거하고 생긴 두 소켓에 
    하나의 bridge_frag를 가져와 양쪽을 동시에 연결합니다.
    """
    target_idx = -1
    for atom in scaffold.GetAtoms():
        if get_actual_label(atom) == target_label:
            target_idx = atom.GetIdx()
            break
    if target_idx == -1: return None

    # 1. Scaffold 분해 및 소켓 생성
    mw = Chem.RWMol(scaffold)
    neighbors = list(scaffold.GetAtomWithIdx(target_idx).GetNeighbors())
    
    # 연결될 이웃 원자들의 인덱스와 결합 타입 저장
    conn_info = []
    for i, nb in enumerate(neighbors):
        btype = scaffold.GetBondBetweenAtoms(target_idx, nb.GetIdx()).GetBondType()
        conn_info.append((nb.GetIdx(), btype, f"{target_label}_{i+1}"))
    
    mw.RemoveAtom(target_idx)
    
    # Scaffold에 소켓 원자(*) 추가
    for nb_idx, btype, label in conn_info:
        new_star = mw.AddAtom(Chem.Atom(0))
        mw.GetAtomWithIdx(new_star).SetProp('atomLabel', label)
        mw.AddBond(nb_idx, new_star, btype)

    # 2. 하나의 bridge_frag를 단 한 번만 합침
    combined = Chem.RWMol(Chem.CombineMols(mw.GetMol(), bridge_frag))
    
    # 3. 양쪽 연결점(Socket <-> Frag) 매칭 및 결합
    # 예: scaffold의 _R5_1과 bridge의 _R5_1을 찾아 연결
    for i in range(len(conn_info)):
        label = f"{target_label}_{i+1}"
        p_idx, c_idx = -1, -1
        p_dummy, c_dummy = -1, -1
        
        for atom in combined.GetAtoms():
            if atom.GetSymbol() == '*':
                lbl = get_actual_label(atom)
                if lbl == label:
                    if atom.GetNeighbors():
                        # 이웃이 있으면 실제 분자 쪽 인덱스, 본인은 dummy 인덱스
                        if p_idx == -1: 
                            p_idx = atom.GetNeighbors()[0].GetIdx()
                            p_dummy = atom.GetIdx()
                        else:
                            c_idx = atom.GetNeighbors()[0].GetIdx()
                            c_dummy = atom.GetIdx()
        
        if p_idx != -1 and c_idx != -1:
            combined.AddBond(p_idx, c_idx, Chem.rdchem.BondType.SINGLE)
            # 연결 후 더미 원자 제거 (인덱스 변화 방지를 위해 역순 제거)
            for d_idx in sorted([p_dummy, c_dummy], reverse=True):
                combined.RemoveAtom(d_idx)

    res = combined.GetMol()
    try:
        Chem.SanitizeMol(res)
        return res
    except:
        return None
    
    # 브릿지(R5) 조각과 연결 (_R5_SOCKET_n <-> _R5_n)
    temp = mw.GetMol()
    for i in range(len(neighbors)):
        label = f"{target_label}_{i+1}"
        temp = join_by_label_debug(temp, bridge_frag, label, label, f"브릿지 {label} 연결")
    return temp

# --- [메인 로직] ---

# 2. R-그룹 분류
groups = {}
for _, row in df[df['Category'] == 'R_GROUP'].iterrows():
    label_key = str(row['Label']).split('(')[0]
    mol = Chem.MolFromSmiles(str(row['CXSMILES']))
    if mol: groups.setdefault(label_key, []).append(mol)

# 3. [1단계] R2 복합체 조립
r2_final_pool = []
r2_candidates = groups.get('R2', [])
r2_scaffolds = [m for m in r2_candidates if any('_R3' in get_actual_label(a) for a in m.GetAtoms() if a.GetSymbol() == '*')]
r2_normals = [m for m in r2_candidates if m not in r2_scaffolds]

print(f"\n--- [1단계] R2 복합체 조립 (R5 분해 로직 적용) ---")
for s_idx, scaffold in enumerate(r2_scaffolds):
    combos = list(itertools.product(groups.get('R3', [None]), 
                                   groups.get('R4', [None]), 
                                   groups.get('R5', [None])))
    for c_idx, (r3, r4, r5) in enumerate(combos):
        temp = scaffold
        # R3, R4 결합
        if r3: temp = join_by_label_debug(temp, r3, '_R3', '_R3_1', f"S{s_idx}: R3")
        if r4 and temp: temp = join_by_label_debug(temp, r4, '_R4', '_R4_1', f"S{s_idx}: R4")
        # R5 분해 및 삽입 (사용자 제안 로직)
        if r5 and temp: temp = decompose_and_insert_bridge(temp, r5, '_R5')
        
        if temp: r2_final_pool.append(temp)

r2_final_pool.extend(r2_normals)
print(f"  => R2 후보군 준비 완료: {len(r2_final_pool)}개")

# 4. [2단계] 최종 CORE 조립 및 저장
core_mol = Chem.MolFromSmiles(str(df[df['Category'] == 'CORE'].iloc[0]['CXSMILES']))
final_smiles = set()
count = 0

print(f"\n--- [2단계] 최종 CORE 조립 시작 ---")
for r1 in groups.get('R1', []):
    inter = join_by_label_debug(core_mol, r1, '_R1', '_R1_1', "CORE+R1")
    if not inter: continue
    
    for r2 in r2_final_pool:
        # 모든 R2 후보는 공통적으로 '_R2_1' 플러그를 통해 연결됨
        final = join_by_label_debug(inter, r2, '_R2', '_R2_1', "CORE+R2")
        if final:
            smi = Chem.MolToSmiles(final, isomericSmiles=True)
            if smi not in final_smiles:
                final_smiles.add(smi)
                count += 1
                # SVG 저장
                d2d = rdMolDraw2D.MolDraw2DSVG(600, 600)
                d2d.DrawMolecule(final)
                d2d.FinishDrawing()
                with open(os.path.join(image_out_dir, f"result_{count:03d}.svg"), "w") as f:
                    f.write(d2d.GetDrawingText())

print(f"\n✅ 완료! 총 {len(final_smiles)}개의 고유 화합물 생성")