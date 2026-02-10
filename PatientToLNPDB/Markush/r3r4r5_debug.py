import pandas as pd
from rdkit import Chem
import os

# 1. 데이터 로드 및 유틸리티
current_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(current_dir, "WO2021021634_Unified_MultiAP.xlsx")
df = pd.read_excel(file_path)

def get_actual_label(atom):
    props = atom.GetPropsAsDict()
    candidates = ['atomLabel', '_label', '_APLabel', 'molFileAlias']
    for key in candidates:
        if key in props: return str(props[key])
    return ""

def print_labels(mol, title):
    labels = [get_actual_label(a) for a in mol.GetAtoms() if a.GetSymbol() == '*']
    print(f"    [{title}] 별표 라벨 목록: {labels}")

# 2. [핵심] 분해 및 삽입 로직
def decompose_and_insert_bridge(scaffold, bridge_frag, target_label):
    """
    Scaffold에서 target_label(_R5) 원자를 제거하고, 
    그 자리에 bridge_frag(R5)를 삽입하여 연결함.
    """
    print(f"\n--- [단계: {target_label} 분해 및 삽입 시작] ---")
    
    # A. 대상 원자 찾기
    target_idx = -1
    for atom in scaffold.GetAtoms():
        if get_actual_label(atom) == target_label:
            target_idx = atom.GetIdx()
            break
            
    if target_idx == -1:
        print(f"    ❌ '{target_label}' 원자를 찾지 못했습니다.")
        return scaffold

    # B. 연결 정보 기록 (이웃 원자와 결합 타입)
    target_atom = scaffold.GetAtomWithIdx(target_idx)
    neighbors = []
    for nb in target_atom.GetNeighbors():
        bond = scaffold.GetBondBetweenAtoms(target_idx, nb.GetIdx())
        neighbors.append((nb.GetIdx(), bond.GetBondType()))
    
    print(f"    - '{target_label}' 원자의 이웃 수: {len(neighbors)}개")

    # C. 원자 제거 및 새로운 소켓 생성
    mw = Chem.RWMol(scaffold)
    # 안전하게 이웃 원자에 임시 태그 부여 (인덱스 변화 대비)
    for i, (nb_idx, _) in enumerate(neighbors):
        mw.GetAtomWithIdx(nb_idx).SetIntProp('temp_tag', i + 1)
    
    mw.RemoveAtom(target_idx)
    
    # 끊어진 자리마다 새로운 별표(*) 소켓 추가
    for i in range(len(neighbors)):
        new_star_idx = mw.AddAtom(Chem.Atom(0)) # dummy atom
        mw.GetAtomWithIdx(new_star_idx).SetProp('atomLabel', f"{target_label}_{i+1}")
        
        # 태그를 통해 이웃 원자 다시 찾아 연결
        for atom in mw.GetAtoms():
            if atom.HasProp('temp_tag') and atom.GetIntProp('temp_tag') == i + 1:
                mw.AddBond(atom.GetIdx(), new_star_idx, neighbors[i][1])
                atom.ClearProp('temp_tag')
                break

    decomposed_scaffold = mw.GetMol()
    print_labels(decomposed_scaffold, "분해 후 Scaffold")

    # D. Bridge Fragment(R5) 연결
    # bridge_frag에는 _R5_1, _R5_2가 있다고 가정함
    temp = decomposed_scaffold
    for i in range(len(neighbors)):
        label = f"{target_label}_{i+1}"
        # bridge_frag의 _R5_1, _R5_2를 순차적으로 매칭
        temp = join_by_label_simple(temp, bridge_frag, label, label)
        if temp:
            print(f"    ✅ {label} 지점 연결 성공")
        else:
            print(f"    ❌ {label} 지점 연결 실패")
            
    return temp

def join_by_label_simple(parent, child, p_label, c_label):
    if not parent or not child: return None
    mw = Chem.RWMol(Chem.CombineMols(parent, child))
    p_idx, c_idx = -1, -1
    p_dummy, c_dummy = -1, -1
    for atom in mw.GetAtoms():
        if atom.GetSymbol() == '*':
            val = get_actual_label(atom)
            if val == p_label and p_idx == -1:
                if atom.GetNeighbors():
                    p_idx = atom.GetNeighbors()[0].GetIdx(); p_dummy = atom.GetIdx()
            elif val == c_label and c_idx == -1:
                if atom.GetNeighbors():
                    c_idx = atom.GetNeighbors()[0].GetIdx(); c_dummy = atom.GetIdx()
    if p_idx != -1 and c_idx != -1:
        mw.AddBond(p_idx, c_idx, Chem.rdchem.BondType.SINGLE)
        for idx in sorted([p_dummy, c_dummy], reverse=True): mw.RemoveAtom(idx)
        res = mw.GetMol()
        Chem.SanitizeMol(res); return res
    return None

# 3. 데이터 추출 및 실행
groups = {}
for _, row in df[df['Category'] == 'R_GROUP'].iterrows():
    label_key = str(row['Label']).split('(')[0]
    mol = Chem.MolFromSmiles(str(row['CXSMILES']))
    if mol: groups.setdefault(label_key, []).append(mol)

# 진단 대상: 20번 행 Scaffold 및 R5 Fragment
r2_scaffold = [m for m in groups.get('R2', []) if any('_R3' in get_actual_label(a) for a in m.GetAtoms())][0]
r5_bridge = groups.get('R5', [])[0]

# 실행
final_complex = decompose_and_insert_bridge(r2_scaffold, r5_bridge, '_R5')

if final_complex:
    print_labels(final_complex, "최종 복합체")
    print(f"SMILES: {Chem.MolToSmiles(final_complex)}")