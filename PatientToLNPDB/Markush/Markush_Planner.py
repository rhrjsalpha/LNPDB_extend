import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import collections
import itertools

class MarkushBFSAssembler:
    def __init__(self, df):
        self.df = df
        self.registry = self._prepare_registry()
        
    def _prepare_registry(self):
        """Label별로 후보 분자들을 그룹화하여 저장"""
        reg = collections.defaultdict(list)
        for _, row in self.df.iterrows():
            mol = Chem.MolFromSmiles(str(row['CXSMILES']))
            if mol:
                reg[str(row['Label'])].append(mol)
        return reg

    def get_actual_label(self, atom):
        """원자의 속성에서 라벨 문자열 추출"""
        props = atom.GetPropsAsDict()
        for key in ['atomLabel', '_label', 'molFileAlias', '_APLabel']:
            if key in props: return str(props[key])
        return ""

    def join_simple(self, parent, child, p_label, c_label):
        """단순 R-group 치환 결합"""
        mw = Chem.RWMol(Chem.CombineMols(parent, child))
        p_idx, c_idx = -1, -1
        p_dummy, c_dummy = -1, -1

        for atom in mw.GetAtoms():
            if atom.GetSymbol() == '*':
                lbl = self.get_actual_label(atom)
                if lbl == p_label and p_idx == -1:
                    if atom.GetNeighbors():
                        p_idx, p_dummy = atom.GetNeighbors()[0].GetIdx(), atom.GetIdx()
                elif lbl == c_label and c_idx == -1:
                    if atom.GetNeighbors():
                        c_idx, c_dummy = atom.GetNeighbors()[0].GetIdx(), atom.GetIdx()

        if p_idx != -1 and c_idx != -1:
            mw.AddBond(p_idx, c_idx, Chem.rdchem.BondType.SINGLE)
            for idx in sorted([p_dummy, c_dummy], reverse=True): mw.RemoveAtom(idx)
            res = mw.GetMol()
            Chem.SanitizeMol(res)
            return res
        return None

    def insert_bridge(self, scaffold, bridge_frag, target_label):
        """R5와 같은 브릿지 삽입 결합 (찌꺼기 제거 로직)"""
        target_idx = -1
        for atom in scaffold.GetAtoms():
            if self.get_actual_label(atom) == target_label:
                target_idx = atom.GetIdx()
                break
        if target_idx == -1: return None

        mw = Chem.RWMol(scaffold)
        neighbors = list(scaffold.GetAtomWithIdx(target_idx).GetNeighbors())
        conn_info = []
        for i, nb in enumerate(neighbors):
            btype = scaffold.GetBondBetweenAtoms(target_idx, nb.GetIdx()).GetBondType()
            conn_info.append((nb.GetIdx(), btype, f"{target_label}_{i+1}"))
        mw.RemoveAtom(target_idx)

        for nb_idx, btype, label in conn_info:
            new_star = mw.AddAtom(Chem.Atom(0))
            mw.GetAtomWithIdx(new_star).SetProp('atomLabel', label)
            mw.AddBond(nb_idx, new_star, btype)

        combined = Chem.RWMol(Chem.CombineMols(mw.GetMol(), bridge_frag))
        for i in range(len(conn_info)):
            label = f"{target_label}_{i+1}"
            p_idx, c_idx, p_dummy, c_dummy = -1, -1, -1, -1
            for atom in combined.GetAtoms():
                if atom.GetSymbol() == '*':
                    if self.get_actual_label(atom) == label:
                        if atom.GetNeighbors():
                            if p_idx == -1: p_idx, p_dummy = atom.GetNeighbors()[0].GetIdx(), atom.GetIdx()
                            else: c_idx, c_dummy = atom.GetNeighbors()[0].GetIdx(), atom.GetIdx()
            if p_idx != -1 and c_idx != -1:
                combined.AddBond(p_idx, c_idx, Chem.rdchem.BondType.SINGLE)
                for d_idx in sorted([p_dummy, c_dummy], reverse=True): combined.RemoveAtom(d_idx)
        
        res = combined.GetMol()
        Chem.SanitizeMol(res)
        return res

    def run(self, core_label):
        """BFS 기반 조립 실행"""
        if core_label not in self.registry:
            print(f"Error: {core_label}이 엑셀에 없습니다.")
            return []

        # 큐 초기화 (CORE 후보들로 시작)
        queue = collections.deque(self.registry[core_label])
        print(f"queue = {queue}")
        final_results = []
        visited_smiles = set()

        print(f"BFS 조립 시작: {core_label}")

        while queue:
            curr = queue.popleft()
            print(f"pop_left {curr}, 큐 길이: {len(queue)} 남음")
            
            # 현재 분자의 소켓 추출
            sockets = [self.get_actual_label(a) for a in curr.GetAtoms() 
                       if a.GetSymbol() == '*' and '_R' in self.get_actual_label(a)
                       and not any(x in self.get_actual_label(a) for x in ['_1', '_2'])]

            if not sockets:
                smi = Chem.MolToSmiles(curr, isomericSmiles=True)
                if smi not in visited_smiles:
                    visited_smiles.add(smi)
                    final_results.append(curr)
                continue

            # 가장 앞에 있는 소켓 하나 해결
            target_socket = sockets[0]
            target_key = target_socket.replace('_', '')
            child_label = next((k for k in self.registry.keys() if target_key in k), None)
            print(f"처리 소켓: {target_socket}, 대응 라벨: {child_label}")

            if not child_label:
                print(f"{target_socket} 매칭 실패, 스킵.")
                continue

            for child_mol in self.registry[child_label]:
                try:
                    # 별표 개수로 Action 결정
                    stars = len([a for a in child_mol.GetAtoms() if a.GetSymbol() == '*'])
                    if stars >= 2:
                        next_mol = self.insert_bridge(curr, child_mol, target_socket)
                    else:
                        next_mol = self.join_simple(curr, child_mol, target_socket, f"{target_socket}_1")
                    
                    if next_mol:
                        queue.append(next_mol)
                except:
                    continue

        print(f"조립 완료: 총 {len(final_results)}개 생성")
        return final_results

# --- 실행부 ---
df = pd.read_excel("/Users/kogeon/python_projects_path/LNPDB_extend/PatientToLNPDB/Markush/WO2021021634_Unified_MultiAP.xlsx")
assembler = MarkushBFSAssembler(df)
results = assembler.run("Main_Skeleton_Original")

import os
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

def save_assembled_images(mol_list, output_dir="assembled_images"):
    """
    조립된 분자 리스트를 받아 SVG 이미지 파일로 저장하는 독립 함수
    """
    # 1. 저장 디렉토리 생성
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"📂 디렉토리 생성됨: {output_dir}")

    print(f"📸 총 {len(mol_list)}개의 이미지 생성을 시작합니다...")

    for i, mol in enumerate(mol_list):
        if mol is None: continue
        
        # 2. 이미지용 2D 좌표 계산
        # 조립 직후의 분자는 좌표가 없으므로 정렬이 필요합니다.
        mol_to_draw = Chem.Mol(mol) # 원본 보호를 위해 복사본 사용
        AllChem.Compute2DCoords(mol_to_draw)
        
        # 3. SVG 드로잉 설정 (600x600 px)
        d2d = rdMolDraw2D.MolDraw2DSVG(600, 600)
        
        # 옵션: 원자 번호나 라벨을 숨기고 구조만 깔끔하게 출력
        dopts = d2d.drawOptions()
        dopts.addStereoAnnotation = True # 입체 정보 표시
        
        d2d.DrawMolecule(mol_to_draw)
        d2d.FinishDrawing()
        
        # 4. 파일 저장
        file_path = os.path.join(output_dir, f"mol_{i+1:03d}.svg")
        with open(file_path, "w") as f:
            f.write(d2d.GetDrawingText())

    print(f"✅ 모든 이미지 저장 완료: {output_dir}")

# --- 실제 사용 예시 ---
# 1. 조립 실행
# assembler = MarkushBFSAssembler(df)
# results = assembler.run("Main_Skeleton_Original")

# 2. 클래스 외부에서 함수 호출
save_assembled_images(results, output_dir="final_output_svgs")