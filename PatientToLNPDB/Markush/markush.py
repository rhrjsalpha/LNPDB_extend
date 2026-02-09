import re
import os
from rdkit import Chem

def extract_detailed_markush(file_path):
    if not os.path.exists(file_path):
        return "파일을 찾을 수 없습니다."

    with open(file_path, 'r', encoding='utf-8') as f:
        text = f.read().strip()

    # --- [추가] Core 구조 추출 및 분석 ---
    # 첫 번째 | 기호 이전이 Core SMILES입니다.
    core_raw = text.split('|')[0].strip()
    # 전체 텍스트에서 Core에 대한 메타데이터(|...|) 부분만 추출
    core_meta_match = re.search(r"^.*?\|(.*?)\|", text)
    
    core_info = {"smiles": core_raw, "labeled_smiles": "변환 실패", "mapping": []}
    
    if core_meta_match:
        # Core의 메타데이터를 포함한 조각을 만들어 RDKit으로 분석
        core_chunk = f"{core_raw} |{core_meta_match.group(1)}|"
        mol_core = Chem.MolFromSmiles(core_chunk)
        if mol_core:
            mapping = []
            for atom in mol_core.GetAtoms():
                if atom.GetSymbol() == '*':
                    if atom.HasProp('_label'):
                        l = atom.GetProp('_label')
                        atom.SetIsotope(int(l.replace('R', '')))
                        mapping.append(f"* -> {l}")
            core_info["labeled_smiles"] = Chem.MolToSmiles(mol_core, isomericSmiles=True)
            core_info["mapping"] = mapping

    # 1. RG: 섹션 추출
    rg_match = re.search(r"RG:(.*)(?=\|atomProp|\||$)", text)
    if not rg_match: return {"Core": core_info, "R_Groups": "RG 데이터 없음"}
    rg_content = rg_match.group(1)

    # 2. R-그룹별 분리
    r_groups = {}
    parts = re.split(r'(_R\d+=)', rg_content)
    
    for i in range(1, len(parts), 2):
        label = parts[i].replace('=', '').strip('_')
        raw_chunks = re.findall(r'\{(.*?)\}', parts[i+1])
        
        candidates = []
        for chunk in raw_chunks:
            mol = Chem.MolFromSmiles(chunk)
            labeled_smi = "변환 실패"
            mapping_info = []

            if mol:
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() == '*':
                        if atom.HasProp('_label'):
                            l = atom.GetProp('_label')
                            atom.SetIsotope(int(l.replace('R', '')))
                            mapping_info.append(f"* -> {l}")
                        for prop in atom.GetPropNames():
                            if 'AP' in prop:
                                atom.SetIsotope(0)
                                mapping_info.append(f"* -> {prop}")
                labeled_smi = Chem.MolToSmiles(mol, isomericSmiles=True)

            candidates.append({
                "original_chunk": f"{{{chunk}}}",
                "labeled_smiles": labeled_smi,
                "connection_map": mapping_info
            })
        r_groups[label] = candidates

    return {"Core": core_info, "R_Groups": r_groups}

# --- 실행 및 결과 출력 ---
file_path = '/Users/kogeon/python_projects_path/LNPDB_extend/PatientToLNPDB/Markush/WO2021021634_formula1.cxsmiles'
results = extract_detailed_markush(file_path)

if isinstance(results, dict):
    # Core 출력
    core = results["Core"]
    print("="*50)
    print(f"CORE STRUCTURE")
    print(f"  원본 SMILES: {core['smiles']}")
    print(f"  번호 매핑: {core['labeled_smiles']}")
    print(f"  연결 지점: {', '.join(core['mapping'])}")
    print("="*50)

    # R-Groups 출력
    for label, items in results["R_Groups"].items():
        print(f"\n[{label}] 정의부")
        for idx, item in enumerate(items, 1):
            print(f"  후보 {idx}: {item['original_chunk']}")
            print(f"    └─ 번호 매핑: {item['labeled_smiles']}")
            print(f"    └─ 연결 상세: {', '.join(item['connection_map'])}")
else:
    print(results)