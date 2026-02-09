import pandas as pd
import re
import os

def extract_detailed_markush(file_path):
    # 1. 파일 존재 여부 확인 및 읽기
    if not os.path.exists(file_path):
        # 파일이 없으면 빈 결과 구조 반환 (에러 메시지 대신)
        return {"Raw_Full": f"Error: File not found at {file_path}", "Core": {"smiles": "", "mapping": []}, "R_Groups": {}}

    with open(file_path, 'r', encoding='utf-8') as f:
        full_text = f.read().strip()

    if not full_text:
        return {"Raw_Full": "Error: Empty file", "Core": {"smiles": "", "mapping": []}, "R_Groups": {}}

    # 2. Core와 Meta 섹션 분리
    main_parts = full_text.split('|')
    core_smiles = main_parts[0].strip()
    full_meta = main_parts[1] if len(main_parts) > 1 else ""

    # 3. Core 연결 정보 (세미콜론 위치 기반)
    core_mappings = []
    core_label_match = re.search(r"\$(.*?)\$", full_meta)
    if core_label_match:
        labels = core_label_match.group(1).split(';')
        for idx, lbl in enumerate(labels):
            if '_R' in lbl:
                clean_lbl = lbl.replace('$', '').strip()
                core_mappings.append(f"AtIdx {idx} -> {clean_lbl}")

    core_info = {"smiles": core_smiles, "mapping": core_mappings}

    # 4. R-Groups 섹션 (RG: 이후 데이터 추출)
    r_groups = {}
    # RG: 시작부터 다음 파이프(|)나 문자열 끝까지 추출
    rg_match = re.search(r"RG:(.*?)(?=\|atomProp|$)", full_text)
    if rg_match:
        rg_content = rg_match.group(1)
        # _R1=, _R2= 등으로 분할
        rg_blocks = re.split(r'(_R\d+=)', rg_content)
        
        for i in range(1, len(rg_blocks), 2):
            label = rg_blocks[i].replace('=', '').strip('_')
            # {} 내부 후보군 추출
            candidates_raw = re.findall(r'\{(.*?)\}', rg_blocks[i+1])
            
            candidate_list = []
            for chunk in candidates_raw:
                sub_parts = chunk.split('|')
                smi = sub_parts[0].strip()
                meta = sub_parts[1] if len(sub_parts) > 1 else ""
                
                # Connection 정보 추출 (AP 및 R-label)
                conn_tags = []
                # _AP1, _AP2 등 추출
                aps = re.findall(r"_AP\d+", meta)
                if aps: conn_tags.extend(aps)
                # 내부 매핑된 _R3, _R4 등 추출
                internal_rs = re.findall(r"_R\d+", meta)
                if internal_rs: conn_tags.extend(internal_rs)
                
                if not conn_tags and '*' in smi:
                    conn_tags.append("*")

                candidate_list.append({
                    "smiles": smi,
                    "map": ", ".join(sorted(list(set(conn_tags)))),
                    "raw": f"{{{chunk}}}"
                })
            r_groups[label] = candidate_list

    return {"Raw_Full": full_text, "Core": core_info, "R_Groups": r_groups}

def save_unified_markush(results, input_file_path):
    # results가 딕셔너리인지 체크 (안전장치)
    if not isinstance(results, dict):
        print(f"❌ 데이터 추출 실패: {results}")
        return

    cx_filename = os.path.splitext(os.path.basename(input_file_path))[0]
    script_dir = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(script_dir, f"{cx_filename}_Unified.xlsx")

    rows = []
    # 1. 원본 소스
    rows.append({
        "Category": "FULL_SOURCE", 
        "Label": "Original CXSMILES", 
        "SMILES": results.get("Raw_Full", ""), 
        "Connection": "", 
        "Original_Chunk": ""
    })
    
    # 2. Core (빈 줄 없이)
    core = results.get("Core", {})
    rows.append({
        "Category": "CORE", 
        "Label": "Main_Skeleton", 
        "SMILES": core.get('smiles', ""), 
        "Connection": ", ".join(core.get('mapping', [])), 
        "Original_Chunk": core.get('smiles', "")
    })

    # 3. R-Groups (빈 줄 없이)
    for label, candidates in results.get("R_Groups", {}).items():
        for idx, item in enumerate(candidates, 1):
            rows.append({
                "Category": "R_GROUP",
                "Label": f"{label} (Candidate {idx})",
                "SMILES": item['smiles'],
                "Connection": item['map'],
                "Original_Chunk": item['raw']
            })

    df = pd.DataFrame(rows)
    df.to_excel(save_path, index=False, sheet_name='Markush_Data')
    print(f"✅ 통합 파일 생성 완료: {save_path}")

# 실행
file_path = '/Users/kogeon/python_projects_path/LNPDB_extend/PatientToLNPDB/Markush/WO2021021634_formula1.cxsmiles'
data_results = extract_detailed_markush(file_path)
save_unified_markush(data_results, file_path)