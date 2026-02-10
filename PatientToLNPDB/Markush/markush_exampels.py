import pandas as pd
import re
import os

def unified_markush_regex_multi_ap(input_path, output_path):
    if not os.path.exists(input_path):
        return print(f"오류: {input_path} 파일을 찾을 수 없습니다.")

    with open(input_path, 'r', encoding='utf-8') as f:
        text = f.read().strip()

    all_data = []

    # --- 1. ,RG: 패턴을 기준으로 CORE와 R-Groups 분리 ---
    if ",RG:" in text:
        core_part, rg_part = text.split(",RG:", 1)
    else:
        core_part = text
        rg_part = ""

    # --- 2. CORE 구조 처리 (원본 유지 및 바 보정) ---
    core_part = core_part.strip()
    if "|" in core_part and not core_part.endswith("|"):
        core_part = core_part + "|"
    
    all_data.append({
        "Category": "CORE",
        "Label": "Main_Skeleton_Original",
        "CXSMILES": core_part
    })

    # --- 3. R-Groups 섹션 처리 (다중 AP 라벨 변환) ---
    if rg_part:
        r_parts = re.split(r'(_R\d+=)', rg_part)
        for i in range(1, len(r_parts), 2):
            r_label = r_parts[i].replace('=', '').strip('_') # 예: R5
            
            candidates = re.findall(r'\{(.*?)\}', r_parts[i+1])
            for chunk in candidates:
                clean_chunk = chunk.strip()
                
                # [핵심 수정] 모든 _APn 패턴을 찾아 _{r_label}_n 형태로 교체
                # 예: R5 후보의 _AP1 -> _R5_1, _AP2 -> _R5_2
                # 이렇게 하면 다중 연결 지점도 고유한 ID를 갖게 되어 조립 시 정확히 매칭됩니다.
                def ap_replacer(match):
                    ap_num = match.group(1)
                    return f"_{r_label}_{ap_num}"
                
                clean_chunk = re.sub(r'_AP(\d+)', ap_replacer, clean_chunk)

                all_data.append({
                    "Category": "R_GROUP",
                    "Label": f"{r_label}(Candidate)",
                    "CXSMILES": clean_chunk
                })

    # --- 4. 엑셀 파일 저장 ---
    df = pd.DataFrame(all_data)
    df.to_excel(output_path, index=False)
    print(f"✅ 다중 라벨 매칭 완료: 총 {len(df)}행 저장")

# --- 실행 ---
input_file = '/Users/kogeon/python_projects_path/LNPDB_extend/PatientToLNPDB/Markush/WO2021021634_formula1.cxsmiles'
output_file = '/Users/kogeon/python_projects_path/LNPDB_extend/PatientToLNPDB/Markush/WO2021021634_Unified_MultiAP.xlsx'

unified_markush_regex_multi_ap(input_file, output_file)