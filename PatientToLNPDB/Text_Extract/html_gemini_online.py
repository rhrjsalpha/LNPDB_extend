import google.generativeai as genai
import json
import pandas as pd
import os
import re
import time
from pathlib import Path

# 1. API Configuration
# Please ensure your API key is correct.
genai.configure(api_key="")

def get_gemini_model(model_name):
    """Returns a GenerativeModel instance with JSON mode enabled."""
    return genai.GenerativeModel(
        model_name=model_name,
        generation_config={"response_mime_type": "application/json"}
    )

def read_patent_text(file_path):
    with open(file_path, 'r', encoding='utf-8') as f:
        return f.read()

# --- [Step 1] Document Indexing (Gemini 3 Pro) ---
def step1_index_sections(text, model_instance, cache_path):
    if os.path.exists(cache_path):
        print(f"-> Step 1 cache loaded: {cache_path}")
        with open(cache_path, 'r', encoding='utf-8') as f: return json.load(f)

    print("Step 1: Indexing document structure (Pro)...")
    prompt = """
    Analyze the following LNP patent and categorize it into 7 categories. 
    Find the [start paragraph number] and [end paragraph number] for each.
    Categories: 1.Background, 2.Markush, 3.Specific Compounds, 4.Definitions, 5.Manufacturing, 6.Synthesis, 7.Biological Examples
    Return Format: [{"category_no": 7, "category_name": "Biological Examples", "start_tag": "[1432]", "end_tag": "[1484]"}]
    Respond strictly in JSON format. Provide summaries in English.
    """
    # Using the provided text for analysis
    response = model_instance.generate_content(f"{prompt}\n\nTEXT: {text[:]}")
    data = json.loads(response.text)
    with open(cache_path, 'w', encoding='utf-8') as f: json.dump(data, f, indent=2, ensure_ascii=False)
    return data

# --- [Step 2] Data Point Identification (Gemini 3 Flash) ---
def step2_identify_data_points(sliced_text, model_instance, cache_path):
    if os.path.exists(cache_path):
        print(f"-> Step 2 cache loaded: {cache_path}")
        with open(cache_path, 'r', encoding='utf-8') as f: return json.load(f)

    print("Step 2: Identifying data points and tables (Flash)...")
    prompt = """
    Analyze the text to identify paragraph numbers ([tag_no]) containing Tables, Figures, or experimental data.
    
    Strict Rules for the "data" field:
    1. If a Table is only mentioned or described in the text (e.g., "Results are shown in Table 1A"), 
       label it as "Table [Name] description" (e.g., "Table 1A description").
    2. Only use the label "Table [Name]" (e.g., "Table 1A") if the actual data block/grid is present in that paragraph.
    3. Include other labels like "Experimental conditions", "Experimental results", or "Fig" where applicable.
    4. Provide the summary strictly in English.

    Return a list of objects with these keys:
    - "tag_no": The paragraph number in brackets like "[1437]".
    - "data": A list of identified components (e.g., ["Table 1A description", "Experimental conditions"]).
    - "summary_english": A brief English summary.

    Return Format: [{"tag_no": "[1436]", "data": ["Table 1A description", "Experimental conditions"], "summary_english": "Description of LNP parameters later shown in Table 1A."}]
    Respond strictly in JSON format.
    """
    response = model_instance.generate_content(f"{prompt}\n\nTEXT: {sliced_text}")
    data = json.loads(response.text)
    with open(cache_path, 'w', encoding='utf-8') as f: json.dump(data, f, indent=2, ensure_ascii=False)
    return data

# --- [Step 2.5] Extract Experimental Conditions (Metadata) in Batches ---
def step2_5_extract_conditions(relevant_text, identified_data_points, model_instance, cache_path):
    """
    Identifies experimental conditions for specific tables found in Step 2.
    Processes tables in batches to ensure exhaustive extraction.
    """
    if os.path.exists(cache_path):
        print(f"-> Step 2.5 cache loaded: {cache_path}")
        with open(cache_path, 'r', encoding='utf-8') as f: return json.load(f)

    # 1. Filter only actual Tables from Step 2 data
    target_tables = []
    for dp in identified_data_points:
        # Extract table names from the "data" list (e.g., "Table 1A")
        # We exclude entries labeled as "description"
        tables_in_tag = [d for d in dp.get('data', []) if "table" in d.lower() and "description" not in d.lower()]
        target_tables.extend(tables_in_tag)
    
    # Remove duplicates while preserving order
    target_tables = list(dict.fromkeys(target_tables))
    
    if not target_tables:
        print("!! No actual tables identified for metadata extraction.")
        return {"table_mapping": []}

    print(f"Step 2.5: Extracting metadata for {len(target_tables)} identified tables in batches...")
    
    lnpdb_cols = ["IL_name", "IL_molratio", "HL_name", "HL_molratio", "CHL_name", "CHL_molratio", "PEG_name", "PEG_molratio", "Aqueous_buffer", "Mixing_method", "Model_type", "Route_of_administration", "Cargo_type"]
    
    master_mapping = []
    batch_size = 5 # Process 5 tables per request to maintain accuracy

    for i in range(0, len(target_tables), batch_size):
        batch = target_tables[i:i + batch_size]
        print(f"  - Processing batch {i//batch_size + 1}: {batch}")
        
        prompt = f"""
        Analyze the provided text to extract 'experimental conditions' for the following SPECIFIC tables: {batch}
        
        Rules:
        1. Only extract data for the tables listed above.
        2. Do NOT include experimental results/values (e.g., skip flux, size, or concentration values).
        3. If a table represents multiple conditions (e.g., different dose levels or different administration routes described in the text), create a separate condition object for each row.
        4. Use the following columns: {lnpdb_cols}
        5. Respond strictly in English in the specified JSON format.
        
        Return Format:
        {{
          "table_mapping": [
            {{
              "table_name": "Table Name from list",
              "conditions": [
                {{ "IL_name": "...", "IL_molratio": "...", ... }},
                {{ "IL_name": "...", "IL_molratio": "...", ... }}
              ]
            }}
          ]
        }}
        """
        
        try:
            response = model_instance.generate_content(f"{prompt}\n\nTEXT: {relevant_text}")
            batch_data = json.loads(response.text)
            if "table_mapping" in batch_data:
                master_mapping.extend(batch_data["table_mapping"])
        except Exception as e:
            print(f"    !! Error processing batch {batch}: {e}")
        
        time.sleep(1) # Respect API rate limits

    final_data = {"table_mapping": master_mapping}
    
    # Save integrated result to cache
    with open(cache_path, 'w', encoding='utf-8') as f: 
        json.dump(final_data, f, indent=2, ensure_ascii=False)
        
    return final_data

# --- [Step 3 Helper] 개별 테이블 구조 추출 ---
def step3_get_raw_table(specific_text, table_name, model_instance):
    """
    단일 테이블 이름(예: 'Table 1A')을 지정하여 해당 표만 정밀하게 2D 리스트로 추출합니다.
    """
    print(f"  - Reconstructing structured data for: {table_name}...")
    prompt = f"""
    Find and extract the data ONLY for '{table_name}' from the provided text.
    Convert it into a structured 2D list (Rows/Columns) including headers as the first row.
    
    Rules:
    1. Focus strictly on {table_name}. Ignore other tables present in the text.
    2. Ensure the output is a valid JSON object.
    3. If multiple instances of the same table name exist in the chunk, merge them into one continuous structure.
    
    Return Format: {{ "raw_data": [["Header1", "Header2"], ["Val1", "Val2"]] }}
    """
    response = model_instance.generate_content(f"{prompt}\n\nTEXT: {specific_text}")
    return json.loads(response.text)

def main_pipeline(file_path):
    file_stem = Path(file_path).stem
    output_dir = Path(f"./output_{file_stem}")
    tables_dir = output_dir / "tables"
    output_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    # Use Gemini 3 model names as per provided image
    m3_pro = get_gemini_model('gemini-3-pro-preview')
    m3_flash = get_gemini_model('gemini-3-flash-preview')

    full_text = read_patent_text(file_path)
    
    # Step 1
    sections = step1_index_sections(full_text, m3_pro, output_dir / "step1_sections.json")
    
    # Slice Manufacturing (5) and Examples (7)
    relevant_slices = ""
    for sec in sections:
        start_tag = sec.get('start_tag')
        end_tag = sec.get('end_tag')
        if start_tag and end_tag and sec.get('category_no') in [5, 7]:
            s_idx = full_text.find(str(start_tag))
            e_idx = full_text.find(str(end_tag))
            if s_idx != -1 and e_idx != -1:
                relevant_slices += full_text[s_idx : e_idx + 200] + "\n"

    # Step 2 & 2.5
    data_points = step2_identify_data_points(relevant_slices, m3_flash, output_dir / "step2_data_points.json")
    # main_pipeline 함수 내부 수정 후 (data_points 인수 추가)
    condition_master = step2_5_extract_conditions(
    relevant_slices, 
    data_points,      # 이 부분이 누락되었습니다.
    m3_flash, 
    output_dir / "step2_5_conditions.json"
    )

    # --- [main_pipeline 내 Step 3 실행 구간] ---
# data_points와 condition_master가 준비된 상태에서 실행됩니다.

    print("Step 3: Creating individual Excel files (Naming: TableName_ParagraphNumber)...")
    for dp in data_points:
       data_list = dp.get('data', [])
       
       # 실제 데이터 블록이 있는 테이블 명칭만 필터링 (description 제외)
       table_names = [str(d) for d in data_list if "table" in str(d).lower() and "description" not in str(d).lower()]
       
       if not table_names:
           continue
       
       tag = dp['tag_no'] # 예: "[1437]"
       clean_tag_num = re.sub(r'\D', '', tag) # 숫자만 추출: "1437"
       t_idx = full_text.find(tag)
       if t_idx == -1: continue

       # 현재 단락부터 다음 단락 직전까지 텍스트 슬라이싱
       match = re.search(r'(\d+)', tag)
       if match:
           curr_num = int(match.group(1))
           next_tag = tag.replace(match.group(1), str(curr_num + 1).zfill(len(match.group(1))))
           next_idx = full_text.find(next_tag, t_idx)
           chunk = full_text[t_idx : next_idx] if next_idx != -1 else full_text[t_idx : t_idx+8000]
       else: 
           chunk = full_text[t_idx : t_idx+8000]

       # 식별된 각 테이블에 대해 개별 추출 및 저장 루프
       for specific_table in table_names:
           try:
               # 1. Gemini를 통한 표 구조 복원
               table_data = step3_get_raw_table(chunk, specific_table, m3_flash)
               
               # 2. 파일명 설정: 테이블명_단락번호.xlsx
               # 특수문자 제거 및 공백을 언더바로 변환
               safe_table_name = re.sub(r'\W+', '_', specific_table).strip('_')
               excel_filename = f"{safe_table_name}_{clean_tag_num}.xlsx"
               excel_path = tables_dir / excel_filename
               
               # 3. 엑셀 파일 생성 (시트 2개 구성)
               with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
                   # 시트 1: Gemini가 복원한 원본 형태의 표
                   pd.DataFrame(table_data['raw_data']).to_excel(writer, sheet_name='Raw_Table', index=False, header=False)
                   
                   # 시트 2: Step 2.5에서 미리 추출된 해당 테이블의 실험 조건 매핑
                   matched_conditions = []
                   for entry in condition_master.get('table_mapping', []):
                       # 테이블명이 부분 일치하거나 포함되는지 확인
                       if specific_table.lower() in entry['table_name'].lower():
                           matched_conditions.extend(entry['conditions'])
                   
                   if matched_conditions:
                       pd.DataFrame(matched_conditions).to_excel(writer, sheet_name='LNPDB_Conditions', index=False)
               
               print(f"    => {excel_path.name} saved successfully.")
               time.sleep(0.5) # API 속도 제한 고려

           except Exception as e:
               print(f"    !! Extraction failed for {specific_table} in {tag}: {e}")

if __name__ == "__main__":
    target_path = "/Users/kogeon/python_projects_path/LNPDB_extend/PatientToLNPDB/Text_Extract/US20170210697A1.txt"
    main_pipeline(target_path)