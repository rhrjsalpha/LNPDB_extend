import pandas as pd
import os
import re
import json
import time
from pathlib import Path
import google.generativeai as genai

# API 설정 (실제 키 입력 필요)
genai.configure(api_key="")

def get_gemini_model(model_name='gemini-2.5-flash-lite'):
    return genai.GenerativeModel(
        model_name=model_name,
        generation_config={"response_mime_type": "application/json"}
    )

def split_by_paragraphs(text):
    """[숫자] 패턴을 기준으로 텍스트를 분할하여 리스트로 반환"""
    # [1], [12], [1234] 등의 패턴을 찾아 분할
    pattern = r'(\[\d+\])'
    parts = re.split(pattern, text)
    
    paragraphs = []
    for i in range(1, len(parts), 2):
        tag = parts[i]
        content = parts[i+1] if i+1 < len(parts) else ""
        paragraphs.append({"tag": tag, "text": content.strip()})
    return paragraphs

def extract_compound_info(paragraph_list, model_instance):
    """단락 리스트를 배치 단위로 Gemini에게 전달하여 화합물 정보 추출"""
    extracted_data = []
    batch_size = 10  # 속도와 정확도를 위해 10단락씩 묶어서 처리
    
    prompt_template = """
    Analyze the following patent paragraphs and extract information about chemical compounds.
    For each compound found, provide:
    1. compound_name (e.g., "Compound 1", "DLin-MC3-DMA")
    2. iupac_name
    3. inchi
    4. smiles
    5. has_structure_fig (True/False - based on mention of 'Figure', 'Fig', 'Structure of')
    6. other_info (molecular weight, description, etc.)

    Return a list of objects. If no compound is found in a paragraph, return an empty list for that tag.
    Strictly follow this JSON format:
    [
      {
        "tag": "[123]",
        "compounds": [
          {"compound_name": "...", "iupac_name": "...", "inchi": "...", "smiles": "...", "has_structure_fig": true, "other_info": "..."}
        ]
      }
    ]
    """

    for i in range(0, len(paragraph_list), batch_size):
        batch = paragraph_list[i:i+batch_size]
        batch_text = "\n\n".join([f"TAG: {p['tag']}\nTEXT: {p['text']}" for p in batch])
        
        print(f"Processing paragraphs {batch[0]['tag']} to {batch[-1]['tag']}...")
        
        try:
            response = model_instance.generate_content(f"{prompt_template}\n\nDATA:\n{batch_text}")
            batch_result = json.loads(response.text)
            extracted_data.extend(batch_result)
        except Exception as e:
            print(f"Error in batch {i}: {e}")
        
        time.sleep(0.5)  # Rate limit 고려

    return extracted_data

def main(file_path):
    file_stem = Path(file_path).stem
    output_dir = Path(f"./compound_extraction_{file_stem}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. 텍스트 읽기 및 분할
    with open(file_path, 'r', encoding='utf-8') as f:
        full_text = f.read()
    
    print("Splitting text by paragraphs...")
    paragraphs = split_by_paragraphs(full_text)
    
    # 2. 정보 추출 (Flash 모델 사용)
    model = get_gemini_model()
    results = extract_compound_info(paragraphs, model)
    
    # 3. 데이터 평탄화 (DataFrame 생성용)
    flat_data = []
    for entry in results:
        tag = entry.get('tag')
        for comp in entry.get('compounds', []):
            comp['paragraph_no'] = tag
            flat_data.append(comp)
    
    # 4. 결과 저장
    df = pd.DataFrame(flat_data)
    # 컬럼 순서 조정
    cols = ['paragraph_no', 'compound_name', 'iupac_name', 'inchi', 'smiles', 'has_structure_fig', 'other_info']
    df = df[cols]
    
    csv_path = output_dir / f"{file_stem}_compounds.csv"
    df.to_csv(csv_path, index=False, encoding='utf-8-sig')
    print(f"Extraction complete. Saved to: {csv_path}")

if __name__ == "__main__":
    target_path = "/Users/kogeon/python_projects_path/LNPDB_extend/PatientToLNPDB/Text_Extract/US20170210697A1.txt"
    main(target_path)