import google.generativeai as genai
import json
import pandas as pd
import os

# 1. Gemini API 설정
genai.configure(api_key="YOUR_GOOGLE_API_KEY")
model = genai.GenerativeModel('gemini-1.5-pro')

def read_patent_text(file_path):
    """텍스트 파일을 읽어오는 함수"""
    with open(file_path, 'r', encoding='utf-8') as f:
        return f.read()

def step1_index_document(text):
    """
    1단계: 문서의 전체 구조를 파악하고 실험 데이터가 있는 범위를 리스트로 받음
    """
    prompt = f"""
    아래는 LNP 관련 특허/논문 텍스트입니다. 
    이 문서를 분석하여 [지질 합성, LNP 제조법, 실험 결과(Test Examples)] 섹션이 각각 어디에 해당하는지 파악하세요.
    
    결과는 반드시 아래의 JSON 리스트 형식으로만 답변하세요:
    [
      {{"section_name": "Example 1", "content_type": "Synthesis", "summary": "IL-1 합성"}},
      {{"section_name": "Test Example 2", "content_type": "Biological Evaluation", "summary": "세포 내 siRNA 전달 효율 측정"}}
    ]
    
    TEXT:
    {text[:8000]}  # 토큰 제한을 고려하여 초기 분석 범위 설정
    """
    
    response = model.generate_content(prompt)
    # 응답에서 JSON 부분만 추출 (Gemini가 마크다운을 섞을 수 있음)
    content = response.text.strip().replace('```json', '').replace('```', '')
    return json.loads(content)

def step2_extract_lnpdb_data(experimental_text):
    """
    2단계: 실제 실험이 언급된 부분에서 LNPDB 컬럼에 맞는 데이터 추출
    SMILES는 제외하고 이름과 수치 위주로 추출
    """
    columns = [
        "IL_name", "IL_molratio", "HL_name", "HL_molratio", "CHL_name", "CHL_molratio", 
        "PEG_name", "PEG_molratio", "IL_to_nucleicacid_massratio", "Aqueous_buffer",
        "Mixing_method", "Model", "Model_type", "Model_target", "Route_of_administration",
        "Cargo", "Cargo_type", "Dose_ug_nucleicacid", "Experiment_method", "Experiment_value"
    ]
    
    prompt = f"""
    당신은 전문 데이터 엔지니어입니다. 제공된 실험 텍스트에서 LNPDB 템플릿 형식에 맞게 데이터를 추출하세요.
    
    규칙:
    1. SMILES 컬럼은 모두 제외하거나 빈칸으로 둡니다.
    2. 지질 조성(Mol ratio)이 '35/16/46.5/2.5'와 같이 제공되면 각 성분에 맞게 분리하세요.
    3. 실험 결과 값(Experiment_value)은 수치 데이터 위주로 추출하세요.
    4. 각 실험구(Preparation 또는 Example)를 하나의 행(Row)으로 만드세요.
    
    대상 컬럼: {columns}
    반환 형식: {{"results": [{{column: value}}, ...]}} (JSON)
    
    TEXT:
    {experimental_text}
    """
    
    response = model.generate_content(prompt)
    content = response.text.strip().replace('```json', '').replace('```', '')
    return json.loads(content)

def main_pipeline(file_path):
    # 파일 읽기
    full_text = read_patent_text(file_path)
    
    # [Step 1] 문서 인덱싱
    print("문서 구조 분석 중...")
    sections = step1_index_document(full_text)
    print(f"분석된 섹션 수: {len(sections)}")
    
    # [Step 2] 실험 섹션 필터링 및 재질문
    # 'Biological Evaluation'이나 'Test'가 포함된 섹션을 위주로 텍스트를 다시 구성하여 질문
    # 여기서는 예시로 전체 텍스트 중 실험 파트라고 판단되는 부분을 재질문
    print("상세 실험 데이터 추출 중 (LNPDB 형식)...")
    
    # 실제로는 sections에서 찾은 키워드로 텍스트를 slice하여 질문하는 것이 효율적
    final_data = []
    # 텍스트가 너무 길면 나눠서 호출하는 루프가 필요함
    extraction_result = step2_extract_lnpdb_data(full_text) # 예시를 위해 전체 텍스트 사용
    final_data.extend(extraction_result['results'])
    
    # [Step 3] 결과 저장
    df = pd.DataFrame(final_data)
    output_file = "lnpdb_extracted_results.csv"
    df.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"추출 완료! 저장된 파일: {output_file}")
    return df

# 실행
if __name__ == "__main__":
     main_pipeline("patent_sample.txt")