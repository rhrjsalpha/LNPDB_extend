import pandas as pd
import json
import os
import sys
import time
from google import genai
from google.genai import types

# 1. API 키 설정
API_KEY = "" # 여기에 본인의 키를 입력하세요
client = genai.Client(api_key=API_KEY)

def extract_lipid_info(claim_text):
    # 2. 프롬프트 수정: Cargo 카테고리 추가
    prompt = f"""
    당신은 특허 분석 전문가입니다. 아래의 특허 청구항 텍스트를 분석하여 정보를 추출하세요.
    
    분류 기준:
    1. Ionizable_Lipid: 이온화 가능 지질, 양이온성 지질
    2. PEG_Lipid: PEG 지질, PEG-지질 컨쥬게이트
    3. Helper_Lipid: 인지질(Phospholipid), 중성 지질
    4. Structural_Lipid: 콜레스테롤 및 기타 구조적 지질
    5. Cargo: 내부에 적재되는 치료/예방 제제 (예: mRNA, siRNA, Corticosteroid, 단백질 등)

    각 카테고리에 대해 다음 정보를 추출하여 JSON 형식으로 응답하세요:
    - names: 문서에 언급된 모든 구체적인 이름이나 화합물 리스트
    - details: 몰 백분율(mol%), 중량비(wt/wt), 혹은 Cargo의 경우 구체적인 특징(예: 'chemically modified', 'logP > 3.0')

    결과 형식:
    {{
      "Ionizable_Lipid": {{ "names": [], "details": "" }},
      "PEG_Lipid": {{ "names": [], "details": "" }},
      "Helper_Lipid": {{ "names": [], "details": "" }},
      "Structural_Lipid": {{ "names": [], "details": "" }},
      "Cargo": {{ "names": [], "details": "" }}
    }}

    청구항 텍스트:
    {claim_text}
    """

    # 3. 모델 호출 및 429 에러 대응
    try:
        response = client.models.generate_content(
            model='gemini-2.0-flash', 
            contents=prompt,
            config=types.GenerateContentConfig(
                response_mime_type='application/json'
            )
        )
        return json.loads(response.text)
    except Exception as e:
        if "429" in str(e):
            print("한도 초과(429) 발생. 10초 후 재시도합니다...")
            time.sleep(10)
            return extract_lipid_info(claim_text)
        print(f"에러 발생: {e}")
        return None

def convert_to_dataframe(extracted_data):
    if not extracted_data:
        return pd.DataFrame()
        
    rows = []
    for category, info in extracted_data.items():
        rows.append({
            "Category": category,
            "Names": ", ".join(info['names']) if info['names'] else "N/A",
            "Details": info['details'] if info['details'] else "N/A"
        })
    return pd.DataFrame(rows)

# --- 실행부 ---
patent_claims = """
What is claimed is:
1. A lipid nanoparticle comprising:
(a) a cationic lipid;
(b) a non-cationic lipid;
(c) a corticosteroid and;
(d) a nucleic acid, wherein the nucleic acid and the corticosteroid are encapsulated within the lipid nanoparticle.
2. A population of lipid nanoparticles comprising the lipid nanoparticle of claim 1.
3. A population of lipid nanoparticles, comprising at least one population of lipid nanoparticles selected from:
(a) a first population of lipid nanoparticles that each comprise a cationic lipid, a non-cationic lipid, and a corticosteroid; and
(b) a second population of lipid nanoparticles that each comprise a cationic lipid, a non-cationic lipid, and a nucleic acid,
wherein the first population of lipid nanoparticles does not comprise a nucleic acid, and wherein the second population of lipid nanoparticles does not comprise a corticosteroid.
4. A population of lipid nanoparticles comprising the first and second populations of lipid nanoparticles of claim 3.
5. The lipid nanoparticle or population thereof of any one of claims 1-4, wherein the nucleic acid is mRNA.
6. The lipid nanoparticle or population thereof of any one of claims 1-5, wherein the corticosteroid has a logP greater than 3.0.
7. The lipid nanoparticle or population thereof of any one of claims 1-6, wherein the corticosteroid is a glucocorticoid.
8. The lipid nanoparticle or population thereof of any one of claims 1-6, wherein the corticosteroid is a mineralocorticoid.
9. The lipid nanoparticle or population thereof of any one of claims 1-6, wherein the corticosteroid is a clobetasol.
10. The lipid nanoparticle or population thereof of any one of claims 1-9, wherein the non-cationic lipid is selected from a PEG-lipid conjugate and a phospholipid.
11. The lipid nanoparticle or population thereof of any one of claims 1-10, wherein the lipid nanoparticle further comprises cholesterol.
12. The lipid nanoparticle or population thereof of any one of claims 10 and 11, wherein the phospholipid comprises dipalmitoylphosphatidylcholine (DPPC), distearoylphosphatidylcholine (DSPC), or a mixture thereof.
13. The lipid nanoparticle or population thereof of any one of claims 10 and 11, wherein the PEG-lipid conjugate is selected from the group consisting of a PEG-diacylglycerol (PEG-DAG) conjugate, a PEG-dialkyloxypropyl (PEG-DAA) conjugate, a PEG-phospholipid conjugate, a PEG-ceramide (PEG-Cer) conjugate, and a mixture thereof.
14. The lipid nanoparticle or population thereof of claim 13, wherein the PEG-lipid conjugate is a PEG-DAA conjugate.
15. The lipid nanoparticle or population thereof of claim 14, wherein the PEG-DAA conjugate is selected from the group consisting of a PEG-didecyloxypropyl (C10) conjugate, a PEG-dilauryloxypropyl (C12) conjugate, a PEG-dimyristyloxypropyl (C14) conjugate, a PEG-dipalmityloxypropyl (C16) conjugate, a PEG-distearyloxypropyl (C18) conjugate, and a mixture thereof.
16. The lipid nanoparticle or population thereof of any one of claims 1-15, wherein the lipid nanoparticle has a lipid:nucleic mass ratio of from about 9:1 to about 20:1.
17. The lipid nanoparticle or population thereof of any one of claims 1-16, wherein the mRNA is chemically modified.
18. The lipid nanoparticle or population thereof of any one of claims 1-17, wherein the lipid nanoparticle comprises an electron dense core.
19. The lipid nanoparticle or population thereof of any one of claims 1-17, wherein the lipid nanoparticle comprises an electron dense core and wherein the mRNA is located within the electron dense core.
20. A population of lipid particles comprising a multiplicity of lipid nanoparticles of claim 1.
21. The lipid nanoparticle or population thereof of any one of claims 1-18, which has an IFIT response that is no more than 30 fold greater than a reference IFIT response of phosphate buffered saline.
22. The lipid nanoparticle or population thereof of any one of claims 1-21, which has a PEG-lipid conjugate present in an amount of at least 3 mole percent; and mRNA encapsulated within the lipid particle; provided that the lipid particle comprises less than 0.5 mole percent phospholipid.
23. A pharmaceutical composition comprising the lipid nanoparticle or population thereof of any one of claims 1-22, and a pharmaceutically acceptable carrier.
24. A method for introducing an mRNA that encodes a protein into a cell, the method comprising contacting the cell with the lipid nanoparticle or population thereof of any one of claims 5-23, under conditions whereby the mRNA is introduced into the cell and expressed therein to produce the protein.
25. A method for treating and/or ameliorating one or more symptoms associated with a disease in a human, caused by impaired expression of a protein in the human, the method comprising administering to the human a therapeutically effective amount of the lipid nanoparticle or population thereof of any one of claims 5-23, wherein the mRNA encapsulated within the lipid nanoparticle encodes the protein.
26. A lipid nanoparticle formulation comprising a multiplicity of lipid nanoparticles, wherein each lipid nanoparticle comprises:
(a) a cationic lipid;
(b) a non-cationic lipid; and
(c) mRNA encapsulated within the lipid particle,
wherein the lipid nanoparticle formulation has an IFIT response that is no more than 30 fold greater than a reference IFIT response of phosphate buffered saline.
27. A method of making a lipid nanoparticle, comprising combining:
(a) a cationic lipid;
(b) a non-cationic lipid; and
(c) purified mRNA so as to form a lipid nanoparticle, wherein the mRNA is encapsulated within the lipid nanoparticle, and wherein the lipid nanoparticle has an IFIT response that is no more than 30 fold greater than a reference IFIT response of phosphate buffered saline.
28. A method of making a lipid nanoparticle formulation comprising a multiplicity of lipid nanoparticles, the method comprising the step of combining:
(a) a cationic lipid;
(b) a non-cationic lipid; and
(c) purified mRNA so as to form a lipid nanoparticle formulation comprising a multiplicity of lipid nanoparticles, wherein the mRNA is encapsulated within the lipid particles in the lipid nanoparticle formulation, and wherein the lipid nanoparticle formulation has an IFIT response that is no more than 30 fold greater than a reference IFIT response of phosphate buffered saline.
29. A lipid nanoparticle formulation comprising a multiplicity of lipid nanoparticles made by a process comprising the steps of combining:
(a) a cationic lipid;
(b) a non-cationic lipid; and
(c) purified mRNA so as to form a lipid nanoparticle formulation comprising a multiplicity of lipid nanoparticles, wherein the mRNA is encapsulated within the lipid particles in the lipid nanoparticle formulation, and wherein the lipid nanoparticle formulation has an IFIT response that is no more than 30 fold greater than a reference IFIT response of phosphate buffered saline.
30. A lipid nanoparticle comprising:
(a) a cationic lipid;
(b) a PEG-lipid conjugate present in an amount of at least 3 mole percent; and
(c) mRNA encapsulated within the lipid particle;
provided that the lipid particle comprises less than 0.5 mole percent phospholipid.
31. A population of lipid nanoparticles wherein each lipid nanoparticle in the population comprises:
(a) a cationic lipid;
(b) a PEG-lipid conjugate present in an amount of at least 3 mole percent; and
(c) mRNA encapsulated within the lipid nanoparticle;
provided that the lipid nanoparticle comprises less than 0.5 mole percent phospholipid.
"""

print("데이터 추출 중 (Cargo 정보 포함)...")
result = extract_lipid_info(patent_claims)
default_name = "US20190240354A1"

if result:
    df = convert_to_dataframe(result)
    print("\n--- 추출 결과 표 ---")
    print(df)

    # 폴더 생성 및 저장
    base_dir = os.path.dirname(os.path.abspath(__file__))
    extract_folder = os.path.join(base_dir, "extracted")
    os.makedirs(extract_folder, exist_ok=True)

    export_name = sys.argv[1] if len(sys.argv) > 1 else default_name
    if not export_name.lower().endswith('.csv'):
        export_name = export_name + '.csv'

    csv_path = os.path.join(extract_folder, export_name)
    df.to_csv(csv_path, index=False, encoding='utf-8-sig') # 한글 깨짐 방지
    print(f"\n성공적으로 저장되었습니다: {csv_path}")