import os
import pandas as pd
import json
import glob
import time

# IDE의 빨간 줄(Import 에러)을 방지하기 위한 방식
try:
    from google import generativeai as genai
except ImportError:
    print("❌ 패키지가 없습니다. 터미널에서 'pip install -U google-generativeai'를 실행하세요.")

# ==========================================
# 1. 설정 (사용자 수정 필요)
# ==========================================
GENAI_API_KEY = "AIzaSyAWzeQ1b4Mrb4XfpJJG8LwtUNqkTBF_2GM"  # 발급받은 API 키 입력
genai.configure(api_key=GENAI_API_KEY)

# 모델 설정 (긴 텍스트 분석에 최적화된 Pro 모델 권장)
model = genai.GenerativeModel('gemini-3-pro-preview')

INPUT_FOLDER = "EXTRACTED_TEXT"  # Selenium으로 뽑은 텍스트 폴더
METADATA_CSV = "20250409_LipidLibrarySample.xlsx"  # 기존 메타데이터
OUTPUT_FILE = "LNP_Experimental_Database_Final.xlsx"

# ==========================================
# 2. 고도화된 LNPDB 추출 프롬프트
# ==========================================
SYSTEM_PROMPT = """
당신은 지질 나노입자(LNP) 특허 문헌에서 실험 데이터를 구조화하여 추출하는 전문 분석가입니다.
제공된 텍스트에서 '실제 실험(Examples)' 및 '결과 데이터(Tables)'를 찾아 아래 스키마에 맞춰 추출하세요.

[핵심 규칙]
1. 로우 분리: 하나의 특허 내에서 이온화 지질(IL), 조성비, 실험 모델 중 하나라도 다르면 별개의 JSON 객체로 생성하세요.
2. 실제 데이터 우선: Claims가 아닌 실제 제조 및 테스트된 'Examples' 수치에 집중하세요.
3. ID 생성: '저자/출원인 이니셜_출원연도_순번' (예: GT_2021_1, GT_2021_2) 형식으로 Experiment_ID를 만드세요.
4. 데이터가 없는 항목은 "N/A" 또는 공백으로 두세요.

반드시 아래 필드를 포함한 JSON 리스트([...])로만 답변하세요:
{
  "Experiment_ID": "ID (예: GT_2021_1)",
  "IL_name": "이온화 지질 명칭",
  "IL_SMILES": "이온화 지질 전체 SMILES",
  "IL_head_name": "Amine head 명칭",
  "IL_head_SMILES": "Amine head SMILES",
  "IL_linker_name": "Linker 명칭",
  "IL_linker_SMILES": "Linker SMILES",
  "IL_tail1_name": "Tail1 명칭",
  "IL_tail1_SMILES": "Tail1 SMILES",
  "IL_tail2_name": "Tail2 명칭",
  "IL_tail2_SMILES": "Tail2 SMILES",
  "IL_molratio": "IL 몰비(0-100)",
  "IL_to_nucleicacid_massratio": "IL/핵산 질량비",
  "IL_to_nucleicacid_chargeratio": "N/P ratio",
  "HL_name": "Helper Lipid 명칭",
  "HL_SMILES": "Helper Lipid SMILES",
  "HL_molratio": "HL 몰비",
  "CHL_name": "Cholesterol 명칭",
  "CHL_SMILES": "Cholesterol SMILES",
  "CHL_molratio": "Cholesterol 몰비",
  "PEG_name": "PEG 지질 명칭",
  "PEG_SMILES": "PEG 지질 SMILES",
  "PEG_molratio": "PEG 몰비",
  "fifthcomponent_name": "제5성분 명칭",
  "fifthcomponent_SMILES": "제5성분 SMILES",
  "fifthcomponent_molratio": "제5성분 몰비",
  "Aqueous_buffer": "수상 버퍼 (acetate/citrate 등)",
  "Dialysis_buffer": "투석 버퍼 (PBS/None 등)",
  "Mixing_method": "혼합방법 (microfluidics 등)",
  "Model": "in_vitro 또는 in_vivo",
  "Model_type": "세포주 또는 동물종",
  "Model_target": "타겟 조직/장기",
  "Route_of_administration": "투여 경로",
  "Cargo": "탑재 단백질/유전자 이름",
  "Cargo_type": "핵산 종류 (mRNA/siRNA 등)",
  "Dose_ug_nucleicacid": "투여량(ug)",
  "Experiment_method": "측정 방법 (luminescence 등)",
  "Experiment_batching": "individual 또는 barcoded",
  "Experiment_value": "실험 결과 수치"
}
"""


# ==========================================
# 3. 분석 및 추출 실행 로직
# ==========================================
def run_lnpdb_extraction():
    # 메타데이터 로드
    try:
        meta_df = pd.read_excel(METADATA_CSV)
        print(f"📊 메타데이터 로드 성공: {len(meta_df)}건")
    except Exception as e:
        print(f"⚠️ 메타데이터 로드 실패: {e}")
        meta_df = pd.DataFrame()

    all_rows = []
    txt_files = glob.glob(os.path.join(INPUT_FOLDER, "*.txt"))

    if not txt_files:
        print(f"❌ '{INPUT_FOLDER}' 폴더에 텍스트 파일이 없습니다.")
        return

    for file_path in txt_files:
        patent_id = os.path.basename(file_path).replace(".txt", "")
        print(f"🔍 분석 중... {patent_id}")

        # 메타데이터 매칭 (Applicant, Link 등 보존)
        current_meta = {}
        if not meta_df.empty and patent_id in meta_df['PATENT_ID'].values:
            current_meta = meta_df[meta_df['PATENT_ID'] == patent_id].iloc[0].to_dict()
        else:
            current_meta = {"PATENT_ID": patent_id}

        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()

        try:
            # Gemini API 호출 (Native JSON 모드)
            response = model.generate_content(
                f"{SYSTEM_PROMPT}\n\n[특허 본문]\n{content[:30000]}",
                generation_config={"response_mime_type": "application/json", "temperature": 0.1}
            )

            extracted_data = json.loads(response.text)

            # 리스트 형태가 아니면 변환
            if isinstance(extracted_data, dict):
                extracted_data = [extracted_data]

            # 각 실험별로 로우 생성 및 메타데이터 결합
            for exp in extracted_data:
                combined_row = {**current_meta, **exp}
                all_rows.append(combined_row)

            print(f"   ✅ {len(extracted_data)}개의 실험 데이터 행을 생성했습니다.")

        except Exception as e:
            print(f"   ❌ {patent_id} 분석 중 오류 발생: {e}")

        # API 할당량 고려
        time.sleep(2)

    # 엑셀 저장
    if all_rows:
        df_final = pd.DataFrame(all_rows)
        # 컬럼 순서 정리 (필요 시 조절)
        df_final.to_excel(OUTPUT_FILE, index=False)
        print(f"\n✨ 모든 작업 완료! 파일 저장됨: {OUTPUT_FILE}")
    else:
        print("❌ 추출된 데이터가 없습니다.")


if __name__ == "__main__":
    run_lnpdb_extraction()