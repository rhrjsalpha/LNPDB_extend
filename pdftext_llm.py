import ollama
import json
import re
import fitz  # PyMuPDF
from typing import List, Dict, Any

# ============================================================
# 유틸 함수
# ============================================================

def remove_think_block(text: str) -> str:
    """<think>...</think> 블록 제거"""
    return re.sub(
        r"<think>.*?</think>",
        "",
        text,
        flags=re.DOTALL | re.IGNORECASE
    ).strip()


def extract_json_only(text: str) -> str:
    """
    LLM 출력에서 ```json ... ``` 또는 ``` ... ``` 제거 후
    JSON 본문만 반환
    """
    text = remove_think_block(text)

    # ```json ... ``` 형태
    m = re.search(
        r"```(?:json)?\s*(\{.*?\})\s*```",
        text,
        re.DOTALL | re.IGNORECASE
    )
    if m:
        return m.group(1).strip()

    # fallback: 첫 { ~ 마지막 }
    m = re.search(r"(\{.*\})", text, re.DOTALL)
    if m:
        return m.group(1).strip()

    raise ValueError("JSON 객체를 찾지 못함")


def iter_pdf_pages(pdf_path: str):
    """PDF를 페이지 단위로 순회"""
    doc = fitz.open(pdf_path)
    for page_idx, page in enumerate(doc, start=1):
        text = page.get_text("text").strip()
        if text:
            yield page_idx, text
    doc.close()

# ============================================================
# 설정
# ============================================================

PDF_FILE = "MP_2025.pdf"   # ← 필요시 MP_2025_SI.pdf 로 변경
THINK_MODEL = "jinbora/deepseek-r1-Bllossom:8b"
JSON_MODEL  = "jinbora/deepseek-r1-Bllossom:8b"

# ============================================================
# 프롬프트
# ============================================================

THINK_PROMPT = """
다음 TEXT를 읽고
LNP (Lipid Nanoparticle) 조성 정보를 중심으로 정리하라.

출력 규칙 (중요):
- TEXT에서 LNP 조성 정보가 하나도 없으면
  반드시 아래 한 줄만 정확히 출력하라.

없음

- 위 경우를 제외하고는 "없음"이라는 단어를 절대 사용하지 마라.
- 불확실하거나 추론이 필요한 경우에도 추론하지 말고 "없음"으로 처리하라.
- 원문에 없는 정보는 절대 추가하지 마라.

정리 기준 (LNP가 존재하는 경우에만 적용):
- LNP 이름 또는 formulation 구분
- 핵산 종류 (mRNA, siRNA, pDNA 등)
- IL to nucleic acid mass ratio
- lipid 이름과 molar ratio
- "respectively", "for B10 formulation" 등의 문맥을 명확히 풀어서 작성

출력 형식 (예시):

LNP1:
- mRNA
- IL 10 : nucleic acid 1
- lipid name 1 : ratio
- lipid name 2 : ratio
- lipid name 3 : ratio
- lipid name 4 : ratio

LNP2:
- sgRNA
- IL 10 : nucleic acid 1
- lipid name 1 : ratio
- lipid name 2 : ratio
- lipid name 3 : ratio
- lipid name 4 : ratio

TEXT:
"""


JSON_PROMPT = """
너는 정보 변환기이다.

아래 SUMMARY에서 LNP 조성 정보를 추출하라.

출력 규칙:
- 반드시 JSON만 출력하라
- 설명, 문장, 주석 출력 금지
- 값이 없으면 null로 채워라
- 여러 LNP가 있으면 배열로 출력하라

JSON 스키마:
{
  "LNPs": [
    {
      "LNP_name": string | null,
      "Nucleic_acid": "mRNA" | "siRNA" | "pDNA" | null,
      "IL_to_Nucleic_Acid_Mass_Ratio": string | null,
      "Lipids": {
        "Ionizable_lipid": string | null,
        "Helper_lipid": string | null,
        "Cholesterol": string | null,
        "PEG_lipid": string | null
      },
      "Molar_ratios": {
        "Ionizable_lipid": number | string | null,
        "Helper_lipid": number | string | null,
        "Cholesterol": number | string | null,
        "PEG_lipid": number | string | null
      }
    }
  ]
}

SUMMARY:
"""

# ============================================================
# 메인 파이프라인 (페이지 단위)
# ============================================================

all_lnps: List[Dict[str, Any]] = []

START_PAGE = 12   # 시작 페이지 (1-indexed)
END_PAGE   = 16  # 끝 페이지 (포함)
for page_idx, page_text in iter_pdf_pages(PDF_FILE):
    if page_idx < START_PAGE:
        continue
    if page_idx > END_PAGE:
        break
    print(f"\n================ PAGE {page_idx} =================")

    try:
        # ---------- STEP 1: THINK ----------
        think_response = ollama.generate(
            model=THINK_MODEL,
            prompt=THINK_PROMPT + page_text,
            options={"temperature": 0}
        )

        summary_text = remove_think_block(think_response["response"])
        print(summary_text)

        if not summary_text.strip():
            print("⚠️ 요약 결과 없음 → 페이지 스킵")
            continue

        # ---------- STEP 2: JSON ----------
        json_response = ollama.generate(
            model=JSON_MODEL,
            prompt=JSON_PROMPT + summary_text,
            options={"temperature": 0}
        )

        json_text = extract_json_only(json_response["response"])

        data = json.loads(json_text)

        # ---------- 결과 누적 ----------
        lnps = data.get("LNPs", [])
        if lnps:
            for lnp in lnps:
                lnp["_source_page"] = page_idx
            all_lnps.extend(lnps)
            print(f"✅ LNP {len(lnps)}개 추출")
        else:
            print("ℹ️ LNP 없음")

    except Exception as e:
        print(f"❌ PAGE {page_idx} 실패: {e}")

# ============================================================
# 최종 결과 출력
# ============================================================

print("\n================ FINAL RESULT ================\n")
print(json.dumps(all_lnps, indent=2, ensure_ascii=False))

# 필요 시 파일 저장
with open("lnp_extracted_from_pdf.json", "w", encoding="utf-8") as f:
    json.dump(all_lnps, f, indent=2, ensure_ascii=False)

print(f"\n💾 저장 완료: lnp_extracted_from_pdf.json (총 {len(all_lnps)} LNP)")

# ============================================================
# CSV 변환
# ============================================================

import csv

csv_rows = []

for lnp in all_lnps:
    row = {
        "source_page": lnp.get("_source_page"),
        "LNP_name": lnp.get("LNP_name"),
        "Nucleic_Acid": lnp.get("Nucleic_acid"),
        "IL_to_Nucleic_Acid_Mass_Ratio": lnp.get("IL_to_Nucleic_Acid_Mass_Ratio"),

        "Ionizable_lipid": lnp.get("Lipids", {}).get("Ionizable_lipid"),
        "Helper_lipid": lnp.get("Lipids", {}).get("Helper_lipid"),
        "Cholesterol": lnp.get("Lipids", {}).get("Cholesterol"),
        "PEG_lipid": lnp.get("Lipids", {}).get("PEG_lipid"),

        "Ionizable_lipid_ratio": lnp.get("Molar_ratios", {}).get("Ionizable_lipid"),
        "Helper_lipid_ratio": lnp.get("Molar_ratios", {}).get("Helper_lipid"),
        "Cholesterol_ratio": lnp.get("Molar_ratios", {}).get("Cholesterol"),
        "PEG_lipid_ratio": lnp.get("Molar_ratios", {}).get("PEG_lipid"),
    }
    csv_rows.append(row)

csv_file = "lnp_extracted_from_pdf.csv"

with open(csv_file, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=csv_rows[0].keys())
    writer.writeheader()
    writer.writerows(csv_rows)

print(f"\n📄 CSV 저장 완료: {csv_file} (총 {len(csv_rows)} rows)")