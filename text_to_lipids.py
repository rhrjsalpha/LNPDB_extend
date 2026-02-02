import ollama
import pdfplumber

MODEL = "jinbora/deepseek-r1-Bllossom:8b"
PDF_PATH = "MP_2025.pdf"   # 분석할 PDF

PROMPT_PREFIX = """
다음 텍스트에
lipid nanoparticle (LNP)
lipid 이름 2개 이상
molar ratio, weight ratio, % 조성 등 비율 정보
이 정보들이 명시적으로 포함되어 있으면 True,
그렇지 않으면 False만 출력하라.
또한 텍스트중 포함되어 있는 부분만 추출하여 출력하라
"""

def page_has_lnp_composition(text: str) -> bool:
    if not text or len(text.strip()) < 50:
        return False

    response = ollama.generate(
        model=MODEL,
        prompt=PROMPT_PREFIX + "\n\n" + text,
        options={"temperature": 0}
    )

    answer = response["response"].strip()
    return answer


results = []

with pdfplumber.open(PDF_PATH) as pdf:
    for i, page in enumerate(pdf.pages):
        text = page.extract_text() or ""
        has_lnp = page_has_lnp_composition(text)

        results.append({
            "page": i + 1,
            "has_lnp_composition": has_lnp
        })

        print(f"Page {i+1}: {has_lnp}")

print("\nSummary:")
print(results)