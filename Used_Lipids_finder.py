import fitz  # PyMuPDF
import ollama
import os
import re
from typing import List, Dict, Tuple

# ============================================================
# 설정
# ============================================================
PDF_FILE = r"C:\Users\kogun\PycharmProjects\LNPDB\pdfs\MP_2025.pdf"
MODEL_SUBSECTION = "jinbora/deepseek-r1-Bllossom:8b"
MODEL_LIPID = "jinbora/deepseek-r1-Bllossom:8b"

MAX_PAGE_CHARS = 3500   # subsection 추출용
MAX_BLOCK_CHARS = 4000  # lipid 추출용

# ============================================================
# 유틸
# ============================================================
def remove_think_block(text: str) -> str:
    return re.sub(
        r"<think>.*?</think>",
        "",
        text,
        flags=re.DOTALL | re.IGNORECASE
    ).strip()


# ============================================================
# 1️⃣ Page 단위 subsection 제목 추출 (LLM)
# ============================================================
def extract_subsection_titles_from_page(
    page_text: str,
    page_index: int
) -> List[str]:

    PROMPT = f"""
You are a scientific paper STRUCTURE extractor.

Your task:
- From the given page text, extract ONLY subsection titles.

Rules:
- Output subsection titles ONLY
- One title per line
- No explanation
- No numbering
- No bullet points
- If no subsection titles exist, output NOTHING

What qualifies as a subsection title:
- Conceptual headings describing results, methods, analysis, comparisons
- Usually noun phrases or short sentences
- NOT reaction conditions
- NOT tables
- NOT figure labels (a, b, c)
- NOT single words like "Normal", "KO", "ClMg"

====================
PAGE {page_index}
TEXT:
"""

    page_text = page_text[:MAX_PAGE_CHARS]

    response = ollama.generate(
        model=MODEL_SUBSECTION,
        prompt=PROMPT + page_text,
        options={"temperature": 0}
    )

    raw = remove_think_block(response["response"])
    print("subsection raw response:", raw)
    titles = [
        line.strip()
        for line in raw.splitlines()
        if line.strip()
    ]

    return titles


# ============================================================
# 2️⃣ 전체 PDF 텍스트 로드
# ============================================================
def extract_full_text(pdf_path: str) -> str:
    doc = fitz.open(pdf_path)
    texts = []
    for page in doc:
        texts.append(page.get_text("text"))
    doc.close()
    return "\n".join(texts)


# ============================================================
# 3️⃣ subsection anchor 기반 텍스트 분할
# ============================================================
def split_text_by_subsections(
    full_text: str,
    subsection_titles: List[str]
) -> Dict[str, str]:

    blocks = {}
    positions: List[Tuple[int, str]] = []

    for title in subsection_titles:
        for m in re.finditer(re.escape(title), full_text):
            positions.append((m.start(), title))

    positions.sort()

    for i, (start, title) in enumerate(positions):
        end = positions[i + 1][0] if i + 1 < len(positions) else len(full_text)
        blocks[title] = full_text[start:end].strip()

    return blocks


# ============================================================
# 4️⃣ Lipid 이름 extractor (subsection 단위)
# ============================================================
def extract_lipids_from_text(
    subsection_title: str,
    text: str
) -> List[str]:

    PROMPT = f"""
You are an INFORMATION EXTRACTOR.

Task:
Extract LIPID NAMES ONLY from the given text.

This is NOT summarization.
This is NOT explanation.
If you output anything other than lipid names, it is WRONG.

Target lipid types:
- Ionizable lipids (e.g. DLin-MC3-DMA, ALC-0315, SM-102)
- Helper lipids (e.g. DOPE, DSPC, DOTAP)
- Cholesterol / sterol derivatives
- PEG-lipids (e.g. *-PEG2000, *-PEG5000)
- Rare formulation components (e.g. 6,6'-trehalose dioleate)

Rules:
- Output lipid names ONLY
- One lipid per line
- No numbering
- No explanation
- If none exist, output NOTHING

Subsection:
{subsection_title}

====================
TEXT:
"""

    text = text[:MAX_BLOCK_CHARS]

    response = ollama.generate(
        model=MODEL_LIPID,
        prompt=PROMPT + text,
        options={"temperature": 0}
    )

    raw = remove_think_block(response["response"])

    return sorted({
        line.strip()
        for line in raw.splitlines()
        if line.strip()
    })


# ============================================================
# MAIN
# ============================================================
def main():

    print("📄 PDF 로딩 중...")
    doc = fitz.open(PDF_FILE)

    all_subsection_titles: List[str] = []

    print("📑 Page 단위 subsection 추출 중...")
    for i, page in enumerate(doc):
        page_text = page.get_text("text")
        titles = extract_subsection_titles_from_page(page_text, i)
        if titles:
            print(f"Page {i} → {titles}")
            all_subsection_titles.extend(titles)

    doc.close()

    # 중복 제거 (순서 유지)
    seen = set()
    subsection_titles = []
    for t in all_subsection_titles:
        if t not in seen:
            seen.add(t)
            subsection_titles.append(t)

    print("\n📄 전체 텍스트 로드...")
    full_text = extract_full_text(PDF_FILE)

    print("✂️ Subsection 기준 텍스트 분할...")
    subsection_blocks = split_text_by_subsections(full_text, subsection_titles)

    print("\n🧪 Lipid 추출 시작...")
    results = {}

    for title, text in subsection_blocks.items():
        lipids = extract_lipids_from_text(title, text)
        if lipids:
            print(f"\n[{title}]")
            print(lipids)
            results[title] = lipids

    print("\n===== FINAL SUMMARY =====")
    all_lipids = sorted({l for v in results.values() for l in v})
    print(all_lipids)


if __name__ == "__main__":
    main()
