from bs4 import BeautifulSoup
import re
from collections import defaultdict
import ollama

HTML_PATH = "downloaded_html/patents.google.com_patent_US20170210697A1_en.html"

# ======================================================
# Pattern definitions
# ======================================================
PATTERNS = {
    "InChIKey": re.compile(r"\b[A-Z]{14}-[A-Z]{10}-[A-Z]\b"),
    "InChI": re.compile(r"InChI=1S?\/"),
    "SMILES": re.compile(r"\b[BCNOPSFIlcnopfsibrc0-9=\-\(\)\\\/]{10,}\b"),
    "IUPAC_like": re.compile(
        r"\b(acid|ester|amine|amide|alkane|alkene|ol|one|ium)\b",
        re.IGNORECASE
    ),
    "Lipid_keywords": re.compile(
        r"\b(lipid|ionizable|nanoparticle|lnp|helper lipid)\b",
        re.IGNORECASE
    ),
    "Composition_keywords": re.compile(
        r"\b(composition|formulation|molar ratio|wt%|mol%)\b",
        re.IGNORECASE
    )
}

# ======================================================
# Load HTML
# ======================================================
with open(HTML_PATH, encoding="utf-8") as f:
    soup = BeautifulSoup(f, "html.parser")

# ======================================================
# Sequential scan
# ======================================================
results = defaultdict(list)

for idx, el in enumerate(soup.find_all(True)):
    text = el.get_text(" ", strip=True)
    if not text:
        continue

    for key, pattern in PATTERNS.items():
        if pattern.search(text):
            results[key].append({
                "element_index": idx,
                "tag": el.name,
                "text_snippet": text[:200]
            })

# ======================================================
# Summary report
# ======================================================
summary = {
    k: len(v) for k, v in results.items()
}

print("=== Presence summary ===")
for k, v in summary.items():
    print(f"{k}: {v}")

MODEL = "jinbora/deepseek-r1-Bllossom:8b"
def ask_ollama_to_verify(text):
    """
    Ollama에게 해당 텍스트가 유효한 조성 정보(Composition/Ratio)인지 확인 요청
    반드시 TRUE 또는 SKIP만 반환하도록 강제
    """
    prompt = f"""
You are a strict classifier.

Decide whether the following patent text contains
EXPLICIT quantitative formulation information such as:
- lipid composition ratios (e.g. IL:HL:CHOL:PEG)
- molar percentage or weight percentage
- N:P ratio or nanoparticle:nucleic acid ratio

Rules:
- If such quantitative composition or ratio information is present → reply ONLY with TRUE
- If not present → reply ONLY with SKIP
- Do NOT explain.
- Do NOT output anything else.

Text:
\"\"\"{text}\"\"\"
"""

    response = ollama.generate(
        model="gemma3:1b",
        prompt=prompt
    )

    # 🔒 응답 정규화
    verdict = response["response"].strip().upper()

    if verdict.startswith("TRUE"):
        return "TRUE"
    else:
        return "SKIP"


# ======================================================
# Sequential scan & AI Verification
# ======================================================
import time
from datetime import datetime

elements = list(soup.find_all(True))
total_len = len(elements)

start_time = time.time()

total_elements = 0
candidate_elements = 0
ollama_calls = 0
accepted = 0
skipped = 0

refined_results = []

def format_time(sec):
    if sec < 60:
        return f"{sec:.1f}s"
    elif sec < 3600:
        return f"{sec/60:.1f}m"
    else:
        return f"{sec/3600:.2f}h"
# 너무 많은 호출을 방지하기 위해 Composition_keywords가 걸린 항목만 정밀 검사
for idx, el in enumerate(elements):
    total_elements += 1

    text = el.get_text(" ", strip=True)
    if not text or len(text) < 10:
        continue

    # 1차 필터
    if PATTERNS["Composition_keywords"].search(text) or \
       PATTERNS["Lipid_keywords"].search(text):

        candidate_elements += 1
        snippet = text[:300]

        print("\n" + "="*80)
        print(f"[{datetime.now().strftime('%H:%M:%S')}]")
        print(f"[*] Candidate #{candidate_elements}")
        print(f"[*] Element index: {idx}/{total_len} "
              f"({idx/total_len*100:.2f}%)")
        print(f"[*] Snippet:\n{snippet}")
        print("-"*80)

        # Ollama 호출
        t0 = time.time()
        ollama_calls += 1

        verdict = ask_ollama_to_verify(snippet)

        dt = time.time() - t0

        print(f"[OLLAMA] response time: {dt:.2f}s")
        print(f"[OLLAMA] verdict: {verdict}")

        if verdict != "TRUE":
            skipped += 1
            print("[RESULT] ❌ SKIPPED")
            continue

        accepted += 1
        print("[RESULT] ✅ ACCEPTED")

        refined_results.append({
            "index": idx,
            "raw_text": snippet,
            "ai_analysis": verdict,
            "ollama_time_sec": round(dt, 2)
        })

        # ===== 진행률 & ETA 출력 =====
        if ollama_calls >= 3:  # 초반 노이즈 방지
            elapsed = time.time() - start_time
            avg_per_call = elapsed / ollama_calls

            remaining_candidates_est = (
                candidate_elements / total_elements
            ) * (total_len - idx)

            eta_sec = avg_per_call * remaining_candidates_est

            print("\n--- Progress report ---")
            print(f"Scanned elements   : {idx}/{total_len}")
            print(f"Candidates found   : {candidate_elements}")
            print(f"Ollama calls       : {ollama_calls}")
            print(f"Accepted / Skipped : {accepted} / {skipped}")
            print(f"Elapsed time      : {format_time(elapsed)}")
            print(f"Avg / Ollama call : {avg_per_call:.2f}s")
            print(f"Estimated ETA     : {format_time(eta_sec)}")
            print("-----------------------\n")



# ======================================================
# 결과 출력
# ======================================================
for item in refined_results:
    print(f"\n[Found at Index {item['index']}]")
    print(f"AI Analysis: {item['ai_analysis']}")
