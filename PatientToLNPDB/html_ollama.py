import os
import re
import json
import requests
from bs4 import BeautifulSoup

# =========================
# Config
# =========================
HTML_PATH = "downloaded_html/patents.google.com_patent_US20170210697A1_en.html"
OLLAMA_URL = "http://localhost:11434/api/generate"
MODEL = "jinbora/deepseek-r1-Bllossom:8b"            # 또는 llama3.1:8b

CHUNK_CHARS = 8000               # 특허 HTML 길면 6k~10k 추천
OVERLAP = 500                    # 청크 경계 정보 누락 방지

# =========================
# Helpers
# =========================
def html_to_clean_text(html: str) -> str:
    soup = BeautifulSoup(html, "html.parser")

    # script/style 제거
    for tag in soup(["script", "style", "noscript"]):
        tag.decompose()

    # 줄 단위로 보기 좋게
    text = soup.get_text("\n")
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()

def chunk_text(text: str, chunk_chars: int, overlap: int):
    chunks = []
    i = 0
    n = len(text)
    while i < n:
        j = min(n, i + chunk_chars)
        chunks.append(text[i:j])
        if j == n:
            break
        i = max(0, j - overlap)
    return chunks

def ollama_json(prompt: str) -> dict:
    payload = {
        "model": MODEL,
        "prompt": prompt,
        "stream": False,
        "format": "json"
    }
    r = requests.post(OLLAMA_URL, json=payload, timeout=300)
    r.raise_for_status()
    out = r.json()
    # out["response"] 자체가 JSON 문자열일 수 있음
    return json.loads(out["response"])

# =========================
# Prompt template
# =========================
SYSTEM = """You are a scientific patent data curator.
Your task: scan the given chunk of a Google Patents HTML text (already cleaned).
Identify whether the chunk contains any of the following signals and return ONLY JSON.

Signals to detect:
- SMILES
- InChI
- InChIKey
- IUPAC chemical name (or IUPAC-like)
- lipid composition/formulation components (ionizable lipid, helper lipid, cholesterol, PEG-lipid, DSPC, DOPE, etc.)
- molar ratio / wt% / mol% / N:P ratio / formulation ratios
- formulation procedure (mixing, microfluidics, ethanol injection, pH, buffer, dialysis, etc.)

For each detected signal, include:
- where: a short local locator (e.g., nearest heading or first 8 words)
- evidence: a short exact snippet (max 160 chars)

Also infer section title if present in chunk (e.g., "Claims", "Description", "Examples", "Links", "Info").
Return JSON with keys exactly as specified.
"""

def make_prompt(chunk_id: int, chunk_text: str) -> str:
    return f"""{SYSTEM}

CHUNK_ID: {chunk_id}

TEXT:
{chunk_text}

Return JSON with schema:
{{
  "chunk_id": {chunk_id},
  "section_guess": "",
  "signals": {{
    "smiles": [],
    "inchi": [],
    "inchikey": [],
    "iupac_name": [],
    "lipid_composition": [],
    "molar_ratio": [],
    "formulation_procedure": []
  }}
}}
"""

# =========================
# Main
# =========================
def main():
    with open(HTML_PATH, "r", encoding="utf-8") as f:
        html = f.read()

    text = html_to_clean_text(html)
    chunks = chunk_text(text, CHUNK_CHARS, OVERLAP)

    all_chunk_results = []
    for cid, ch in enumerate(chunks):
        prompt = make_prompt(cid, ch)
        res = ollama_json(prompt)
        all_chunk_results.append(res)
        print(f"[chunk {cid+1}/{len(chunks)}] section={res.get('section_guess','')}")

    # ---- merge into a single map ----
    merged = {
        "signals": {
            "smiles": [],
            "inchi": [],
            "inchikey": [],
            "iupac_name": [],
            "lipid_composition": [],
            "molar_ratio": [],
            "formulation_procedure": []
        },
        "section_map": [],
        "overall_summary": {}
    }

    # collect signals with chunk id
    for r in all_chunk_results:
        cid = r["chunk_id"]
        sec = r.get("section_guess", "")
        for k in merged["signals"].keys():
            for item in r["signals"].get(k, []):
                item2 = dict(item)
                item2["chunk_id"] = cid
                item2["section_guess"] = sec
                merged["signals"][k].append(item2)

    # section map: contiguous chunk ranges with same section_guess
    cur = None
    for r in all_chunk_results:
        sec = r.get("section_guess", "") or "UNKNOWN"
        cid = r["chunk_id"]
        if cur is None or cur["section_title"] != sec:
            if cur:
                merged["section_map"].append(cur)
            cur = {
                "section_title": sec,
                "start_chunk": cid,
                "end_chunk": cid,
                "what_is_here": []
            }
        else:
            cur["end_chunk"] = cid
    if cur:
        merged["section_map"].append(cur)

    # overall summary
    merged["overall_summary"] = {
        "has_chem_structures": any(len(merged["signals"][k]) > 0 for k in ["smiles","inchi","inchikey","iupac_name"]),
        "has_formulation_info": any(len(merged["signals"][k]) > 0 for k in ["lipid_composition","molar_ratio","formulation_procedure"]),
        "notes": "Signals are heuristic; review evidence snippets for confirmation."
    }

    out_path = os.path.splitext(HTML_PATH)[0] + "_llm_map.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(merged, f, ensure_ascii=False, indent=2)

    print(f"\nSaved: {out_path}")
    print(json.dumps(merged["overall_summary"], indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()
