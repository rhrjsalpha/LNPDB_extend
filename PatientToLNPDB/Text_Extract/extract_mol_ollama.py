import pandas as pd
import re
import json
import time
from pathlib import Path
import ollama


OLLAMA_MODEL = "jinbora/deepseek-r1-Bllossom:8b"   # 필요시 변경

# ==============================
# 단락 분할
# ==============================
def split_by_paragraphs(text):
    pattern = r'(\[\d+\])'
    parts = re.split(pattern, text)

    paragraphs = []
    for i in range(1, len(parts), 2):
        tag = parts[i]
        content = parts[i+1] if i+1 < len(parts) else ""
        paragraphs.append({"tag": tag, "text": content.strip()})

    return paragraphs


# ==============================
# Ollama 호출 (로컬 SDK)
# ==============================
def call_ollama(prompt):
    response = ollama.chat(
        model=OLLAMA_MODEL,
        messages=[
            {"role": "system", "content": "Return STRICT JSON only."},
            {"role": "user", "content": prompt}
        ],
        format="json"   # JSON 강제
    )

    return response["message"]["content"]


# ==============================
# 단락 1개씩 처리
# ==============================
def extract_compound_info(paragraph_list):

    extracted_data = []

    prompt_template = """
Extract chemical compound information from the paragraph below.

Return STRICT JSON ONLY in this format:

[
  {
    "tag": "[123]",
    "compounds": [
      {
        "compound_name": "...",
        "iupac_name": "...",
        "inchi": "...",
        "smiles": "...",
        "has_structure_fig": true,
        "other_info": "..."
      }
    ]
  }
]

If no compound exists, return:

[
  {
    "tag": "[123]",
    "compounds": []
  }
]
"""

    total = len(paragraph_list)
    start_time = time.time()

    for idx, p in enumerate(paragraph_list):

        tag = p["tag"]
        text = p["text"]

        current_start = time.time()

        # 진행률
        progress = (idx + 1) / total * 100

        # 평균 처리 시간
        elapsed = time.time() - start_time
        avg_time = elapsed / (idx + 1)

        # 남은 시간 계산
        remaining = avg_time * (total - (idx + 1))

        # 시간 포맷 변환
        remaining_min = int(remaining // 60)
        remaining_sec = int(remaining % 60)

        print(
            f"[{idx+1}/{total}] {tag} | "
            f"{progress:.2f}% | "
            f"ETA: {remaining_min}m {remaining_sec}s"
        )

        full_prompt = f"""
TAG: {tag}
TEXT:
{text}

{prompt_template}
"""

        try:
            response_text = call_ollama(full_prompt)
            result = json.loads(response_text)
            extracted_data.extend(result)

        except Exception as e:
            print(f"Error at {tag}: {e}")
            continue

        # 실제 처리 시간
        current_elapsed = time.time() - current_start
        print(f"   → Took {current_elapsed:.2f}s")

    total_elapsed = time.time() - start_time
    print("\n==============================")
    print(f"Total Time: {total_elapsed/60:.2f} minutes")
    print("==============================\n")

    return extracted_data



# ==============================
# 메인 실행
# ==============================
def main(file_path):

    file_stem = Path(file_path).stem
    output_dir = Path(f"./compound_extraction_{file_stem}")
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(file_path, 'r', encoding='utf-8') as f:
        full_text = f.read()

    print("Splitting text...")
    paragraphs = split_by_paragraphs(full_text)

    results = extract_compound_info(paragraphs)

    flat_data = []
    for entry in results:
        tag = entry.get('tag')
        for comp in entry.get('compounds', []):
            comp['paragraph_no'] = tag
            flat_data.append(comp)

    df = pd.DataFrame(flat_data)

    if not df.empty:
        cols = [
            'paragraph_no',
            'compound_name',
            'iupac_name',
            'inchi',
            'smiles',
            'has_structure_fig',
            'other_info'
        ]
        df = df.reindex(columns=cols)

    csv_path = output_dir / f"{file_stem}_compounds_ollama.csv"
    df.to_csv(csv_path, index=False, encoding='utf-8-sig')

    print(f"Extraction complete → {csv_path}")


if __name__ == "__main__":
    target_path = "US20170210697A1.txt"
    main(target_path)
