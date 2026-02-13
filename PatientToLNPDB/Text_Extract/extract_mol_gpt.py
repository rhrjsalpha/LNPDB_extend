import pandas as pd
import re
import json
import time
from pathlib import Path
from openai import OpenAI

# ==============================
# OpenAI API 설정
# ==============================
client = OpenAI(api_key="")

GPT_MODEL = "gpt-4o-mini"   # 빠르고 JSON 안정적
# GPT_MODEL = "gpt-4.1-mini"  # 더 정확 (느림)


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
# GPT 호출
# ==============================
def call_gpt(prompt):

    response = client.chat.completions.create(
        model=GPT_MODEL,
        messages=[
            {"role": "system", "content": "Return STRICT JSON only."},
            {"role": "user", "content": prompt}
        ],
        response_format={"type": "json_object"},  # JSON 강제
        temperature=0
    )

    return response.choices[0].message.content


# ==============================
# 화합물 정보 추출
# ==============================
def extract_compound_info(paragraph_list):

    extracted_data = []
    batch_size = 5

    prompt_template = """
Analyze the following patent paragraphs and extract information about chemical compounds.

For each compound found, provide:
1. compound_name
2. iupac_name
3. inchi
4. smiles
5. has_structure_fig (True/False)
6. other_info

Return STRICT JSON in this format:

{
  "results": [
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
}
"""

    for i in range(0, len(paragraph_list), batch_size):

        batch = paragraph_list[i:i+batch_size]
        batch_text = "\n\n".join(
            [f"TAG: {p['tag']}\nTEXT:\n{p['text']}" for p in batch]
        )

        print(f"Processing {batch[0]['tag']} → {batch[-1]['tag']}")

        full_prompt = f"{prompt_template}\n\nDATA:\n{batch_text}"

        try:
            response_text = call_gpt(full_prompt)
            parsed = json.loads(response_text)
            extracted_data.extend(parsed["results"])

        except Exception as e:
            print(f"Error in batch {i}: {e}")

        time.sleep(0.3)  # Rate control

    return extracted_data


# ==============================
# 메인
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

    csv_path = output_dir / f"{file_stem}_compounds_gpt.csv"
    df.to_csv(csv_path, index=False, encoding='utf-8-sig')

    print(f"Extraction complete → {csv_path}")


if __name__ == "__main__":
    target_path = "US20170210697A1.txt"
    main(target_path)
