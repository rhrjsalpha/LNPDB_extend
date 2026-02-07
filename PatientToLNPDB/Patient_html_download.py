import os
import re
import pandas as pd
import requests
from urllib.parse import urlparse

# ======================================================
# Config
# ======================================================
EXCEL_PATH = "20250409_LipidLibrarySample.xlsx"
URL_COLUMN = "Link"   # 🔴 실제 컬럼명에 맞게 수정
OUTPUT_DIR = "downloaded_html"

HEADERS = {
    "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36"
}

# ======================================================
# Helper functions
# ======================================================
def safe_filename(text, max_len=120):
    """파일명으로 안전하게 변환"""
    text = re.sub(r"[^\w\-_.]", "_", text)
    return text[:max_len]

def fetch_html(url):
    r = requests.get(url, headers=HEADERS, timeout=20)
    r.raise_for_status()
    return r.text

# ======================================================
# Main
# ======================================================
def main():

    # Excel 로드
    df = pd.read_excel(EXCEL_PATH)

    if URL_COLUMN not in df.columns:
        raise ValueError(f"❌ '{URL_COLUMN}' 컬럼이 Excel에 없습니다")

    # 출력 폴더 생성
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for idx, url in enumerate(df[URL_COLUMN].dropna(), start=1):
        try:
            print(f"[{idx}] Fetching: {url}")

            html = fetch_html(url)

            parsed = urlparse(url)
            base_name = safe_filename(parsed.netloc + parsed.path)

            if not base_name.endswith(".html"):
                base_name += ".html"

            save_path = os.path.join(OUTPUT_DIR, base_name)

            with open(save_path, "w", encoding="utf-8") as f:
                f.write(html)

            print(f"    ✅ Saved → {save_path}")

        except Exception as e:
            print(f"    ❌ Failed: {e}")

if __name__ == "__main__":
    main()