import pandas as pd
import re
import requests
import time

# -----------------------------
# PMID → DOI
# -----------------------------
def pmid_to_doi(pmid):
    if pd.isna(pmid):
        return None

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "pubmed",
        "id": str(pmid),
        "retmode": "xml"
    }

    try:
        # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi 가 params에 의해 아래와 같이 변하여 request를 보냄
        # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=(pmid)&retmode=xml
        r = requests.get(url, params=params, timeout=10)
        if r.status_code != 200:
            return None
        print(r.text)

        m = re.search(r'<ArticleId IdType="doi">(.*?)</ArticleId>', r.text)
        return m.group(1) if m else None

    except Exception:
        return None


# -----------------------------
# LNPDB 로드
# -----------------------------
df = pd.read_csv("../LNPDB.csv")

# -----------------------------
# 논문 단위 테이블 생성
# -----------------------------
paper_df = (
    df
    .groupby(["Publication_PMID", "Publication_link"])
    .size()
    .reset_index(name="count")
)

print(f"총 고유 논문 수: {len(paper_df)}")

# -----------------------------
# PMID → DOI (중복 호출 방지)
# -----------------------------
paper_df["DOI_from_pmid"] = None

mask = paper_df["Publication_PMID"].notna()

for idx, pmid in paper_df.loc[mask, "Publication_PMID"].items():
    doi = pmid_to_doi(pmid)
    paper_df.at[idx, "DOI_from_pmid"] = doi
    time.sleep(0.34)  # NCBI rate limit 보호

# -----------------------------
# 결과 확인
# -----------------------------
print("\n=== DOI 수집 요약 ===")
print("PMID 있는 논문 수:", mask.sum())
print("DOI 확보 논문 수:", paper_df["DOI_from_pmid"].notna().sum())
print("DOI 미확보 논문 수:", paper_df["DOI_from_pmid"].isna().sum())

# -----------------------------
# 저장
# -----------------------------
ROTATION_ANGLE = 80   # x축 DOI 기울기
TOP_N = None          # 상위 N개만 표시 (None이면 전체)
FIG_WIDTH = 14
FIG_HEIGHT = 6

paper_df.to_csv("LNPDB_paper_table_with_DOI.csv", index=False)

print("\n저장 완료: LNPDB_paper_table_with_DOI.csv")

plot_df = paper_df.copy()

# DOI 정규화
plot_df["DOI_Norm"] = plot_df["DOI_from_pmid"]

# DOI 없는 논문 제거 (선택)
plot_df = plot_df[plot_df["DOI_Norm"].notna()]

# count 기준 정렬
plot_df = plot_df.sort_values("count", ascending=False)

# 상위 N개 제한
if TOP_N is not None:
    plot_df = plot_df.head(TOP_N)

import matplotlib.pyplot as plt

plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))

plt.bar(
    plot_df["DOI_Norm"],
    plot_df["count"]
)

plt.xlabel("DOI")
plt.ylabel("Number of LNP entries")
plt.title("LNPDB: Number of entries per DOI")

plt.xticks(
    rotation=ROTATION_ANGLE,
    ha="right"
)

plt.tight_layout()
plt.show()