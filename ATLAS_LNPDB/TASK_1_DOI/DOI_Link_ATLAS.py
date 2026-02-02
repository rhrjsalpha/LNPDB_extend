import pandas as pd
import re
import matplotlib.pyplot as plt

# DOI 정규화 : https 부터 doi.org/ 있으면 제거함
def normalize_doi(doi):
    if pd.isna(doi):
        return None

    doi = str(doi).strip()

    # https://doi.org/ 제거
    doi = re.sub(r"^https?://(dx\.)?doi\.org/", "", doi, flags=re.IGNORECASE)

    # doi: prefix 제거
    doi = re.sub(r"^doi:\s*", "", doi, flags=re.IGNORECASE)

    return doi if doi else None


# ATLAS
atlas_df = pd.read_csv("../lnp_atlas_export_20260124_0129.csv", encoding="cp949")

print(f"ATLAS 전체 row 수: {len(atlas_df)}")

# paper_doi 단독 groupby
doi_only_df = (
    atlas_df
    .groupby("paper_doi", dropna=False, sort=False)
    .size()
    .reset_index(name="count")
)

doi_only_df["DOI_ATLAS_Norm"] = doi_only_df["paper_doi"].apply(normalize_doi)

print(f"\n[paper_doi 기준 고유 논문 수]: {len(doi_only_df)}")

# paper_title 단독 groupby
title_only_df = (
    atlas_df
    .groupby("paper_title", dropna=False, sort=False)
    .size()
    .reset_index(name="count")
)

print(f"[paper_title 기준 고유 논문 수]: {len(title_only_df)}")


# DOI 와 title 불일치 진단


# 3-1. 하나의 DOI에 여러 title이 매칭되는 경우
doi_title_mismatch = (
    atlas_df
    .groupby("paper_doi")["paper_title"]
    .nunique(dropna=False)
    .reset_index(name="n_unique_titles")
)

doi_title_mismatch = doi_title_mismatch[
    doi_title_mismatch["n_unique_titles"] > 1
]

print(f"\n[주의 DOI 하나에 여러 title 존재]: {len(doi_title_mismatch)}")

# 3-2. 하나의 title에 여러 DOI가 매칭되는 경우
title_doi_mismatch = (
    atlas_df
    .groupby("paper_title")["paper_doi"]
    .nunique(dropna=False)
    .reset_index(name="n_unique_dois")
)

title_doi_mismatch = title_doi_mismatch[
    title_doi_mismatch["n_unique_dois"] > 1
]

print(f"[주의 title 하나에 여러 DOI 존재]: {len(title_doi_mismatch)}")


# 최종: paper_doi + paper_title 동시 groupby
final_paper_df = (
    atlas_df
    .groupby(["paper_doi", "paper_title"], dropna=False, sort=False)
    .size()
    .reset_index(name="count")
)

final_paper_df["DOI_ATLAS_Norm"] = final_paper_df["paper_doi"].apply(normalize_doi)

print(f"\n[paper_doi + paper_title 기준 고유 논문 수]: {len(final_paper_df)}")


# 저장
#doi_only_df.to_csv("ATLAS_paper_doi_only.csv", index=False)
#title_only_df.to_csv("ATLAS_paper_title_only.csv", index=False)
final_paper_df.to_csv("ATLAS_paper_doi_title_final.csv", index=False)

print("\n저장 완료:")
print("- ATLAS_paper_doi_only.csv")
print("- ATLAS_paper_title_only.csv")
print("- ATLAS_paper_doi_title_final.csv")

# ===============================
# 시각화 파라미터 (여기만 조절)
# ===============================
ROTATION_ANGLE = 80   # x축 DOI 글자 기울기 (0~90 권장)
TOP_N = None            # 상위 N개 DOI만 표시 (None이면 전체)
FIG_WIDTH = 14
FIG_HEIGHT = 6

# ===============================
# 데이터 준비
# ===============================
plot_df = doi_only_df.copy()

# 정규화 DOI 없는 경우 제외 (선택)
plot_df = plot_df[plot_df["DOI_ATLAS_Norm"].notna()]

# count 기준 정렬
plot_df = plot_df.sort_values("count", ascending=False)

# 상위 N개 제한
if TOP_N is not None:
    plot_df = plot_df.head(TOP_N)

# ===============================
# Bar plot
# ===============================
plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))

plt.bar(
    plot_df["DOI_ATLAS_Norm"],
    plot_df["count"]
)

plt.xlabel("DOI")
plt.ylabel("Number of entries")
plt.title("ATLAS: Number of entries per DOI")

plt.xticks(
    rotation=ROTATION_ANGLE,
    ha="right"   # 글자가 겹치지 않게 우측 정렬
)

plt.tight_layout()
plt.show()