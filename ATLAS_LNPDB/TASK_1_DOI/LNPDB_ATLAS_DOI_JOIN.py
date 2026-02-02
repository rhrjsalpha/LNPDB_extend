import pandas as pd
import re

# -----------------------------
# DOI 정규화 함수 (공통)
# -----------------------------
def normalize_doi(doi):
    if pd.isna(doi):
        return None

    doi = str(doi).strip()
    doi = re.sub(r"^https?://(dx\.)?doi\.org/", "", doi, flags=re.IGNORECASE)
    doi = re.sub(r"^doi:\s*", "", doi, flags=re.IGNORECASE)

    return doi if doi else None


# =============================
# 1. CSV 로드
# =============================
atlas_df = pd.read_csv("ATLAS_paper_doi_title_final.csv")
lnpdb_df = pd.read_csv("LNPDB_paper_table_with_DOI.csv")

# =============================
# 2. DOI 정규화
# =============================
atlas_df["DOI_norm"] = atlas_df["DOI_ATLAS_Norm"].apply(normalize_doi)
lnpdb_df["DOI_norm"] = lnpdb_df["DOI_from_pmid"].apply(normalize_doi)

# DOI 없는 행 제거
atlas_df = atlas_df.dropna(subset=["DOI_norm"])
lnpdb_df = lnpdb_df.dropna(subset=["DOI_norm"])

# =============================
# 3. DOI 집합
# =============================
atlas_dois = set(atlas_df["DOI_norm"].unique())
lnpdb_dois = set(lnpdb_df["DOI_norm"].unique())

atlas_only = atlas_dois - lnpdb_dois        # XX
overlap = atlas_dois & lnpdb_dois            # YY
lnpdb_only = lnpdb_dois - atlas_dois          # ZZ

# =============================
# 4. ATLAS only (XX)
# =============================
atlas_only_df = (
    atlas_df[atlas_df["DOI_norm"].isin(atlas_only)]
    [["DOI_norm", "paper_doi", "paper_title", "count"]]
    .drop_duplicates()
    .rename(columns={"DOI_norm": "DOI"})
)

# =============================
# 5. LNPDB only (ZZ)
# =============================
lnpdb_only_df = (
    lnpdb_df[lnpdb_df["DOI_norm"].isin(lnpdb_only)]
    [["DOI_norm", "Publication_PMID", "Publication_link", "count"]]
    .drop_duplicates()
    .rename(columns={"DOI_norm": "DOI"})
)

# =============================
# 6. Overlap (YY) — 핵심
# =============================
overlap_atlas = atlas_df[atlas_df["DOI_norm"].isin(overlap)][
    ["DOI_norm", "paper_doi", "paper_title", "count"]
].drop_duplicates()

overlap_lnpdb = lnpdb_df[lnpdb_df["DOI_norm"].isin(overlap)][
    ["DOI_norm", "Publication_PMID", "Publication_link", "count"]
].drop_duplicates()

overlap_df = (
    overlap_atlas
    .merge(
        overlap_lnpdb,
        on="DOI_norm",
        how="inner",
        suffixes=("_atlas", "_lnpdb")
    )
    .rename(columns={"DOI_norm": "DOI"})
)

# =============================
# 7. 결과 요약
# =============================
print("\n=== DOI Set Summary ===")
print(f"ATLAS only (XX): {len(atlas_only_df)}")
print(f"Overlap (YY):    {len(overlap_df)}")
print(f"LNPDB only (ZZ): {len(lnpdb_only_df)}")

# =============================
# 8. CSV 저장
# =============================
atlas_only_df.to_csv("DOI_ATLAS_only_with_title.csv", index=False)
overlap_df.to_csv("DOI_overlap_ATLAS_LNPDB_with_meta.csv", index=False)
lnpdb_only_df.to_csv("DOI_LNPDB_only_with_pmid.csv", index=False)

print("\nCSV 저장 완료:")
print("- DOI_ATLAS_only_with_title.csv")
print("- DOI_overlap_ATLAS_LNPDB_with_meta.csv")
print("- DOI_LNPDB_only_with_pmid.csv")

