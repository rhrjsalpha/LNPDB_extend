import pandas as pd
import os

lnpdb = pd.read_csv("../LNPDB.csv")
atlas = pd.read_csv("../lnp_atlas_export_20260124_0129.csv", encoding="cp949")
print(len(atlas))
joined_doi = pd.read_csv("../TASK_1_DOI/DOI_overlap_ATLAS_LNPDB_with_meta.csv")

print(joined_doi.head())
print(len(lnpdb))

lnpdb_joined =pd.merge(joined_doi, lnpdb, on="Publication_PMID")
print(len(lnpdb_joined))

atlas_joined = pd.merge(joined_doi, atlas, on="paper_doi")
print(len(atlas_joined))
print(atlas_joined)

BASE_DIR = "./YY_BY_DOI"
os.makedirs(BASE_DIR, exist_ok=True)
shape_records = []

for _, doi_row in joined_doi.iterrows():
    DOI = doi_row["DOI"]
    lnpdb_pmid = int(doi_row["Publication_PMID"])
    atlas_paper_doi = doi_row["paper_doi"]

    # -------------------------
    # DOI 폴더 생성
    # -------------------------
    safe_doi = DOI.replace("/", "_")
    doi_dir = os.path.join(BASE_DIR, safe_doi)
    os.makedirs(doi_dir, exist_ok=True)

    # -------------------------
    # LNPDB: PMID 기준 필터
    # -------------------------
    lnpdb_each_doi_df = lnpdb_joined[
        lnpdb_joined["Publication_PMID"] == lnpdb_pmid
    ]

    # -------------------------
    # ATLAS: paper_doi 기준 필터
    # -------------------------
    atlas_each_doi_df = atlas_joined[
        atlas_joined["paper_doi"] == atlas_paper_doi
    ]

    # -------------------------
    # CSV 저장
    # -------------------------
    lnpdb_path = os.path.join(doi_dir, "lnpdb.csv")
    atlas_path = os.path.join(doi_dir, "atlas.csv")

    lnpdb_each_doi_df.to_csv(lnpdb_path, index=False)
    atlas_each_doi_df.to_csv(atlas_path, index=False)

    shape_records.append({
        "DOI": DOI,
        "LNPDB_rows": lnpdb_each_doi_df.shape[0],
        "ATLAS_rows": atlas_each_doi_df.shape[0],
    })

    print(
        f"[{DOI}] "
        f"LNPDB: {lnpdb_each_doi_df.shape} | "
        f"ATLAS: {atlas_each_doi_df.shape}"
    )

shape_df = pd.DataFrame(shape_records)
shape_df.to_csv("DOI_row_counts_LNPDB_ATLAS.csv", index=False)

print(shape_df)

plot_df = shape_df.melt(
    id_vars="DOI",
    value_vars=["LNPDB_rows", "ATLAS_rows"],
    var_name="Database",
    value_name="Rows"
)

plot_df["Database"] = plot_df["Database"].replace({
    "LNPDB_rows": "LNPDB",
    "ATLAS_rows": "ATLAS"
})

import matplotlib.pyplot as plt
import numpy as np

ROTATION_ANGLE = 80
FIG_WIDTH = 14
FIG_HEIGHT = 6

dois = shape_df["DOI"].values
x = np.arange(len(dois))
width = 0.38

plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))

plt.bar(
    x - width/2,
    shape_df["LNPDB_rows"],
    width,
    label="LNPDB",
    color="#1f77b4"
)

plt.bar(
    x + width/2,
    shape_df["ATLAS_rows"],
    width,
    label="ATLAS",
    color="#ff7f0e"
)

plt.xticks(
    ticks=x,
    labels=dois,
    rotation=ROTATION_ANGLE,
    ha="right"
)

plt.ylabel("Number of rows")
plt.xlabel("DOI")
plt.title("Number of entries per DOI (LNPDB vs ATLAS)")
plt.legend()
plt.tight_layout()
plt.show()
