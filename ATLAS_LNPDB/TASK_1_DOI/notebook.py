import pandas as pd

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)
pd.set_option("display.width", None)
pd.set_option("display.max_colwidth", None)

lnpdb_df = pd.read_csv("../LNPDB.csv", encoding="utf-8")

print("=== LNPDB columns ===")
print(lnpdb_df.columns.tolist())

pmid_counts = (
    lnpdb_df
    .groupby("Publication_PMID", dropna=False, sort=False)
    .size()
    .reset_index(name="count")
)

print("=== LNPDB: publication_pubmed_id 기준 count ===")
print(f"\n고유 publication_pubmed_id 수: {len(pmid_counts)}")
print(pmid_counts)

link_counts = (
    lnpdb_df
    .groupby("Publication_link", dropna=False, sort=False)
    .size()
    .reset_index(name="count")
)

print("\n=== LNPDB: publication_link 기준 count ===")
print(f"\n고유 publication_link 수: {len(link_counts)}")
print(link_counts)
