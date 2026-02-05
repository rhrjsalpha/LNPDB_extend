import pandas as pd
from collections import Counter

# ======================================================
# Load data
# ======================================================
path_to_csv = "/Users/kogeon/python_projects_path/LNPDB_extend/ATLAS_LNPDB/lnp_atlas_export_20260124_0129.csv"
atlas = pd.read_csv(path_to_csv, encoding="latin1")

# ======================================================
# Extract keys from bioactivity_profile
# ======================================================
key_counter = Counter()

for text in atlas["bioactivity_profile"].dropna():

    # 1) split by ';'
    fields = text.split(";")

    for field in fields:
        field = field.strip()
        if not field:
            continue

        # 2) split key and value by first ':'
        if ":" in field:
            key = field.split(":", 1)[0].strip().lower()
            key_counter[key] += 1

# ======================================================
# Convert to DataFrame
# ======================================================
key_df = (
    pd.DataFrame(key_counter.items(), columns=["field_name", "count"])
      .sort_values("count", ascending=False)
      .reset_index(drop=True)
)

print(key_df)
