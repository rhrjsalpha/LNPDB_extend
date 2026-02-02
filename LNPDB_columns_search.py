import pandas as pd

LNPDB = pd.read_csv("LNPDB.csv")
print(LNPDB.columns)

def inspect_unique_values(df, col, top_n=50):
    values = (
        df[col]
        .dropna()
        .astype(str)
        .str.strip()
    )
    uniq = values.value_counts()
    print(f"[{col}] unique count =", uniq.shape[0])
    print(uniq.head(top_n))

inspect_unique_values(LNPDB, "IL_name")
inspect_unique_values(LNPDB, "HL_name")
inspect_unique_values(LNPDB, "CHL_name")
inspect_unique_values(LNPDB, "PEG_name")
inspect_unique_values(LNPDB, "fifthcomponent_name")