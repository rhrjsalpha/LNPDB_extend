import pandas as pd
import os
from collections import Counter


def select_inchi(row):

    sources = {
        "OPSIN": row["opsin_inchi"],
        "PubChem": row["pubchem_inchi"],
        "CIR": row["cir_inchi"]
    }

    # None 제거
    valid = {k: v for k, v in sources.items() if pd.notna(v) and v}

    if len(valid) == 0:
        return None, 0, None

    # InChI 값별로 그룹화
    counter = Counter(valid.values())

    # 가장 많이 나온 InChI
    most_common_inchi, count = counter.most_common(1)[0]

    # 동일 InChI 가진 source 찾기
    matched_sources = [k for k, v in valid.items() if v == most_common_inchi]

    # === 규칙 적용 ===

    # 3개 존재
    if len(valid) == 3:
        if count >= 2:
            return most_common_inchi, count, ";".join(matched_sources)
        else:
            return "CHECK_REQUIRED", 1, ";".join(valid.keys())

    # 2개 존재
    if len(valid) == 2:
        if count == 2:
            return most_common_inchi, 2, ";".join(matched_sources)
        else:
            return "CHECK_REQUIRED", 1, ";".join(valid.keys())

    # 1개 존재
    if len(valid) == 1:
        only_source = list(valid.keys())[0]
        return list(valid.values())[0], 1, only_source


def process_file(csv_path):

    df = pd.read_csv(csv_path)

    selected_inchi = []
    selected_count = []
    selected_sources = []

    for _, row in df.iterrows():
        s_inchi, s_count, s_sources = select_inchi(row)

        selected_inchi.append(s_inchi)
        selected_count.append(s_count)
        selected_sources.append(s_sources)

    df["selected_inchi"] = selected_inchi
    df["selected_count"] = selected_count
    df["selected_sources"] = selected_sources

    output_path = os.path.splitext(csv_path)[0] + "_with_selection.csv"
    df.to_csv(output_path, index=False)

    print(f"완료 → {output_path}")


if __name__ == "__main__":
    process_file(
        "compound_extraction_US20170210697A1/resolved_output_all_sources.csv"
    )
