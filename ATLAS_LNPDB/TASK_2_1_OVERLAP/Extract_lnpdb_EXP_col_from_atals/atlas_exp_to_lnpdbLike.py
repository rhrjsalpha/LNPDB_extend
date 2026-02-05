import pandas as pd
from collections import Counter

# ======================================================
# Load data
# ======================================================
path_to_csv = "/Users/kogeon/python_projects_path/LNPDB_extend/ATLAS_LNPDB/TASK_2_1_OVERLAP/YY_BY_DOI/10.1016_j.jconrel.2025.01.071/matched_atlas.csv"
atlas = pd.read_csv(path_to_csv)

# ======================================================
# 1. ATLAS bioactivity_profile 파싱
# ======================================================
def parse_atlas_profile(text):
    """
    'key: value; key2: value2' -> dict
    """
    result = {}
    if not isinstance(text, str):
        return result

    for field in text.split(";"):
        field = field.strip()
        if ":" not in field:
            continue
        key, value = field.split(":", 1)
        result[key.strip().lower()] = value.strip()

    return result


# ======================================================
# 2. LNPDB 컬럼별 추출 함수
# ======================================================
def extract_model_type(d):
    vals = []
    if "animal_model" in d:
        vals.append(d["animal_model"])
    if "cell_line" in d:
        vals.append(d["cell_line"])
    return vals


def extract_route(d):
    if "administration_route" not in d:
        return []

    val = d["administration_route"].lower()
    route_map = {
        "intramuscular": "intramuscular",
        "i.m.": "intramuscular",
        "intravenous": "intravenous",
        "i.v.": "intravenous",
        "subcutaneous": "subcutaneous",
        "s.c.": "subcutaneous"
    }

    routes = []
    for k, v in route_map.items():
        if k in val:
            routes.append(v)
    return list(set(routes))


def extract_model_target(d):
    if "biodistribution_result" not in d:
        return []

    val = d["biodistribution_result"].lower()
    targets = []

    for organ in [
        "liver", "muscle", "lung", "kidney",
        "spleen", "heart", "skin", "lymph"
    ]:
        if organ in val:
            targets.append(organ)

    return targets


def extract_experiment_method(d):
    if "gene_expression_result" not in d:
        return []

    val = d["gene_expression_result"].lower()
    methods = []

    if "ivis" in val:
        methods.append("IVIS imaging")
    if "luciferase" in val:
        methods.append("luciferase expression")

    return methods


def extract_dose_raw(d):
    return d.get("dose", None)


# ======================================================
# 3. 전체 row 처리
# ======================================================
rows = []

for text in atlas["bioactivity_profile"]:
    print(text)
    d = parse_atlas_profile(text)

    row = {
        "Model_type": extract_model_type(d),
        "Route_of_administration": extract_route(d),
        "Model_target": extract_model_target(d),
        "Experiment_method": extract_experiment_method(d),
        "Dose_raw": extract_dose_raw(d)
    }
    rows.append(row)

lnpdb_df = pd.DataFrame(rows)

print(lnpdb_df.head())
