import pandas as pd
import re
import os

# ======================================================
# Config (DOI_DIR만 바꾸면 됨)
# ======================================================
DOI_DIR = "/Users/kogeon/python_projects_path/LNPDB_extend/ATLAS_LNPDB/TASK_2_1_OVERLAP/YY_BY_DOI/10.1016_j.jconrel.2025.01.071"
ATLAS_CSV = os.path.join(DOI_DIR, "matched_atlas.csv")
LNPDB_CSV = os.path.join(DOI_DIR, "matched_lnpdb.csv")
OUTPUT_DIR = os.path.join(DOI_DIR, "ATLAS_LNPDB_by_matched_id")
os.makedirs(OUTPUT_DIR, exist_ok=True)

assert os.path.exists(ATLAS_CSV), f"Not found: {ATLAS_CSV}"
assert os.path.exists(LNPDB_CSV), f"Not found: {LNPDB_CSV}"

# ======================================================
# Helpers
# ======================================================
def split_bioactivity(text):
    d = {}
    if pd.isna(text):
        return d
    for part in str(text).split(";"):
        part = part.strip()
        if ":" in part:
            k, v = part.split(":", 1)
            d[k.strip().lower()] = v.strip()
    return d

def normalize_list(x):
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return []
    s = str(x).strip()
    if not s:
        return []
    return [i.strip() for i in re.split(r",| and ", s) if i.strip()]

def find_targets(text):
    text = (text or "").lower()
    hits = []
    for t in ["liver","muscle","spleen","lung","kidney","lymph node","skin","heart","bone marrow","whole body","multiorgan","ear"]:
        if t in text:
            hits.append(t)
    hits = [
        h.replace("lymph node","lymph_node")
         .replace("bone marrow","bone_marrow")
         .replace("whole body","whole_body")
        for h in hits
    ]
    return sorted(set(hits))

def guess_experiment_method(gene_expr_text):
    t = (gene_expr_text or "").lower()
    if "ivis" in t:
        return ["IVIS"]
    if "luciferase" in t:
        return ["luciferase"]
    if "flow cytometry" in t or "facs" in t:
        return ["flow_cytometry"]
    if "elisa" in t:
        return ["ELISA"]
    return []

def safe_get(row, col):
    return row[col] if col in row.index else None

def lower_or_none(x):
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return None
    s = str(x).strip()
    return s.lower() if s else None

# ======================================================
# ATLAS_EXP (bioactivity + target_type/nucleic_acid_sequence + atlas_ratio/ratio_order + loading_capacity_std)
# ======================================================
def build_atlas_exp(atlas_df):
    rows = []

    for _, r in atlas_df.iterrows():
        parsed = split_bioactivity(r.get("bioactivity_profile"))

        # (그림1) 매핑
        cargo = safe_get(r, "target_type")                 # LNPDB Cargo
        cargo_type = safe_get(r, "nucleic_acid_sequence")  # LNPDB Cargo_type

        # (그림3) 조성 관련: atlas_ratio / ratio_order
        atlas_ratio = safe_get(r, "atlas_ratio")
        ratio_order = safe_get(r, "ratio_order")

        # (그림4) 추가 컬럼
        loading_capacity_std = safe_get(r, "loading_capacity_std")

        # Model 분기: animal_model → in_vivo, cell_line → in_vitro
        models = []
        for a in normalize_list(parsed.get("animal_model")):
            models.append(("in_vivo", a))
        for c in normalize_list(parsed.get("cell_line")):
            models.append(("in_vitro", c))
        if not models:
            models = [(None, None)]

        for model, model_type in models:
            rows.append({
                "matched_id": safe_get(r, "matched_id"),
                "atlas_row_id": safe_get(r, "atlas_row_id"),

                # ===== EXP 요약 =====
                "Cargo": cargo,
                "Cargo_type": cargo_type,
                "Model": model,
                "Model_type": model_type,
                "Route_of_administration": lower_or_none(parsed.get("administration_route")),
                "Model_target": find_targets(parsed.get("biodistribution_result")),
                "Experiment_method": guess_experiment_method(parsed.get("gene_expression_result")),
                "Dose_raw": parsed.get("dose"),
                "Toxicity_profile": parsed.get("toxicity_profile"),
                "Cytokine_response": parsed.get("cytokine_response"),

                # ===== 조성/로딩 =====
                "atlas_ratio": atlas_ratio,
                "ratio_order": ratio_order,
                "loading_capacity_std": loading_capacity_std,
            })

    return pd.DataFrame(rows)

# ======================================================
# LNPDB_EXP (EXP + 조성정보(그림2) + IL_to_nucleicacid_massratio)
# ======================================================
def build_lnpdb_exp(lnpdb_df):
    rows = []
    for _, r in lnpdb_df.iterrows():
        rows.append({
            "matched_id": safe_get(r, "matched_id"),
            "lnpdb_row_id": safe_get(r, "lnpdb_row_id"),

            # ===== EXP =====
            "Cargo": safe_get(r, "Cargo"),
            "Cargo_type": safe_get(r, "Cargo_type"),
            "Model": safe_get(r, "Model"),
            "Model_type": safe_get(r, "Model_type"),
            "Model_target": safe_get(r, "Model_target"),
            "Route_of_administration": lower_or_none(
                safe_get(r, "Route_of_administration") or safe_get(r, "Route_of_adm") or safe_get(r, "Route_of_admin")
            ),
            "Experiment_method": safe_get(r, "Experiment_method"),
            "Dose_ug_nucleicacid": safe_get(r, "Dose_ug_nucleicacid"),
            "Dose_raw": safe_get(r, "Dose_raw"),

            # ===== (그림2) 조성 =====
            "IL_name": safe_get(r, "IL_name"),
            "IL_SMILES": safe_get(r, "IL_SMILES"),
            "HL_name": safe_get(r, "HL_name"),
            "HL_SMILES": safe_get(r, "HL_SMILES"),
            "CHL_name": safe_get(r, "CHL_name"),
            "CHL_SMILES": safe_get(r, "CHL_SMILES"),
            "PEG_name": safe_get(r, "PEG_name"),
            "PEG_SMILES": safe_get(r, "PEG_SMILES"),

            "IL_molratio": safe_get(r, "IL_molratio"),
            "HL_molratio": safe_get(r, "HL_molratio"),
            "CHL_molratio": safe_get(r, "CHL_molratio"),
            "PEG_molratio": safe_get(r, "PEG_molratio"),

            # ===== (그림4) LNPDB 측 =====
            "IL_to_nucleicacid_massratio": safe_get(r, "IL_to_nucleicacid_massratio"),
        })
    return pd.DataFrame(rows)

# ======================================================
# Load
# ======================================================
atlas = pd.read_csv(ATLAS_CSV)
lnpdb = pd.read_csv(LNPDB_CSV)

atlas_exp = build_atlas_exp(atlas)
lnpdb_exp = build_lnpdb_exp(lnpdb)

# ======================================================
# Save per matched_id (RAW + EXP + COMPARE)
# ======================================================
common_ids = sorted(set(atlas["matched_id"]) & set(lnpdb["matched_id"]))

for mid in common_ids:
    atlas_raw_sub = atlas[atlas["matched_id"] == mid]
    lnpdb_raw_sub = lnpdb[lnpdb["matched_id"] == mid]

    atlas_exp_sub = atlas_exp[atlas_exp["matched_id"] == mid]
    lnpdb_exp_sub = lnpdb_exp[lnpdb_exp["matched_id"] == mid]

    # 비교: ATLAS_EXP 행 x LNPDB_EXP 행 (suffix로 구분)
    compare = atlas_exp_sub.merge(
        lnpdb_exp_sub,
        on="matched_id",
        how="outer",
        suffixes=("_ATLAS", "_LNPDB")
    )

    out_xlsx = os.path.join(OUTPUT_DIR, f"{mid}_ATLAS_LNPDB.xlsx")

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        atlas_raw_sub.to_excel(writer, sheet_name=f"{mid}_ATLAS_RAW", index=False)
        lnpdb_raw_sub.to_excel(writer, sheet_name=f"{mid}_LNPDB_RAW", index=False)

        atlas_exp_sub.to_excel(writer, sheet_name=f"{mid}_ATLAS_EXP", index=False)
        lnpdb_exp_sub.to_excel(writer, sheet_name=f"{mid}_LNPDB_EXP", index=False)

        compare.to_excel(writer, sheet_name=f"{mid}_COMPARE", index=False)

    print(f"✅ saved: {out_xlsx}")

print(f"\nDone. Total files created: {len(common_ids)}")

