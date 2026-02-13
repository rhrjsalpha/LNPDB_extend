import pandas as pd
import requests
import urllib.parse
import time
import os
from collections import deque


# ===============================
# 1️⃣ 필터링
# ===============================
def filter_dataframe(df):

    df = df[
        df["compound_name"].str.contains("COMPOUND", case=False, na=False)
    ]

    df = df[df["iupac_name"].notna()]
    df = df[df["iupac_name"].str.strip() != ""]

    return df


# ===============================
# 2️⃣ OPSIN
# ===============================
def query_opsin(name):
    try:
        url = f"https://opsin.ch.cam.ac.uk/opsin/{urllib.parse.quote(name)}.json"
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            data = r.json()
            return data.get("smiles"), data.get("inchi")
    except:
        pass
    return None, None


# ===============================
# 3️⃣ PubChem
# ===============================
def query_pubchem(name):
    try:
        base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name"
        url = f"{base}/{urllib.parse.quote(name)}/property/IsomericSMILES,InChI/JSON"
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            data = r.json()
            props = data["PropertyTable"]["Properties"][0]
            return props.get("IsomericSMILES"), props.get("InChI")
    except:
        pass
    return None, None


# ===============================
# 4️⃣ CIR
# ===============================
def query_cir(name):
    try:
        base = "https://cactus.nci.nih.gov/chemical/structure"
        smiles_url = f"{base}/{urllib.parse.quote(name)}/smiles"
        inchi_url = f"{base}/{urllib.parse.quote(name)}/inchi"

        smiles = requests.get(smiles_url, timeout=10).text.strip()
        inchi = requests.get(inchi_url, timeout=10).text.strip()

        if "Page not found" not in smiles:
            return smiles, inchi
    except:
        pass
    return None, None


# ===============================
# 5️⃣ 메인 처리 (3군데 전부 조회 + ETA)
# ===============================
def enrich_smiles_inchi(csv_path):

    csv_path = os.path.abspath(csv_path)
    base_dir = os.path.dirname(csv_path)

    print(f"CSV 위치: {csv_path}")

    df = pd.read_csv(csv_path)
    df = filter_dataframe(df)
    df = df.reset_index(drop=True)

    total = len(df)
    start_time = time.time()
    recent_times = deque(maxlen=20)

    # 결과 저장 리스트
    opsin_smi, opsin_inchi = [], []
    pubchem_smi, pubchem_inchi = [], []
    cir_smi, cir_inchi = [], []

    for idx, row in df.iterrows():

        iter_start = time.time()
        name = row["iupac_name"]

        # ===== ETA 계산 =====
        progress = (idx + 1) / total * 100
        if len(recent_times) > 0:
            avg_recent = sum(recent_times) / len(recent_times)
        else:
            avg_recent = (time.time() - start_time) / (idx + 1)

        remaining = total - (idx + 1)
        eta_seconds = remaining * avg_recent
        eta_min = int(eta_seconds // 60)
        eta_sec = int(eta_seconds % 60)

        print(f"[{idx+1}/{total}] {progress:.2f}% | ETA: {eta_min}m {eta_sec}s")
        print(f"   Resolving: {name}")

        # ===== 3군데 모두 조회 =====
        s1, i1 = query_opsin(name)
        s2, i2 = query_pubchem(name)
        s3, i3 = query_cir(name)

        opsin_smi.append(s1)
        opsin_inchi.append(i1)

        pubchem_smi.append(s2)
        pubchem_inchi.append(i2)

        cir_smi.append(s3)
        cir_inchi.append(i3)

        iter_time = time.time() - iter_start
        recent_times.append(iter_time)

        print(f"   → Took {iter_time:.2f}s\n")

        time.sleep(0.3)

    # ===== 컬럼 추가 =====
    df["opsin_smiles"] = opsin_smi
    df["opsin_inchi"] = opsin_inchi

    df["pubchem_smiles"] = pubchem_smi
    df["pubchem_inchi"] = pubchem_inchi

    df["cir_smiles"] = cir_smi
    df["cir_inchi"] = cir_inchi

    output_path = os.path.join(base_dir, "resolved_output_all_sources.csv")
    df.to_csv(output_path, index=False)

    total_time = time.time() - start_time
    print(f"\n완료 → {output_path}")
    print(f"총 소요 시간: {int(total_time//60)}m {int(total_time%60)}s")


if __name__ == "__main__":
    enrich_smiles_inchi(
        "compound_extraction_US20170210697A1/US20170210697A1_compounds_gemini_flash-light_5batch.csv"
    )
