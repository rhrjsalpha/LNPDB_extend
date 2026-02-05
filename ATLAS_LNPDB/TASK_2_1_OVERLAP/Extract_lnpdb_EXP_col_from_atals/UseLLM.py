from openai import OpenAI
import pandas as pd
import json

my_api_key = ""
client = OpenAI(api_key=my_api_key)

path_to_csv = "/Users/kogeon/python_projects_path/LNPDB_extend/ATLAS_LNPDB/TASK_2_1_OVERLAP/YY_BY_DOI/10.1016_j.jconrel.2025.01.071/matched_atlas.csv"
atlas = pd.read_csv(path_to_csv)

bio_text = atlas.loc[0, "bioactivity_profile"]

prompt = f"""
You are an expert in lipid nanoparticle (LNP) experimental data curation.

From the following ATLAS bioactivity description, extract values corresponding
to the LNPDB schema.

Rules:
- Only use information explicitly stated or clearly implied.
- If a value is not mentioned, return null.
- Normalize terminology to LNPDB-style labels.
- Return strictly valid JSON.
- Do not include explanations.

LNPDB schema:
- Cargo
- Cargo_type (one of: mRNA, siRNA, pDNA, unknown)
- Model_target
- Model (in vitro or in vivo)
- Route_of_administration
- Experiment_method
- Dose_ug_nucleicacid

ATLAS bioactivity profile:
\"\"\"
{bio_text}
\"\"\"
"""

response = client.chat.completions.create(
    model="gpt-4.1",
    messages=[{"role": "user", "content": prompt}],
    temperature=0
)

result = json.loads(response.choices[0].message.content)
print(result)

