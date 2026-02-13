from rdkit import Chem
from rdkit.Chem import Draw
import os

# ===============================
# 설정
# ===============================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))  # 현재 스크립트 폴더
OUTPUT_DIR = os.path.join(BASE_DIR, "output_imgs")
IMG_SIZE = (300, 300)

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ===============================
# .smi / .smiles 파일 탐색
# ===============================
smi_files = [
    f for f in os.listdir(BASE_DIR)
    if f.lower().endswith((".smi", ".smiles"))
]

if not smi_files:
    raise RuntimeError("❌ 현재 폴더에 .smi / .smiles 파일이 없습니다.")

print(f"📂 Found {len(smi_files)} SMILES files")

# ===============================
# 변환 실행
# ===============================
for smi_file in smi_files:
    input_path = os.path.join(BASE_DIR, smi_file)
    file_prefix = os.path.splitext(smi_file)[0]

    file_out_dir = os.path.join(OUTPUT_DIR, file_prefix)
    os.makedirs(file_out_dir, exist_ok=True)

    print(f"\n▶ Processing {smi_file}")

    with open(input_path, "r", encoding="utf-8") as f:
        for idx, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            smiles = parts[0]
            name = parts[1] if len(parts) > 1 else f"mol_{idx}"

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"  [SKIP] Invalid SMILES: {smiles}")
                continue

            out_path = os.path.join(file_out_dir, f"{name}.png")
            Draw.MolToFile(mol, out_path, size=IMG_SIZE)

    print(f"  ✅ Saved to: {file_out_dir}")

print("\n🎉 All SMILES converted to images")
