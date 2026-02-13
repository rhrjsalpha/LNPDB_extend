from rdkit import Chem
from rdkit.Chem import rdChemReactions
from itertools import product

def attach(core_smiles, frag_smiles, mapnum=1):
    """
    core_smiles: 코어(결합점 [*:mapnum] 포함)
    frag_smiles: 치환기(결합점 [*:mapnum] 포함)
    """
    core = Chem.MolFromSmiles(core_smiles)
    frag = Chem.MolFromSmiles(frag_smiles)
    if core is None or frag is None:
        return None

    # core 쪽: [X]-[*:n] 형태, frag 쪽: [*:n]-[Y] 형태라고 가정하고 결합
    rxn = rdChemReactions.ReactionFromSmarts(
        f"[*:a]-[*:{mapnum}].[*:{mapnum}]-[*:b]>>[*:a]-[*:b]"
    )
    ps = rxn.RunReactants((core, frag))
    if not ps:
        return None

    # 첫 생성물 사용 (여러 개면 규칙으로 필터링 가능)
    prod = ps[0][0]
    Chem.SanitizeMol(prod)
    return Chem.MolToSmiles(prod)

# -----------------------------
# 1) 코어에 결합점 넣기 (예시)
#   실제로는 그림의 (I) 골격을 SMILES로 만든 뒤,
#   R1 자리와 R2 자리에 [*:1], [*:2]를 둡니다.
# -----------------------------
core = "OCC(OC(=O)[*:2])CO[*:1]"   # (예시용) R1, R2 자리에 더미

# 2) R1 후보들 (결합점 [*:1] 포함)
R1_list = [
    "[*:1]C(=O)CCCC",              # 예시: acyl
    "[*:1]C(=O)CCCCCCCCC",         # 더 긴 acyl
]

# 3) R2 후보들 (결합점 [*:2] 포함)
R2_list = [
    "[*:2]CCN1CCN(C)CC1",           # morpholine-like 예시(실제는 구조에 맞게)
    "[*:2]CCN1CCCC1",               # pyrrolidine-like
]

# 4) 전개
products = []
for r1, r2 in product(R1_list, R2_list):
    step1 = attach(core, r1, mapnum=1)
    if step1 is None:
        continue
    step2 = attach(step1, r2, mapnum=2)
    if step2 is None:
        continue
    products.append(step2)

print("N products:", len(products))
print(products[:5])
