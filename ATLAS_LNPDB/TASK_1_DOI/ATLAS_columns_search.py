import pandas as pd
from collections import Counter, defaultdict


# ============================================================
# 1. CSV 로드
# ============================================================
df = pd.read_csv(
    "../lnp_atlas_export_20260124_0129.csv",
    encoding="cp949"
)

bio_series = df["bioactivity_profile"].dropna()

# ============================================================
# 2. key 조합 추출 (순서 유지)
# ============================================================
combo_counter = Counter()
combo_order_repr = {}   # frozenset -> 원본 순서 리스트

for text in bio_series:
    keys_in_order = []
    key_set = set()

    for part in text.split(";"):
        part = part.strip()
        if not part or ":" not in part:
            continue

        key = part.split(":", 1)[0].strip()

        if key not in key_set:
            keys_in_order.append(key)
            key_set.add(key)

    if key_set:
        combo = frozenset(key_set)
        combo_counter[combo] += 1

        # 최초 등장한 조합의 순서를 대표로 저장
        if combo not in combo_order_repr:
            combo_order_repr[combo] = keys_in_order

# ============================================================
# 3. 출력 (정렬 ❌, 원본 순서 유지)
# ============================================================
print(f"\n총 bioactivity_profile row 수: {len(bio_series)}")
print(f"서로 다른 key 조합 수: {len(combo_counter)}\n")

for combo, count in combo_counter.most_common():
    keys_ordered = combo_order_repr[combo]
    combo_str = ", ".join(keys_ordered)
    print(f"[{count:4d} rows]  {combo_str}")
