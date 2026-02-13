import pandas as pd

# 1. 파일 읽기
lnp_db = pd.read_excel('LNPDB_1.xlsx')
smiles_db = pd.read_csv('compounds_with_smiles.csv')
table1 = pd.read_excel('Table1.xlsx')
table2 = pd.read_excel('Table2.xlsx')

# --- 데이터 전처리 함수: 숫자면 추출하고, MC3 같은 문자열은 그대로 유지 ---
def extract_cpd_id(series):
    extracted = series.astype(str).str.extract('(\d+)')[0]
    return extracted.fillna(series.astype(str))

# 공통 키 생성
smiles_db['cpd_key'] = extract_cpd_id(smiles_db['Compound'])
lnp_db['cpd_key'] = extract_cpd_id(lnp_db['IL_name'])
table1['cpd_key'] = extract_cpd_id(table1['화합물 (Compound)'])
table2['cpd_key'] = extract_cpd_id(table2['화합물 (Compound)'])

# 2. SMILES 추가 (LNPDB에 원래 SMILES 컬럼이 있다면 덮어쓰기 위해 컬럼 제거 후 병합)
if 'SMILES' in lnp_db.columns:
    lnp_db = lnp_db.drop(columns=['IL_SMILES'])
lnp_db = pd.merge(lnp_db, smiles_db[['cpd_key', 'SMILES']], on='cpd_key', how='left')

# 3. Table 1 데이터 추가 (물리화학적 특성)
# 기존 LNPDB에 해당 컬럼들이 있다면 중복을 피하기 위해 제거 후 병합
cols_to_add = ["크기 (Size, nm)", "다분산 지수 (PDI)", "봉입 효율 (EE, %)", "pKa"]
lnp_db = lnp_db.drop(columns=[c for c in cols_to_add if c in lnp_db.columns])
lnp_db = pd.merge(lnp_db, table1[['cpd_key'] + cols_to_add], on='cpd_key', how='left')
# 병합된 'SMILES' 컬럼명을 'IL_SMILES'로 변경
lnp_db = lnp_db.rename(columns={'SMILES': 'IL_SMILES'})

# 4. Table 2 실험 데이터 통합 및 Row 확장 (Melt)
table2_melted = table2.melt(
    id_vars=['cpd_key'],
    value_vars=['3시간 후 (Total Flux)', '6시간 후 (Total Flux)', '24시간 후 (Total Flux)'],
    var_name='time_point',
    value_name='New_Experiment_value'
)

# New_Experiment_method 생성 (n=6 mice, 0.5 mg/kg i.v.)
table2_melted['New_Experiment_method'] = (
    table2_melted['time_point'] +
    ", Intravenous administration to mice (n=6), 0.5 mg/kg, Luciferase expression (Total Flux)"
)

# 5. 최종 결합 및 컬럼 치환
# 기존의 빈 Experiment_value, Experiment_method 컬럼을 제거하고 새로 만든 값으로 대체
final_db = pd.merge(lnp_db, table2_melted, on='cpd_key', how='inner')

# 기존 컬럼이 있다면 제거하고 새로운 값으로 이름 변경
for old_col, new_col in [('Experiment_value', 'New_Experiment_value'),
                         ('Experiment_method', 'New_Experiment_method')]:
    if old_col in final_db.columns:
        final_db = final_db.drop(columns=[old_col])
    final_db = final_db.rename(columns={new_col: old_col})

# 불필요한 보조 컬럼 정리
final_db = final_db.drop(columns=['cpd_key', 'time_point'])

# 결과 저장
final_db.to_csv('LNP_Integrated_Database.csv', index=False, encoding='utf-8-sig')

print("데이터 통합 완료: 'LNP_Integrated_Database.csv' 파일이 생성되었습니다.")