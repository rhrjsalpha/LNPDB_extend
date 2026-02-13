import pandas as pd
import os
import re
from rdkit import Chem
from collections import Counter

def get_canonical_smiles(smiles):
    if pd.isna(smiles) or str(smiles).strip() == "": return None
    try:
        mol = Chem.MolFromSmiles(str(smiles).strip())
        return Chem.MolToSmiles(mol, canonical=True) if mol else None
    except:
        return None

def get_inchi_from_smiles(smiles):
    if pd.isna(smiles) or str(smiles).strip() == "": return None
    try:
        mol = Chem.MolFromSmiles(str(smiles).strip())
        return Chem.MolToInchi(mol) if mol else None
    except:
        return None

def extract_number(name):
    if pd.isna(name): return None
    match = re.search(r'\d+', str(name))
    return match.group() if match else None

def process_chemical_data(file_path):
    base_dir = os.path.dirname(os.path.abspath(file_path))
    output_dir = os.path.join(base_dir, "analysis_results")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df = pd.read_csv(file_path)
    df['InChI'] = df['SMILES'].apply(get_inchi_from_smiles)
    
    # 1. 비교를 위한 데이터 분리 (Compound vs XX)
    # 이름에 'Compound'가 포함된 행과 'XX'가 포함된 행을 구분
    df['Name_Upper'] = df['Chemical_Name'].fillna('').str.upper()
    df['Num_ID'] = df['Chemical_Name'].apply(extract_number)

    comp_df = df[df['Name_Upper'].str.contains('COMPOUND')].dropna(subset=['Num_ID'])
    xx_df = df[df['Name_Upper'].str.contains('XX')].dropna(subset=['Num_ID'])

    # 각 ID별로 가장 빈번한 InChI 하나만 선택 (대표값)
    comp_ref = comp_df.groupby('Num_ID')['InChI'].agg(lambda x: Counter(x.dropna()).most_common(1)[0][0] if not x.dropna().empty else None).reset_index()
    xx_ref = xx_df.groupby('Num_ID')['InChI'].agg(lambda x: Counter(x.dropna()).most_common(1)[0][0] if not x.dropna().empty else None).reset_index()

    # 2. 숫자(ID) 기준으로 머지
    comparison = pd.merge(comp_ref, xx_ref, on='Num_ID', suffixes=('_comp', '_xx'), how='inner')

    # 3. 일치 여부 판별 (True/False)
    # 둘 다 None이 아니고 값이 같을 때 True
    comparison['Match_Result'] = (comparison['InChI_comp'] == comparison['InChI_xx']) & (comparison['InChI_comp'].notna())

    # 4. 컬럼 순서 재배치 및 이름 변경
    # 요청: compound number - compound inchi - xx inchi - xx number - true/false
    comparison['xx_number'] = comparison['Num_ID']
    comparison = comparison.rename(columns={'Num_ID': 'compound_number', 'InChI_comp': 'compound_inchi', 'InChI_xx': 'xx_inchi', 'Match_Result': 'Is_Identical'})
    
    final_comparison = comparison[['compound_number', 'compound_inchi', 'xx_inchi', 'xx_number', 'Is_Identical']]

    # 5. 파일 저장
    final_comparison.to_csv(os.path.join(output_dir, "compound_xx_comparison.csv"), index=False, encoding='utf-8-sig')

    # 기존 랭킹 코드 등은 그대로 유지 (생략 가능하나 구조상 유지)
    print(f"비교 분석 완료! 파일이 저장되었습니다:\n{os.path.join(output_dir, 'compound_xx_comparison.csv')}")

# 실행
target_file = '/Users/kogeon/python_projects_path/LNPDB_extend/PatientToLNPDB/Text_Extract/260209_001/US20170210697A1/03.Analysis_Results/US20170210697A1_Integrated_Analysis.csv'
process_chemical_data(target_file)