import re
import requests
import pandas as pd
import time

# 설정: 입력 파일명과 출력 CSV 파일명
input_file_name = 'US_20170210697_A1_263_Cpd_61.txt'
output_csv_name = 'compounds_with_smiles.csv'


def get_smiles_from_cir(compound_name):
    """CIR API를 사용하여 화합물 이름으로부터 SMILES를 가져옵니다."""
    try:
        # CIR API URL (화합물 이름을 SMILES로 변환)
        url = f"https://cactus.nci.nih.gov/chemical/structure/{compound_name}/smiles"
        response = requests.get(url, timeout=10)

        if response.status_code == 200:
            return response.text.strip()
        else:
            return "Not Found"
    except Exception as e:
        print(f"Error fetching SMILES for {compound_name[:30]}...: {e}")
        return "Error"


def process_compounds_to_csv(file_path):
    # 'Compound 숫자: IUPAC이름' 패턴 추출
    # 예: Compound 169: Heptadecan-9-yl ...
    pattern = re.compile(r'Compound\s+(\d+):\s+(.*?)(?=\n|\[|\Z)', re.DOTALL)

    data_list = []

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
            matches = pattern.findall(content)

            print(f"총 {len(matches)}개의 화합물을 발견했습니다. SMILES 변환을 시작합니다...")

            for i, (cpd_num, iupac_name) in enumerate(matches):
                # 데이터 정제
                cpd_label = f"Compound {cpd_num}"
                clean_iupac = " ".join(iupac_name.split())

                # CIR API 호출 (서버 부하를 고려하여 약간의 지연을 둘 수 있음)
                print(f"[{i + 1}/{len(matches)}] {cpd_label} 처리 중...")
                smiles = get_smiles_from_cir(clean_iupac)

                data_list.append({
                    "Compound": cpd_label,
                    "IUPAC": clean_iupac,
                    "SMILES": smiles
                })

                # API 매너: 연속 호출 시 지연 시간 추가 (필요 시)
                time.sleep(0.5)

    except FileNotFoundError:
        print(f"오류: {file_path} 파일을 찾을 수 없습니다.")
        return

    # 데이터프레임 생성 및 CSV 저장
    if data_list:
        df = pd.DataFrame(data_list)
        df.to_csv(output_csv_name, index=False, encoding='utf-8-sig')
        print(f"\n저장 완료! '{output_csv_name}' 파일을 확인하세요.")
        print(df.head())
    else:
        print("추출된 데이터가 없습니다.")


if __name__ == "__main__":
    process_compounds_to_csv(input_file_name)