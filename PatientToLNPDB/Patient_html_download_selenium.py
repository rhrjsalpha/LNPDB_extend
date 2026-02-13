import os
import pandas as pd
import time
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager

# 1. 브라우저 설정 (창이 보이도록 설정)
options = webdriver.ChromeOptions()
# options.add_argument('--headless')  # 필요 시 주석 해제
options.add_argument('--disable-gpu')
options.add_argument('--no-sandbox')
options.add_argument('--window-size=1200,1000')
options.add_argument(
    "user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/119.0.0.0 Safari/537.36")

driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=options)

# 2. 경로 및 파일 설정
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_FILE = os.path.join(BASE_DIR, "20250409_LipidLibrarySample.xlsx")
OUTPUT_FOLDER = os.path.join(BASE_DIR, "EXTRACTED_TEXT")

if not os.path.exists(OUTPUT_FOLDER):
    os.makedirs(OUTPUT_FOLDER)
    print(f"📂 폴더 생성 완료: {OUTPUT_FOLDER}")

# 3. 데이터 로드
try:
    if INPUT_FILE.endswith('.csv'):
        df = pd.read_csv(INPUT_FILE)
    else:
        df = pd.read_excel(INPUT_FILE)
    print(f"📊 데이터 로드 완료: {len(df)}개의 링크를 처리합니다.")
except Exception as e:
    print(f"❌ 파일을 읽을 수 없습니다: {e}")
    driver.quit()
    exit()

# 4. 추출할 XPath 리스트 설정
# XPATH1: 이전에 요청하신 본문/설명 영역
# XPATH2: 이번에 새로 추가하신 섹션 영역
XPATH1 = "/html/body/search-app/search-result/search-ui/div/div/div/div/div/result-container/patent-result/div/div/div/div[2]"
XPATH2 = "/html/body/search-app/search-result/search-ui/div/div/div/div/div/result-container/patent-result/div/div/div/div[1]/div[2]/section"

for index, row in df.iterrows():
    patent_id = str(row['PATENT_ID']).strip()
    url = str(row['Link']).strip()

    if not url.startswith("http"):
        continue

    print(f"🌐 [{index + 1}/{len(df)}] {patent_id} 접속 및 추출 중...")

    try:
        driver.get(url)
        wait = WebDriverWait(driver, 10)  # 각 요소를 기다리는 시간

        extracted_texts = []

        # 영역 1 추출
        try:
            el1 = wait.until(EC.presence_of_element_located((By.XPATH, XPATH1)))
            extracted_texts.append(el1.text)
        except:
            print(f"  - 영역 1을 찾을 수 없음")

        # 영역 2 추출 (추가된 부분)
        try:
            el2 = driver.find_element(By.XPATH, XPATH2)
            extracted_texts.append(el2.text)
        except:
            print(f"  - 영역 2를 찾을 수 없음")

        # 두 영역의 텍스트를 합치기 (구분선 포함)
        final_content = "\n\n" + "=" * 50 + "\n[SECTION 1]\n" + "=" * 50 + "\n\n"
        final_content += extracted_texts[0] if len(extracted_texts) > 0 else "N/A"

        if len(extracted_texts) > 1:
            final_content += "\n\n" + "=" * 50 + "\n[SECTION 2]\n" + "=" * 50 + "\n\n"
            final_content += extracted_texts[1]

        # 파일 저장
        output_filename = f"{patent_id}.txt".replace("/", "_").replace("\\", "_")
        output_path = os.path.join(OUTPUT_FOLDER, output_filename)

        with open(output_path, "w", encoding="utf-8") as f:
            f.write(final_content)

        print(f"✅ 저장 완료: {output_filename}")
        time.sleep(1.5)  # 부하 방지를 위한 대기

    except Exception as e:
        print(f"❌ {patent_id} 처리 실패: {e}")

# 5. 종료
driver.quit()
print("\n✨ 모든 작업이 완료되었습니다.")