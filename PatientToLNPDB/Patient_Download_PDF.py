import requests

patent_id = "WO2021021634A1"
pdf_url = f"https://patents.google.com/patent/{patent_id}.pdf"

r = requests.get(pdf_url, timeout=30)

if r.status_code == 200 and r.headers["Content-Type"].startswith("application/pdf"):
    with open(f"{patent_id}.pdf", "wb") as f:
        f.write(r.content)
    print("PDF 다운로드 성공")
else:
    print("PDF 없음 또는 접근 차단")