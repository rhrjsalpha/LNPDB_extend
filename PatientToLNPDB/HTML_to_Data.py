from bs4 import BeautifulSoup
import json

HTML_PATH = "downloaded_html/patents.google.com_patent_US20170210697A1_en.html"   # 경로 맞게 수정

def clean(x):
    return " ".join(x.split())

with open(HTML_PATH, encoding="utf-8") as f:
    soup = BeautifulSoup(f, "html.parser")

# ==================================================
# 1. Info 범위 결정 (h2 Info ~ 다음 h2)
# ==================================================
info_h2 = soup.find("h2", string="Info")
if not info_h2:
    raise RuntimeError("Info section not found")

info_nodes = []
for el in info_h2.find_all_next():
    if el.name == "h2" and el is not info_h2:
        break
    info_nodes.append(el)

# ==================================================
# 2. dt / dd 순차 파싱
# ==================================================
info = {}
i = 0

while i < len(info_nodes):
    el = info_nodes[i]

    if el.name == "dt":
        key = clean(el.get_text())
        values = []
        j = i + 1

        while j < len(info_nodes) and info_nodes[j].name == "dd":
            dd = info_nodes[j]

            # ---- time 태그 ----
            time_tag = dd.find("time")
            if time_tag and time_tag.has_attr("datetime"):
                values.append(time_tag["datetime"])

            # ---- 일반 itemprop ----
            elif dd.has_attr("itemprop") and dd["itemprop"] != "events":
                values.append(clean(dd.get_text()))

            j += 1

        if values:
            info[key] = values if len(values) > 1 else values[0]

        i = j
    else:
        i += 1

# ==================================================
# 3. Events는 전역 수집 (핵심)
# ==================================================
events = []

for dd in soup.find_all("dd", itemprop="events"):
    event = {}

    t = dd.find("time", itemprop="date")
    if t:
        event["date"] = t.get("datetime", clean(t.get_text()))

    for span in dd.find_all("span"):
        prop = span.get("itemprop")
        if not prop:
            continue
        event.setdefault(prop, [])
        event[prop].append(clean(span.get_text()))

    # list → 단일값 정리
    for k, v in event.items():
        if len(v) == 1:
            event[k] = v[0]

    events.append(event)

info["Events"] = events

# ==================================================
# 4. 출력
# ==================================================
print(json.dumps(info, indent=2, ensure_ascii=False))


