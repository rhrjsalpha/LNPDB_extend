from vllm import LLM, SamplingParams
from vllm.model_executor.models.deepseek_ocr import (
    NGramPerReqLogitsProcessor
)
from PIL import Image

# 1. 모델 로드 (한 번만)
llm = LLM(
    model="deepseek-ai/DeepSeek-OCR",
    enable_prefix_caching=False,
    mm_processor_cache_gb=0,
    logits_processors=[NGramPerReqLogitsProcessor],
    gpu_memory_utilization=0.70,  #
)

# 2. 이미지 로드
image = Image.open(r"C:\Users\kogun\PycharmProjects\LNPDB\pdf_pages\MP_2025_page_000.png").convert("RGB")
image = image.resize((640, 640))

# 3. 프롬프트 (README 권장)
prompt = "<image>\n<|grounding|>Convert the document to markdown."

# 4. 입력 포맷
inputs = [{
    "prompt": prompt,
    "multi_modal_data": {"image": image}
}]

# 5. 샘플링 파라미터 (README 기준)
params = SamplingParams(
    temperature=0.0,
    max_tokens=1024,
    extra_args=dict(
        ngram_size=30,
        window_size=90,
        whitelist_token_ids={128821, 128822},  # <td>, </td>
    ),
    skip_special_tokens=False,
)

# 6. 추론
outputs = llm.generate(inputs, params)

# 7. 결과 출력
print(outputs[0].outputs[0].text)
