import os
import fitz  # pymupdf

def pdf_to_images(pdf_path: str, out_dir: str = "pdf_pages", dpi: int = 300):
    os.makedirs(out_dir, exist_ok=True)
    doc = fitz.open(pdf_path)
    base = os.path.splitext(os.path.basename(pdf_path))[0]

    img_paths = []
    for i, page in enumerate(doc):
        pix = page.get_pixmap(dpi=dpi)
        img_path = os.path.join(out_dir, f"{base}_page_{i:03d}.png")
        pix.save(img_path)
        img_paths.append(img_path)

    doc.close()
    return img_paths

import os
import torch
from transformers import AutoModel, AutoTokenizer

MODEL_ID = "deepseek-ai/DeepSeek-OCR"

def load_deepseek_ocr():
    device = "cuda" if torch.cuda.is_available() else "cpu"

    tokenizer = AutoTokenizer.from_pretrained(MODEL_ID, trust_remote_code=True)

    # Windows에서는 flash-attn 설치가 어려울 수 있어 sdpa/eager 권장
    model = AutoModel.from_pretrained(
        MODEL_ID,
        trust_remote_code=True,
        use_safetensors=True,
        _attn_implementation="sdpa",   # "eager"도 가능 / flash_attention_2는 flash-attn 설치 필요
    )

    model = model.eval()

    if device == "cuda":
        # 모델 카드 기준 BF16 권장 :contentReference[oaicite:3]{index=3}
        model = model.cuda().to(torch.bfloat16)
    else:
        model = model.to(torch.float32)

    return model, tokenizer

def ocr_one_image(model, tokenizer, image_path: str, output_dir: str = "ocr_out"):
    os.makedirs(output_dir, exist_ok=True)

    # 모델 카드 예시 프롬프트 :contentReference[oaicite:4]{index=4}
    prompt = "<image>\n<|grounding|>Convert the document to markdown. "

    res = model.infer(
        tokenizer,
        prompt=prompt,
        image_file=image_path,
        output_path=output_dir,
        base_size=1024,
        image_size=640,
        crop_mode=True,
        save_results=True,
        test_compress=True,
    )
    return res

def ocr_pdf(pdf_path: str, pages_dir="pdf_pages", out_dir="ocr_out"):
    img_paths = pdf_to_images(pdf_path, out_dir=pages_dir, dpi=300)

    model, tokenizer = load_deepseek_ocr()

    for i, img in enumerate(img_paths):
        print(f"OCR {i+1}/{len(img_paths)}: {img}")
        ocr_one_image(model, tokenizer, img, output_dir=out_dir)

if __name__ == "__main__":
    ocr_pdf("MP_2025.pdf", pages_dir="MP_2025_pages", out_dir="MP_2025_ocr")