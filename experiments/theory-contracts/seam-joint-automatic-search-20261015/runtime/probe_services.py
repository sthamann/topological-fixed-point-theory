#!/usr/bin/env python3
"""Run a small real inference on each dedicated endpoint and save API receipts."""
import concurrent.futures
import json
import pathlib
import time
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parent

def probe(port):
    base = f"http://127.0.0.1:{port}"
    with urllib.request.urlopen(base + "/v1/models", timeout=10) as response:
        models = json.load(response)
    payload = {"model": "tfpt-qwen3-32b", "messages": [
        {"role": "user", "content": "Let A=diag(1,-1). Can an invertible complex 2 by 2 matrix U satisfy U A U^{-1}=2 A? Give the decisive invariant and answer in under 80 words."}],
        "max_tokens": 2048, "temperature": 0.6, "top_p": 0.95,
        "top_k": 20, "chat_template_kwargs": {"enable_thinking": True}}
    req = urllib.request.Request(base + "/v1/chat/completions",
        data=json.dumps(payload).encode(), headers={"Content-Type": "application/json"})
    started = time.monotonic()
    with urllib.request.urlopen(req, timeout=300) as response:
        result = json.load(response)
    elapsed = time.monotonic() - started
    message = result["choices"][0]["message"]
    assert message.get("content"), "No final answer; thinking budget may be insufficient"
    receipt = {"endpoint": base, "models": models, "request": payload,
        "response": result, "elapsed_seconds": elapsed,
        "completion_tokens_per_second": result["usage"]["completion_tokens"] / elapsed}
    (ROOT / f"probe-{port}.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return {"port": port, "seconds": round(elapsed, 2), "tokens": result["usage"],
        "tokens_per_second": round(receipt["completion_tokens_per_second"], 2),
        "final_answer": message["content"]}

if __name__ == "__main__":
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for result in pool.map(probe, [18931, 18932]):
            print(json.dumps(result, ensure_ascii=False), flush=True)
