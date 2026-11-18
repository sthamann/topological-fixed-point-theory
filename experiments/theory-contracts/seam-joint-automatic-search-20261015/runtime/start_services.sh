#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PYTHON=${TFPT_VLLM_PYTHON:-/home/aime/vllm/.venv/bin/python}
mkdir -p "$ROOT/logs" "$ROOT/pids" "$ROOT/compile-cache"
test -f "$ROOT/model/model.safetensors.index.json"
for spec in '1:18931' '2:18932'; do
  gpu=${spec%:*}
  port=${spec#*:}
  gpu_uuid=$(nvidia-smi -i "$gpu" --query-gpu=uuid --format=csv,noheader)
  gpu_name=$(nvidia-smi -i "$gpu" --query-gpu=name --format=csv,noheader)
  [[ "$gpu_name" == *H200* ]] || { echo "GPU $gpu is not an H200" >&2; exit 1; }
  pidfile="$ROOT/pids/vllm-$port.pid"
  if [[ -f "$pidfile" ]] && kill -0 "$(cat "$pidfile")" 2>/dev/null; then
    echo "Existing dedicated service $port PID $(cat "$pidfile")"
    continue
  fi
  if ss -ltnH | awk '{print $4}' | grep -Eq ":$port$"; then
    echo "Port $port already occupied; refusing replacement" >&2
    exit 1
  fi
  nohup env CUDA_VISIBLE_DEVICES="$gpu_uuid" CUDA_DEVICE_ORDER=PCI_BUS_ID \
    HF_HOME="$ROOT/hf-home" XDG_CACHE_HOME="$ROOT/cache/gpu-$gpu" \
    TORCH_HOME="$ROOT/compile-cache/gpu-$gpu/torch" \
    TORCHINDUCTOR_CACHE_DIR="$ROOT/compile-cache/gpu-$gpu/inductor" \
    TRITON_CACHE_DIR="$ROOT/compile-cache/gpu-$gpu/triton" \
    VLLM_CACHE_ROOT="$ROOT/compile-cache/gpu-$gpu" \
    "$PYTHON" -m vllm.entrypoints.openai.api_server \
      --model "$ROOT/model" --served-model-name tfpt-qwen3-32b \
      --host 127.0.0.1 --port "$port" --dtype bfloat16 \
      --max-model-len 32768 --gpu-memory-utilization 0.88 \
      --max-num-seqs 4 --max-num-batched-tokens 8192 \
      --reasoning-parser qwen3 --disable-log-requests \
      > "$ROOT/logs/vllm-$port.log" 2>&1 < /dev/null &
  echo "$!" > "$pidfile"
  echo "Started H200 GPU $gpu, loopback port $port, PID $!"
done
