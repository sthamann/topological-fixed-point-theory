"""Submit a reviewed, revision-bound record through the canonical TFPT Atlas."""
from pathlib import Path
import argparse
import json
import subprocess

from verify import digest


def atlas_call(root,tool,arguments):
    command = [str(root/"semantic_explorer/run.sh"),"mcp"]
    requests = [
        {"jsonrpc":"2.0","id":1,"method":"initialize",
         "params":{"protocolVersion":"2024-11-05","capabilities":{},
                   "clientInfo":{"name":"OriginLab","version":"1"}}},
        {"jsonrpc":"2.0","method":"notifications/initialized"},
        {"jsonrpc":"2.0","id":2,"method":"tools/call",
         "params":{"name":tool,"arguments":arguments}},
    ]
    output = subprocess.run(command,cwd=root,text=True,check=True,capture_output=True,
                            timeout=90,input="".join(json.dumps(x)+"\n" for x in requests)).stdout
    for line in output.splitlines():
        message = json.loads(line)
        if message.get("id") == 2:
            if "error" in message:
                raise RuntimeError(message["error"]["message"])
            result = message["result"]
            if result.get("isError"):
                raise RuntimeError(str(result["content"]))
            return json.loads(result["content"][0]["text"])
    raise RuntimeError("Atlas did not return the requested result")


def publish(root,review_path):
    review = json.loads(review_path.read_text())
    if set(review) != {"reviewer","outbox_sha256","outbox","update"} or not review["reviewer"].strip():
        raise ValueError("a source reviewer and an exact outbox binding are required")
    outbox = (review_path.parent/review["outbox"]).resolve()
    if digest(outbox)!=review["outbox_sha256"]:
        raise ValueError("reviewed outbox has changed")
    draft = json.loads(outbox.read_text())
    for name,sha in draft["candidate_artifacts"].items():
        if name not in {"results","summary","frontier","engine","sources"} or digest(outbox.with_name(name+".json"))!=sha:
            raise ValueError("reviewed campaign artifacts have changed")
    update = review["update"]
    if update["expected_revision"]!=draft["atlas_revision"]:
        raise ValueError("review and campaign bind different Atlas revisions")
    return atlas_call(root,"record_research_update",update)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",type=Path,default=Path(__file__).resolve().parents[2])
    parser.add_argument("--review",type=Path)
    parser.add_argument("--status",action="store_true")
    args = parser.parse_args()
    if args.status == bool(args.review):
        parser.error("choose --status or supply --review")
    result = atlas_call(args.repo,"graph_status",{})["meta"]["research_program"] if args.status else publish(args.repo,args.review)
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
