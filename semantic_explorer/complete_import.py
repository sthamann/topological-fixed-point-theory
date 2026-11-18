#!/usr/bin/env python3
"""Continue one existing import, then retry its errors with current extractor.

Run at most two repair passes; preserve successes and don't start a duplicate writer.
No scheduler, daemon or modification of research sources.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser();p.add_argument('--key-file',required=True);p.add_argument('--workers',type=int,default=16)
    a=p.parse_args();lock=ROOT/'data/.extract.lock'
    while lock.exists():
        pid=int(lock.read_text())
        try:os.kill(pid,0)
        except ProcessLookupError:raise SystemExit('Stale import lock; inspect before removing it')
        time.sleep(10)
    for attempt in range(2):
        cmd=[sys.executable,str(ROOT/'extract.py'),'--key-file',a.key_file,'--workers',str(a.workers),'--resume-errors','--embeddings']
        result=subprocess.run(cmd)
        if result.returncode:raise SystemExit(result.returncode)
        subprocess.run([sys.executable,str(ROOT/'build_context.py')],check=True)
        g=json.loads((ROOT/'data/knowledge.json').read_text());coverage=g['meta']['coverage']
        print(json.dumps({'repair_pass':attempt+1,'coverage':coverage}),flush=True)
        if coverage['claims_pending']==0 and coverage['claims_error']==0:return
    raise SystemExit('Import retained remaining errors explicitly; inspect data/knowledge.json meta.errors')

if __name__=='__main__':main()
