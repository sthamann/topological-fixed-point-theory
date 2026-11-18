"""Run only the three decision-bearing checks, sequentially per fixed output."""
from pathlib import Path
import hashlib,json,subprocess,sys
r=Path(__file__).resolve().parent
checks=[('check_covariant_rigidity.py','covariant-rigidity-results.json',[]),('admissible_pair_source.py','result.json',['--out',str(r/'result.json')]),('check_source_clock_wedge_cycle.py','source-clock-wedge-cycle.json',[])]
rows=[]
for checker,result,args in checks:
 old=None; old_stdout=None
 for flags in [[],['-OO']]:
  p=subprocess.run([sys.executable,*flags,str(r/checker),*args],capture_output=True)
  if p.returncode: raise RuntimeError(checker+' failed: '+p.stderr.decode())
  data=(r/result).read_bytes()
  if old is not None and (old!=data or old_stdout!=p.stdout): raise RuntimeError(checker+' optimized mismatch')
  old=data;old_stdout=p.stdout
 rows.append({'checker':checker,'normal_optimized_equal':True,'checker_sha256':hashlib.sha256((r/checker).read_bytes()).hexdigest(),'result':result,'result_sha256':hashlib.sha256(old).hexdigest(),'stdout_sha256':hashlib.sha256(old_stdout).hexdigest()})
(r/'VERIFICATION.json').write_text(json.dumps({'scope':'Three decisive checks only; not a full publication or physical TFPT closure.','checks':rows},indent=2)+'\n')
print(json.dumps(rows,indent=2))
