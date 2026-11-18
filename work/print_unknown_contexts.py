import ast,json,sys
from pathlib import Path
repo=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4');sys.path.insert(0,str(repo/'semantic_explorer'));import extract
corp=extract.Corpus(repo); D=json.load(open('/Users/stefanhamann/Documents/Codex/2026-10-13/wie-x20/work/unknown89-baseline.json'))
for j,g in enumerate(D['groups'],1):
 p=repo/g['script_path'];text,_=corp.read(p);ls=text.splitlines();tree=ast.parse(text);line=g['candidates'][0]['line']
 scopes=[n for n in ast.walk(tree) if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef)) and n.lineno<=line<=n.end_lineno]
 owner=min(scopes,key=lambda n:n.end_lineno-n.lineno) if scopes else None
 if owner: a=max(owner.lineno,line-4);b=min(owner.end_lineno,line+4); name=owner.name
 else:a=max(1,line-4);b=min(len(ls),line+4);name='<module>'
 print(f"\n===== {j:02d} {g['claim_id']} {g['script_path']}:{line} scope={name} =====")
 print('\n'.join(f'{k}: {ls[k-1]}' for k in range(a,b+1)))
