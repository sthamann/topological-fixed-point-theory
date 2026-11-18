"""Sixteen actual whole-program examples from the first finite frontier."""
from pathlib import Path
import hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from catalogue import inputs,operations
from completion import OriginalCompletion
from construction_batches import make_batch,batch_summary
from search.synthesis import Synthesizer,ast_hash
from source_recovery.native_constructors import NativeInputMissing
repo=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
frontier=make_batch()
requested=['action_even_degree_2','action_even_degree_4','action_even_degree_6','action_even_degree_8',
 'action_UY_even_degree_4','action_UY_even_degree_8','kernel_unmarked_degree_2','kernel_UY_degree_8',
 'two_even_vertices_2_4','two_even_vertices_4_6','two_even_vertices_6_8',
 'ordered_AB','ordered_BA','reflected_A','sandwiched_B','reflected_AB']
examples=[next(r for r in frontier if r['structure_family']==name)for name in requested]
validator=OriginalCompletion(repo,'tfpt_complete')
deps=[ROOT/'construction_batches.py',ROOT/'catalogue.py',ROOT/'completion.py',
 ROOT/'source_recovery/native_constructors.py',validator.q_path,
 validator.q_path.parents[1]/'checks/original_checks.py',validator.q_path.parent/'original_source_adapter.py']
identity=hashlib.sha256(b''.join(p.read_bytes()for p in deps)).hexdigest()
engine=Synthesizer(operations(),inputs(),validator,validator.required_obligation_ids,validator.target_scope,
 identity,ROOT/'review/construction-batch-native-run',verifier_dependencies=deps,
 missing_input_exceptions=(NativeInputMissing,))
# Static checks cover the entire finite frontier; execution below covers16.
all_types=all(engine.type_of(r['program'])=='WholeSourceCandidate'for r in frontier)
rows=[]
for index,row in enumerate(examples):
 start=time.monotonic();record=engine.evaluate_program(row['program'],max_cost=64)
 value=engine.values.get(ast_hash(row['program']))
 goal=record.get('goal',{});obligations={x['id']:x['status']for x in goal.get('obligations',[])}
 values={k:str(value[k])for k in ('Z0','ZU','ZY','ZUY')}if isinstance(value,dict)and all(k in value for k in ('Z0','ZU','ZY','ZUY'))else None
 out={'id':row['id'],'structure_family':row['structure_family'],'ast_cost':row['ast_cost'],
 'evaluation_status':record['evaluation_status'],'constructed_common_values':obligations.get('constructed_common_values'),
 'target_status':goal.get('status'),'values':values,'connected_UY':str(value.get('normalized',{}).get('connected_UY'))if isinstance(value,dict)else None,
 'hypothesis_count':len(value.get('hypotheses',[]))if isinstance(value,dict)else None,
 'original_missing_obligations':[k for k,v in obligations.items()if v=='MISSING_INPUT'],
 'elapsed_seconds':time.monotonic()-start}
 rows.append(out); print(json.dumps({'example':index+1,**out}),flush=True)
result={'scope':'16 executed whole candidate programs on the stated hypothetical charts; no primitive physical source claim',
 'frontier':batch_summary(frontier),'all_256_ASTs_typechecked':all_types,
 'executed_examples':len(rows),'all_16_values_built_and_recomputed':all(r['evaluation_status']=='VALUE_BUILT'and r['constructed_common_values']=='PASS'and r['values']is not None for r in rows),
 'original_completion_claimed':False,'rows':rows,
 'source_sha256':{str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in deps}}
(ROOT/'review/construction-batch-native-controls.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items()if k not in ('rows','source_sha256')},indent=2),flush=True)
if not result['all_16_values_built_and_recomputed']:raise SystemExit(1)
