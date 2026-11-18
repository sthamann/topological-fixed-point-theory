"""Finite typed program synthesis over an injected, trusted operation registry.

JSON supplies ASTs only, never Python/import/shell commands. Missing source data
are local findings. The injected verifier owns the goal scope and obligations.
"""
from dataclasses import dataclass
from pathlib import Path
import hashlib,inspect,itertools,json,marshal,time,types

class MissingInput(Exception):pass
class SearchBudgetStop(OverflowError):pass
@dataclass(frozen=True)
class Unavailable:
    reason:str
@dataclass(frozen=True)
class InputBinding:
    type_name:str
    value:object
    fingerprint:str
@dataclass(frozen=True)
class Operation:
    id:str
    input_types:tuple
    output_type:str
    execute:object
    version:str
    @classmethod
    def from_value(cls,value):
        if isinstance(value,cls):return value
        return cls(value['id'],tuple(value['input_types']),value['output_type'],value['execute'],value['version'])

def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)
def digest(value):return hashlib.sha256(canonical(value).encode()).hexdigest()
def ast_hash(ast):return digest(ast)
def cost(ast):return 1+sum(cost(a)for a in ast['args'])
def callable_code_hash(fn):
    target=inspect.unwrap(fn)
    code=getattr(target,'__code__',None)
    if code is None:code=getattr(getattr(target,'__call__',None),'__code__',None)
    if code is None:raise ValueError('A trusted executable Python callable with inspectable verifier/operation code is required')
    h=hashlib.sha256(marshal.dumps(code))
    try:
        path=inspect.getsourcefile(target)
        if path and Path(path).is_file():h.update(Path(path).read_bytes())
    except (TypeError,OSError):pass
    return h.hexdigest()
def semantic_fingerprint(value):
    """Hash actual bound data/code, never a truncated display or self-label.

    Stateful verifiers, closures, defaults and their referenced trusted helper
    functions participate. Unsupported opaque objects must be made explicit.
    """
    seen={}
    def form(obj):
        if obj is None or isinstance(obj,(str,bool,int)):return obj
        if isinstance(obj,float):return {'float_hex':obj.hex()}
        if isinstance(obj,complex):return {'complex':[obj.real.hex(),obj.imag.hex()]}
        if isinstance(obj,bytes):return {'bytes_sha256':hashlib.sha256(obj).hexdigest(),'length':len(obj)}
        if isinstance(obj,Path):
            path=obj.resolve();return {'path':str(path),'file_sha256':hashlib.sha256(path.read_bytes()).hexdigest()if path.is_file()else None}
        module=type(obj).__module__
        if module.startswith('sympy'):
            import sympy
            return {'sympy':sympy.srepr(obj)}
        if hasattr(obj,'dtype')and hasattr(obj,'shape')and callable(getattr(obj,'tobytes',None)):
            if getattr(obj.dtype,'hasobject',False):return {'array_shape':list(obj.shape),'object_values':form(obj.tolist())}
            return {'array_shape':list(obj.shape),'dtype':str(obj.dtype),'data_sha256':hashlib.sha256(obj.tobytes(order='C')).hexdigest()}
        if isinstance(obj,types.ModuleType):
            path=getattr(obj,'__file__',None)
            return {'module':obj.__name__,'version':getattr(obj,'__version__',None),
                    'file_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest()if path and Path(path).is_file()else None}
        ident=id(obj)
        if ident in seen:return {'reference':seen[ident]}
        seen[ident]=len(seen)
        if isinstance(obj,dict):
            keys=sorted(obj,key=semantic_fingerprint)
            return {'mapping':[[form(k),form(obj[k])]for k in keys]}
        if isinstance(obj,(list,tuple)):return {'sequence_type':type(obj).__name__,'items':[form(x)for x in obj]}
        if isinstance(obj,(set,frozenset)):return {'set':[form(x) for x in sorted(obj,key=semantic_fingerprint)]}
        if inspect.ismethod(obj):return {'method':form(obj.__func__),'bound_self':form(obj.__self__)}
        if isinstance(obj,type):
            methods={k:callable_code_hash(v)for k,v in vars(obj).items()if inspect.isfunction(v)}
            return {'class':obj.__module__+'.'+obj.__qualname__,'methods':methods}
        if callable(obj):
            target=inspect.unwrap(obj)
            if target is not obj:return {'unwrapped_callable':form(target)}
            if inspect.isfunction(obj):
                globals_used={k:obj.__globals__[k] for k in sorted(set(obj.__code__.co_names)) if k in obj.__globals__ and k!='__builtins__'}
                cells=[]
                for cell in obj.__closure__ or ():
                    try:cells.append(form(cell.cell_contents))
                    except ValueError:cells.append({'empty_closure_cell':True})
                return {'function_code_sha256':callable_code_hash(obj),'closure':cells,
                        'defaults':form(obj.__defaults__),'keyword_defaults':form(obj.__kwdefaults__),
                        'attributes':form(vars(obj)),'referenced_globals':form(globals_used)}
            if inspect.isbuiltin(obj):return {'builtin':obj.__module__+'.'+obj.__qualname__}
            state=getattr(obj,'__dict__',None)
            if state is not None:return {'callable_class':form(type(obj)),'actual_instance_state':form(state),'call_code':callable_code_hash(obj)}
        state=getattr(obj,'__dict__',None)
        if state is not None:return {'object_class':form(type(obj)),'actual_instance_state':form(state)}
        raise TypeError('No safe semantic fingerprint for '+module+'.'+type(obj).__name__+'; expose its actual data explicitly')
    return digest(form(value))
def evidence_json(value):
    if value is None or isinstance(value,(str,bool,int)):return value
    if isinstance(value,float):
        # Do not write NaN/Infinity as valid evidence.
        return value if value==value and abs(value)!=float('inf')else str(value)
    if isinstance(value,complex):return {'complex':[value.real,value.imag]}
    if isinstance(value,dict):return {str(k):evidence_json(v)for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [evidence_json(v)for v in value]
    # Evidence is display data, never deserialized/evaluated as Python.
    return {'python_type':type(value).__module__+'.'+type(value).__name__,'display':str(value)}

class Synthesizer:
    def __init__(self,operations,inputs,validator,required_obligation_ids,target_scope,
                 input_fingerprint,out_dir,target_type='WholeSourceCandidate',
                 verifier_dependencies=(),type_validators=None,context=None,
                 missing_input_exceptions=()):
        ops=[Operation.from_value(o)for o in operations]
        if len({o.id for o in ops})!=len(ops):raise ValueError('Duplicate operation id')
        if not required_obligation_ids or len(set(required_obligation_ids))!=len(required_obligation_ids):raise ValueError('Explicit nonempty unique goal obligations required')
        if not target_scope or not input_fingerprint:raise ValueError('Explicit trusted goal scope and actual input fingerprint required')
        for o in ops:
            if not o.id or not o.output_type or not o.version or not callable(o.execute):raise ValueError('Malformed trusted operation')
        for key,binding in inputs.items():
            if not isinstance(binding,InputBinding)or not binding.fingerprint:raise ValueError('Typed original input binding and fingerprint required: '+key)
        self.ops={o.id:o for o in sorted(ops,key=lambda o:o.id)};self.inputs=dict(inputs)
        self.validator=validator;self.required=tuple(required_obligation_ids);self.scope=target_scope;self.target=target_type
        self.type_validators=type_validators or {};self.context=context;self.out=Path(out_dir);self.out.mkdir(parents=True,exist_ok=True)
        self.missing_exceptions=(MissingInput,*tuple(missing_input_exceptions))
        dependency_hashes={str(Path(p).resolve()):hashlib.sha256(Path(p).read_bytes()).hexdigest()for p in verifier_dependencies}
        self.fingerprint_components={'search_core_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'missing_input_exception_types':[e.__module__+'.'+e.__qualname__ for e in self.missing_exceptions],
            'input_fingerprint':input_fingerprint,'inputs':{k:{'type':v.type_name,'fingerprint':v.fingerprint,'actual_value_sha256':semantic_fingerprint(v.value)}for k,v in sorted(inputs.items())},
            'operations':[{'id':o.id,'input_types':o.input_types,'output_type':o.output_type,'version':o.version,'code_sha256':callable_code_hash(o.execute),'bound_state_sha256':semantic_fingerprint(o.execute)}for o in self.ops.values()],
            'verifier_code_sha256':callable_code_hash(validator),'verifier_dependency_sha256':dependency_hashes,
            'verifier_bound_state_sha256':semantic_fingerprint(validator),'actual_context_sha256':semantic_fingerprint(context),
            'required_obligation_ids':self.required,'target_scope':self.scope,'target_type':self.target,
            'type_validator_code_sha256':{k:callable_code_hash(v)for k,v in sorted(self.type_validators.items())},
            'type_validator_bound_state_sha256':{k:semantic_fingerprint(v)for k,v in sorted(self.type_validators.items())}}
        self.fingerprint=digest(self.fingerprint_components);self.records={};self.values={};self.new_evaluations=0;self.dependency_replays=0;self.queue_offset=0
        self._evaluation_limit=None
        self.archive=self.out/'results.jsonl'
        if self.archive.exists():
            for line in self.archive.read_text().splitlines():
                record=json.loads(line)
                if record['search_fingerprint']==self.fingerprint:self.records[record['ast_hash']]=record
        self.manifest_path=self.out/('manifest-'+self.fingerprint+'.json')
        self._write(self.manifest_path,{'search_fingerprint':self.fingerprint,**self.fingerprint_components})
    def _write(self,path,value):
        temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(evidence_json(value),indent=2,allow_nan=False)+'\n');temp.replace(path)
    def type_of(self,ast):
        if not isinstance(ast,dict)or set(ast)!={'op','id','args'}or not isinstance(ast['id'],str)or not isinstance(ast['args'],list):raise ValueError('AST must contain only op, id and args')
        if ast['op']=='input':
            if ast['id']not in self.inputs or ast['args']:raise ValueError('Unknown or malformed typed input')
            return self.inputs[ast['id']].type_name
        if ast['op']!='apply'or ast['id']not in self.ops:raise ValueError('Operation not in trusted registry')
        operation=self.ops[ast['id']]
        if tuple(self.type_of(a)for a in ast['args'])!=operation.input_types:raise ValueError('Operation argument type mismatch')
        return operation.output_type
    def enqueue(self,ast):
        self.type_of(ast)
        with (self.out/'queue.jsonl').open('a')as f:f.write(canonical({'ast':ast})+'\n')
    def _queued(self,max_cost):
        path=self.out/'queue.jsonl'
        if not path.exists():return []
        found=[]
        with path.open()as f:
            f.seek(self.queue_offset)
            while True:
                start=f.tell();line=f.readline()
                if not line:break
                if not line.endswith('\n'):f.seek(start);break
                try:
                    value=json.loads(line);ast=value['ast']if isinstance(value,dict)and set(value)=={'ast'}else value
                    self.type_of(ast)
                    if cost(ast)>max_cost:raise ValueError('Queued AST outside declared finite cost bound')
                    found.append(ast)
                except (ValueError,TypeError,KeyError,RecursionError)as e:
                    with (self.out/'queue-rejections.jsonl').open('a')as rejected:rejected.write(canonical({'error':str(e),'line':line.rstrip(),'search_fingerprint':self.fingerprint})+'\n')
            self.queue_offset=f.tell()
        return found
    def _validation(self,value):
        try:answer=self.validator(value,self.context)
        except self.missing_exceptions as e:answer={'missing_input':str(e)}
        if not isinstance(answer,dict):answer={'missing_input':'Verifier did not return an obligation dictionary'}
        if answer.get('scope')not in (None,self.scope):raise ValueError('Verifier scope differs from the trusted goal scope')
        supplied=answer.get('obligations',[])
        if not isinstance(supplied,list):supplied=[]
        rows={}
        for r in supplied:
            if not isinstance(r,dict)or r.get('id')not in self.required:continue
            if r['id']in rows:raise ValueError('Duplicate verifier obligation')
            if r.get('status')not in ('PASS','FAIL','MISSING_INPUT','NOT_PROVED','ERROR'):raise ValueError('Invalid verifier obligation status')
            rows[r['id']]=r
        actual=[rows.get(k,{'id':k,'status':'MISSING_INPUT','reason':answer.get('missing_input','Verifier omitted the actual required obligation')})for k in self.required]
        complete=all(r['status']=='PASS'for r in actual)
        status='SCOPED_HIT'if complete else ('CANDIDATE_FAIL'if any(r['status']=='FAIL'for r in actual)else'PARTIAL')
        return {'status':status,'complete_for_target_scope':complete,'target_scope':self.scope,'obligations':actual,
                'evidence':answer.get('evidence'),'candidate_self_attestation_used':False}
    def _value(self,ast):
        key=ast_hash(ast)
        if key in self.values:return self.values[key]
        previous=self.records.get(key)
        if previous and previous['evaluation_status']=='MISSING_INPUT':
            value=Unavailable(previous['missing_input']);self.values[key]=value;return value
        if previous and previous['evaluation_status']in ('EXECUTION_ERROR','TYPE_VALUE_ERROR'):
            value=Unavailable('Previously rejected dependency: '+previous.get('error',''));self.values[key]=value;return value
        typ=self.type_of(ast);record={'search_fingerprint':self.fingerprint,'ast_hash':key,'ast':ast,'cost':cost(ast),'output_type':typ,'target_scope':self.scope}
        try:
            if previous is None and self._evaluation_limit is not None and self.new_evaluations>=self._evaluation_limit:raise SearchBudgetStop('new_evaluation_budget')
            if ast['op']=='input':value=self.inputs[ast['id']].value
            else:
                args=[self._value(a)for a in ast['args']]
                missing=[a.reason for a in args if isinstance(a,Unavailable)]
                if missing:raise MissingInput('; '.join(missing))
                if previous is None and self._evaluation_limit is not None and self.new_evaluations>=self._evaluation_limit:raise SearchBudgetStop('new_evaluation_budget')
                value=self.ops[ast['id']].execute(*args)
            if isinstance(value,Unavailable):raise MissingInput(value.reason)
            check=self.type_validators.get(typ)
            if check is not None and check(value)is not True:
                record.update(evaluation_status='TYPE_VALUE_ERROR',error='Trusted runtime type validator rejected produced value');value=Unavailable(record['error'])
            else:
                record.update(evaluation_status='VALUE_BUILT',value_python_type=type(value).__module__+'.'+type(value).__name__)
                if typ==self.target and previous is None:record['goal']=self._validation(value)
        except SearchBudgetStop:raise
        except self.missing_exceptions as e:
            value=Unavailable(str(e));record.update(evaluation_status='MISSING_INPUT',missing_input=str(e))
            if typ==self.target:record['goal']={'status':'PARTIAL','complete_for_target_scope':False,'target_scope':self.scope,
                'obligations':[{'id':k,'status':'MISSING_INPUT','reason':str(e)}for k in self.required]}
        except Exception as e:
            value=Unavailable(type(e).__name__+': '+str(e));record.update(evaluation_status='EXECUTION_ERROR',error=value.reason)
            if typ==self.target:record['goal']={'status':'PARTIAL','complete_for_target_scope':False,'target_scope':self.scope,
                'obligations':[{'id':k,'status':'ERROR','reason':value.reason}for k in self.required]}
        self.values[key]=value
        if previous is None:
            record=evidence_json(record)
            with self.archive.open('a')as f:f.write(canonical(record)+'\n')
            self.records[key]=record;self.new_evaluations+=1
        else:self.dependency_replays+=1
        return value
    def publish_status(self,updates=None):
        """Publish current results after either enumeration OR an AI proposal."""
        path=self.out/'status.json';status={}
        if path.exists():
            previous=json.loads(path.read_text())
            if previous.get('search_fingerprint')==self.fingerprint:status=previous
        goals=[r for r in self.records.values()if r['output_type']==self.target]
        hits=[r for r in goals if r.get('goal',{}).get('complete_for_target_scope')is True]
        unresolved=[r for r in goals if r.get('goal',{}).get('status')=='PARTIAL']
        def order(r):return (sum(o['status']!='PASS'for o in r['goal']['obligations']),r['cost'],r['ast_hash'])
        status.setdefault('coverage','no_enumeration_claim_from_direct_proposals')
        status.update({'search_fingerprint':self.fingerprint,'target_type':self.target,'target_scope':self.scope,
            'verifier_code_sha256':self.fingerprint_components['verifier_code_sha256'],
            'required_obligation_ids':self.required,'coverage_is_all_physical_sources':False,
            'goal_candidates_in_current_fingerprint_archive':len(goals),'scoped_hits':len(hits),
            'partial_candidates':len(unresolved),'failed_candidates':sum(r.get('goal',{}).get('status')=='CANDIDATE_FAIL'for r in goals),
            'missing_input_candidates':sum(r['evaluation_status']=='MISSING_INPUT'for r in goals),
            'AST_execution_cost_is_physical_weight':False,'updated_unix_time':time.time()})
        if updates:status.update(updates)
        self._write(path,status);self._write(self.out/'best-unresolved.json',min(unresolved,key=order)if unresolved else None)
        self._write(self.out/'scoped-hits.json',hits)
        return status
    def evaluate_program(self,ast,max_cost=64):
        """Execute one typed proposal directly on the SAME constructor/verifier.

        A proposal can be above the current enumeration bound. That does not
        assert enumeration coverage of its cost, or change physical weights.
        """
        try:
            self.type_of(ast)
            if type(max_cost)is not int or max_cost<1 or cost(ast)>max_cost:raise ValueError('Proposal outside explicit AST cost bound')
        except (ValueError,TypeError,KeyError,RecursionError)as e:
            rejected={'search_fingerprint':self.fingerprint,'evaluation_status':'INVALID_AST',
                      'target_scope':self.scope,'error':str(e),'ast':evidence_json(ast)}
            with (self.out/'proposal-rejections.jsonl').open('a')as f:f.write(canonical(rejected)+'\n')
            self.publish_status({'last_activity':'direct_program_rejected','last_proposal_error':str(e)})
            return rejected
        started=time.monotonic();before=self.new_evaluations
        self._evaluation_limit=None
        key=ast_hash(ast)
        if key not in self.records:self._value(ast)
        record=self.records[key]
        self.publish_status({'last_activity':'direct_program_evaluation','last_proposal_ast_hash':key,
            'last_proposal_max_cost':max_cost,'last_proposal_cost':cost(ast),
            'new_evaluations':self.new_evaluations-before,'elapsed_seconds':time.monotonic()-started})
        return record
    def _partitions(self,total,count):
        if count==0:
            if total==0:yield ()
        elif count==1:
            if total>=1:yield (total,)
        else:
            for first in range(1,total-count+2):
                for rest in self._partitions(total-first,count-1):yield (first,*rest)
    def run(self,max_cost,max_asts=10000,max_new_evaluations=None):
        if type(max_cost)is not int or max_cost<1 or max_asts<1:raise ValueError('Positive finite search bounds required')
        started=time.monotonic();start_evaluations=self.new_evaluations;start_replays=self.dependency_replays
        self._evaluation_limit=None if max_new_evaluations is None else start_evaluations+max_new_evaluations
        bins={};seen=set();pending={};grammar_asts=0;processed=0;last_complete_cost=0;stop=None
        def stage(ast,grammar=False):
            nonlocal grammar_asts
            typ=self.type_of(ast);key=ast_hash(ast)
            if key in seen:return
            if len(seen)>=max_asts:raise OverflowError('finite_ast_budget')
            seen.add(key);c=cost(ast);bins.setdefault((typ,c),[]).append(ast);pending.setdefault(c,{})[key]=ast
            if grammar:grammar_asts+=1
        def drain(through):
            nonlocal processed
            for ast in self._queued(max_cost):stage(ast)
            for c in sorted(k for k in pending if k<=through):
                for key,ast in sorted(pending.pop(c).items()):
                    if max_new_evaluations is not None and self.new_evaluations-start_evaluations>=max_new_evaluations:raise OverflowError('new_evaluation_budget')
                    if key not in self.records:self._value(ast)
                    processed+=1
        try:
            for key in sorted(self.inputs):stage({'op':'input','id':key,'args':[]},True)
            for n in range(1,max_cost+1):
                for op in self.ops.values():
                    for partition in self._partitions(n-1,len(op.input_types)):
                        lists=[bins.get((typ,c),[])for typ,c in zip(op.input_types,partition)]
                        if any(not pool for pool in lists):continue
                        for args in itertools.product(*lists):stage({'op':'apply','id':op.id,'args':list(args)},True)
                drain(n);last_complete_cost=n
            drain(max_cost)
        except OverflowError as e:stop=str(e)
        self._evaluation_limit=None
        status={'max_cost':max_cost,'max_asts':max_asts,'last_activity':'enumeration',
                'coverage':'complete_for_declared_registry_inputs_and_cost_bound'if stop is None else'partial_finite_search_budget',
                'coverage_is_all_physical_sources':False,'last_complete_cost':last_complete_cost,'stop_reason':stop,
                'grammar_asts_discovered':grammar_asts,'unique_asts_discovered':len(seen),'processed_this_invocation':processed,
                'new_evaluations':self.new_evaluations-start_evaluations,'dependency_replays':self.dependency_replays-start_replays,
                'elapsed_seconds':time.monotonic()-started,'AST_execution_cost_is_physical_weight':False}
        return self.publish_status(status)
