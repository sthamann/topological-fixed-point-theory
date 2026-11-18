"""Versioned spatial field plugin for the EXISTING candidate search engine."""
from pathlib import Path
import hashlib,importlib.util,json,os,sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from assemble import Assembler,V252,HIGGS
DOMAIN_PATH=Path(os.environ['TFPT_FIELD_DOMAIN']).resolve()
A=Assembler(DOMAIN_PATH);DOMAIN=A.domain
SCHEMA='tfpt.field-configuration.v1'
RUN_LABEL='full96_spatial_field_candidates_NOT_complete_physical_source'
SOURCE_PATH=HERE.parent/'source/evaluator.py'
SOURCE_FILES=[str(DOMAIN_PATH),str(HERE/'assemble.py'),str(HERE/'launch.py'),str(HERE/'repo_paths.py'),*DOMAIN['_resolved_source_originals']]
if SOURCE_PATH.exists():SOURCE_FILES.append(str(SOURCE_PATH))
if isinstance(DOMAIN.get('whole_source_bundle'),str):
    bundle=Path(DOMAIN['whole_source_bundle'])
    if not bundle.is_absolute():bundle=DOMAIN_PATH.parent/bundle
    SOURCE_FILES.append(str(bundle.resolve()))
    # Actual referenced inputs participate in the engine's existing resume key.
    decoded=json.loads(bundle.read_text())
    for ref in decoded.get('source_refs',[]):
        if isinstance(ref,dict) and ref.get('path'):
            SOURCE_FILES.append(str(Path(ref['path']).expanduser().resolve()))
elif isinstance(DOMAIN.get('whole_source_bundle'),dict):
    for ref in DOMAIN['whole_source_bundle'].get('source_refs',[]):
        if isinstance(ref,dict) and ref.get('path'):
            SOURCE_FILES.append(str(Path(ref['path']).expanduser().resolve()))
_SPEC=importlib.util.spec_from_file_location('tfpt_field_source_adapter',SOURCE_PATH) if SOURCE_PATH.exists() else None
ADAPTER=importlib.util.module_from_spec(_SPEC)if _SPEC else None
if ADAPTER is not None:_SPEC.loader.exec_module(ADAPTER)
COUNT=len(DOMAIN['field_proposal']['slots']);BASE=len(DOMAIN['field_proposal']['allowed_multipliers'])
ENUMERATION_SIZE=BASE**COUNT

def sample(rng,n):return [{'schema':SCHEMA,'domain_sha256':DOMAIN['_sha256'],'genes':rng.integers(0,BASE,COUNT).tolist()}for _ in range(n)]
def canonicalize(g):
    if not valid(g):raise ValueError('Malformed field genotype or changed domain')
    return {'schema':SCHEMA,'domain_sha256':DOMAIN['_sha256'],'genes':list(g['genes'])}
def valid(g):return isinstance(g,dict)and set(g)=={'schema','domain_sha256','genes'}and g['schema']==SCHEMA and g['domain_sha256']==DOMAIN['_sha256']and isinstance(g['genes'],list)and len(g['genes'])==COUNT and all(type(x)is int and 0<=x<BASE for x in g['genes'])
def mutate(rng,parent):
    g=canonicalize(parent);j=int(rng.integers(COUNT));g['genes'][j]=int(rng.integers(BASE));return g
def crossover(rng,a,b):
    a=canonicalize(a);b=canonicalize(b);a['genes']=[x if rng.random()<.5 else y for x,y in zip(a['genes'],b['genes'])];return a
def enumerate_candidates(shard=0,shards=1):
    if shards<1 or not 0<=shard<shards:raise ValueError('Invalid shard')
    for index in range(shard,ENUMERATION_SIZE,shards):
        value=index;genes=[]
        for _ in range(COUNT):genes.append(value%BASE);value//=BASE
        yield {'schema':SCHEMA,'domain_sha256':DOMAIN['_sha256'],'genes':genes}
def capability_report():
    if ADAPTER is None:return {'ready':False,'missing_operations':['original whole source adapter source/evaluator.py'],'physical_source_complete':False}
    return ADAPTER.contribution_capability_report(DOMAIN)
def ensemble_registration():
    required=DOMAIN.get('ensemble_requirements')
    if not isinstance(required,list) or not required or any(not isinstance(key,str) or not key.strip() for key in required) or len(set(required))!=len(required):
        return {'required_ids':[], 'registered_ids':[],
                'missing_implementation_ids':['explicit_nonempty_unique_ensemble_requirements_declaration'],
                'status':'MISSING_REQUIREMENT_DECLARATION'}
    registry=ADAPTER.ensemble_test_registry(DOMAIN) if ADAPTER is not None and hasattr(ADAPTER,'ensemble_test_registry') else {}
    registered=[key for key,value in registry.items() if callable(value)]
    return {'required_ids':required,'registered_ids':registered,
            'missing_implementation_ids':[key for key in required if key not in registered],
            'status':'IMPLEMENTED_NOT_YET_MEASURED' if all(key in registered for key in required) else 'MISSING_IMPLEMENTATION'}
def evaluate(candidates,device):
    if device!='cpu':raise ValueError('This first full96 numerical backend is CPU; no unimplemented CUDA path is claimed')
    results=[]
    for c in candidates:
        c=canonicalize(c);assembled=A.assemble(c)
        tol=float(DOMAIN['matrix_residual_tolerance']);rr=assembled['kinematic_residuals']
        hard=[{'id':key,'kind':'per_configuration_operator_numeric','passed':value<=tol,'residual':value,'tolerance':tol}for key,value in rr.items()]
        hard.append({'id':'integer_genotype_original_domain','kind':'exact_schema','passed':True})
        if ADAPTER is None:
            source={'status':'NOT_EVALUATED_MISSING_PRODUCER','contributions':{key:None for key in ['W0','WU','WY','WUY']},
             'physical_source_complete':False,'ensemble_checks':[],'missing_operations':['original whole source adapter']}
        else:source=ADAPTER.evaluate_contribution(c,assembled,DOMAIN)
        # Genotype fitness ranks only implemented hard violations and coverage.
        # It is NEVER an action, proposal probability, or physical source weight.
        failed=sum(not x['passed']for x in hard)
        values_complete=all(source.get('contributions',{}).get(k)is not None for k in ['W0','WU','WY','WUY'])
        results.append({'objectives':[failed,0 if values_complete else 1],
          'feasible':bool(failed==0),'descriptor':[c['genes']],
          'details':{'candidate_id':assembled['candidate_id'],'configuration_hard_pass':failed==0,
          'hard_checks':hard,'matrix':assembled['diagnostics'],'source':source,
          'feasible_scope':'per_configuration_kinematic_conditions_only',
          'ensemble_registration':ensemble_registration(),
          'physical_solution':False,'full_search_completion_claimed':False,
          'fitness_is_physical_weight':False,'numerical_matrix_checks_are_exact_proof':False,
          'native_inputs_diagnostic':DOMAIN.get('input_role')!='physical_original_input',
          'h_N_role':assembled['h_N_role']}})
    return results
