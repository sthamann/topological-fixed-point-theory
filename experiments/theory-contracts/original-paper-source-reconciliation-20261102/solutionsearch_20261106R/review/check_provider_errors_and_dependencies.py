"""Targeted provider-stage and existing dependency binding checks, no models."""
from pathlib import Path
import hashlib
import json
import sys
import tempfile
import types
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from source_identity import SourceProvider, SourceTimeFamily, SourceFamilyMissing, ORIGINAL_SCOPE, IDS
from completion import OriginalCompletion
from catalogue import inputs, operations
from source_recovery import native_constructors as native
from search.synthesis import semantic_fingerprint

REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
bindings=inputs()
operator=native.nil_action(native.quartic_source_program(bindings['chart_2'].value,
    native.scalar_c3(),native.scalar_c6_ratio(),bindings['source_symbols'].value))
candidate=native.source_jets(operator,bindings['dictionary_2'].value)


def call(provider, value=candidate):
    verifier=OriginalCompletion(REPO,original_providers={'original.common_source':provider})
    return {r['id']:r for r in verifier(value)['obligations'] if r['id'] in IDS}


checks={}
after={}
for exception in (ValueError,KeyError,TypeError):
    def broken(context,exception=exception):
        raise exception('independent upstream producer regression')
    result=call(SourceProvider(broken,ORIGINAL_SCOPE,'host-bound producer regression only'))
    checks[exception.__name__+'_provider_load_is_ERROR']=all(r['status']=='ERROR' and r['error_stage']=='original_provider_load' for r in result.values())
    after[exception.__name__]=result

def missing(context):
    raise SourceFamilyMissing('explicit missing whole-source input')
checks['explicit_missing_input_stays_MISSING_INPUT']=all(r['status']=='MISSING_INPUT' for r in call(SourceProvider(missing,ORIGINAL_SCOPE,'explicit missing input control')).values())
ju,jy=candidate['source_symbols'];t=s.Symbol('t',positive=True)
d=candidate['dictionary']
family=SourceTimeFamily(candidate['matrix'],{'U':ju,'Y':jy},{'interval':t},d['P'],d['cap'],d['reference_scalar'],d['outer_scalar'])
bad={**candidate,'original_packet':{'provider_id':'original.common_source','source_symbols':{'U':ju,'Y':jy},'time_symbols':{'interval':t},'intertwiner':s.zeros(1)}}
checks['candidate_comparison_defect_stays_FAIL']=all(r['status']=='FAIL' for r in call(SourceProvider(lambda ctx:family,ORIGINAL_SCOPE,'host-bound comparison control'),bad).values())

# Exercise the actual semantic_fingerprint on the actual verifier/provider
# closure, replacing only its domain_path with a temporary copy. A byte
# change at the same captured path must alter verifier_bound_state_sha256.
verifier=OriginalCompletion(REPO)
original=verifier.original_providers['original.common_source']
cells=dict(zip(original.load.__code__.co_freevars,original.load.__closure__))
actual_domain=cells['domain_path'].cell_contents
def cell(value):
    return (lambda:value).__closure__[0]
with tempfile.TemporaryDirectory(prefix='tfpt-domain-dependency-') as td:
    path=Path(td)/'native-diagnostic-v1.json'
    path.write_bytes(actual_domain.read_bytes())
    closure=tuple(cell(path) if name=='domain_path' else cells[name] for name in original.load.__code__.co_freevars)
    load=types.FunctionType(original.load.__code__,original.load.__globals__,original.load.__name__,original.load.__defaults__,closure)
    verifier.original_providers={'original.common_source':SourceProvider(load,original.acceptance_scope,original.location)}
    before=semantic_fingerprint(verifier)
    path.write_bytes(path.read_bytes()+b'\n')
    after_fingerprint=semantic_fingerprint(verifier)
checks['actual_domain_is_captured_as_Path']=isinstance(actual_domain,Path)
checks['domain_byte_change_changes_verifier_bound_fingerprint']=before!=after_fingerprint
files=[native.NATIVE/'kernel_functional.py',native.NATIVE/'joint_source.py',native.NATIVE/'source_sewing.py',
       native.ROOT/'verification/v330_qgamb_admissible_measure.py',
       native.ROOT/'documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/002-section-intrinsic-closed-branch-and-canonical-rigid-obje.tex']
expected=hashlib.sha256(b''.join(p.read_bytes() for p in files)).hexdigest()
checks['native_version_binds_actual_kernel_bytes']=native.version()==expected
registry=operations()
checks['registered_native_operations_carry_that_version']=all(o['version']==expected for o in registry if o['id'].startswith('native.'))
manifests=sorted((ROOT/'runs/grounded-context-preview').glob('manifest-*.json'))
recorded=[json.loads(p.read_text()) for p in manifests]
checks['actual_recorded_fingerprint_contains_kernel_bound_versions']=any(any(o['id']=='native.even_kernel' and o['version']==expected for o in m['operations']) for m in recorded)
result={'all_pass':all(checks.values()),'checks':checks,'actual_completion_provider_errors_after_fix':after,
    'observed_before_fix':{name:{key:'FAIL' for key in IDS} for name in ('ValueError','KeyError','TypeError')},
    'domain_binding':'actual provider.load closure Path -> semantic_fingerprint -> verifier_bound_state_sha256',
    'kernel_binding':'native.version concatenated original bytes -> registry Operation.version -> fingerprint_components.operations',
    'native_version':expected,'source_identity_sha256':hashlib.sha256((ROOT/'source_identity.py').read_bytes()).hexdigest(),
    'completion_unchanged_sha256':hashlib.sha256((ROOT/'completion.py').read_bytes()).hexdigest(),
    'canonical_mutations':0,'model_calls':0}
(ROOT/'review/provider-errors-dependency-controls.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
raise SystemExit(0 if result['all_pass'] else 1)
