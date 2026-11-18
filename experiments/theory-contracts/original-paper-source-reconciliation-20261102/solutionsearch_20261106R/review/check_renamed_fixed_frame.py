"""Targeted actual OriginalCompletion regression for renamed moving frames.

The host-injected two-dimensional family is arithmetic control data, never an
original physical provider. This tests the fixed-frame comparator's premise.
"""
from pathlib import Path
import copy
import hashlib
import json
import sys
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from source_identity import SourceTimeFamily, SourceProvider, ORIGINAL_SCOPE, IDS
from completion import OriginalCompletion

u,v,a,b=s.symbols('u v a b',real=True)
t,tau=s.symbols('t tau',positive=True)
H=s.diag(1,2)
O=s.diag(1+u+t,2+v+2*t)
family=SourceTimeFamily(O,{'U':u,'Y':v},{'interval':t},s.eye(2),s.eye(2),1,1,
                        carried_operators={'H':H})
R=lambda z:s.Matrix([[s.cos(z),-s.sin(z)],[s.sin(z),s.cos(z)]])
verifier=OriginalCompletion('/Users/stefanhamann/Projekte/tfpt-theoryv4',
    original_providers={'review':SourceProvider(lambda ctx:family,ORIGINAL_SCOPE,
        'HOST-INJECTED ARITHMETIC REGRESSION ONLY; not the default physical provider')})


def candidate(original_z=None,candidate_z=None):
    C=s.eye(2) if candidate_z is None else R(candidate_z)
    W=s.eye(2) if original_z is None else R(original_z)
    K=C*O.subs({u:a,v:b,t:tau})*C.T
    forms={'Z0':K.subs({a:0,b:0}),'ZU':K.diff(a).subs({a:0,b:0}),
           'ZY':K.diff(b).subs({a:0,b:0}),'ZUY':K.diff(a,b).subs({a:0,b:0})}
    values={k:s.simplify(s.trace(m)) for k,m in forms.items()}
    z=values['Z0']
    return {'matrix':K,'source_symbols':(a,b),
        'dictionary':{'P':s.eye(2),'cap':s.eye(2),'reference_scalar':1,'outer_scalar':1},
        'ordered_source_forms':forms,**values,
        'normalized':{'connected_UY':s.simplify(values['ZUY']/z-values['ZU']*values['ZY']/z**2)},
        'original_packet':{'provider_id':'review','source_symbols':{'U':a,'Y':b},
            'time_symbols':{'interval':tau},'intertwiner':W,
            'carried_operators':{'H':C*H*C.T}}}


def rows(c):
    return {x['id']:x for x in verifier(c)['obligations']}


good=candidate()
good_rows=rows(good)
time_rows=rows(candidate(t,tau))
source_rows=rows(candidate(u,a))
high=copy.deepcopy(good)
high['matrix']+=a**2*b**3*s.eye(2)
high_rows=rows(high)
controls={
    'fixed_frame_still_passes_entire_family':all(good_rows[k]['status']=='PASS' for k in IDS),
    'actual_completion_computes_four_values':good_rows['constructed_common_values']['status']=='PASS',
    'renamed_original_time_dependent_frame_requires_contacts':all(time_rows[k]['status']=='MISSING_INPUT' for k in IDS),
    'renamed_original_source_dependent_frame_requires_contacts':all(source_rows[k]['status']=='MISSING_INPUT' for k in IDS),
    'higher_order_mutation_preserves_four_jets':all(high['matrix'].diff(a,i).diff(b,j).subs({a:0,b:0})==good['matrix'].diff(a,i).diff(b,j).subs({a:0,b:0}) for i,j in [(0,0),(1,0),(0,1),(1,1)]),
    'same_four_jets_pass_only_common_values':high_rows['constructed_common_values']['status']=='PASS',
    'whole_family_rejects_higher_source_orders':all(high_rows[k]['status']=='FAIL' for k in IDS),
}
result={'scope':'Independent fixed-frame software regression, not physical source identity',
    'checks':controls,'all_pass':all(controls.values()),
    'observed_before_fix':{'original_time_t_after_tau_to_t_renaming':{k:'PASS' for k in IDS},
                         'original_source_u_after_a_to_u_renaming':{k:'PASS' for k in IDS}},
    'after_fix_time':{k:time_rows[k] for k in IDS},'after_fix_source':{k:source_rows[k] for k in IDS},
    'source_identity_sha256':hashlib.sha256((ROOT/'source_identity.py').read_bytes()).hexdigest(),
    'completion_sha256':hashlib.sha256((ROOT/'completion.py').read_bytes()).hexdigest()}
(ROOT/'review/renamed-fixed-frame-controls.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
raise SystemExit(0 if result['all_pass'] else 1)
