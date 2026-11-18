"""Whole original v330 label-word controls; never primitive physical credit."""
from pathlib import Path
import copy,json,sys,hashlib
import sympy as s
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from source_identity import *
from completion import OriginalCompletion
from catalogue import inputs
from source_recovery import native_constructors as native
n,m=s.symbols('n m',integer=True,nonnegative=True);j0,j1=s.symbols('j0 j1',real=True)
u,y=s.symbols('jU jY',real=True);a,b=s.symbols('left right',integer=True,nonnegative=True)
T=s.Matrix([[1651,1267,1456],[1267,1651,1456],[1456,1456,1462]])/4374
Pi=s.ones(3)/3;v=s.Matrix([1,-1,0]);w=s.Matrix([1,1,-2]);P1=v*v.T/2;P2=w*w.T/6
power=lambda k:Pi+s.Rational(64,729)**k*P1+s.Rational(1,729)**k*P2
assert power(1)==T
E=[s.diag(*(int(i==k)for i in range(3)))for k in range(3)]
# Actual complete finite label word, with both original block counts. Its
# fixed middle insertion means no K(n+m)=K(n)K(m) assumption is imposed.
K=s.diag(1,s.exp(j1),1)*power(n)*E[2]*power(m)*s.diag(s.exp(j0),1,1)
family=SourceTimeFamily(K,{'initial_label0':j0,'final_label1':j1},
 {'left_original_block_count':n,'right_original_block_count':m},s.eye(3),Pi,1,1,
 carried_operators={'original_transfer':T,'label0':E[0],'label1':E[1],'middle_label2':E[2]},
 antilinear_operators={'original_real_conjugation':s.eye(3)})
# A complex unitary frame also checks the distinct antilinear rule.
W=s.diag(s.I,1,1)*s.Matrix([[0,1,0],[0,0,1],[1,0,0]])
moved=transport_family(family,W);rename={j0:u,j1:y,n:a,m:b}
change=lambda A:s.Matrix(A).subs(rename,simultaneous=True)
candidate={'matrix':change(moved.kernel),'dictionary':{'P':change(moved.P),'cap':change(moved.cap),
 'reference_scalar':moved.reference_scalar,'outer_scalar':moved.outer_scalar}}
packet={'provider_id':'finite.original_v330_label_words','intertwiner':W,
 'source_symbols':{'initial_label0':u,'final_label1':y},
 'time_symbols':{'left_original_block_count':a,'right_original_block_count':b},
 'carried_operators':{k:change(v)for k,v in moved.carried_operators.items()},
 'antilinear_operators':{k:change(v)for k,v in moved.antilinear_operators.items()}}
candidate['original_packet']=packet
scope='original_v330_finite_label_word_NOT_primitive_TFPT'
providers={packet['provider_id']:SourceProvider(lambda ctx:family,scope,'verification/v330_qgamb_admissible_measure.py; exact original finite label-word construction')}
check=lambda c:original_identity_results(c,providers,required_scope=scope)
controls={};good=check(candidate)
controls['entire_original_finite_label_source_and_time_match']=all(r['status']=='PASS'for r in good.values())
controls['no_semigroup_identity_required_for_marked_interval_family']=controls['entire_original_finite_label_source_and_time_match']
high=copy.deepcopy(candidate);high['matrix']+=u**2*y**3*s.eye(3)
controls['mutation_preserves_all_four_zero_source_jets']=all(high['matrix'].diff(u,i).diff(y,j).subs({u:0,y:0})==candidate['matrix'].diff(u,i).diff(y,j).subs({u:0,y:0})for i,j in [(0,0),(1,0),(0,1),(1,1)])
controls['same_four_jets_do_not_pass_whole_source_identity']=all(r['status']=='FAIL'for r in check(high).values())
bad=copy.deepcopy(candidate);bad['original_packet']['time_symbols']['left_original_block_count']=2*a
controls['new_time_scale_rejected']=all(r['status']=='FAIL'for r in check(bad).values())
bad=copy.deepcopy(candidate);bad['original_packet']['carried_operators']['label0']=change(moved.carried_operators['label1'])
controls['changed_original_mark_rejected']=all(r['status']=='FAIL'for r in check(bad).values())
bad=copy.deepcopy(candidate);bad['dictionary']['reference_scalar']=2
controls['changed_reference_rejected']=all(r['status']=='FAIL'for r in check(bad).values())
bad=copy.deepcopy(candidate);bad['dictionary']['cap']+=s.eye(3)
controls['changed_original_cap_rejected']=all(r['status']=='FAIL'for r in check(bad).values())
bad=copy.deepcopy(candidate);bad['original_packet']['antilinear_operators']['original_real_conjugation']=W*s.eye(3)*W.H
controls['antilinear_not_misread_as_linear_transport']=all(r['status']=='FAIL'for r in check(bad).values())
controls['finite_original_provider_never_earns_primitive_scope']=all(r['status']=='MISSING_INPUT'for r in original_identity_results(candidate,providers).values())
self_attested=copy.deepcopy(candidate);self_attested.update(physical_source_complete=True,source_id='primitive_original',proof_status='PASS')
controls['source_id_and_self_attestation_not_authority']=all(r['status']=='MISSING_INPUT'for r in original_identity_results(self_attested,{}).values())
empty=replace(family,source_symbols={},time_symbols={});ep=copy.deepcopy(packet);ep['source_symbols']={};ep['time_symbols']={}
controls['empty_family_does_not_vacuously_pass']=all(r['status']=='MISSING_INPUT'for r in original_identity_results({**candidate,'original_packet':ep},{ep['provider_id']:SourceProvider(lambda ctx:empty,scope,'empty control')},required_scope=scope).values())
# Pure typed-form transport arithmetic, not an extra v330 physical claim.
form=s.Matrix([[0,1,0],[-1,0,2],[0,-2,0]])
typed=replace(family,bilinear_forms={'arithmetic_skew_form':form},domain_relations={'arithmetic_domain_row':s.Matrix([[1,2,0]])})
typed_moved=transport_family(typed,W)
controls['bilinear_congruence_and_domain_pullback_rules']=W.T*typed_moved.bilinear_forms['arithmetic_skew_form']*W==form and typed_moved.domain_relations['arithmetic_domain_row']*W==typed.domain_relations['arithmetic_domain_row']
# The actual R candidate and actual installed default original contribution
# operation are exercised; no model/backend or private network is called.
ins=inputs();chart=ins['chart_2'].value;symbols=ins['source_symbols'].value
op=native.nil_action(native.quartic_source_program(chart,native.scalar_c3(),native.scalar_c6_ratio(),symbols))
real_candidate=native.source_jets(op,ins['dictionary_2'].value)
completion=OriginalCompletion('/Users/stefanhamann/Projekte/tfpt-theoryv4','tfpt_complete')(real_candidate)
rows={r['id']:r for r in completion['obligations']}
controls['default_calls_actual_installed_original_operation']=all('actual_native_call'in rows[k]for k in IDS)
controls['default_missing_original_body_remains_precise']=all(rows[k]['status']=='MISSING_INPUT'and'upstream' in rows[k]['missing_inputs'][0]for k in IDS)
controls['actual_candidate_values_still_computed']=rows['constructed_common_values']['status']=='PASS'
controls['T1_through_T8_remain_missing_not_unsupported']=all(rows['TFPT.TOE.COMPLETE.01.T'+str(i)]['status']=='MISSING_INPUT'for i in range(1,9))
result={'scope':'finite original v330 word identity plus software mutations; no primitive or T1-T8 closure',
 'controls':controls,'count':len(controls),'passed':all(controls.values()),'whole_family_comparison':good,
 'default_original_attempt':{k:rows[k]for k in IDS},
 'source_sha256':hashlib.sha256((ROOT/'source_identity.py').read_bytes()).hexdigest(),
 'completion_sha256':hashlib.sha256((ROOT/'completion.py').read_bytes()).hexdigest()}
(ROOT/'review/source-identity-controls.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2));assert result['passed']
