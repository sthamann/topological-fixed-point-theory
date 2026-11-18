"""Small arithmetic/mutation controls on existing source branches, not physics."""
from pathlib import Path
import json,copy,hashlib
import sympy as s
from original_checks import *
HERE=Path(__file__).resolve().parent
OLD=HERE/'arithmetic-control-bundle.json'
b=json.loads(OLD.read_text())
# Recover the intended exact rational arithmetic-control matrices from their
# previously generated float serialization. This is NOT physical native input.
def exact(x):
 if isinstance(x,list):return [exact(y)for y in x]
 if isinstance(x,dict):return {k:exact(v)for k,v in x.items()}
 if isinstance(x,float):return s.nsimplify(x,tolerance=1e-12,rational=True)
 return x
b=exact(b);d={key:b[key]for key in ['gram_form','physical_selector_form','boundary_pairing_form','reference_scalar','scalar_factor']}
d.update(quotient_lift=s.Matrix([[s.Rational(1,2),0],[0,1],[0,0]]),null_basis=s.Matrix([[0],[0],[1]]),atoms=b['fiber_records'])
controls={};out=common_four(d);controls['common_four_exact_control']=out['status']=='PASS'
wrong=copy.deepcopy(d);wrong['observed_values']=dict(out['values']);wrong['observed_values']['ZU']+=1
controls['changed_integrated_Z_detected']=common_four(wrong)['status']=='FAIL'
wrong=copy.deepcopy(d);wrong['atoms'][0]['ordered_source_forms']['ZUY'][2][2]=[1,0]
controls['marked_null_class_detected']=common_four(wrong)['status']=='FAIL'
controls['missing_actual_mass_rejected']=common_four({**d,'atoms':[{'original_configuration_id':'missing'}]})['status']=='MISSING_INPUT'
# Actual old v330 transfer, not a map of its labels to dynamical Higgs fields.
T=s.Matrix([[1651,1267,1456],[1267,1651,1456],[1456,1456,1462]])/4374;Pi=s.ones(3)/3
amps=[T**3];Q=s.Matrix([[s.trace((T**3).H*T**3)]])
rp={'basis_words':[[]],'actual_Q0':Q,'ordered_word_amplitudes':amps,'reflection_normalization':1,
 'gram_form':Q,'quotient_lift':s.Matrix([[1/s.sqrt(Q[0,0])]]),'null_basis':s.zeros(1,0)}
controls['native_transfer_reflected_amplitude_positive']=reflected_rp(rp)['status']=='PASS'
controls['unrelated_input_Gram_not_enough']=reflected_rp({k:v for k,v in rp.items()if k!='ordered_word_amplitudes'})['status']=='MISSING_INPUT'
controls['changed_actual_reflected_Q0_detected']=reflected_rp({**rp,'actual_Q0':-Q})['status']=='FAIL'
controls['empty_reflected_family_rejected']=reflected_rp({**rp,'basis_words':[]})['status']=='MISSING_INPUT'
# R048 literal fixed-R quadratic branch: symbolic positive precision, no new
# numerical physical parameter, and no full-source/whole-Higgs interpretation.
x=s.Symbol('eta',real=True);k=s.Symbol('K',positive=True);S=k*x*x/2
sd={'coordinates':[x],'action':S,'coordinate_density':1,'variation_vector':[1],
 'whole_marked_insertion':x,'original_domain_bounds':[(-s.oo,s.oo)],
 'integrate_against_exp_minus_action':lambda f:s.integrate(f*s.exp(-S),(x,-s.oo,s.oo)),
 'integrate_boundary_face':lambda i,e,f:s.limit(f,x,e)}
controls['original_quadratic_SD_integrated_contact']=ward_sd(sd)['status']=='PASS'
controls['SD_boundary_flux_not_unit_default']=ward_sd({key:value for key,value in sd.items()if key!='integrate_boundary_face'})['status']=='MISSING_INPUT'
# Test density derivative algebra independently of a physical integration claim.
nonflat={**sd,'coordinate_density':1+x*x}
sdr=ward_sd(nonflat);controls['nonconstant_measure_Jacobian_retained']=sdr['status']=='PASS'and equal(sdr['measure_jacobian_contact'],2*x*x)is True
# Whole-word branch keeps each original arithmetic mark, no compressed evolution.
V=d['quotient_lift'];P=clean(V.H*matrix(d['physical_selector_form'])*V);C=clean(V.H*matrix(d['boundary_pairing_form'])*V)
forms={key:clean(V.H*out['unnormalized_forms'][key]*V)for key in MARKS}
placement={'retained_embedding':s.eye(2),'whole_P':P,'whole_cap':C,'outer_scalar':scalar(d['scalar_factor']),
 'reference_scalar':scalar(d['reference_scalar']),'whole_word_source':lambda recipe:forms[recipe],
 'original_whole_word_recipes':{key:key for key in MARKS},'observed_full_ordered_forms':forms,
 'observed_retained_source_forms':forms,'original_domain_transport_relations':[{'actual':s.eye(2),'required':s.eye(2)}],
 'original_factor_transport_relations':[{'actual':scalar(d['reference_scalar']),'required':scalar(d['reference_scalar'])}]}
controls['whole_word_placement_without_semigroup']=whole_placement(placement)['status']=='PASS'
wrong={**placement,'observed_retained_source_forms':{**forms,'ZUY':P*forms['ZUY']*P}}
controls['premature_internal_projection_detected']=whole_placement(wrong)['status']=='FAIL'
fine=[{'id':(i,j),'i':i,'unnormalized_mass':T[i,j]/3}for i in range(3)for j in range(3)]
coarse=[{'id':i,'unnormalized_mass':s.Rational(1,3)}for i in range(3)]
proj={'kind':'exact_original_weighted_pushforward','fine_atoms':fine,'coarse_atoms':coarse,
 'original_retained_map':lambda a:a['i'],'required_limit_scope':'finite_weighted_blocking'}
controls['native_transfer_all_label_pushforward']=regulator_limit(proj)['status']=='PASS'
controls['empty_pushforward_not_passed']=regulator_limit({**proj,'fine_atoms':[],'coarse_atoms':[]})['status']=='MISSING_INPUT'
tail={'kind':'original_finite_transfer_cauchy_tail','actual_transfer':T,'actual_stationary_projector':Pi,'n':6,'required_limit_scope':'finite_transfer_iterations'}
tr=regulator_limit(tail);controls['native_transfer_analytic_entire_tail']=tr['status']=='PASS'and tr['q']==s.Rational(64,729)
controls['finite_tail_not_UV_proof']=regulator_limit({**tail,'required_limit_scope':'physical_field_UV'})['status']=='MISSING_INPUT'
controls['missing_native_check_inputs_stay_missing']=all(f({})['status']=='MISSING_INPUT'for f in FUNCTIONS.values())
# Exercise the actual registered consumer on exact represented arithmetic
# source entries. Only this explicitly arithmetic fixture is rationalized.
def rational_json(x):
 if isinstance(x,list):return [rational_json(y)for y in x]
 if isinstance(x,dict):return {k:rational_json(v)for k,v in x.items()}
 return str(x)if isinstance(x,s.Basic)else x
bound_bundle=rational_json(b)
digest=hashlib.sha256(json.dumps(bound_bundle,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
entries=[]
for atom in bound_bundle['fiber_records']:
 w=scalar(atom['unnormalized_scalar_weight']);forms={}
 for z,wkey in zip(MARKS,('W0','WU','WY','WUY')):
  A=w*matrix(atom['ordered_source_forms'][z])
  forms[wkey]={'value':A}
 entries.append({'configuration_id':atom['original_configuration_id'],'domain_id':bound_bundle['domain_id'],
  'native_inputs_sha256':bound_bundle['native_inputs_sha256'],'source_bundle_sha256':digest,
  'original_measure_mass':atom['original_measure_mass'],'contributions':forms})
domain={'whole_source_bundle':bound_bundle,'domain_id':bound_bundle['domain_id'],'_sha256':bound_bundle['native_inputs_sha256'],
 'original_check_inputs':{IDS[0]:{'source_bundle_sha256':digest,'quotient_lift':d['quotient_lift'],'null_basis':d['null_basis']}}}
current=out['values'];callback=ensemble_test_registry(domain)[IDS[0]]
controls['registered_callback_actual_current_entries_and_result']=callback(current,entries,domain)['status']=='PASS'
bad_result={**current,'ZU':current['ZU']+1}
controls['registered_callback_changed_current_Z_detected']=callback(bad_result,entries,domain)['status']=='FAIL'
# Follow the real Q plugin -> adapter -> registry path, rather than a separate
# registry lookalike. The data remain the SAME stated arithmetic control.
import importlib.util,os
os.environ['TFPT_FIELD_DOMAIN']=str(HERE.parent/'domain/native-diagnostic-v1.json')
runtime=HERE.parent/'runtime/field_problem.py'
spec=importlib.util.spec_from_file_location('tfpt_q_actual_field_plugin_control',runtime)
plugin=importlib.util.module_from_spec(spec);spec.loader.exec_module(plugin)
actual_registry=plugin.ADAPTER.ensemble_test_registry(domain)
controls['actual_runtime_all_five_math_callbacks_registered']=all(cid in actual_registry and callable(actual_registry[cid])for cid in IDS)
controls['actual_runtime_positive_current_callback']=actual_registry[IDS[0]](current,entries,domain)['status']=='PASS'
controls['actual_runtime_changed_current_Z_rejected']=actual_registry[IDS[0]](bad_result,entries,domain)['status']=='FAIL'
result={'scope':'existing arithmetic-control matrices and old v330 transfer; symbolic R048 quadratic contact. No primitive physical source values.',
 'passed':all(controls.values()),'controls':controls,'count':len(controls),'actual_transfer_q':str(tr['q']),'actual_transfer_tail_n6':str(tr['tail_bound']),
 'native_missing_check_results':run_ensemble_checks({},[],{}),'source_sha256':hashlib.sha256((HERE/'original_checks.py').read_bytes()).hexdigest()}
(HERE/'controls.json').write_text(json.dumps(jsonable(result),indent=2)+'\n')
print(json.dumps(jsonable(result),indent=2));assert result['passed']
