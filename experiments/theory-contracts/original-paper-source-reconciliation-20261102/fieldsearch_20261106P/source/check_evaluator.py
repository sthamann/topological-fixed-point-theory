"""End-to-end arithmetic controls, explicitly NOT TFPT physical source inputs."""
import copy
import json
import tempfile
from pathlib import Path
import numpy as np
from evaluator import (evaluate_source, capability_report, integration_plan,
                       evaluate_contribution, aggregate_contributions,
                       run_ensemble_checks, ensemble_test_registry, MARKS)


def enc(a):
    a=np.asarray(a,dtype=complex)
    return np.stack([a.real,a.imag],axis=-1).tolist()


def pair(z):
    return [complex(z).real,complex(z).imag]


def form(a):
    B=np.array([[2,0,0],[0,1,0]],dtype=complex)
    return enc(B.conj().T@a@B)


def main():
    assembled={'candidate_id':'ARITHMETIC_ONLY','native_inputs_sha256':'ARITHMETIC_ONLY',
               'domain_id':'ARITHMETIC_ONLY',
               'h_N_role':'original_tangential_one_particle_NOT_joint_OS_generator'}
    candidate={'candidate_id':'ARITHMETIC_ONLY'}
    v=np.array([1,1j])/np.sqrt(2)
    P=np.outer(v,v.conj())
    cap=np.array([[2,1+1j],[3-1j,4]],dtype=complex)
    q=[np.array(x,dtype=complex) for x in
       ([[2,3j],[-1j,5]],[[3,1],[-2j,4]],[[1,2j],[3,-2]],[[4,-1j],[2,7]])]
    rows=[]
    factors=[(2+1j,3-2j),(-1+2j,1+3j)]
    for i,(mass,weight) in enumerate(factors):
        rows.append({'original_configuration_id':f'ARITHMETIC_ATOM_{i}',
                     'native_candidate_id':'SAME_VISIBLE_CONFIGURATION',
                     'original_measure_binding':'arithmetic control, not physical',
                     'whole_fiber_binding':'arithmetic control, not physical',
                     'original_measure_mass':pair(mass),
                     'unnormalized_scalar_weight':pair(weight),
                     'ordered_source_forms':{k:form((i+1)*a) for k,a in zip(MARKS,q)}})
    b={'schema':'tfpt.original-blocked-source.v1','input_role':'arithmetic_control_only',
       'representation':'original_exact_finite_fiber_sum',
       'matrix_convention':'sesquilinear_forms_in_original_word_basis',
       'ensemble_scope':'whole_original_ensemble',
       'native_inputs_sha256':'ARITHMETIC_ONLY','domain_id':'ARITHMETIC_ONLY',
       'source_refs':[{'id':'arithmetic-control-only-not-a-theory-source'}],
       'component_bindings':{k:'arithmetic control only' for k in
           ('all_bosonic_fibres','whole_physical_selector','boson_generator_or_history',
            'fermion_pairing_and_higher_vertices','reference_and_eta','absolute_rest',
            'same_original_time_and_source_marks')},
       'scalar_factor':pair(2-3j),'reference_scalar':pair(3+4j),
       'configuration_unit':{'kind':'original_retained_field_fiber_with_whole_history_endpoint_packet',
                             'native_map_binding':'arithmetic control only'},
       'measure_rule':{'kind':'original_finite_atomic_measure','source_binding':'arithmetic control only'},
       'retained_words':[[],['ARITHMETIC_A'],['NULL_CONTROL']],
       'time':{'id':'ARITHMETIC_TIME','kind':'original_euclidean_ordered','T':3,'t_U':1,'t_Y':2},
       'gram_form':enc(np.diag([4,1,0])), 'physical_selector_form':form(P),
       'boundary_pairing_form':form(cap), 'fiber_records':rows}
    domain={'whole_source_bundle':b}
    Path(__file__).with_name('arithmetic-control-bundle.json').write_text(json.dumps(b,indent=2)+'\n')
    got=evaluate_source(candidate,assembled,domain)
    multiplier=sum((i+1)*mass*weight for i,(mass,weight) in enumerate(factors))
    expected={k:(2-3j)/(3+4j)*multiplier*np.trace(cap@P@a@P) for k,a in zip(MARKS,q)}
    controls={'complex_unnormalized_fiber_sum_and_nonorthogonal_null_quotient':
        all(np.allclose(complex(*got[k]),value) for k,value in expected.items()),
        'arithmetic_control_never_physical':not got['physical_source_complete']}
    z,zu,zy,zuy=[expected[k] for k in MARKS]
    controls['connected_subtraction_after_whole_integration']=np.allclose(
        complex(*got['normalized']['connected_UY']),zuy/z-zu*zy/z**2)
    # Same result from the already aggregated whole endpoint matrices.
    direct=copy.deepcopy(b)
    direct['representation']='original_exact_blocked_endpoint_source'
    direct['ordered_source_forms']={k:form(multiplier*a) for k,a in zip(MARKS,q)}
    res=evaluate_source(candidate,assembled,{'whole_source_bundle':direct})
    controls['tower_aggregation_matches_whole_kernel_contraction']=all(np.allclose(res[k],got[k]) for k in MARKS)
    bad=copy.deepcopy(b)
    bad['fiber_records'][0]['ordered_source_forms']['ZUY'][2][0]=[1,0]
    controls['null_class_cannot_acquire_marked_value']=evaluate_source(candidate,assembled,{'whole_source_bundle':bad})['status']=='original_source_input_inconsistent'
    missing=evaluate_source(candidate,assembled,{})
    controls['no_h_N_to_joint_source_fallback']=all(missing[k] is None for k in MARKS) and not missing['physical_source_complete']
    bad=copy.deepcopy(b)
    del bad['fiber_records'][0]['unnormalized_scalar_weight']
    controls['no_unit_weight_fallback']=evaluate_source(candidate,assembled,{'whole_source_bundle':bad})['status']=='original_source_input_inconsistent'
    controls['empty_bundle_cannot_pass_preflight']=capability_report({'whole_source_bundle':{}})['numeric_bundle_ready'] is False
    relabelled=copy.deepcopy(b)
    relabelled['input_role']='source_reviewed_original'
    controls['role_string_cannot_assert_physical_completion']=evaluate_source(candidate,assembled,{'whole_source_bundle':relabelled})['physical_source_complete'] is False
    controls['valid_numeric_consumer_does_not_certify_physical_source']=capability_report(domain)['physical_source_complete'] is False and capability_report(domain)['numeric_bundle_ready'] is True
    plan=integration_plan(domain)
    unit_assembled={k:v for k,v in assembled.items() if k!='candidate_id'}
    entries=[evaluate_contribution(item['configuration'],unit_assembled,domain) for item in plan['items']]
    aggregated=aggregate_contributions(entries,domain)
    controls['contribution_driver_matches_original_finite_sum']=all(np.allclose(aggregated[k],got[k]) for k in MARKS)
    controls['same_visible_configuration_not_deduplicated']=len(plan['items'])==2 and entries[0]['native_candidate_id']==entries[1]['native_candidate_id'] and entries[0]['configuration_id']!=entries[1]['configuration_id']
    controls['partial_archive_not_promoted_to_whole_Z']=all(aggregate_contributions(entries[:1],domain)[k] is None for k in MARKS)
    cross=copy.deepcopy(entries)
    cross[0]['domain_id']='OTHER_DOMAIN'
    controls['cross_domain_packets_cannot_be_aggregated']=aggregate_contributions(cross,domain)['status']=='original_contribution_binding_mismatch'
    controls['explicit_source_unit_cannot_override_native_id']=evaluate_contribution(plan['items'][0]['configuration'],assembled,domain)['status']=='original_configuration_binding_mismatch'
    controls['only_implemented_ensemble_callback_registered']=set(ensemble_test_registry(domain))=={'original_domain_coverage'} and run_ensemble_checks(aggregated,entries,domain)[0]['status']=='PASS'
    missing_weight=copy.deepcopy(b)
    del missing_weight['fiber_records'][0]['unnormalized_scalar_weight']
    controls['bad_fiber_is_typed_error_not_crash']=evaluate_contribution(plan['items'][0]['configuration'],unit_assembled,{'whole_source_bundle':missing_weight})['status']=='original_fiber_input_inconsistent'
    with tempfile.TemporaryDirectory() as td:
        fp=Path(td)/'invalid.json'
        fp.write_text('[]')
        controls['json_list_is_reported_invalid_not_crash']=capability_report({'whole_source_bundle':str(fp)})['numeric_bundle_ready'] is False
    result={'scope':'arithmetic consumer controls only; no physical TFPT values',
            'controls':{k:bool(v) for k,v in controls.items()},
            'pass':all(controls.values()),'arithmetic_output':got}
    path=Path(__file__).with_name('evaluator-controls.json')
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'pass':result['pass'],'controls':result['controls']},indent=2))
    if not result['pass']: raise SystemExit(1)

if __name__=='__main__': main()
