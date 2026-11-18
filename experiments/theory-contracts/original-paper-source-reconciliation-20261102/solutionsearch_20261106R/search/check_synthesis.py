"""Method controls using the ALREADY GIVEN v330 finite transfer, not TFPT completion."""
from pathlib import Path
import hashlib,json,shutil
import sympy as s
from synthesis import *
HERE=Path(__file__).resolve().parent
SOURCE=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4/verification/v330_qgamb_admissible_measure.py')
# Exact matrix already proved in original BLOCKING-PATH.txt section1.
T=s.Matrix([[1651,1267,1456],[1267,1651,1456],[1456,1456,1462]])/4374
Pi=s.ones(3)/3

def actual_candidate(A):return {'matrix':A,'complete':True}
def corrupted_candidate(A):return {'matrix':A+s.eye(3),'complete':True}
def unproved_candidate(A):return {'matrix':A+s.Symbol('unproved_constructor_coefficient')*s.eye(3),'complete':True}
def missing_primitive():return Unavailable('The actual primitive source operation has not been evaluated')

def finite_original_checker(candidate,context):
    A=candidate['matrix'];conditions={'original_transfer_match':A-context['T'],
        'actual_row_sum':A*s.ones(3,1)-s.ones(3,1),'actual_stationary_contact':Pi*A-Pi}
    rows=[]
    for name,delta in conditions.items():
        signs=[s.simplify(v).is_zero for v in delta]
        status='FAIL'if False in signs else('PASS'if all(v is True for v in signs)else'NOT_PROVED')
        rows.append({'id':name,'status':status,'actual_residual':str(delta)})
    return {'obligations':rows,'evidence':{'self_claim_complete_ignored':True,'class':'given v330 finite transfer only'}}

OPS=[Operation('wrap_original_transfer',('OriginalTransfer',),'WholeSourceCandidate',actual_candidate,'existing finite transfer'),
     Operation('corrupt_transfer_mutation',('OriginalTransfer',),'WholeSourceCandidate',corrupted_candidate,'method mutation only'),
     Operation('unproved_constructor',('OriginalTransfer',),'WholeSourceCandidate',unproved_candidate,'unproved symbolic coefficient, no chosen physical parameter'),
     Operation('missing_actual_primitive',(),'WholeSourceCandidate',missing_primitive,'explicit missing original operation')]
REQUIRED=('original_transfer_match','actual_row_sum','actual_stationary_contact')

def make(out,value=T,version='actual-original-v330'):
    return Synthesizer(OPS,{'given_transfer':InputBinding('OriginalTransfer',value,version)},finite_original_checker,
        REQUIRED,'finite_original_v330_constructor_method_control_NOT_TFPT_F1',hashlib.sha256(SOURCE.read_bytes()).hexdigest(),out,
        context={'T':T},verifier_dependencies=[SOURCE])

def fault_dictionary_checker(candidate,context):return {'error':'Actual primitive matching information not supplied'}
def checker_for_reference(reference):
    def checker(candidate,context):return finite_original_checker(candidate,{'T':reference})
    return checker
class BoundReferenceChecker:
    def __init__(self,reference):self.reference=reference;self.trusted={'actual':finite_original_checker}
    def __call__(self,candidate,context):return self.trusted['actual'](candidate,{'T':self.reference})
def type_guard_for_reference(reference):
    def guard(candidate):return candidate['matrix']==reference
    return guard

def main(out):
    if out.exists():raise RuntimeError('Use a fresh method-control output; preserve previous native runs')
    engine=make(out)
    goal={'op':'apply','id':'wrap_original_transfer','args':[{'op':'input','id':'given_transfer','args':[]}]}
    engine.enqueue(goal);engine.enqueue(goal)
    # Model-supplied code/import commands are neither registered nor executable.
    with (out/'queue.jsonl').open('a')as f:f.write(canonical({'op':'apply','id':'os.system','args':[]})+'\n')
    first=engine.run(max_cost=2,max_asts=32)
    records=list(engine.records.values());goals=[r for r in records if r['output_type']=='WholeSourceCandidate']
    by_op={r['ast']['id']:r for r in goals}
    resume=make(out);second=resume.run(max_cost=2,max_asts=32)
    changed=make(out,value=T+s.eye(3),version='actual-original-v330')
    third=changed.run(max_cost=2,max_asts=32)
    direct=make(out/'direct-path')
    limited=direct.run(max_cost=1,max_asts=1)
    proposal=direct.evaluate_program(goal,max_cost=64)
    invalid=direct.evaluate_program({'op':'apply','id':'os.system','args':[]},max_cost=64)
    direct_status=json.loads((direct.out/'status.json').read_text())
    fault=Synthesizer(OPS,{'given_transfer':InputBinding('OriginalTransfer',T,'actual-original-v330')},fault_dictionary_checker,
        REQUIRED,'finite_original_v330_constructor_method_control_NOT_TFPT_F1','actual original finite transfer',out/'fault-dictionary',context={'T':T})
    fault_record=fault.evaluate_program(goal)
    closure_out=out/'bound-verifier-regression'
    def bound(checker):
        return Synthesizer(OPS,{'given_transfer':InputBinding('OriginalTransfer',T,'unchanged-declared-original-input')},checker,
            REQUIRED,'finite_original_v330_constructor_method_control_NOT_TFPT_F1','same-declared-fingerprint',closure_out,context=None)
    old_checker=bound(checker_for_reference(T));old_goal=old_checker.evaluate_program(goal)
    new_checker=bound(checker_for_reference(T+s.eye(3)));new_goal=new_checker.evaluate_program(goal)
    bound_identity=(old_checker.fingerprint!=new_checker.fingerprint and old_goal['goal']['status']=='SCOPED_HIT'and new_goal['goal']['status']=='CANDIDATE_FAIL'
      and semantic_fingerprint(BoundReferenceChecker(T))!=semantic_fingerprint(BoundReferenceChecker(T+s.eye(3)))
      and semantic_fingerprint({'T':T})!=semantic_fingerprint({'T':T+s.eye(3)})
      and semantic_fingerprint(type_guard_for_reference(T))!=semantic_fingerprint(type_guard_for_reference(T+s.eye(3))))
    guards={
      'actual_original_finite_checker_finds_correct_program':by_op['wrap_original_transfer']['goal']['status']=='SCOPED_HIT',
      'executable_unproved_program_remains_partial':by_op['unproved_constructor']['goal']['status']=='PARTIAL',
      'corrupted_program_fails_actual_obligations':by_op['corrupt_transfer_mutation']['goal']['status']=='CANDIDATE_FAIL',
      'missing_primitive_is_local_not_global_stop':by_op['missing_actual_primitive']['evaluation_status']=='MISSING_INPUT'and first['scoped_hits']==1,
      'self_attestation_has_no_effect':all(r['goal'].get('candidate_self_attestation_used')is False for r in goals if r['evaluation_status']=='VALUE_BUILT'),
      'cost_fair_order':all(a['cost']<=b['cost']for a,b in zip(records,records[1:])),
      'queued_duplicates_not_reevaluated':len(records)==5 and len({r['ast_hash']for r in records})==5,
      'unsafe_operation_rejected_without_execution':(out/'queue-rejections.jsonl').exists(),
      'resume_skips_completed_AST_fingerprint':second['new_evaluations']==0,
      'changed_native_input_changes_fingerprint':third['search_fingerprint']!=first['search_fingerprint']and third['new_evaluations']==5,
      'changed_input_does_not_reuse_old_success':third['scoped_hits']==0,
      'finite_grammar_coverage_labeled':first['coverage']=='complete_for_declared_registry_inputs_and_cost_bound'and first['coverage_is_all_physical_sources']is False,
      'goal_scope_and_verifier_code_hash_persisted':first['target_scope'].startswith('finite_original_v330')and len(first['verifier_code_sha256'])==64,
      'direct_proposal_works_after_enumeration_budget':limited['stop_reason']=='finite_ast_budget'and proposal['goal']['status']=='SCOPED_HIT',
      'direct_proposal_uses_actual_same_verifier_path':proposal['goal']['candidate_self_attestation_used']is False,
      'direct_proposal_does_not_claim_enum_coverage':direct_status['coverage']=='partial_finite_search_budget'and direct_status['max_cost']==1,
      'unsafe_direct_proposal_is_local_rejection':invalid['evaluation_status']=='INVALID_AST'and direct_status['scoped_hits']==1,
      'missing_fault_dictionary_stays_partial':fault_record['goal']['status']=='PARTIAL'and fault_record['goal']['complete_for_target_scope']is False,
      'actual_closure_bound_checker_context_and_type_guard_cache_identity':bound_identity,
    }
    result={'passed':all(guards.values()),'guards':guards,'first':first,'resume':second,'changed_input':third,
            'control_scope':'Software method controls on the given old v330 finite transfer. No F1/global TFPT success fixture.',
            'native_source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            'search_core_sha256':hashlib.sha256((HERE/'synthesis.py').read_bytes()).hexdigest()}
    (HERE/'method-controls.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2));assert result['passed']
if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=HERE/'control-run')
    main(parser.parse_args().out)
