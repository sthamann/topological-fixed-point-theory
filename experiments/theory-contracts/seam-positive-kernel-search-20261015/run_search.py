#!/usr/bin/env python3
"""Construct, certify and distinguish positive kernels in a declared search class.

This bounded campaign ends with exact finite candidates and a discriminating
source question. It never feeds a candidate prediction back as a source fact.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
from fractions import Fraction
from kernel_search import MomentSearch, OBJECTIVES, TRANSFER_KINDS, normalize_fixed_moments, load_source_control
from finite_completion import complete
from word_algebra import WordAlgebra
from linear_relations import normalize_linear_relations


def dump(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')


def fingerprint():
    folder=Path(__file__).parent
    names=['run_search.py','kernel_search.py','finite_completion.py','gns_completion.py','exact_certificate.py','word_algebra.py','linear_relations.py','response_search.py']
    return {name:hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in names}


def discriminate(first, second, payload):
    if not first.get('ok') or not second.get('ok'):
        return {'status':'NO_TWO_CERTIFIED_COMPLETIONS'}
    if first['input_sha256'] != second['input_sha256']:
        raise ValueError('Cannot compare candidates fitted to different supplied data')
    a={tuple(p['word']):p['value'] for p in first['candidate_predictions']}
    b={tuple(p['word']):p['value'] for p in second['candidate_predictions']}
    differences=[]
    for word in sorted(a.keys() & b.keys(),key=lambda w:(len(w),w)):
        av=tuple(Fraction(a[word][p]) for p in ('re','im'))
        bv=tuple(Fraction(b[word][p]) for p in ('re','im'))
        if av!=bv:
            differences.append({'word':list(word),'first':a[word],'second':b[word]})
    complete_prefix=set(map(tuple,(r['word'] for r in payload['moments']))) == set(WordAlgebra(payload['m']).words(2*payload['degree']))
    return {'status':'EXACT_DISTINGUISHING_WORDS_FOUND' if differences else 'NO_DIFFERENCE_IN_TESTED_PROBE_FAMILY',
            'same_input_moments_replayed':len(payload['moments']),
            'all_words_through_degree_replayed':2*payload['degree'] if complete_prefix else None,
            'shortest_found_degree':len(differences[0]['word']) if differences else None,
            'shortest_certified':bool(differences and complete_prefix and len(differences[0]['word'])==2*payload['degree']+1),
            'differences':differences,
            'question_de':'Welche originale Nahtrelation legt die Transferantwort nach zwei Feldoperationen fest und koppelt sie an die Einfeld-Antwort?',
            'question_en':'Which original seam relation fixes transfer after two field operations and couples it to the single-field response?',
            'action':'Check a source-derived composite-field/transfer relation. Neither prediction becomes a new premise merely because a solver produced it.'}


def campaign(out, m=5, degree=2, seed=7, transfer_floor=0.0, cluster_gap=0.0,*,objective='feasibility',transfer_ceiling=1.0,transfer_kind='vacuum_effect',reference_eigenvalue=1.0,fixed_moments=None,linear_relations=None,source_observables=None,completion_representation='minimal_fock',withheld_response=None):
    out=Path(out)
    start=time.monotonic()
    code=fingerprint()
    if not 0 <= cluster_gap < 1:
        raise ValueError('cluster_gap must be in [0,1)')
    if objective not in OBJECTIVES: raise ValueError('Unknown search objective')
    if completion_representation not in ('minimal_fock','faithful_car_gns'):
        raise ValueError('Unknown completion representation')
    fixed_moments=normalize_fixed_moments(m,degree,fixed_moments)
    linear_relations=normalize_linear_relations(m,degree,linear_relations)
    config={'m':m,'degree':degree,'seed':seed,'transfer_floor':transfer_floor,'cluster_gap':cluster_gap,
            'transfer_ceiling':transfer_ceiling,'transfer_kind':transfer_kind,
            'reference_eigenvalue':reference_eigenvalue,
            'objective':objective,'solver':'SCS','fixed_moments':fixed_moments,'linear_relations':linear_relations,'source_observables':source_observables or {},
            'completion_representation':completion_representation}
    if withheld_response is not None:
        config['withheld_response']=withheld_response
    identity=hashlib.sha256(json.dumps({'config':config,'code':code},sort_keys=True).encode()).hexdigest()
    reserved=('moments.json','min-trace.json','max-trace.json','min-probe.json','max-probe.json','feasibility.json','next-source-question.json','withheld-response.json')
    if not (out/'summary.json').exists() and any((out/name).exists() for name in reserved):
        raise ValueError('Existing campaign artifacts without a matching summary; choose a fresh output directory')
    if (out/'summary.json').exists():
        previous=json.loads((out/'summary.json').read_text())
        if previous.get('campaign_identity')!=identity:
            raise ValueError('Existing output belongs to different code/configuration; choose a fresh output directory')
        # Reuse only if every output hash still matches; no stale or corrupted cache.
        if all((out/name).exists() and hashlib.sha256((out/name).read_bytes()).hexdigest()==digest
               for name,digest in previous.get('artifacts',{}).items()):
            return previous
        raise ValueError('Stored campaign artifacts changed; choose a fresh output directory')
    search=MomentSearch(m,degree,transfer_floor,fixed_moments=fixed_moments,linear_relations=linear_relations,
                         transfer_ceiling=transfer_ceiling,transfer_kind=transfer_kind,
                         reference_eigenvalue=reference_eigenvalue,
                         source_observables=source_observables)
    if withheld_response is not None:
        from response_search import normalize_response
        normalize_response(search,withheld_response)
    payload=search.solve(objective=objective,seed=seed,solver='SCS')
    payload['cluster_gap']=cluster_gap
    payload['cluster_scope']='Finite transfer clustering is source-motivated; a positive numerical gap is a declared search bound, not a physical time scale. It is imposed during completion, not the moment relaxation.'
    dump(out/'moments.json',payload)
    certified=[]
    if payload.get('moments'):
        probe=None
        completion_objectives=('min-trace','max-trace')
        if completion_representation=='faithful_car_gns':
            from gns_completion import first_missing_probe
            probe=first_missing_probe(payload)
            completion_objectives=('min-probe','max-probe') if probe is not None else ('feasibility',)
        for completion_objective in completion_objectives:
            if completion_representation=='faithful_car_gns':
                result=complete(payload,completion_objective,representation=completion_representation,probe_word=probe)
            elif transfer_kind=='vacuum_effect' and transfer_ceiling==1 and reference_eigenvalue==1:
                result=complete(payload,completion_objective)
            else:
                result={'objective':completion_objective,'ok':False,
                        'verdict':'INCONCLUSIVE_UNSUPPORTED_TRANSFER_TYPE',
                        'reason':'Minimal pure-Fock completion handles vacuum_effect with ceiling1 and reference_eigenvalue1 only; scaled-reference, thermal GNS/bounded positive transfers need a compatible representation and recovery, without silently rescaling, forcing a pure state or contraction.'}
            dump(out/(completion_objective+'.json'),result)
            certified.append(result)
    discriminator=discriminate(*certified,payload) if len(certified)==2 else {'status':'NO_TWO_COMPLETIONS_FOR_COMPARISON'}
    dump(out/'next-source-question.json',discriminator)
    accepted_count=sum(x.get('ok',False) for x in certified)
    verdict=('CONSTRUCTIVE_FINITE_KERNEL_SEARCH_WITH_EXACT_DISCRIMINATOR'
             if accepted_count==2 and discriminator['status']=='EXACT_DISTINGUISHING_WORDS_FOUND'
             else 'CERTIFIED_FINITE_COMPLETION_WITHOUT_DISCRIMINATOR' if accepted_count
             else 'SEARCH_INCONCLUSIVE_NO_CERTIFIED_COMPLETION')
    summary={'schema':'tfpt.positive-kernel-campaign.v1','campaign_identity':identity,
             'configuration':config,'code_sha256':code,
             'search_restrictions':{'objective':objective,'objective_role':payload['objective_role'],
                 'objective_definition':payload['objective_definition'],
                 'global_transfer_floor':transfer_floor,'completion_cluster_gap':cluster_gap,
                 'global_transfer_ceiling':transfer_ceiling,'transfer_kind':transfer_kind,
                 'reference_eigenvalue':reference_eigenvalue,
                 'reference_eigenvalue_scope':payload['reference_eigenvalue_scope'],
                 'transfer_scope':payload['transfer_scope'],
                 'completion_objectives':'min/max trace explore compatible finite completions; they do not select physical sources',
                 'scope':'Explicit search/calibration restrictions, not a physical source selector or native scale'},
             'completion_ansatz':{'representation':'minimal irreducible CAR Fock space of dimension2^m',
                 'transfer_type_requirement':'vacuum_effect with ceiling1 and reference_eigenvalue1 only; no scaled-reference or thermal GNS completion claimed',
                 'state_requirement':'all inferred occupations within1e-6 of0 or1 in the given mode basis',
                 'status':'additional realization ansatz, not selected by feasibility or inferred from nonflat moments',
                 'mixed_state_handling':'inconclusive if this completion cannot handle the state; no fallback forcing pure occupations'},
             'inconclusive_completion_is_not_source_exclusion':True,
             'conditional_normalized_state':True,'absolute_partition_function_recovered':False,
             'state_normalization_scope':payload['state_normalization_scope'],
             'solver_status':payload['solver_status'],
             'certified_finite_completions':accepted_count,
             'completions':[{'objective':x['objective'],'ok':x.get('ok',False),
                             'verdict':x['verdict'],'reason':x.get('reason'),
                             'dimension':x.get('dimension'),
                             'finite_clustering':x.get('finite_clustering'),
                             'finite_unitary_field_time':x.get('finite_selfadjoint_H_minus_log_T_exists',False)} for x in certified],
             'next_question':discriminator,
             'original_source_derivation':False,'physical_marking_selected':False,
             'physical_quarter_holonomy_identified':False,'event_instrument_derived':False,
             'verdict':verdict,
             'scope':'Finite CAR search and an explicitly additional minimal Fock completion; not a physical seam construction.',
             'elapsed_seconds':time.monotonic()-start}
    names=['moments.json','next-source-question.json']+[x['objective']+'.json' for x in certified]
    if withheld_response is not None:
        from response_search import evaluate
        response_result=evaluate(search,withheld_response,solver='SCS')
        response_result['completion_cluster_gap_not_in_response_relaxation']=cluster_gap
        dump(out/'withheld-response.json',response_result)
        names.append('withheld-response.json')
        summary['withheld_response']={
            'artifact':'withheld-response.json',
            'affine_status':response_result['affine_certificate']['status'],
            'identity_proved':response_result['affine_certificate']['identity_proved'],
            'exact_value':response_result['affine_certificate'].get('value'),
            'numerical_endpoints':response_result['numerical_endpoints'],
            'source_derivation_or_physical_uniqueness_proved':False}
    summary['artifacts']={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in names}
    if completion_representation=='faithful_car_gns':
        summary['search_restrictions']['completion_objectives']='min/max of a generated missing Hermitian response explore ambiguity; not physical source selection'
        summary['scope']='Finite faithful CAR-cyclic GNS completion; not a physical seam construction or original source selection.'
        summary['completion_ansatz']={
            'representation':'faithful CAR cyclic GNS of dimension4^m, currently at most two modes',
            'state_requirement':'complete positive definite CAR-only moment Gram; no pure state or Wick rule',
            'transfer_requirement':'T preserves the stipulated CAR cyclic space and satisfies declared bounds',
            'proof':'exact GNS Gram, bounds, adjoints, CAR, quarter, reference and all supplied words replayed',
            'status':'additional finite realization ansatz; not derived from truncated flatness or P1',
            'source_query':'first missing Hermitian one-transfer word, followed by exact candidate separation'}
    summary['elapsed_seconds']=time.monotonic()-start
    dump(out/'summary.json',summary)
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--m',type=int,default=5)
    p.add_argument('--degree',type=int,default=2)
    p.add_argument('--seed',type=int,default=7)
    p.add_argument('--transfer-floor',type=float,default=0)
    p.add_argument('--transfer-ceiling',type=float,default=1)
    p.add_argument('--reference-eigenvalue',type=float,default=1)
    p.add_argument('--transfer-kind',choices=TRANSFER_KINDS,default='vacuum_effect')
    p.add_argument('--cluster-gap',type=float,default=0)
    p.add_argument('--objective',choices=OBJECTIVES,default='feasibility')
    p.add_argument('--source-control',type=Path,help='JSON source fixed_moments/linear_relations/source_observables')
    p.add_argument('--withheld-response',type=Path,help='Hermitian word polynomial JSON with terms only; evaluated after fitting, never imposed as source data')
    p.add_argument('--completion-representation',choices=('minimal_fock','faithful_car_gns'),default='minimal_fock')
    args=p.parse_args()
    result=campaign(args.out,args.m,args.degree,args.seed,args.transfer_floor,args.cluster_gap,
                    objective=args.objective,transfer_ceiling=args.transfer_ceiling,
                    transfer_kind=args.transfer_kind,reference_eigenvalue=args.reference_eigenvalue,
                    completion_representation=args.completion_representation,
                    withheld_response=json.loads(args.withheld_response.read_text()) if args.withheld_response else None,
                    **load_source_control(args.source_control))
    display={k:result[k] for k in ('verdict','configuration','search_restrictions','completion_ansatz','solver_status','certified_finite_completions','completions','elapsed_seconds')}
    display['configuration']={k:v for k,v in result['configuration'].items() if k not in ('fixed_moments','linear_relations','source_observables')}
    display['configuration']['fixed_moment_count']=len(result['configuration']['fixed_moments'])
    display['configuration']['linear_relation_count']=len(result['configuration']['linear_relations'])
    if 'withheld_response' in result: display['withheld_response']=result['withheld_response']
    print(json.dumps(display,indent=2))
    print(json.dumps(result['next_question'],indent=2))


if __name__=='__main__': main()
