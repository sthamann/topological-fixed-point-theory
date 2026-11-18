#!/usr/bin/env python3
"""Existing search plugin plus an original-mass integration entry point."""
from pathlib import Path
import argparse,importlib.util,json,os,sys,hashlib
from repo_paths import repo_root,ENGINE_RELATIVE
ROOT=repo_root()
ENGINE=ROOT/ENGINE_RELATIVE

def write_json(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n')

def integrate(plugin,args,preflight):
    adapter=plugin.ADAPTER
    registration=plugin.ensemble_registration()
    cap=plugin.capability_report()
    plan=adapter.integration_plan(plugin.DOMAIN) if adapter is not None else {'ready':False,'items':[],'missing_operations':['original source adapter']}
    # Call the real contribution path even when the native producer is absent.
    # Its typed missing operation is part of the actual launch result.
    probe=plan['items'][0]['configuration'] if plan.get('items') else {'schema':'tfpt.original-source-unit.v1','original_configuration_id':None}
    assembled={'candidate_id':probe.get('native_candidate_id'),
               'native_inputs_sha256':plugin.DOMAIN['_sha256'],'domain_id':plugin.DOMAIN['domain_id'],
               'h_N_role':'not_supplied_not_aliased_to_joint_generator'}
    attempt=adapter.evaluate_contribution(probe,assembled,plugin.DOMAIN) if adapter is not None else {'status':'MISSING_IMPLEMENTATION','missing_operations':['original source adapter']}
    values_ready=all(attempt.get('contributions',{}).get(k)is not None for k in ('W0','WU','WY','WUY'))
    ready=cap.get('ready') is True and plan.get('ready') is True and values_ready and not registration['missing_implementation_ids']
    preflight.update(integration_startable=ready,ensemble_registration=registration,
                     actual_contribution_attempt=attempt,integration_plan={key:value for key,value in plan.items() if key!='items'},
                     original_atom_count=len(plan.get('items',[])),
                     GA_selection_used_for_integration=False)
    write_json(args.out/'preflight.json',preflight)
    if not ready:
        result={'status':'ORIGINAL_SOURCE_INTEGRATION_NOT_STARTED','preflight':preflight,
                'whole_physical_search_completed':False}
        write_json(args.out/'result.json',result);print(json.dumps(result,indent=2));return 2
    if plan['kind']!='exact_original_finite_atoms':
        raise ValueError('No implicit quadrature/importance rule: this implementation handles supplied original finite atoms only')
    entries=[]
    for index,item in enumerate(plan['items']):
        # Finite original atoms already carry their complete original mass.
        # Distinct atoms can have identical visible candidates and are NOT deduplicated.
        if item.get('multiplicity')!=1:
            raise ValueError('This exact atom plan enumerates each original atom once; multiplicity must be 1')
        candidate=item['configuration']
        if item.get('assembler_required'):
            candidate=plugin.canonicalize(candidate);assembled=plugin.A.assemble(candidate)
        else:
            assembled={'candidate_id':candidate.get('native_candidate_id'),
                       'native_inputs_sha256':plugin.DOMAIN['_sha256'],'domain_id':plugin.DOMAIN['domain_id'],
                       'h_N_role':'not_required_for_supplied_whole_history_endpoint_packet'}
        entry=adapter.evaluate_contribution(candidate,assembled,plugin.DOMAIN)
        if entry.get('missing_operations') or any(entry.get('contributions',{}).get(k)is None for k in ('W0','WU','WY','WUY')):
            write_json(args.out/'failed-contribution.json',{'index':index,'candidate':candidate,'result':entry})
            raise ValueError('A required original contribution is not evaluated; no partial integral is accepted')
        entries.append(entry)
    result=adapter.aggregate_contributions(entries,plugin.DOMAIN)
    checks=adapter.run_ensemble_checks(result,entries,plugin.DOMAIN)
    required=registration['required_ids'];by_id={row['id']:row for row in checks}
    complete=all(key in by_id and by_id[key].get('status')=='PASS' for key in required)
    result.update(ensemble_checks=checks,required_ensemble_ids=required,required_ensemble_checks_passed=complete,
                  integration_mode=plan['kind'],evaluated_original_atoms=len(entries),
                  GA_selection_used_for_integration=False,unique_only_mean_used=False,
                  whole_physical_search_completed=False,
                  primitive_source_identity_proved=False)
    write_json(args.out/'contributions.json',entries);write_json(args.out/'result.json',result)
    print(json.dumps({'status':'SUPPLIED_ORIGINAL_ATOMS_NUMERICALLY_AGGREGATED',
                      'input_role':cap.get('input_role'),'result':result,
                      'whole_physical_search_completed':False},indent=2))
    return 0 if all(result.get(k)is not None for k in ('Z0','ZU','ZY','ZUY')) else 2

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--domain',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--mode',choices=['integrate','diagnostic_search','full','diagnostic'],default='integrate')
    p.add_argument('--generations',type=int,default=2);p.add_argument('--population',type=int,default=4);p.add_argument('--batch',type=int,default=2)
    p.add_argument('--seed',type=int,default=731);p.add_argument('--exhaustive-after',action='store_true')
    p.add_argument('--hours',type=float,default=None);p.add_argument('--island',type=int,default=0)
    p.add_argument('--exhaustive-shards',type=int,default=1)
    args=p.parse_args();args.mode={'full':'integrate','diagnostic':'diagnostic_search'}.get(args.mode,args.mode)
    os.environ['TFPT_FIELD_DOMAIN']=str(args.domain.resolve());os.environ['PYTHONDONTWRITEBYTECODE']='1'
    spec=importlib.util.spec_from_file_location('tfpt_existing_search_engine',ENGINE);engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
    plugin=engine.load_problem(Path(__file__).with_name('field_problem.py'))
    preflight={'mode':args.mode,'domain_id':plugin.DOMAIN['domain_id'],'domain_sha256':plugin.DOMAIN['_sha256'],
      'capabilities':plugin.capability_report(),'finite_diagnostic_candidate_count':plugin.ENUMERATION_SIZE,
      'spatial_coordinate_is_OS_time':False,
      'sources':{str(path):hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in plugin.SOURCE_FILES}}
    args.out.mkdir(parents=True,exist_ok=True)
    if args.mode=='integrate':return integrate(plugin,args,preflight)
    write_json(args.out/'preflight.json',preflight)
    args.output=str(args.out.resolve());args.device='cpu';args.migration_every=5;args.source_manifest=str(args.domain.resolve());args.problem=str(Path(__file__).with_name('field_problem.py'))
    result=engine.run(args,problem=plugin)
    print(json.dumps({'status':'DIAGNOSTIC_FIELD_CANDIDATES_EVALUATED','result':result,
                      'feasible_means_kinematic_only':True,'whole_physical_search_completed':False},indent=2));return 0
if __name__=='__main__':sys.exit(main())
