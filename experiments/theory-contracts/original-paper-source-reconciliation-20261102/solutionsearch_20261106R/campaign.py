"""Target-directed batches of whole source programs, with OpenAI Sol only."""
from pathlib import Path
import argparse, datetime as dt, json, os, signal, sys, time
import sympy as sp
from catalogue import inputs, operations, starting_programs
from completion import OriginalCompletion
from construction_batches import make_batch
from search.synthesis import Synthesizer, ast_hash, evidence_json
from source_recovery.native_constructors import NativeInputMissing
from sol_backend import SolBackend, MODEL
from proposal_context import all_trial_history, build_messages, validate_grounding, require_new_program, physical_signal_audit
from run_search import evidence, write, flat_program, decode_program, source_signature


def now(): return dt.datetime.now(dt.timezone.utc).isoformat()
def append(path, value):
    with path.open('a') as stream:
        stream.write(json.dumps(evidence_json(value), ensure_ascii=False)+'\n')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--model',choices=[MODEL],default=MODEL)
    p.add_argument('--rounds',type=int,default=24)
    p.add_argument('--batch-size',type=int,default=64)
    p.add_argument('--proposal-batch',type=int,default=8)
    p.add_argument('--proposal-cost',type=int,default=64)
    p.add_argument('--target',choices=['joint_source','tfpt_complete','tfpt_validated'],default='tfpt_complete')
    p.add_argument('--without-model',action='store_true')
    p.add_argument('--hypothesis-only',action='store_true',help='Explicit finite hypothesis study; does not claim a fundamental discovery route')
    p.add_argument('--prepare-sol-payload',action='store_true',help='Write the exact first proposal messages locally; do not call a model or run the campaign')
    args=p.parse_args(argv)
    if args.rounds<1 or not 1<=args.batch_size<=512: p.error('positive rounds and a whole-program batch1..512 required')
    repo,out=args.repo.resolve(),args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    pack=evidence(repo);write(out/'original-evidence.json',pack)
    history=all_trial_history(out.parent,exclude=out);write(out/'all-prior-trial-feedback.json',history)
    prior_ids={x['id']for x in history['whole_programs']}
    verifier=OriginalCompletion(repo,args.target)
    bindings,registry=inputs(),operations()
    dependencies=[Path(__file__),Path(__file__).with_name('run_search.py'),Path(__file__).with_name('catalogue.py'),
                  Path(__file__).with_name('completion.py'),Path(__file__).with_name('construction_batches.py'),
                  Path(__file__).with_name('sol_backend.py'),Path(__file__).with_name('knowledge_context.py'),
                  Path(__file__).with_name('proposal_context.py'),Path(__file__).with_name('context_codec.py'),
                  Path(__file__).with_name('source_identity.py'),Path(__file__).with_name('source_recovery')/'native_constructors.py',
                  verifier.q_path,verifier.q_path.parent/'original_source_adapter.py',
                  verifier.q_path.parents[1]/'checks/original_checks.py',repo/'semantic_explorer/research/program.json']
    engine=Synthesizer(registry,bindings,verifier,verifier.required_obligation_ids,verifier.target_scope,
                       pack['manifest']['atlas_sha256'],out,verifier_dependencies=dependencies,
                       missing_input_exceptions=(NativeInputMissing,))
    totals={'generated':0,'evaluated':0,'new_source_signatures':0,'duplicate_calculations':0,
            'repeated_programs':0,'structure_families':0}
    families=set();signatures=set();stopping=False
    worker={'pid':os.getpid(),'status':'RUNNING','started_utc':now(),'target':args.target,
            'requested_rounds':args.rounds,'finished_rounds':0,'search_fingerprint':engine.fingerprint,
            'search_unit':'whole executable source construction, with independent calculation and open original proof tests',
            'automatic_full_TFPT_proof_discovery_implemented':False,
            'automatic_whole_source_identity_comparison_implemented':True,
            'capability_limit':'Ausführbare Quellenhypothesen werden konstruiert und ausgerechnet. Der ganze ursprüngliche Quellen/Zeit/Markierungsvergleich ist angeschlossen. Seine ursprüngliche Quellenfamilie muss tatsächlich ausgewertet vorliegen; die vollständigen fundamentalen Beweise erzeugt diese Programmsprache noch nicht selbst.',
            'backend_identity':{'provider':'openai_codex','model':MODEL,'receipt_verified':False},
            'trial_summary':totals,'settings':vars(args),
            'knowledge_context':pack['manifest'],
            'prior_trial_count':len(prior_ids)}
    def persist():
        worker['last_update_utc']=now();worker['trial_summary']=dict(totals);write(out/'worker.json',worker)
    def phase(name,label,description,current=None,total=None):
        worker['phase']={'id':name,'label':label,'description':description,'current':current,'total':total,'started_at':now()};persist()
    def stop(signum,frame):
        nonlocal stopping
        stopping=True
    signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
    persist()

    def evaluate_trial(program,meta):
        totals['generated']+=1;key=ast_hash(program)
        if key in engine.records:
            totals['repeated_programs']+=1;persist()
            return {'ast_hash':key,'status':'REPEATED_PROGRAM','reason':'Dieses vollständige Versuchsprogramm wurde bereits ausgewertet.'}
        begun=now();start=time.monotonic()
        record=engine.evaluate_program(program,max_cost=args.proposal_cost)
        elapsed=time.monotonic()-start
        value=engine.values.get(key);signature=source_signature(value)
        totals['evaluated']+=1
        if signature:
            if signature in signatures:totals['duplicate_calculations']+=1
            else:signatures.add(signature);totals['new_source_signatures']+=1
        family=meta.get('structure_family')
        if family:families.add(family);totals['structure_families']=len(families)
        difference={}
        baseline=engine.values.get(meta.get('parent_id'))
        if isinstance(value,dict) and isinstance(baseline,dict):
            for k in ('Z0','ZU','ZY','ZUY'):
                if k in value and k in baseline:
                    difference[k]={'before':str(baseline[k]),'after':str(value[k]),
                                   'delta':str(sp.simplify(value[k]-baseline[k]))}
        metadata={'id':key,'generated_at':begun,'finished_at':now(),'elapsed_seconds':elapsed,
                  **{k:v for k,v in meta.items() if k not in ('program','id')},
                  'tested_difference':difference,'source_signature':signature,
                  'source_signature_schema':'whole-source-family-v2',
                  'calculation_status':record.get('evaluation_status')}
        append(out/'trial-metadata.jsonl',metadata)
        worker['last_status']=engine.publish_status()
        persist()
        return {'ast_hash':key,'evaluation_status':record.get('evaluation_status'),
                'goal':record.get('goal'),'error':record.get('error'),
                'missing_input':record.get('missing_input'),'elapsed_seconds':elapsed,
                'trial':metadata}

    phase('initialization','Vergleichsprogramme werden berechnet','Vier explizite Startprogramme dienen zum Vergleich der tatsächlich berechneten Antworten.')
    for i,program in enumerate(starting_programs()):
        if stopping:break
        evaluate_trial(program,{'origin':'systematic','change_type':'starting_program',
                               'structure_family':'declared_quartic_start',
                               'what_changes':'Explizite Startvariante der endlichen Quellenkonstruktion.',
                               'what_is_held_fixed':'Die ursprüngliche Kernelalgebra bleibt erhalten; die physische Zuordnung der gewählten Karte ist eine offene Hypothese.'})

    backend=None
    physical_types={'OriginalRuntimeBinding','SourceTimeFamily','OriginalFiniteIntegrand','WholeOriginalSourcePacket'}
    producers=[op['id'] for op in registry if op['output_type']in physical_types or op.get('physical_data_producer')is True]
    producers += [k for k,b in bindings.items()if b.type_name in physical_types and type(b.value).__name__!='Unavailable']
    route=physical_signal_audit(engine.records,verifier.required_obligation_ids,verifier.trusted,producers)
    write(out/'completion-route-audit.json',route);worker['completion_route_audit']=route
    base_messages=build_messages(pack,bindings,registry,history,verifier.required_obligation_ids,
                                 args.proposal_cost,args.proposal_batch)
    write(out/'model-context-coverage.json',{'knowledge_manifest':pack['manifest'],
         'prior_whole_programs':len(prior_ids),'untruncated_inputs':len(bindings),
         'prompt_utf8_bytes':sum(len(m['content'].encode())for m in base_messages),
         'current_run_feedback_policy':'all stored round feedback, no last-N window',
         'all_known_successes_are_scoped':True,'all_known_exclusions_keep_premises':True})
    if args.prepare_sol_payload:
        write(out/'sol-request-preview.json',{'destination':{'provider':'OpenAI via existing Codex session','model':MODEL},
                                             'messages':base_messages,
                                             'followup_scope':'The same selected TFPT material plus actual candidate calculations and counterfeedback for these search rounds.'})
        worker['status']='LOCAL_PREVIEW_FINISHED'
        phase('local_preview','Lokale Übergabevorschau erstellt',
              'Das vollständige Paket wurde lokal gespeichert. Dieser Vorschaumodus sendet keine Daten an ein Modell.')
        return 0
    if not route['completion_route_available'] and not args.hypothesis_only:
        worker['status']='NO_COMPLETION_ROUTE'
        worker['stop_reason']='In dieser Konfiguration fehlen auswertbare physikalische Prüfwerte und ein registrierter Produzent der ursprünglichen Quelldaten. Weitere Varianten würden dieselben fehlenden Eingaben wiederholen.'
        phase('no_completion_route','Der aktuelle Suchraum erreicht das Ziel noch nicht',worker['stop_reason'])
        worker['finished_utc']=now();persist();return 2
    if not args.without_model:
        backend=SolBackend(log_dir=out/'sol_calls',model=args.model,batch_size=args.proposal_batch)
    feedback=[];offset=0
    for iteration in range(args.rounds):
        if stopping:break
        row={'round':iteration+1,'started_utc':now(),'search_fingerprint':engine.fingerprint,'proposals':[]}
        if backend:
            phase('openai_proposal','OpenAI Sol schlägt neue Versuchsprogramme vor',
                  'Ein Aufruf erzeugt mehrere konkrete Programme. Danach werden sie durch die vorhandenen Originalfunktionen ausgeführt.')
            messages=base_messages+([{'role':'user','content':json.dumps({'actual_previous_results':feedback,
                'instruction':'Use actual counterfeedback. Do not repeat an evaluated construction or change only explanatory prose.'},ensure_ascii=False)}]if feedback else[])
            try:
                raw,metrics=backend(messages);batch=json.loads(raw)['proposals']
                worker['backend_identity']={**backend.identity,'receipt_verified':True}
                worker['backend_receipt']={'verified':True,'at':now(),'call_dir':metrics['call_dir'],
                                           'exit_code':metrics['exit_code'],'tool_events':metrics['tool_events']}
                row.update(metrics=metrics,raw=raw)
                for index,proposal in enumerate(batch):
                    if stopping:break
                    phase('candidate_evaluation','Ein Sol Vorschlag wird ausgerechnet',
                          'Der gesamte Quellenkern und seine vier gemeinsamen Antworten werden neu berechnet.',index+1,len(batch))
                    try:
                        validation=validate_grounding(proposal,pack,prior_ids|set(engine.records))
                        program=decode_program(proposal['program'])
                        novelty=require_new_program(program,history['whole_programs'],engine.fingerprint)
                        result=evaluate_trial(program,{'origin':'openai_sol','change_type':'model_construction',
                           'structure_family':'model_program_proposal','what_changes':proposal['changed_premise'],
                           'hypothesis':proposal['hypothesis'],'grounding':proposal['grounding'],
                           'grounding_validation':validation,'novelty_assessment':novelty,'proposed_origin_argument':proposal['proposed_origin_argument'],
                           'missing_derivation':proposal['missing_derivation']})
                        item={'proposal':proposal,'evaluation':result};feedback.append(item)
                    except Exception as exc:item={'proposal':proposal,'proposal_error':type(exc).__name__+': '+str(exc)};feedback.append(item)
                    row['proposals'].append(item)
            except Exception as exc:
                row['provider_error']=type(exc).__name__+': '+str(exc)
                worker['provider_error']=row['provider_error']
                # Never replace the explicitly requested OpenAI model by a local model.
                phase('provider_error','Sol Aufruf konnte nicht abgeschlossen werden',
                      'Der Fehler wird gespeichert. Die direkte Programmsuche bleibt ausführbar; es gibt keinen Ersatz durch Ollama.')
        batch=make_batch(size=args.batch_size,max_cost=args.proposal_cost,offset=offset)
        offset+=len(batch)
        for index,trial in enumerate(batch):
            if stopping:break
            phase('systematic_batch','Ganze Versuchsprogramme werden systematisch geprüft',
                  'Direkt erzeugte vollständige Programme werden ausgeführt, ohne zuerst sämtliche kleinen Zwischenprogramme aufzuzählen.',index+1,len(batch))
            evaluate_trial(trial['program'],{'origin':'systematic',**{k:v for k,v in trial.items()if k not in ('program','id')}})
        row['batch']={'generated':len(batch),'offset_after':offset,'totals':dict(totals)}
        row['finished_utc']=now();append(out/'model-rounds.jsonl',row)
        worker['finished_rounds']=iteration+1;worker['last_status']=engine.publish_status();persist()
        print(json.dumps({'round':iteration+1,'trial_summary':totals,'sol_proposals':len(row['proposals']),
                          'provider_error':row.get('provider_error')},ensure_ascii=False),flush=True)
        if worker['last_status']['scoped_hits']>0:break
    worker['status']='INTERRUPTED'if stopping else'DECLARED_RUN_FINISHED'
    worker['finished_utc']=now();worker['phase']={'id':'finished','label':'Lauf angehalten'if stopping else'Geplanter Suchlauf beendet',
        'description':'Alle tatsächlichen Ergebnisse und offenen ursprünglichen Nachweise bleiben gespeichert.'};persist()
    return 0


if __name__=='__main__':raise SystemExit(main())
