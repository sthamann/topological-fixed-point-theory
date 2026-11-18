"""Actual prior trials and grounded proposal references; never physical proofs."""
from pathlib import Path
import json
from search.synthesis import ast_hash, evidence_json
from run_search import flat_program


def all_trial_history(runs_root, exclude=None):
    records={}; trials={}; origins={}; bad_lines=[]; evaluation_history={}
    for run in sorted(Path(runs_root).iterdir()):
        if not run.is_dir() or (exclude is not None and run.resolve()==Path(exclude).resolve()):continue
        for filename,target,keyname in [('results.jsonl',records,'ast_hash'),('trial-metadata.jsonl',trials,'id')]:
            path=run/filename
            if not path.is_file():continue
            for n,line in enumerate(path.read_text().splitlines(),1):
                if not line.strip():continue
                try: row=json.loads(line)
                except json.JSONDecodeError:
                    bad_lines.append({'run':run.name,'file':filename,'line':n});continue
                key=row.get(keyname)
                if key:
                    target[key]=row;origins.setdefault(key,set()).add(run.name)
                    if filename=='results.jsonl' and row.get('output_type')=='WholeSourceCandidate':
                        ev={'run':run.name,'search_fingerprint':row.get('search_fingerprint'),
                            'evaluation_status':row.get('evaluation_status'),'goal':row.get('goal')}
                        evaluation_history.setdefault(key,[]).append(ev)
    values=[]
    for key,record in records.items():
        if record.get('output_type')!='WholeSourceCandidate':continue
        goal=record.get('goal') or {};proofs=goal.get('obligations') or {}
        if isinstance(proofs,list):proofs={x['id']:x for x in proofs}
        checks={k:v for k,v in proofs.items()if not isinstance(v,dict) or v.get('status')not in ('MISSING_INPUT','NOT_PROVED')}
        meta=trials.get(key,{})
        values.append({'id':key,'program':flat_program(record['ast']),
            'evaluation_status':record.get('evaluation_status'),'search_fingerprint':record.get('search_fingerprint'),'scope':record.get('target_scope'),
            'goal_status':goal.get('status'),'complete_for_target_scope':goal.get('complete_for_target_scope',False),
            'actual_counterexamples_and_passed_scoped_checks':checks,
            'open_obligations':{k:v.get('status')for k,v in proofs.items()if isinstance(v,dict)and v.get('status')in ('MISSING_INPUT','NOT_PROVED')},
            'computed_values':goal.get('evidence',{}).get('computed_values'),
            'source_signature':meta.get('source_signature'),'what_changes':meta.get('what_changes'),
            'proposed_origin_argument':meta.get('proposed_origin_argument'),'observed_in':sorted(origins[key]),
            'all_recorded_evaluation_versions':evaluation_history[key]})
    return {'scope':'Executed whole-source hypothesis programs from prior R runs. Finite scoped checks and source signatures prove no primitive physical identity.',
            'whole_programs':values,'incomplete_jsonl_lines':bad_lines,
            'interpretation':'All stored whole programs are represented without a last-N window. A different explanation for the same program is not a new construction.'}


def grounding_catalog(pack):
    problems=pack['problems']
    if isinstance(problems,dict):problems=problems.values()
    result={'finding_ids':set(),'premise_ids':set(),'strategy_ids':set()}
    premises=pack['premises']
    if isinstance(premises,dict):result['premise_ids'].update(premises)
    else:result['premise_ids'].update(x['id']for x in premises)
    for problem in problems:
        result['finding_ids'].update(x['id']for x in problem.get('findings',[]))
        result['strategy_ids'].update(x['id']for x in problem.get('strategies',[]))
    return result


def validate_grounding(proposal, pack, prior_ids):
    refs=proposal['grounding'];known=grounding_catalog(pack)
    known['prior_program_ids']=set(prior_ids)
    unknown={k:[x for x in refs[k]if x not in known[k]]for k in known}
    unknown={k:v for k,v in unknown.items()if v}
    if unknown:raise ValueError('Unknown knowledge or prior-trial reference: '+json.dumps(unknown))
    if not refs['finding_ids'] or not refs['premise_ids']:
        raise ValueError('A proposal must identify actual known findings AND explicit premises it uses')
    if not refs['why_not_excluded'].strip():raise ValueError('Actual exclusion comparison is required')
    return {'references_exist':True,'semantic_derivation_verified':False,
            'scope':'Reference validity and explicit comparison only. Original mathematical identity is independently checked, never inherited from citations.'}


def build_messages(pack, bindings, registry, history, required_obligations, max_cost=64, batch_size=8):
    import inspect
    from context_codec import model_knowledge_packet, compact_trial_history
    descriptions=[]
    for op in registry:
        row={k:op[k]for k in ('id','input_types','output_type')}
        try:row['actual_operation_definition']=inspect.getsource(op['execute'])
        except (OSError,TypeError):row['definition_state']='source not available; do not infer additional semantics'
        descriptions.append(row)
    values={k:{'type':binding.type_name,'value':evidence_json(binding.value),
               'input_fingerprint':binding.fingerprint}for k,binding in bindings.items()}
    system=(
        'Du entwickelst konkrete neue ausführbare Programme für dieselbe ursprüngliche gemeinsame TFPT Quelle. '
        'Alle wissenschaftlichen Befunde im vollständigen deutschen Atlaspaket behalten ihren Status, ihre '
        'Voraussetzungen und ihren genauen Geltungsbereich. established_scoped bedeutet nur den erklärten Scope; '
        'ein Quellenhash bedeutet Dateifrische, keinen Beweis. Keine bedingte Rechnung als ursprüngliche Physik ausgeben. '
        'Nutze ALLE bereits gelösten Anschlüsse, c3/P1/P2, Ursprung, E8, Clocks, Flavor/Yukawa/Alpha, die gesamte '
        'geschlossene Quellen/Zustands/Zeit/Referenzverkettung und die passenden Lean/Verifikationsbefunde. '
        'Prüfe die bekannten ausgeschlossenen Ansätze und alle gespeicherten Versuche bevor du neue vorschlägst. '
        'Wiederhole keinen ausgeschlossenen Ansatz unter unveränderten Voraussetzungen. Nenne bei einem anderen '
        'Geltungsbereich genau die geänderte Prämisse und warum der Ausschluss dafür nicht gilt. '
        'Deine grounding Felder nennen existierende finding_ids, premise_ids, strategy_ids und einschlägige prior_program_ids. '
        'Diese Referenzen werden geprüft, gelten aber nicht selbst als mathematischer Nachweis. '
        'Nutze ausschließlich die aufgeführten Operationen und Eingaben; deren vollständige Werte und Rollen sind gegeben. '
        'Unavailable bedeutet tatsächlich fehlende Originaleingabe, keine Erlaubnis zu einem Einheitsersatz. '
        'Die endlichen chart_2/chart_4 und dictionary_2/dictionary_4 sind ausdrücklich Hypothesen. Sie ersetzen nicht '
        'das ursprüngliche Feldintegral und liefern allein keinen fundamentalen Treffer. Wenn der fehlende primitive '
        'Produzent aus dieser Grammatik nicht folgt, benenne die genaue fehlende Rechenregel statt einen weiteren '
        'endlichen Parameterfit als Lösung auszugeben. Bewahre das ganze äußere P, Cap, Eta und die absolute Referenz. '
        'Kein inneres P, keine Semigruppenannahme, kein aus KMS erfundenes Profil, keine Wahl eines neuen H_OS. '
        'Für jedes Programm: konkrete zusätzliche Annahme, noch fehlende Herleitung und substanzielles Ursprungsargument. '
        'Antwort nur als SchemaJSON; program ist eine Liste mit Verweisen auf frühere Instruktionen, letzterOutput '
        'WholeSourceCandidate. Keine Befehle, kein Python ausführen, keine Tools. Vollständige Originalabnahme '
        'erfolgt durch unabhängige codegebundene Originalprüfer, niemals durch deine PASS Angaben. '
        'Erklärungen auf Deutsch. Ein bereits ausgeführtes Programm mit anderer Beschreibung ist kein neuer Versuch.')
    request={'operations':descriptions,'inputs':values,'original_context':model_knowledge_packet(pack),
             'all_prior_trial_feedback':compact_trial_history(history),'required_obligations':list(required_obligations),
             'max_ast_cost':max_cost,'desired_proposals':batch_size}
    return [{'role':'system','content':system},
            {'role':'user','content':json.dumps(request,ensure_ascii=False,separators=(',',':'))}]


def require_new_program(program, history_rows, current_fingerprint):
    key=ast_hash(program)
    previous=[r for r in history_rows if r['id']==key]
    if not previous:return {'archived_program':False,'recheck_under_changed_context':False}
    versions=[e for row in previous for e in row.get('all_recorded_evaluation_versions',[])]
    if any(e.get('search_fingerprint')==current_fingerprint for e in versions):
        raise ValueError('REPEATED_ARCHIVED_PROGRAM_IN_SAME_CONTEXT: same construction, inputs and verifier were already evaluated')
    return {'archived_program':True,'recheck_under_changed_context':True,
            'scope':'A changed verifier/input/premise context permits rechecking. This is not a newly discovered source.'}


def physical_signal_audit(records, required_ids, trusted_checker_ids=(), declared_data_producers=()):
    import math
    def bool_test(value):
        if type(value)is bool:return True
        if isinstance(value,dict):return any(bool_test(v)for v in value.values())
        if isinstance(value,(list,tuple)):return any(bool_test(v)for v in value)
        return False
    def numeric_bounds(value):
        if isinstance(value,dict):
            if set(('lower','upper'))<=set(value):value=(value['lower'],value['upper'])
            else:return False
        if not isinstance(value,(list,tuple))or len(value)!=2 or any(type(x)is bool for x in value):return False
        try:
            lo,hi=map(float,value)
            return math.isfinite(lo)and math.isfinite(hi)and lo<=hi
        except (TypeError,ValueError,OverflowError):return False
    rows=[r for r in records.values()if r.get('output_type')=='WholeSourceCandidate'and r.get('goal')]
    vectors=set();executed=0;missing=set()
    for row in rows:
        obs=[o for o in row['goal']['obligations']if o['id']!='constructed_common_values']
        vectors.add(tuple((o['id'],o['status'])for o in obs))
        for o in obs:
            if o['status']not in ('MISSING_INPUT','ERROR')and (o['status']in ('PASS','FAIL')or bool_test(o.get('tests'))or numeric_bounds(o.get('residual_bounds'))):
                executed+=1
            if o['status']=='MISSING_INPUT':missing.add(o['id'])
    missing_complete=[x for x in required_ids if x.startswith('TFPT.TOE.')and x not in trusted_checker_ids]
    producers=list(declared_data_producers)
    return {'whole_programs_observed':len(rows),'physical_status_vector_classes':len(vectors),
            'actual_physical_checks_evaluated':executed,
            'missing_registered_completion_checkers':missing_complete,
            'registered_physical_data_producers':producers,
            'observed_missing_obligations':sorted(missing),
            'completion_route_available':executed>0 or bool(producers),
            'all_requested_completion_checkers_registered':not missing_complete,
            'scope':'Actual registered configuration and observed probe results. Positive route means a physical-data investigation can proceed, not full target acceptance. No general TFPT impossibility or ban on other original-source constructions.'}
