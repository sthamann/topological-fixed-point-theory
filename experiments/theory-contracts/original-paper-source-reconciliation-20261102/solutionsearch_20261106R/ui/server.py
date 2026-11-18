#!/usr/bin/env python3
"""Read-only loopback UI snapshot of actual local source-search run files."""
from pathlib import Path
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from urllib.parse import urlsplit,parse_qs
from datetime import datetime,timezone
import argparse,errno,json,os,re,threading,time

MAX_ITEMS=100
RUN_ID=re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]{0,119}$')
STATIC={'/':'index.html','/index.html':'index.html','/style.css':'style.css','/app.js':'app.js',
        '/fonts/manrope.ttf':'fonts/manrope.ttf','/fonts/dm-sans.ttf':'fonts/dm-sans.ttf'}
MIME={'.html':'text/html; charset=utf-8','.css':'text/css; charset=utf-8','.js':'text/javascript; charset=utf-8','.ttf':'font/ttf'}
CHECK_LABELS={
 'constructed_common_values':'Vier Markierungen gemeinsam berechnet',
 'original_joint_weighted_Z0_ZU_ZY_ZUY':'Abgleich mit der ursprünglichen Quelle',
 'original_OS_history_RP':'Positivität der reflektierten Quelle',
 'actual_quantum_Ward_SD_contacts':'Variationsidentitäten und Kontakte',
 'whole_physical_P_Ref_Eta_domains':'Ganze Wörter, Referenz und Domänen',
 'original_regulator_limit_order':'Kontrollierter Regulatorübergang',
 'primitive_original_source_identity':'Identität mit der ursprünglichen Quelle',
 'original_marked_time_word_match':'Ursprüngliche Zeit und Markierungen',
 'original_independent_holdouts':'Unabhängige Vorhersageprüfung'}
STEP_LABELS={
 'native.c3':'Ursprüngliche Konstante c₃ verwenden','native.c6_ratio':'Ursprüngliches Transferverhältnis verwenden',
 'native.quartic_source_program':'Kern mit höheren Wechselwirkungstermen aufbauen',
 'native.nil_action':'Den vollständigen geraden Kern berechnen','native.even_kernel':'Kern in einen Operator übersetzen',
 'native.sew':'Zwei ganze Kerne geordnet verbinden','native.source_jets':'Vier gemeinsame Markierungen auswerten',
 'native.endpoint_marks':'Tatsächliche Endpunktmarken einsetzen','native.monomial':'Einen geraden Kernbaustein hinzufügen',
 'construct.kernel_term':'Einen gewichteten Kernbaustein bilden','construct.kernel_add':'Kernbausteine gemeinsam verwenden',
 'native.kernel_add':'Kernbausteine gemeinsam verwenden','native.adjoint':'Den adjungierten Kern bilden',
 'native.action_from_polynomial':'Die ganze polynomiale Wirkung einsetzen'}

def utc(value=None):
    return datetime.fromtimestamp(time.time()if value is None else value,timezone.utc).isoformat()
def clipped(value,limit=2000):
    if value is None:return None
    text=value if isinstance(value,str)else json.dumps(value,ensure_ascii=False,default=str)
    return text if len(text)<=limit else text[:limit]+' … [gekürzt]'
def run_label(name):
    if name=='sol-original-source-research':return 'Sol 6.1 · Originalquellenforschung'
    if name=='openai-sol-source-search':return 'OpenAI Sol · Quellensuche'
    if name.startswith('openai-sol-source-search-'):return 'OpenAI Sol · Quellensuche '+name[len('openai-sol-source-search-'):]
    return name.replace('-',' ')
def actual_count(value):
    return value if type(value)is int and value>=0 else None
def family_data(value):
    if isinstance(value,dict):return {clipped(k,140):actual_count(v)for k,v in list(value.items())[:100]}
    if isinstance(value,list):return [clipped(x,140)for x in value[:100]]
    return actual_count(value)
def backend_view(worker):
    identity=worker.get('backend_identity')or{}
    receipt=worker.get('backend_receipt')or{}
    if not isinstance(identity,dict):identity={}
    if not isinstance(receipt,dict):receipt={}
    provider=str(identity.get('provider','')).lower()
    requested=(worker.get('settings')or{}).get('model',worker.get('requested_model')or identity.get('model'))
    declared=str(worker.get('status','')).upper()
    verified=identity.get('receipt_verified')is True or receipt.get('verified')is True
    requires_receipt=('openai'in provider or'codex'in provider or'sol'in str(identity.get('model',requested)).lower()
                      or declared.startswith(('PENDING','WAITING')))
    if requires_receipt and not verified:
        return {'model':None,'model_label':'Modellwechsel noch nicht bestätigt','model_state':'pending',
                'model_receipt_verified':False,'requested_model':clipped(requested,160),'model_identity':{}}
    model=identity.get('model',receipt.get('model'))
    if model is None:
        return {'model':None,'model_label':'Modell noch nicht aufgezeichnet','model_state':'unknown',
                'model_receipt_verified':None,'requested_model':clipped(requested,160),'model_identity':{}}
    return {'model':clipped(model,160),'model_label':clipped(identity.get('display_name',identity.get('label',model)),160),
            'model_state':'confirmed'if verified else'recorded','model_receipt_verified':True if verified else None,'requested_model':clipped(requested,160),
            'model_identity':{key:identity[key]for key in ('provider','model','label','display_name','digest','reasoning_effort','receipt_verified')if key in identity}}
def safe_research(value,depth=0):
    """Bound recorded JSON for display; never follow source or artifact paths."""
    if depth>12:return '[verschachtelte Details gekürzt]'
    if isinstance(value,dict):return {str(k)[:200]:safe_research(v,depth+1)for k,v in list(value.items())[:200]}
    if isinstance(value,list):return [safe_research(v,depth+1)for v in value[:200]]
    if isinstance(value,str):return clipped(value,24000)
    if value is None or isinstance(value,(bool,int,float)):return value
    return clipped(value,2000)

def research_round(row,index):
    proposal=row.get('proposal')if isinstance(row.get('proposal'),dict)else{}
    review=row.get('review')if isinstance(row.get('review'),dict)else{}
    return {'id':clipped(row.get('round_id')or'research-'+str(index),160),'ordinal':index,
            'question':clipped(row.get('question'),12000),'started_at':row.get('started_utc'),
            'completed_at':row.get('completed_utc'),'status':'CRITIC_REVIEWED'if review else'PROPOSED',
            'proposal':safe_research(proposal),'review':safe_research(review),
            'proposal_call':safe_research(row.get('proposal_call')or{}),
            'review_call':safe_research(row.get('review_call')or{}),
            'actual_result':safe_research(row.get('actual_result')or{}),
            'proposed_status':clipped(row.get('proposed_status'),1000),
            'critic_verdict':clipped(row.get('critic_verdict'),3000),
            'next_task':clipped(row.get('next_task'),12000),
            'physics_admission':'NOT_ADMITTED_BY_RUNTIME'}

def trial_metrics(worker,status,rounds):
    fields=('generated','evaluated','new_source_signatures','duplicate_calculations','repeated_programs')
    batches=[r['batch']for r in rounds if isinstance(r.get('batch'),dict)]
    explicit=None;source=None
    for owner,label in ((worker,'worker'),(status,'status')):
        for key in ('trial_summary','batch_totals'):
            if isinstance(owner.get(key),dict):explicit=owner[key];source=label+'.'+key;break
        if explicit is not None:break
    if explicit is None:
        for batch in reversed(batches):
            if isinstance(batch.get('totals'),dict):explicit=batch['totals'];source='round.batch.totals';break
    if explicit is not None:
        totals={key:actual_count(explicit.get(key))for key in fields}
        families=explicit.get('structure_families')
    else:
        # Only sum a field when every recorded per-round batch reports it.
        # Unknown fields are not silently converted to zero or partial totals.
        per_round=[b for b in batches if b.get('scope')not in ('cumulative','total')and b.get('cumulative')is not True]
        totals={key:(sum(b[key]for b in per_round)if per_round and all(actual_count(b.get(key))is not None for b in per_round)else None)for key in fields}
        families=None;source='round.batch'if per_round else'not_recorded'
        if per_round and all(isinstance(b.get('structure_families'),dict)for b in per_round):
            names=set().union(*(b['structure_families']for b in per_round))
            families={key:(sum(b['structure_families'][key]for b in per_round)
                           if all(actual_count(b['structure_families'].get(key))is not None for b in per_round)else None)
                      for key in sorted(names)}
    families=family_data(families)
    latest=worker.get('batch')if isinstance(worker.get('batch'),dict)else(batches[-1]if batches else None)
    if latest is not None:latest={**{k:actual_count(latest.get(k))for k in fields},'structure_families':family_data(latest.get('structure_families'))}
    family_count=len(families)if isinstance(families,(dict,list))else actual_count(families)
    return {**totals,'structure_families':families,'structure_family_count':family_count,'latest_batch':latest,'metrics_source':source,
            'recorded_batches':len(batches),'reported_rounds_with_batches':len(batches)}
def trial_metadata(record,round_row,external):
    keys=('origin','change_type','structure_family','baseline_id','changed_operations','tested_difference',
          'parent_id','what_changes','what_is_held_fixed','generated_at','finished_at','elapsed_seconds',
          'hypothesis','proposed_origin_argument','missing_derivation')
    sources=[external,record.get('trial'),record.get('metadata'),record,
             ((round_row or{}).get('evaluation')or{}).get('trial'),
             (round_row or{}).get('trial'),(round_row or{}).get('candidate_metadata'),round_row]
    result={}
    for key in keys:
        value=next((row[key]for row in sources if isinstance(row,dict)and row.get(key)is not None),None)
        if key=='tested_difference'and isinstance(value,dict):
            value={clipped(k,100):{part:clipped(change.get(part),1500)for part in ('before','after','delta')if part in change}
                   for k,change in list(value.items())[:8]if isinstance(change,dict)}
        elif key=='elapsed_seconds'and isinstance(value,(int,float))and not isinstance(value,bool):pass
        elif isinstance(value,list):value=[clipped(x,500)for x in value[:30]]
        elif value is not None:value=clipped(value,2000)
        result[key]=value
    return result
def process_alive(pid):
    if type(pid)is not int or pid<=0:return None
    try:os.kill(pid,0);return True
    except OSError as exc:
        if exc.errno==errno.ESRCH:return False
        # EPERM means unknown to this server, never evidence that it exited.
        return None

def process_state(worker,launch):
    pid=worker.get('pid',launch.get('pid'));alive=process_alive(pid)
    declared=str(worker.get('status','')).upper()
    terminal=declared in ('DECLARED_RUN_FINISHED','FINISHED','COMPLETED','DONE','LOCAL_PREVIEW_FINISHED')
    interrupted=declared in ('INTERRUPTED','STOPPED','FAILED','ERROR','NO_COMPLETION_ROUTE','NO_SELECTION_SIGNAL')
    if terminal:status='finished';message='Der protokollierte Lauf ist beendet.'
    elif interrupted:status='stopped';message=worker.get('stop_reason')or'Der protokollierte Lauf wurde angehalten.'
    elif alive is False:status='stopped';message='Der Suchprozess läuft nicht mehr. Der letzte gespeicherte Stand bleibt sichtbar.'
    elif alive is True:status='running';message='Der Suchprozess ist erreichbar.'
    else:status='unknown';message='Der Prozessstatus lässt sich gerade nicht bestätigen.'
    return {'pid':pid,'alive':alive,'status':status,'message':message,'declared_status':declared or None}

class FileCache:
    def __init__(self):self.cache={};self.lock=threading.RLock()
    def read(self,path,jsonl=False):
        with self.lock:
            try:st=path.stat()
            except FileNotFoundError:return {'data':[]if jsonl else{},'exists':False,'mtime':None,'warnings':[],'trailing_ignored':False}
            key=(st.st_mtime_ns,st.st_size,st.st_ino)
            old=self.cache.get((str(path),jsonl))
            if old and old['key']==key:return old
            warnings=[];trailing=False
            try:
                raw=path.read_bytes()
                if jsonl:
                    lines=raw.splitlines(keepends=True);data=[]
                    for index,line in enumerate(lines):
                        if not line.endswith((b'\n',b'\r')) and index==len(lines)-1:
                            trailing=True;continue
                        if not line.strip():continue
                        try:
                            item=json.loads(line.decode('utf-8'))
                            if isinstance(item,dict):data.append(item)
                            else:warnings.append('Ein Protokolleintrag enthält kein Objekt.')
                        except (UnicodeError,json.JSONDecodeError):warnings.append('Ein vollständiger Protokolleintrag ist nicht lesbar.')
                else:
                    data=json.loads(raw.decode('utf-8'))
                    if not isinstance(data,dict):data={};warnings.append('Eine Statusdatei enthält kein Objekt.')
            except (UnicodeError,json.JSONDecodeError,OSError):
                data=old['data']if old else([]if jsonl else{})
                warnings.append('Eine Datei wird gerade geschrieben oder ist vorübergehend nicht lesbar. Der letzte lesbare Stand bleibt erhalten.')
            result={'key':key,'data':data,'exists':True,'mtime':st.st_mtime,'warnings':warnings[:8],'trailing_ignored':trailing}
            self.cache[(str(path),jsonl)]=result
            return result

def ast_summary(ast):
    steps=[];remaining=[64]
    def walk(node):
        if not isinstance(node,dict)or remaining[0]<=0:return None
        remaining[0]-=1
        args=[walk(a)for a in node.get('args',[])if isinstance(a,dict)]
        op=clipped(node.get('op'),20);name=clipped(node.get('id'),140)
        if op=='apply':
            label=STEP_LABELS.get(name,'Programmoperation ausführen')
            if label not in steps:steps.append(label)
        return {'op':op,'id':name,'args':[a for a in args if a is not None]}
    normalized=walk(ast)
    return normalized,steps[-12:]

def detail_of(check):
    missing=check.get('missing_inputs')
    if missing:
        if 'actual_original_runtime_binding'in missing:return 'Der tatsächliche Abgleich mit der ursprünglichen gemeinsamen Quelle fehlt.'
        return 'Benötigter Originalbeleg noch offen: '+', '.join(str(x)for x in missing[:4])
    for key in ('counterexample','reason','error','scope'):
        if check.get(key):return clipped(check[key],1000)
    return clipped(check.get('tests'),1000)

def computed_fields(record):
    goal=record.get('goal')or{};evidence=goal.get('evidence')or{}
    values=evidence.get('computed_values')or{}
    for check in goal.get('obligations')or[]:
        if check.get('id')=='constructed_common_values'and check.get('unnormalized_values'):
            values=check['unnormalized_values'];break
    return {key:values.get(key)for key in ('Z0','ZU','ZY','ZUY')}

def live_model_proposal_count(rounds,trials,has_round_file,has_trial_file):
    if not has_round_file and not has_trial_file:return None
    completed=0;represented_ids=set()
    for row in rounds:
        for proposal in proposal_rows(row):
            if not proposal.get('proposal'):continue
            completed+=1
            evaluation=proposal.get('evaluation')or{}
            cid=evaluation.get('ast_hash')or(evaluation.get('trial')or{}).get('id')
            if isinstance(cid,str):represented_ids.add(cid)
    # The model round is persisted only after its entire calculation batch.
    # Trial metadata appears immediately after each actual model candidate.
    live_ids={cid for cid,trial in trials.items()if trial.get('origin')in ('openai','openai_sol')}
    return completed+len(live_ids-represented_ids)

def candidate(record,round_row=None,external_trial=None,ordinal=None):
    goal=record.get('goal')or{};evidence=goal.get('evidence')or{}
    ast,steps=ast_summary(record.get('ast'))
    checks=[]
    for c in (goal.get('obligations')or[])[:40]:
        if not isinstance(c,dict):continue
        cid=str(c.get('id',''))
        label=CHECK_LABELS.get(cid,'Gesamtanforderung '+cid.rsplit('.',1)[-1]if cid.startswith('TFPT.TOE.')else cid)
        checks.append({'id':cid,'status':c.get('status','NOT_PROVED'),'label':label,'detail':detail_of(c)})
    unresolved=[c for c in checks if c['status']!='PASS'];failed=[c for c in checks if c['status']=='FAIL']
    proposal=(round_row or{}).get('proposal')or{}
    values=computed_fields(record)
    normalized={key:clipped(values.get(key),2048)for key in ('Z0','ZU','ZY','ZUY')}
    normalized['connected']=clipped(evidence.get('connected_UY'),2048)
    state=goal.get('status')or record.get('evaluation_status','NOT_PROVED')
    labels={'PARTIAL':'Berechnet, Originalbelege offen','SCOPED_HIT':'Prüfungen im angegebenen Zielumfang bestanden',
            'CANDIDATE_FAIL':'Eine Zielprüfung ist gescheitert','MISSING_INPUT':'Ursprünglicher Eingang fehlt',
            'EXECUTION_ERROR':'Programm konnte nicht ausgeführt werden','TYPE_VALUE_ERROR':'Programmausgabe passt nicht zum verlangten Typ'}
    trial=trial_metadata(record,round_row,external_trial)
    origin={'systematic':'systematic','openai':'openai','openai_sol':'openai'}.get(trial['origin'])
    return {'id':record.get('ast_hash'),'ordinal':ordinal,'ast':ast,'human_steps':steps,'cost':record.get('cost'),
            'trial':trial,'origin':origin,'origin_label':{'systematic':'Systematische Suche','openai':'OpenAI Sol'}.get(origin),
            **{key:trial[key]for key in ('change_type','structure_family','baseline_id','changed_operations','tested_difference','parent_id','what_changes','what_is_held_fixed')},
            'status':state,'status_label':labels.get(state,state),'evaluation_status':record.get('evaluation_status'),
            'computed':normalized,'has_computed_values':all(normalized[k]is not None for k in ('Z0','ZU','ZY','ZUY')),
            'checks':checks,'checks_passed':sum(c['status']=='PASS'for c in checks),'checks_total':len(checks),
            'scope':goal.get('target_scope',record.get('target_scope')),'scope_label':'Die Zahlen gelten für dieses Versuchsprogramm und seine ausdrücklich angegebenen Voraussetzungen.',
            'first_failure':(failed or unresolved or[None])[0],
            'hypotheses':[clipped(x,1000)for x in evidence.get('hypotheses',[])[:20]],
            'hypothesis':clipped(proposal.get('hypothesis',trial.get('hypothesis')),2000),
            'proposed_origin_argument':clipped(proposal.get('proposed_origin_argument',trial.get('proposed_origin_argument')),3000),
            'missing_derivation':clipped(proposal.get('missing_derivation',trial.get('missing_derivation')),2000),
            'changed_premise':clipped(proposal.get('changed_premise'),1500),
            'round':(round_row or{}).get('round'),'time':(round_row or{}).get('finished_utc',trial.get('generated_at'))}

def round_event(row,index):
    error=str(row.get('proposal_error')or row.get('backend_error')or row.get('provider_error')or'')
    evaluation=row.get('evaluation')or{};proposal=row.get('proposal')or{}
    if evaluation.get('status')=='REPEATED_PROGRAM'or(error and ('repeat'in error.lower()or'same four whole source forms'in error.lower()or'already evaluated'in error.lower())):
        kind='repeat';title='Bereits geprüftes Versuchsprogramm wiederholt';description='Der Vorschlag bringt keinen neuen geprüften Quellenkern.'
    elif error:
        kind='invalid';title='Vorschlag konnte nicht verwendet werden'
        description='Das Programm ist noch nicht ausführbar.'
        if 'earlier instructions'in error:description='Ein Schritt verweist auf einen noch nicht vorhandenen Programmschritt.'
    elif evaluation.get('evaluation_status')=='VALUE_BUILT':
        goal=evaluation.get('goal')or{}
        if goal.get('complete_for_target_scope')is True:kind='scoped_hit';title='Prüfer bestätigt einen Treffer';description='Alle verlangten Prüfungen des angegebenen Zielumfangs sind bestanden.'
        elif goal.get('status')=='CANDIDATE_FAIL':kind='candidatefail';title='Versuchsprogramm scheitert an einer Prüfung';description='Das Versuchsprogramm wurde berechnet und hat eine verlangte Prüfung nicht bestanden.'
        else:kind='kernelcomputed';title='Versuchsprogramm berechnet';description='Das Versuchsprogramm wurde ausgeführt. Seine ursprüngliche physikalische Zuordnung ist weiter zu prüfen.'
    elif evaluation.get('evaluation_status')=='MISSING_INPUT':
        kind='missinginput';title='Originaleingang für das Versuchsprogramm fehlt';description='Dieser Vorschlag bleibt offen. Andere Versuchsprogramme können weiter geprüft werden.'
        error=str(evaluation.get('missing_input')or'')
    elif evaluation.get('evaluation_status')in ('INVALID_AST','EXECUTION_ERROR','TYPE_VALUE_ERROR'):
        kind='invalid';title='Versuchsprogramm konnte nicht ausgeführt werden';description='Der tatsächliche Auswertungsfehler wurde gespeichert.'
        error=str(evaluation.get('error')or'')
    else:kind='proposed';title='Neuer Programmvorschlag';description='Der protokollierte Vorschlag wird mit den vorhandenen Operationen geprüft.'
    return {'id':'round-'+str(index),'time':row.get('finished_utc',row.get('started_utc')),'kind':kind,
            'title':title,'description':description,'candidate_id':evaluation.get('ast_hash'),'round':row.get('round'),
            'error':clipped(error,1500)or None,'hypothesis':clipped(proposal.get('hypothesis'),2000),
            'proposed_origin_argument':clipped(proposal.get('proposed_origin_argument'),3000),
            'model':(row.get('metrics')or{}).get('model'),'backend_change':clipped((row.get('backend_change')or{}).get('identity',{}).get('model'),100)}

def proposal_rows(row):
    """Normalize legacy single proposals and actual new per-round batches."""
    if isinstance(row.get('proposals'),list):
        result=[]
        for index,item in enumerate(row['proposals']):
            if not isinstance(item,dict):continue
            evaluation=item.get('evaluation')if isinstance(item.get('evaluation'),dict)else{}
            trial=evaluation.get('trial')if isinstance(evaluation.get('trial'),dict)else{}
            result.append({'round':row.get('round'),'proposal_index':index+1,'started_utc':row.get('started_utc'),
                'finished_utc':trial.get('finished_at',trial.get('generated_at',row.get('finished_utc'))),
                'metrics':row.get('metrics')or{},'proposal':item.get('proposal')or{},'evaluation':evaluation,
                'proposal_error':item.get('proposal_error')})
        return result
    return [row]

def activity_of(worker):
    phase=worker.get('phase')
    if not isinstance(phase,dict):return None
    return {'id':clipped(phase.get('id'),140),'label':clipped(phase.get('label'),300),
            'description':clipped(phase.get('description'),1500),'current':actual_count(phase.get('current')),
            'total':actual_count(phase.get('total')),'started_at':clipped(phase.get('started_at'),100)}

def candidate_selection(candidates):
    def important(c):
        return (c.get('status')in ('CANDIDATE_FAIL','SCOPED_HIT','EXECUTION_ERROR','TYPE_VALUE_ERROR','ERROR')
                or c.get('evaluation_status')in ('EXECUTION_ERROR','TYPE_VALUE_ERROR','ERROR')
                or any(check.get('status')in ('FAIL','ERROR')for check in c.get('checks',[])))
    cutoff=max(0,len(candidates)-MAX_ITEMS)
    # Preserve the original full-dataset order and ordinal. A historical
    # counterexample must remain reachable by the client's existing filters.
    selected=[c for index,c in enumerate(candidates)if index>=cutoff or important(c)]
    return selected,{'policy':'latest_plus_all_important_findings','recent_limit':MAX_ITEMS,
                     'important_retained':sum(important(c)for c in selected),
                     'historical_important':sum(important(c)for c in candidates[:cutoff]),
                     'returned':len(selected),'total':len(candidates)}

class Snapshots:
    def __init__(self,runs_root,default_run):
        self.root=Path(runs_root).resolve();self.default=default_run;self.files=FileCache()
    def path(self,name):
        if not RUN_ID.fullmatch(name)or name in ('.','..'):raise ValueError('Ungültiger Laufname')
        path=(self.root/name).resolve()
        if path.parent!=self.root or not path.is_dir():raise FileNotFoundError('Lauf nicht vorhanden')
        return path
    def read_file(self,run,name):
        path=run/name
        if path.resolve().parent!=run.resolve():
            return {'data':[]if name.endswith('.jsonl')else{},'exists':False,'mtime':None,
                    'warnings':['Eine verknüpfte Datei außerhalb dieses Laufs wird nicht gelesen.'],'trailing_ignored':False}
        return self.files.read(path,jsonl=name.endswith('.jsonl'))
    def runs(self):
        rows=[]
        if self.root.is_dir():
            for path in self.root.iterdir():
                if not RUN_ID.fullmatch(path.name)or not path.is_dir()or path.resolve().parent!=self.root:continue
                if not any((path/n).is_file()for n in ('worker.json','status.json','results.jsonl','launch.json')):continue
                w=self.read_file(path,'worker.json');launch=self.read_file(path,'launch.json');status=process_state(w['data'],launch['data'])
                mt=max([x for x in (w['mtime'],launch['mtime'])if x is not None],default=None)
                rows.append({'id':path.name,'label':run_label(path.name),'is_default':path.name==self.default,
                             'kind':w['data'].get('kind'),'historical':path.name!=self.default,
                             'status':status['status'],'updated_at':utc(mt)if mt is not None else None})
        rows.sort(key=lambda r:(not r['is_default'],r['updated_at']or''),reverse=False)
        return {'default_run':self.default,'runs':rows}
    def state(self,name=None):
        name=name or self.default;path=self.path(name)
        reads={n:self.read_file(path,n)for n in ('worker.json','status.json','launch.json','results.jsonl','model-rounds.jsonl','trial-metadata.jsonl')}
        worker=reads['worker.json']['data']
        if worker.get('kind')=='original_source_research':return self.research_state(name,path,worker,reads)
        status=reads['status.json']['data']or worker.get('last_status')or{}
        launch=reads['launch.json']['data'];proc=process_state(worker,launch)
        fingerprint=status.get('search_fingerprint',worker.get('search_fingerprint'))
        records=reads['results.jsonl']['data'];records=[r for r in records if not fingerprint or r.get('search_fingerprint')==fingerprint]
        unique={r.get('ast_hash',json.dumps(r.get('ast'),sort_keys=True)):r for r in records}
        whole=[r for r in unique.values()if r.get('output_type')=='WholeSourceCandidate']
        rounds=reads['model-rounds.jsonl']['data'];round_by_candidate={}
        for r in rounds:
            for proposal in proposal_rows(r):
                evaluation=proposal.get('evaluation')or{};cid=evaluation.get('ast_hash')
                if cid and evaluation.get('status')!='REPEATED_PROGRAM':round_by_candidate[cid]=proposal
        trials={}
        for row in reads['trial-metadata.jsonl']['data']:
            cid=row.get('id',row.get('ast_hash'))
            if isinstance(cid,str):trials[cid]=row
        model_proposals=live_model_proposal_count(rounds,trials,reads['model-rounds.jsonl']['exists'],reads['trial-metadata.jsonl']['exists'])
        all_candidates=[candidate(r,round_by_candidate.get(r.get('ast_hash')),trials.get(r.get('ast_hash')),ordinal=index)
                        for index,r in enumerate(whole,start=1)]
        candidates,selection=candidate_selection(all_candidates)
        events=[]
        if worker.get('started_utc'):events.append({'id':'run-start','time':worker['started_utc'],'kind':'started','title':'Suche gestartet','description':'Programme werden konstruiert, berechnet und an den ursprünglichen Anforderungen geprüft.'})
        for i,row in enumerate(rounds):
            if isinstance(row.get('proposals'),list)and(row.get('metrics')or row.get('provider_error')or row['proposals']):
                error=clipped(row.get('provider_error'),1500)
                events.append({'id':'round-'+str(i)+'-model-call','time':row.get('started_utc',row.get('finished_utc')),
                    'kind':'model_error'if error else'model_call','title':'Sol Aufruf fehlgeschlagen'if error else'Sol Modellrunde gestartet',
                    'description':'Der fehlgeschlagene Aufruf wurde protokolliert.'if error else str(len(row['proposals']))+' einzelne Vorschläge sind für diese Runde protokolliert.',
                    'round':row.get('round'),'error':error,'time_source':'round_start'})
            for j,proposal in enumerate(proposal_rows(row)):
                event=round_event(proposal,str(i)+'-proposal-'+str(j));events.append(event)
                if event['kind']=='kernelcomputed':
                    events.append({'id':event['id']+'-proofs-open','time':event['time'],'kind':'originalproofsopen',
                        'title':'Ursprüngliche Belege bleiben offen','description':'Berechnete Markierungen und die Identität mit der ursprünglichen Quelle werden getrennt geprüft.',
                        'candidate_id':event['candidate_id'],'round':event['round']})
            if isinstance(row.get('batch'),dict)and actual_count(row['batch'].get('generated'))is not None:
                events.append({'id':'round-'+str(i)+'-systematic-batch','time':row.get('finished_utc'),
                    'kind':'batch','title':'Systematischer Batch protokolliert',
                    'description':str(row['batch']['generated'])+' ganze Versuchsprogramme wurden für diesen Batch erstellt. Die tatsächlich ausgewertete Menge wird getrennt gezählt.',
                    'round':row.get('round')})
            changed=(row.get('backend_change')or{}).get('identity')or{}
            if changed.get('model'):
                events.append({'id':event['id']+'-model-change','time':event['time'],'kind':'backend_change',
                    'title':'Lokales Modell gewechselt','description':'Das Protokoll meldet den Wechsel zu '+clipped(changed['model'],100)+'.',
                    'round':event['round']})
        if worker.get('finished_utc'):events.append({'id':'run-finish','time':worker['finished_utc'],'kind':'finished','title':'Protokollierter Lauf beendet','description':proc['message']})
        events.sort(key=lambda r:(r.get('time')or'',r['id']));events=events[-MAX_ITEMS:]
        counts={'unique_asts':len(unique)if reads['results.jsonl']['exists']else None,
                'program_fragments':sum(r.get('output_type')!='WholeSourceCandidate'for r in unique.values())if reads['results.jsonl']['exists']else None,
                'whole_candidates':len(whole)if reads['results.jsonl']['exists']else status.get('goal_candidates_in_current_fingerprint_archive'),
                'model_rounds':len(rounds)if reads['model-rounds.jsonl']['exists']else None,
                'model_proposals':model_proposals,
                'scoped_hits':sum((r.get('goal')or{}).get('complete_for_target_scope')is True for r in whole),
                'partial_candidates':sum((r.get('goal')or{}).get('status')=='PARTIAL'for r in whole),
                'failed_candidates':sum((r.get('goal')or{}).get('status')=='CANDIDATE_FAIL'for r in whole),
                'computed_candidates':sum(all(value is not None for value in computed_fields(r).values())for r in whole)}
        metrics=trial_metrics(worker,status,rounds)
        metrics['model_proposals']=model_proposals
        counts.update(generated_programs=metrics['generated'],evaluated_programs=metrics['evaluated'],
                      new_source_signatures=metrics['new_source_signatures'],duplicate_calculations=metrics['duplicate_calculations'],
                      repeated_programs=metrics['repeated_programs'],
                      structure_families=metrics['structure_families'],structure_family_count=metrics['structure_family_count'])
        open_proofs=any(c['status']=='PARTIAL'or any(k['status']!='PASS'for k in c['checks'])for c in candidates)if candidates else None
        if counts['scoped_hits']:
            headline='Ein Treffer ist im angegebenen Zielumfang bestätigt.';boundary='Die tatsächlichen Prüfungen des angegebenen Zielumfangs sind bestanden.'
        elif counts['whole_candidates']:
            headline='Versuchsprogramme berechnet, ursprüngliche Zuordnung offen.';boundary='Berechnete Markierungen bestätigen noch nicht die ursprüngliche gemeinsame Feldquelle.'
        else:
            headline='Die Suche prüft Programmbausteine und Vorschläge.';boundary='Ein bestätigter ganzer Quellenkandidat liegt in diesem Stand noch nicht vor.'
        updated=max([r['mtime']for r in reads.values()if r['mtime']is not None],default=None)
        warnings=[]
        for r in reads.values():warnings.extend(r['warnings'])
        trailing=any(r['trailing_ignored']for r in reads.values())
        backend=backend_view(worker)
        settings=worker.get('settings')or{}
        return {'run':{'id':name,'label':run_label(name),'kind':worker.get('kind','finite_source_search'),'historical':name!=self.default,'started_at':worker.get('started_utc'),
                       'target':worker.get('target'),'scope':status.get('target_scope'),
                       'batch_size':actual_count(settings.get('batch_size')),
                       'proposal_batch':actual_count(settings.get('proposal_batch')),
                       'knowledge_context':worker.get('knowledge_context')if isinstance(worker.get('knowledge_context'),dict)else None,
                       'prior_trial_count':actual_count(worker.get('prior_trial_count')),
                       'whole_source_identity_check':worker.get('automatic_whole_source_identity_comparison_implemented')is True,
                       'completion_route_audit':worker.get('completion_route_audit'),**backend},
                'process':proc,'activity':activity_of(worker),'summary':{'headline':headline,'explanation':proc['message'],'counts':counts,
                   'requested_rounds':worker.get('requested_rounds'),'finished_rounds':worker.get('finished_rounds')},
                'trial_summary':metrics,
                'candidates':candidates,'candidate_selection':selection,'events':events,'freshness':{'snapshot_at':utc(),'last_file_update':utc(updated)if updated else None,
                   'worker_alive':proc['alive'],'trailing_jsonl_ignored':trailing},
                'coverage':{'label':'Begrenzte Programmsprache','detail':'Die angegebene Suchklasse und Kostenbegrenzung sind keine vollständige Abdeckung aller physikalischen Quellen.',
                   'status':status.get('coverage'),'last_complete_cost':status.get('last_complete_cost'),'max_cost':status.get('max_cost')},
                'boundary':{'headline':boundary,'detail':clipped(worker.get('capability_limit'),3000),
                   'original_proofs_open':open_proofs,'target_scope':status.get('target_scope')},
                'warnings':list(dict.fromkeys(warnings))[:8]}


    def research_state(self,name,path,worker,reads):
        reads['research-rounds.jsonl']=self.read_file(path,'research-rounds.jsonl')
        rows=reads['research-rounds.jsonl']['data']
        rounds=[research_round(row,index)for index,row in enumerate(rows,start=1)]
        proc=process_state(worker,reads['launch.json']['data'])
        metrics=worker.get('trial_summary')if isinstance(worker.get('trial_summary'),dict)else{}
        counts={key:actual_count(metrics.get(key))for key in ('rounds','model_calls','tool_events','completed_commands','failed_commands')}
        counts['recorded_rounds']=len(rounds)
        counts['reviewed_rounds']=sum(row['status']=='CRITIC_REVIEWED'for row in rounds)
        events=[]
        if worker.get('started_utc'):events.append({'id':'research-start','time':worker['started_utc'],'kind':'started',
            'title':'Originalquellenforschung gestartet','description':'Sol liest Originalquellen, prüft die tragende Annahme und lässt den Ansatz unabhängig kritisieren.'})
        for row in rounds:
            events.append({'id':row['id']+'-proposal','time':row['started_at'],'kind':'research_proposal',
                'title':'Forschungsansatz protokolliert','description':row['question'],'candidate_id':row['id']})
            if row['status']=='CRITIC_REVIEWED':events.append({'id':row['id']+'-review','time':row['completed_at'],'kind':'research_review',
                'title':'Unabhängige Kritik gespeichert','description':row['critic_verdict']or'Bewertung und offene Beweisschritte sind in der Runde dokumentiert.','candidate_id':row['id']})
        active=worker.get('active_call')if isinstance(worker.get('active_call'),dict)else{}
        if active:events.append({'id':'active-call','time':worker.get('last_update_utc'),'kind':'research_activity',
            'title':'Unabhängige Prüfung läuft'if active.get('role')=='critic'else'Originalquellen werden bearbeitet',
            'description':clipped(active.get('error')or active.get('latest_event_type')or'Ein echter Modellaufruf ist aktiv.',1500)})
        if proc['status']in ('stopped','finished'):events.append({'id':'research-stop','time':worker.get('finished_utc',worker.get('last_update_utc')),
            'kind':'error'if proc['declared_status']in ('ERROR','FAILED')else'finished','title':'Forschung '+('beendet'if proc['status']=='finished'else'angehalten'),'description':proc['message']})
        events.sort(key=lambda row:(row.get('time')or'',row['id']))
        updated=max([r['mtime']for r in reads.values()if r['mtime']is not None],default=None)
        backend=backend_view(worker)
        return {'run':{'id':name,'label':run_label(name),'kind':'original_source_research','historical':name!=self.default,
                       'started_at':worker.get('started_utc'),'knowledge_context':safe_research(worker.get('knowledge_context')),
                       'backend_receipt':safe_research(worker.get('backend_receipt')),**backend},
            'process':proc,'activity':activity_of(worker),'active_call':safe_research(active),
            'research':{'question':clipped(worker.get('active_question'),12000),'rounds':rounds,
                        'physics_admission':'NOT_ADMITTED_BY_RUNTIME','latest_next_task':rounds[-1]['next_task']if rounds else None},
            'summary':{'headline':'Die ursprüngliche Integrationsvorschrift wird untersucht.','explanation':proc['message'],'counts':counts,
                       'requested_rounds':worker.get('requested_rounds'),'finished_rounds':worker.get('finished_rounds')},
            'trial_summary':counts,'candidates':[], 'events':events[-MAX_ITEMS:],
            'freshness':{'snapshot_at':utc(),'last_file_update':utc(updated)if updated else None,'worker_alive':proc['alive'],
                         'trailing_jsonl_ignored':any(r['trailing_ignored']for r in reads.values())},
            'coverage':{'detail':'Alle gespeicherten Forschungsrunden bleiben auswählbar. Jede Aussage gilt nur unter ihren dokumentierten Voraussetzungen.'},
            'warnings':list(dict.fromkeys(w for r in reads.values()for w in r['warnings']))[:8]}

class ReadOnlyServer(ThreadingHTTPServer):
    daemon_threads=True
    def __init__(self,address,runs_root,default_run,static_root=None):
        self.snapshots=Snapshots(runs_root,default_run);self.static_root=Path(static_root or Path(__file__).parent).resolve()
        super().__init__(address,Handler)

class Handler(BaseHTTPRequestHandler):
    def _respond(self,status,body,content_type='application/json; charset=utf-8',head=False):
        self.send_response(status);self.send_header('Content-Type',content_type);self.send_header('Content-Length',str(len(body)))
        self.send_header('Cache-Control','no-store');self.send_header('X-Content-Type-Options','nosniff')
        self.send_header('Referrer-Policy','no-referrer')
        self.send_header('Content-Security-Policy',"default-src 'self'; connect-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; font-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'")
        self.end_headers()
        if not head:self.wfile.write(body)
    def _json(self,status,value,head=False):self._respond(status,json.dumps(value,ensure_ascii=False,allow_nan=False).encode(),head=head)
    def _route(self,head=False):
        allowed={'127.0.0.1:'+str(self.server.server_port),'localhost:'+str(self.server.server_port)}
        if self.headers.get('Host')not in allowed:return self._json(403,{'error':'Nur der lokale UI Aufruf ist erlaubt.'},head)
        request=urlsplit(self.path)
        try:
            if request.path=='/api/runs':return self._json(200,self.server.snapshots.runs(),head)
            if request.path=='/api/state':
                name=parse_qs(request.query,keep_blank_values=True).get('run',[None])[0]
                return self._json(200,self.server.snapshots.state(name),head)
            if request.path not in STATIC:return self._json(404,{'error':'Nicht vorhanden.'},head)
            path=(self.server.static_root/STATIC[request.path]).resolve()
            if not path.is_relative_to(self.server.static_root)or not path.is_file():return self._json(404,{'error':'Nicht vorhanden.'},head)
            return self._respond(200,path.read_bytes(),MIME[path.suffix],head)
        except ValueError:return self._json(400,{'error':'Ungültiger Laufname.'},head)
        except FileNotFoundError:return self._json(404,{'error':'Lauf nicht vorhanden.'},head)
        except (OSError,TypeError,KeyError,json.JSONDecodeError):return self._json(503,{'error':'Der gespeicherte Stand ist vorübergehend nicht lesbar.'},head)
    def do_GET(self):self._route()
    def do_HEAD(self):self._route(head=True)
    def do_POST(self):self._json(405,{'error':'Diese Ansicht liest ausschließlich Daten.'})
    do_PUT=do_POST;do_DELETE=do_POST;do_PATCH=do_POST

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--runs-root',type=Path,required=True)
    parser.add_argument('--port',type=int,default=8768);parser.add_argument('--default-run',default='constructive-local-search-v3')
    args=parser.parse_args()
    server=ReadOnlyServer(('127.0.0.1',args.port),args.runs_root,args.default_run)
    print('Readonly live view on http://127.0.0.1:'+str(server.server_port),flush=True)
    try:server.serve_forever(poll_interval=0.5)
    except KeyboardInterrupt:pass
    finally:server.server_close()
if __name__=='__main__':main()
