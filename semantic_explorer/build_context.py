#!/usr/bin/env python3
"""Authoritative ledger assessment and prior-attempt overlay, no LLM calls.

This separates reported scoped evidence from open physical goals. It never
upgrades a ledger claim because a script is present or a contract says PASS.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
from datetime import datetime, timezone
try:
    from . import source_inventory
    from .json_transport import json_transport_value
except ImportError:  # direct `python semantic_explorer/build_context.py`
    import source_inventory
    from json_transport import json_transport_value

HERE=Path(__file__).resolve().parent
REPO=HERE.parent
LABELS={
    'recorded_established':'Belegt im angegebenen Rahmen',
    'conditional':'Bedingter Anschluss',
    'open':'Offene Frage',
    'mixed':'Belegte Teilresultate · getrennte Geltungsbereiche',
    'excluded':'Eingegrenzter Ausschluss / Gegenbefund',
    'assumed':'Gesetzte Voraussetzung',
    'historical':'Historischer, abgelöster Eintrag',
    'methodological_rule':'Geltende Prüfregel',
    'unknown':'Status nicht automatisch zugeordnet',
}

def assessment_fingerprint(row):
    """Bind a reviewed display correction to the exact authoritative row."""
    payload={k:row.get(k,'') for k in ('claim_id','claim','status','canonical_status','active')}
    return hashlib.sha256(json.dumps(payload,sort_keys=True,ensure_ascii=False).encode()).hexdigest()


def assessment_reviews():
    path=HERE/'classification_reviews.json'
    if not path.exists(): return {}
    document=json.loads(path.read_text(encoding='utf-8'))
    result={}
    for item in document.get('reviews',[]):
        key=item['claim_id']
        if key in result or item['category'] not in LABELS or not re.fullmatch(r'[0-9a-f]{64}',item['row_sha256']):
            raise ValueError(f'invalid classification review: {key}')
        result[key]=item
    return result


CLASSIFICATION_REVIEWS=assessment_reviews()


def assess(row, reviews=None):
    status=row['status']; canon=row.get('canonical_status',''); text=status.lower()
    if row.get('active')=='false':category='historical'
    elif re.match(r'^firewall\b',text) and ('standing discipline rule' in text or 'in force' in text):
        category='methodological_rule'
    else:
        markers=re.findall(r'\[([^\]]+)\]',text)
        marker_codes={code for marker in markers for code in re.findall(r'(?<![a-z0-9])([eilfnocpa])(?=$|[/|,+\s()\-])',marker)}
        positive_text=re.sub(r'\b(?:not|never|no)\s+(?:a\s+)?(?:proven|proved|verified|exact|theorem)\b','',text)
        established=bool(marker_codes & set('eilfn')) or bool(re.search(r'\b(identity|formal|lattice|numerical|proven|proved|exact|measured|verified)\b|kernel.checked|theorem.grade|abstract master theorem',positive_text))
        # An open/conditional clause may coexist with a narrow established result.
        opened='o' in marker_codes or bool(re.search(r'\[o\]|\bopen\b|remain(?:s)? open|not (?:derived|proven|proved|verified|established|forced)|un(?:constructed|proven|proved|verified)|unresolved',text))
        conditional=bool(marker_codes & set('cp')) or bool(re.search(r'\bconditional\b|\bphysical\b|\bconjecture\b|\bhypothes[ei]s\b|\bhypotheses\b',text))
        assumed=bool(re.search(r'\baxiom\b|\[a\]|declared input',text))
        excluded=bool(re.search(r'\bno.go\b|\brefuted\b|\bcounterexample\b|\bkilled\b|\bfalsified\b|\bdead verdict\b',text)) or any(re.search(r'\bkill\b',m) for m in markers)
        recorded_negative=bool(re.search(r'\bfails?\b|\bnegatives?\b|\bsplit verdict',text))
        proof_near=bool(re.search(r'\bproof[- ]near\b|\bproven[- ]modulo\b',text+' '+canon.lower()))
        if excluded:category='mixed' if re.search(r'\bside finding\b|\bside result\b',text+' '+canon.lower()) else 'excluded'
        elif established and (opened or conditional or assumed or recorded_negative or proof_near):category='mixed'
        elif opened:category='open'
        elif conditional:category='conditional'
        elif assumed:category='assumed'
        elif excluded:category='excluded'
        elif established:category='recorded_established'
        else:category='unknown'
        # A leading canonical type describes this result. A later clause such as
        # "T1-T8 remain [O]" describes its parent and must not erase a narrow theorem.
        leading=re.match(r'^\s*\[([^\]]+)\]',canon.lower())
        if leading and category not in ('excluded','assumed'):
            codes=set(re.findall(r'(?<![a-z0-9])([eilfnocpa])(?=$|[/|,+\s()\-])',leading.group(1)))
            if codes & set('cp'):
                category='mixed' if codes & set('eilfn') else 'conditional'
            elif 'o' in codes:
                category='mixed' if established else 'open'
    review=(CLASSIFICATION_REVIEWS if reviews is None else reviews).get(row.get('claim_id'))
    current_review=bool(review and review.get('row_sha256')==assessment_fingerprint(row))
    if current_review and row.get('active')!='false': category=review['category']
    result={'category':category,'label':LABELS[category],'scope':canon or status,'status_raw':status,
        'classification_uncertainty':category=='unknown',
        'authority':'ledger','classification_method':'status_text_rules','notice':'Berichteter Ledgerstatus im genannten Geltungsbereich; keine neue fachliche Abnahme oder Ausführung.',
        'notice_en':'Ledger-reported status within the stated scope; no new expert review or execution.'}
    if review:
        result['classification_review']={'current':current_review,'row_sha256':review['row_sha256'],
            'reason_de':review.get('reason_de',''),'reason_en':review.get('reason_en','')}
        if current_review: result['classification_method']='source_reviewed_row_mapping'
        else: result['classification_uncertainty']=True
    return result

def merge_original_contract(indexed, original):
    """Current explicit source fields override lossy/stale indexed summaries.

    Preserve absent fields and indexed provenance; empty original lists are
    intentional corrections and must override old values as well.
    """
    result=dict(indexed)
    result['original_contract']=original
    for key in ['question','date','claims_tested','scope','claim_scope','summary_de','summary_en',
                'leaves_open','missing','not_derived','closes','resolved',
                'positive_results','analytic_conditional_results','not_claimed',
                'assumptions','physical_gates_closed','complete_TFPT_solution']:
        if key in original:result[key]=original[key]
    if 'verdict' in original:result['verdict_raw']=original['verdict']
    if 'kills' in original:result['kills_free']=original['kills']
    return result

def sha(data):return hashlib.sha256(data).hexdigest()

# Navigation focus, not evidence classification. All broad links remain in the graph.
PRIMARY_PREFIXES={
 'origin':r'^(AX\.|BOOT\.|SEAM\.(?:ORIGIN|P1|P2|MARK|MILNOR)|FOUND\.)',
 'e8':r'^(E8\.|CARR\.|FORM\.LADDER|SEAM\.(?:EQUIV|GLUE|D5|COVER)|SM\.REP)',
 'clocks':r'^(CLOCK\.|CLK\.|TIME\.|COX\.|SEAM\.CLOCK)',
 'flavor':r'^(FLAV\.|FLAVOR\.|MASS\.|YUK\.|NEUT\.|SM\.)',
 'alpha':r'^(EM\.|ALPHA\.|U1\.)',
 'quantum':r'^(QG\.|QGEO\.|QGAMB\.|QFT\.|QFT4D\.|CODE\.|CLIFF\.|PAULI\.|SEAM\.EQUIV)',
 'gravity':r'^(GRAV\.|QGEO\.|QGAMB\.|QG\.|SEAM\.KHALF)',
 'cosmology':r'^(COSMO\.|COSM\.|INFL\.|BARY\.|VAC\.|DE\.)',
 'rh':r'^(RH\.|PRIME\.|HECKE\.|ZETA\.|TFPT\.POSITIVE)',
 'validation':r'^(REVIEW\.|VALID\.|PRED\.|TEST\.|GATE\.|OBS\.)',
 'universalraum':r'^(UR\.|UNIV\.|MSC\.|PARENT\.|SOURCE\.)',
 'metrology':r'^(METRO\.|CALIB\.|SCALE\.|UNITS\.)',
 'spacetime':r'^(QFT4D\.|CHIRAL4D\.|SPACETIME\.|DIM\.|LORENTZ\.|CONTINUUM\.|RG\.)',
 'process':r'^(PROCESS\.|DYN\.|EVENT\.|PARENT\.|SOURCE\.|INSTR\.|RECORD\.)'}
PRIMARY_TERMS={
 'origin':r'\borigin\b|\bdeclared input\b|\bursprung\b',
 'e8':r'\be8\s+(?:lattice|root|glue|closure)|\bd5\s*(?:plus|\+)|\bhalf.spinor\b',
 'clocks':r'\bclock\b|\bcoxeter\b',
 'flavor':r'\bflavou?r\b|\byukawa\b|\bneutrino\b',
 'alpha':r'\bfine.structure\b|\bthomson\b|\balpha fix',
 'quantum':r'\bquantum (?:code|state|instrument|source)|\bpauli\b|\bclifford\b',
 'gravity':r'\bgravit|\beinstein\b|\barea.law\b|\bmodular boost\b',
 'cosmology':r'\bcosmolog|\binflation\b|\bbaryon\b|\bdark energy\b',
 'rh':r'\briemann\b|\bprime\b|\bhecke\b',
 'validation':r'\bholdout\b|\bfalsifi|\breproduction\b|\breproducib',
 'universalraum':r'\buniversalraum\b|\boperational equivalence\b|\bmsc\b',
 'metrology':r'\bmetrolog|\bcalibrat|\babsolute scale\b|\bv_geo\b',
 'spacetime':r'\b3\+1|\bspacetime\b|\blorentz|\bcontinuum limit\b',
 'process':r'\bevent law\b|\bprocess law\b|\bsource process\b|\bereignis'}

def topic_role(row,topic):
    key=topic.removeprefix('topic:')
    # Restrict lexical focus to the principal claim sentence, not historical tails.
    principal=row.get('claim','').split(' | ')[0].split(';')[0][:260]
    primary=bool(re.search(PRIMARY_PREFIXES.get(key,r'(?!)'),row['claim_id']) or re.search(PRIMARY_TERMS.get(key,r'(?!)'),principal,re.I))
    return 'primary' if primary else 'related'
def atomic(path,data):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(json_transport_value(data),ensure_ascii=False,separators=(',',':'),allow_nan=False));tmp.replace(path)

def build(repo=REPO):
    with (repo/'verification/status_ledger.csv').open(newline='') as f:rows=list(csv.DictReader(f))
    by_id={r['claim_id']:r for r in rows};nodes={};edges=[];sources={};seen=set()
    main=json.loads((HERE/'data/knowledge.json').read_text()) if (HERE/'data/knowledge.json').is_file() else {'nodes':[],'edges':[]}
    existing={n['id']:n for n in main['nodes']}
    orientation=json.loads((HERE/'orientation.json').read_text())
    topics={t['id']:t for t in orientation['topics']}
    memberships={r['claim_id']:set(existing.get('claim:'+r['claim_id'],{}).get('attrs',{}).get('semantic',{}).get('topic_ids',[])) for r in rows}
    for e in main['edges']:
        if e['type']=='about' and e['src'].startswith('claim:'):memberships.setdefault(e['src'].removeprefix('claim:'),set()).add(e['dst'])
    for t in topics.values():
        keywords=[k.lower() for k in t.get('keywords',[]) if len(k)>2]
        for r in rows:
            text=(r['claim_id']+' '+r['claim']).lower()
            if any(k in text for k in keywords):memberships[r['claim_id']].add(t['id'])
    memberships={cid:{t for t in ts if t in topics} for cid,ts in memberships.items()}
    # The primary navigation set comes from the stable ledger, so it must not
    # disappear when an AI annotation is pending or being replaced.
    for row in rows:
        for topic in topics:
            if topic_role(row,topic)=='primary':memberships[row['claim_id']].add(topic)
    def edge(src,dst,kind,attrs=None):
        key=(src,dst,kind)
        if key in seen:return
        seen.add(key);edges.append({'id':'edge:'+sha('|'.join(key).encode())[:20],'src':src,'dst':dst,'type':kind,'attrs':attrs or {}})
    def source(path):
        b=path.read_bytes();rel=str(path.relative_to(repo));digest=sha(b);sid='source:'+sha(f'{rel}|{digest}'.encode())[:20]
        sources[sid]={'id':sid,'path':rel,'sha256':sha(b),'line_start':1,'line_end':len(b.decode('utf-8').splitlines())}
        return sid
    for r in rows:
        cid=r['claim_id'];nid='claim:'+cid
        nodes[nid]={'id':nid,'type':'claim','label':cid,'summary':'','attrs':{'assessment':assess(r),'topic_roles':{t:topic_role(r,t) for t in topics},'topic_role_method':'navigation_focus_rules_not_evidence'}}
        for t in memberships[cid]:edge(nid,t,'about',{'method':'topic_mapping','notice':'Thematische Zuordnung; keine Herleitung'})
    for t in topics.values():
        nodes[t['id']]={'id':t['id'],'type':'topic','label':t['label'],'summary':t['summary'],'attrs':{k:t[k] for k in ['label_en','summary_en','topic_orientation'] if k in t}}
    # Read source-native records directly. The Atlas no longer depends on a
    # generated verification/theory_graph.json snapshot.
    inventory=source_inventory.build(repo,set(by_id))
    for s in inventory['sources']:
        sources[s['id']]=s
    records={n['id']:n for n in inventory['nodes'] if n['type'] in ['contract','rh','experiment','scorecard','failure']}
    for nid,record in records.items():
        a=dict(record.get('attrs',{}));refs=list(a.get('source_refs',[]));path=None
        def add_ref(source_path):
            source_id=source(source_path)
            if source_id and source_id not in refs:refs.append(source_id)
        if record['type']=='contract':
            directory=repo/'experiments/theory-contracts'/nid.removeprefix('contract:')
            p=directory/'contract_index.json'
            if p.is_file():
                try:original=json.loads(p.read_text());a=merge_original_contract(a,original)
                except ValueError:original={}
                a['verdict_raw']=original.get('verdict',a.get('verdict_raw'));a['question']=original.get('question',a.get('question'))
                a['claims_tested']=original.get('claims_tested',a.get('claims_tested',[]));a['date']=original.get('date',a.get('date'))
                add_ref(p);path=p
            for filename in ['PROOF.txt','RESULTS.md','README.md']:
                p=directory/filename
                if p.is_file() and p.stat().st_size<100000:add_ref(p)
        elif record['type']=='rh':
            p=repo/nid.removeprefix('rh:')
            if p.is_file() and p.resolve().is_relative_to(repo.resolve()) and p.suffix in ['.py','.lean','.md','.tex','.json'] and p.stat().st_size<100000:add_ref(p);path=p
        a.update({'source_refs':refs,'inventory_source':'source_native_inventory','attempt_kind':record['type'],
            'notice':'Versuch/Befund mit dokumentiertem Geltungsbereich. PASS allein schließt keine physikalische Gesamtfrage.',
            'reviewed':False if a.get('parse_source') in ['contract_index_llm','llm_backfill'] else None})
        label=a.get('question') or record.get('label') or nid
        summary=a.get('summary_de') or a.get('verdict_raw') or a.get('outcome') or a.get('mechanism') or ''
        nodes[nid]={'id':nid,'type':'attempt','label':str(label),'summary':str(summary),'attrs':a}
        tested=set(x for x in a.get('claims_tested',[]) if isinstance(x,str) and x in by_id)
        for cid in tested:edge(nid,'claim:'+cid,'tests',{'source':'contract_claims_tested'})
        for t in set().union(*(memberships[c] for c in tested)) if tested else set():edge(nid,t,'about',{'method':'linked_claim_topics'})
    allowed={'tests','kills','excludes','leaves_open','contradicts','supports','followup_of','references','constrains_claim','sourced_from','about'}
    for e in inventory['edges']:
        if e['type'] not in allowed:continue
        if e['src'] not in records and e['dst'] not in records:continue
        if e['src'] not in nodes and e['src'] not in existing:continue
        if e['dst'] not in nodes and e['dst'] not in existing:continue
        edge(e['src'],e['dst'],e['type'],{**e.get('attrs',{}),'source':e.get('attrs',{}).get('source','legacy_index'),'notice':'Übernommene Quellenrelation; keine neue Beweisprüfung'})
        for src,dst in [(e['src'],e['dst']),(e['dst'],e['src'])]:
            if src in records and dst.startswith('claim:'):
                for t in memberships.get(dst.removeprefix('claim:'),[]):edge(src,t,'about',{'method':'linked_claim_topics'})
    # Current contracts absent from an older generated index must remain visible.
    for p in (repo/'experiments/theory-contracts').glob('*/contract_index.json'):
        nid='contract:'+p.parent.name
        if nid in nodes:continue
        try:a=json.loads(p.read_text())
        except ValueError:continue
        refs=[source(p)]
        for filename in ['PROOF.txt','RESULTS.md','README.md']:
            q=p.parent/filename
            if q.is_file() and q.stat().st_size<100000:refs.append(source(q))
        nodes[nid]={'id':nid,'type':'attempt','label':a.get('question',p.parent.name),'summary':str(a.get('verdict','unknown')),
            'attrs':{**a,'original_contract':a,'source_refs':refs,'attempt_kind':'contract','notice':'Originalvertrag; keine automatische Claim-Promotion'}}
        for cid in a.get('claims_tested',[]):
            if cid not in by_id:continue
            edge(nid,'claim:'+cid,'tests')
            for t in memberships[cid]:edge(nid,t,'about',{'method':'linked_claim_topics'})
    return {'meta':{'schema_version':1,'generated_at':datetime.now(timezone.utc).isoformat(),'kind':'ledger_assessment_and_attempts','notice':'Status aus Ledger; Versuchsergebnisse werden nicht als Gesamtabschluss gezählt.'},
        'nodes':list(nodes.values()),'edges':edges,'sources':list(sources.values())}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,default=REPO);a=p.parse_args()
    g=build(a.repo);atomic(HERE/'data/context.json',g)
    print(json.dumps({'nodes':len(g['nodes']),'edges':len(g['edges']),'attempts':sum(n['type']=='attempt' for n in g['nodes'])}))
