'use strict';

const $ = id => document.getElementById(id);
const ui = { data:null, run:null, defaultRun:null, selected:null, follow:true, filter:'all', step:0,
  frozen:false, fetching:false, listKey:'', detailKey:'', eventKey:'', lastCount:null };
const number = v => new Intl.NumberFormat('de-DE').format(Number(v || 0));
function el(tag, cls, text) {
  const node = document.createElement(tag);
  if (cls) node.className = cls;
  if (text !== undefined && text !== null) node.textContent = String(text);
  return node;
}
function icon(kind) {
  const svg = document.createElementNS('http://www.w3.org/2000/svg','svg');
  svg.setAttribute('viewBox','0 0 20 20'); svg.setAttribute('aria-hidden','true');
  const p = document.createElementNS(svg.namespaceURI,'path');
  p.setAttribute('d', kind==='pass' ? 'm4 10 4 4 8-8' : kind==='failed' ? 'm5 5 10 10M15 5 5 15' : 'M10 4v6m0 3v1');
  svg.append(p); return svg;
}
function date(value) {
  if (!value) return null;
  const d = new Date(typeof value==='number' && value<1e12 ? value*1000 : value);
  return Number.isNaN(d.valueOf()) ? null : d;
}
function clock(value) { const d=date(value); return d ? d.toLocaleTimeString('de-DE',{hour:'2-digit',minute:'2-digit'}) : '—'; }
function ago(value) {
  const d=date(value); if(!d) return 'Zeitpunkt unbekannt';
  const seconds=Math.max(0,Math.round((Date.now()-d)/1000));
  if(seconds<5)return 'gerade aktualisiert'; if(seconds<60)return 'vor '+seconds+' s aktualisiert';
  if(seconds<3600)return 'vor '+Math.floor(seconds/60)+' Min. aktualisiert';
  return 'aktualisiert um '+clock(value);
}
function calculated(c) { return c.has_computed_values === true || (c.computed && c.computed.Z0 !== null && c.computed.Z0 !== undefined); }
function failed(c) { return c.status==='CANDIDATE_FAIL' || c.checks?.some(x=>x.status==='FAIL'); }
function verified(c) { return c.status==='SCOPED_HIT' && c.checks?.length>0 && c.checks.every(x=>x.status==='PASS'); }
function allNodes(ast) { if(!ast||typeof ast!=='object')return []; return [ast,...(ast.args||[]).flatMap(allNodes)]; }
function trialInfo(c) {
  const t=c.trial||{}, type=c.change_type||t.change_type, family=c.structure_family||t.structure_family||'';
  const modes=(JSON.stringify(c.ast||{}).match(/(?:modes|chart)[^0-9]{0,8}([24])/i)||[])[1];
  const names={starting_program:'Vergleichsprogramm',action_vertex:'Zusätzlicher Wechselwirkungsterm',source_contact:'Geänderte Probeeinfügung',kernel_correction:'Zusatz im Quellenkern',coupled_vertices:'Gekoppelte Wechselwirkungen',whole_sewing:'Ganze Quellenstücke verbinden',whole_word_sewing:'Ganze Quellenstücke verbinden',coupled_action_vertices:'Gekoppelte Wechselwirkungen',source_insertion:'Geänderte Probeeinfügung',marked_action_vertex:'Geänderte Probeeinfügung',whole_ordered_sewing:'Ganze Quellenstücke verbinden',whole_kernel_correction:'Zusatz im Quellenkern',model_construction:'Sol Vorschlag'};
  let name=names[type];
  if(!name && /sew|word/.test(family))name='Ganze Quellenstücke verbinden';
  if(!name && /contact/.test(family))name='Geänderte Probeeinfügung';
  if(!name && /kernel/.test(family))name='Zusatz im Quellenkern';
  const degree=family.match(/degree_(\d+)/)?.[1];
  return {name,degree,modes,type,family,t};
}
function title(c) {
  const v=trialInfo(c);
  if(v.name)return v.name+(v.type==='starting_program'&&v.modes?' · '+v.modes+' Moden':v.degree?' · Grad '+v.degree:'');
  const ids=allNodes(c.ast).map(x=>x.id);
  if(ids.includes('native.sew'))return 'Ganze Quellenstücke verbinden';
  if(ids.includes('construct.kernel_add'))return 'Quellenkern aus mehreren Termen';
  return 'Rechenvorschlag zur gemeinsamen Quelle';
}
function changeDescription(c) {
  const v=trialInfo(c), raw=c.what_changes||v.t.what_changes||'';
  if(v.type==='action_vertex') {
    const mask=raw.match(/mask (\d+)/)?.[1];
    return 'Ein zusätzlicher Term mit '+(v.degree||'mehreren')+' internen Rechenvariablen verändert die Wechselwirkung'+(mask?' (Variablenauswahl '+mask+')':'')+'. Die vollständige endliche Exponentialfunktion wird neu ausgewertet.';
  }
  if(v.type==='starting_program'){const ids=allNodes(c.ast).map(x=>x.id);return 'Eine gewählte endliche Rechenkarte mit '+(v.modes||'wenigen')+' Moden dient als Vergleich. '+(ids.includes('native.c3')?'Der Vorfaktor wird aus c3, der vierfache Wechselwirkungsterm aus dem C6 Verhältnis angesetzt.':'Der Vorfaktor wird auf 1 gesetzt und das negative C6 Verhältnis als Wechselwirkungsterm angesetzt.')+' Diese Zuordnung ist eine Testannahme; die ursprüngliche Feldintegration folgt daraus noch nicht.';}
  if(v.type==='model_construction')return raw||c.hypothesis||'Sol kombiniert vorhandene Rechenregeln zu einem neuen ausführbaren Vorschlag.';
  if(/sew|word|sewing|composition/.test(v.type+' '+v.family))return 'Ganze Quellenstücke werden in einer anderen Reihenfolge verbunden. Die gemeinsame Antwort wird erst nach der vollständigen Verbindung berechnet.';
  if(/contact|insertion|marked_action/.test(v.type+' '+v.family))return 'Die Probeeinfügung in dieselbe Quelle wird verändert. Danach wird geprüft, wie sich die Antworten mit einer oder beiden Einfügungen ändern.';
  if(/coupled|combined/.test(v.type+' '+v.family))return 'Zwei Wechselwirkungsterme wirken gemeinsam in derselben Quelle. Ihre Beiträge werden vor der Auswertung kombiniert; das gesamte Programm wird danach neu berechnet.';
  if(/kernel/.test(v.family))return 'Ein zusätzlicher Term verändert den gesamten Quellenkern. Die vier Antworten werden gemeinsam daraus neu berechnet.';
  return raw||'Vorhandene Rechenregeln werden zu einer auswertbaren Quelle kombiniert.';
}
function candidateSummary(c) {
  if(failed(c))return 'Eine berechnete Gegenprobe widerspricht diesem Versuchsprogramm. Der genaue Befund steht in den Details.';
  if(verified(c))return 'Der Kandidat erfüllt die für diesen Zielbereich verlangten Prüfungen. Der genaue Geltungsbereich bleibt unten sichtbar.';
  if(calculated(c))return 'Dieses Versuchsprogramm liefert vier Antworten aus derselben konstruierten Quelle. Die Zuordnung zur ursprünglichen TFPT Quelle ist noch nicht belegt.';
  return 'Dieser Vorschlag kann noch nicht vollständig ausgewertet werden. Es fehlen Eingaben oder ein Rechenschritt ist fehlgeschlagen.';
}
function pretty(value) {
  if(value===null||value===undefined)return 'noch offen';
  if(typeof value==='object')return value.display||JSON.stringify(value);
  return String(value).replaceAll('pi','π').replaceAll('**','^').replaceAll('*',' · ').replaceAll('/',' / ');
}
function candidateById(id) { return ui.data?.candidates?.find(c=>c.id===id); }
const stepCopy=[
  ['Ein Versuch ist eine konkrete Rechenanweisung.','Zum Beispiel: einen Wechselwirkungsterm hinzufügen und dieselbe Quelle erneut ausrechnen. OpenAI Sol schlägt Programme vor. Daneben variiert die direkte Suche die vorhandenen Regeln systematisch. Wiederholungen zählen gesondert.'],
  ['Ein Versuchsprogramm muss selbst Antworten liefern.','Das Programm baut den gewichteten Quellenkern und berechnet vier zusammengehörige Markierungen. Auch höhere Wechselwirkungen und vollständig verbundene Quellenstücke sind möglich.'],
  ['Gegenproben suchen den ersten Widerspruch.','Die Rechnung wird unabhängig nachgerechnet. Ein Fehler verwirft den Kandidaten. Ein noch fehlender Originalbeleg wird als offen gespeichert und mit einem echten Widerspruch unterschieden.'],
  ['Hier entscheidet sich die ursprüngliche Herkunft.','Ein berechneter Quellenkern muss dieselbe Quelle, Zeit und Markierung wie die ursprüngliche TFPT beschreiben. Im damals aufgezeichneten Lauf war der automatische Ganzquellenvergleich noch nicht angeschlossen.']
];
function renderStep() {
  if(isResearch())return renderResearchStep();
  document.querySelectorAll('[data-step]').forEach(b=>b.setAttribute('aria-pressed',String(Number(b.dataset.step)===ui.step)));
  let [head,copy]=stepCopy[ui.step];
  if(ui.step===3 && ui.data?.run?.whole_source_identity_check)copy='Die ganze parametrisierte Quelle wird einschließlich Zeit, Markoperatoren und äußerer Referenz mit einer unabhängig bereitgestellten Originalfamilie verglichen. Vier passende Zahlen genügen dafür nicht. Fehlt die ursprüngliche Feldintegration, bleibt der Nachweis als fehlende Eingabe offen.';
  if(ui.step===0 && ui.data && ui.data.run?.model!=='gpt-6.1-sol')copy='Die aufgezeichneten Rechenregeln dieses ausgewählten Laufs werden zu Programmen kombiniert. Das verwendete Modell und die tatsächlich gespeicherten Ergebnisse stehen daneben. Die Einstellungen früherer Läufe können abweichen.'; $('stepExplanation').replaceChildren(el('strong','',head),el('p','',copy));
}
function renderHeader(data) {
  renderRunMode(data);
  const process=data.process||{}, counts=data.summary?.counts||{};
  const states={running:'Suche läuft',finished:'Lauf abgeschlossen',stopped:'Lauf angehalten',unknown:'Prozessstatus nicht bestätigt'};
  $('runStatus').textContent=states[process.status]||'Status wird gelesen';
  $('runDot').className='status-dot '+(process.status||'unknown');
  $('runAge').textContent=data.run?.started_at ? (process.status==='running'?'seit ':'gestartet um ')+clock(data.run.started_at) : '';
  $('syncTime').textContent=ago(data.freshness?.last_file_update);
  $('countPrograms').textContent=number(counts.whole_candidates);
  $('countRounds').textContent=number(counts.model_rounds);
  $('countNovel').textContent=number(counts.new_source_signatures ?? data.trial_summary?.new_source_signatures);
  $('countRepeats').textContent=number(counts.duplicate_calculations ?? data.trial_summary?.duplicate_calculations);
  $('flowProposals').textContent=number(counts.whole_candidates)+' Programme';
  $('flowCalculated').textContent=number(counts.computed_candidates ?? data.candidates.filter(calculated).length)+' berechnet';
  $('flowChecks').textContent=number(counts.failed_candidates)+' mit Widerspruch';
  $('flowProofs').textContent=number(counts.scoped_hits)+' bestätigte Zieltreffer';
  const latest=[...(data.events||[])].reverse().find(x=>x.kind!=='started');
  const blocked=['NO_COMPLETION_ROUTE','NO_SELECTION_SIGNAL'].includes(process.declared_status);
  const activity=(process.status==='running'||blocked) ? data.activity : null;
  $('nowTitle').textContent=activity?.label||latest?.title||data.summary?.headline||'Die Suche wartet auf den nächsten Schritt.';
  $('nowDescription').textContent=(activity?.description ? activity.description+(activity.total?' Versuch '+number(activity.current)+' von '+number(activity.total)+' in diesem Paket.':'') : null)||latest?.description||data.summary?.explanation||'Sobald neue Daten vorliegen, erscheinen sie hier.';
  $('liveLabel').textContent=blocked?'Ergebnis der Suchprüfung':process.status==='running'?'Jetzt im Suchlauf':'Zuletzt im Lauf';
  const chosenModel=data.run?.model||data.run?.requested_model;
  $('modelName').textContent=(chosenModel==='gpt-6.1-sol'?'OpenAI Sol':chosenModel||'Nicht aufgezeichnet')+(data.run?.model_receipt_verified?' · bestätigt':blocked?' · Aufruf ausgesetzt':' · noch kein bestätigter Aufruf');
  if(data.run?.historical){$('runStatus').textContent='Historischer Lauf · '+(states[process.status]||'Status unbekannt');$('modelName').textContent='Historischer Modelllauf';$('modelLabel').textContent='Archivansicht';}else{$('modelLabel').textContent='Modell';}
  const finished=data.summary?.finished_rounds ?? counts.model_rounds ?? 0;
  const requested=data.summary?.requested_rounds;
  $('roundValue').textContent=requested ? number(finished)+' von '+number(requested) : number(finished);
  $('roundLabel').textContent='Abgeschlossene Suchrunden';
  $('proofStatus').textContent=counts.scoped_hits>0?'Zielbereich prüfen':'Herkunft offen';
  $('boundaryTitle').textContent=counts.scoped_hits>0 ? 'Ein Treffer braucht seinen genauen Geltungsbereich.' : 'Berechneter Kandidat und fundamentale Lösung bleiben getrennt.';
  $('boundaryDescription').textContent='Der Lauf erzeugt auswertbare Quellenhypothesen. Die Herkunft genau dieser Quelle und die vollständige physikalische Herleitung müssen noch bewiesen werden. Die Programmsprache erzeugt diese Beweise bisher nicht selbst.';
  const manifest=data.run?.knowledge_context, coverage=manifest?.counts;
  $('knowledgeCoverage').textContent=coverage
    ? 'Dieser Lauf verwendet '+number(coverage.findings)+' Atlasbefunde, '+number(coverage.premises)+' Voraussetzungen, '+number(coverage.strategies)+' Strategien und '+number(coverage.completed_plan_subtasks)+' abgeschlossene Planunteraufgaben mit ihrem Geltungsbereich. Dazu kommen '+number(data.run.prior_trial_count)+' verschiedene frühere Versuchsprogramme. Ausschlüsse und fehlgeschlagene Versuche bleiben sichtbar; Quellenhashes bestätigen nur die Dateifrische.'
    : 'Für diesen früheren Lauf ist kein vollständiges Wissensmanifest aufgezeichnet. Sein damaliges Paket kann nur einen Teil des bekannten Stands enthalten. Der neue Lauf überträgt das vollständige geprüfte Atlaswissen.';
  const trials=data.trial_summary||{};
  const directBatch=data.run?.batch_size, proposalBatch=data.run?.proposal_batch;
  $('scaleTitle').textContent=blocked?'Warum die Suche hier anhält':'Warum jetzt mehr Versuche möglich sind';
  $('scaleExplanation').textContent=blocked
    ? 'Die Vergleichsprogramme liefern noch keinen auswertbaren physikalischen Prüfwert. Die aktuelle Programmsprache erzeugt die dafür fehlenden Originaldaten nicht. Weitere Varianten würden denselben Engpass wiederholen. Als Nächstes muss die ursprüngliche Integrationsvorschrift selbst berechenbar werden.'
    : directBatch
    ? 'Dieser Lauf erzeugt direkt ganze Rechenversuche: bis zu '+number(directBatch)+' systematische Programme pro Runde'+(proposalBatch?' und bis zu '+number(proposalBatch)+' Vorschläge je Modellaufruf':'')+'. Früher gingen viele Schritte für kleine Zwischenprogramme und wiederholte Modellideen drauf.'
    : 'Für diesen früheren Lauf sind keine direkten Programmpakete dieser Art aufgezeichnet. Die Zahlen zeigen seine tatsächlich gespeicherten Ergebnisse. Modellaufrufe und ausgewertete Programme sind unterschiedliche Mengen.';
  $('scaleDetail').textContent=number(trials.generated ?? counts.whole_candidates)+' Versuche erzeugt · '+number(trials.evaluated ?? counts.whole_candidates)+' ausgewertet · '+number(trials.model_proposals ?? counts.model_proposals)+' Modellvorschläge. Verschiedene Quellenrechnungen zählen den gesamten berechneten Kern; die vier hier gezeigten Antworten können trotzdem gleich sein. Mehr Varianten liefern erst dann eine physikalische Lösung, wenn auch ihre ursprüngliche Herkunft bewiesen ist.';
  if(data.run?.whole_source_identity_check)$('boundaryDescription').textContent='Der Ganzquellenvergleich ist angeschlossen. Ein Zieltreffer verlangt seine bestandene Quellen-, Zeit- und Markierungsprüfung sowie die weiteren physikalischen Nachweise. Fehlende Originaldaten bleiben offen.';
  $('coverageText').textContent=data.coverage?.detail||'Die Suche erfasst die deklarierte Programmsprache und ihren begrenzten Suchraum.';
  if(ui.lastCount!==null && counts.whole_candidates>ui.lastCount) {
    $('liveAnnounce').textContent=number(counts.whole_candidates-ui.lastCount)+' neue Versuchsprogramme wurden gespeichert.';
  }
  ui.lastCount=counts.whole_candidates;
}
function countsTotal(){return ui.data?.summary?.counts?.whole_candidates||0;}
function selectCandidate(id, fromUser=false) {
  ui.selected=id;
  if(fromUser) ui.follow=false;
  $('followButton').setAttribute('aria-pressed',String(ui.follow));
  renderCandidates(); renderDetail();
}
function renderCandidates() {
  if(isResearch())return renderResearchRounds();
  const candidates=ui.data?.candidates||[];
  const visible=candidates.filter(c=>ui.filter==='all'||ui.filter==='calculated'&&calculated(c)||ui.filter==='failed'&&failed(c)||ui.filter==='open'&&!failed(c)&&!verified(c));
  const key=JSON.stringify([visible.map(c=>[c.id,c.status,c.computed]),ui.selected,ui.filter]);
  $('candidateCount').textContent=number(visible.length)+' angezeigt'+(countsTotal()>candidates.length?' · Auswahl von '+number(countsTotal()):'');
  if(key===ui.listKey)return;
  ui.listKey=key;
  const list=$('candidateList');
  if(!visible.length) {list.replaceChildren(el('div','empty-state',ui.filter==='all'?'Noch kein ganzes Versuchsprogramm gespeichert. Die Suche baut zunächst Rechenbausteine auf.':'In diesem Filter liegen noch keine Versuchsprogramme vor.'));return;}
  const fragment=document.createDocumentFragment();
  [...visible].reverse().forEach(c=>{
    const index=c.ordinal||candidates.findIndex(x=>x.id===c.id)+1;
    const b=el('button','candidate-row'+(ui.selected===c.id?' selected':''));b.type='button';b.setAttribute('aria-pressed',String(ui.selected===c.id));
    const copy=el('div'), name=el('div','candidate-name');name.append(el('span','candidate-number',String(index).padStart(2,'0')),el('span','',title(c)));
    const description=changeDescription(c);
    copy.append(name,el('div','candidate-description',description));
    const status=el('div','candidate-row-status');
    status.append(el('span','status-label '+(calculated(c)?'calculated':'open'),calculated(c)?'Berechnet':'Rechnung offen'));
    status.append(el('span','status-label '+(failed(c)?'failed':verified(c)?'calculated':'open'),failed(c)?'Gegenprobe gescheitert':verified(c)?'Zielbereich erfüllt':'Herkunft offen'));
    b.append(copy,status);b.addEventListener('click',()=>selectCandidate(c.id,true));fragment.append(b);
  });
  list.replaceChildren(fragment);
}
function proofItem(kind,copy) {const li=el('li',kind);li.append(icon(kind),el('span','',copy));return li;}
function disclosure(titleText, content, id) {const d=el('details','detail-disclosure');d.dataset.disclosure=id;d.append(el('summary','',titleText),...content);return d;}
function renderDetail() {
  if(isResearch())return renderResearchDetail();
  const c=candidateById(ui.selected), target=$('candidateDetail');
  if(!c){target.replaceChildren(el('div','empty-state','Für diesen Lauf ist noch kein Versuchsprogramm gespeichert.'));return;}
  const key=JSON.stringify([c.id,c.status,c.computed,c.checks,c.trial,c.tested_difference]);
  if(key===ui.detailKey)return;
  const opened=new Set([...target.querySelectorAll('details[open]')].map(d=>d.dataset.disclosure));
  ui.detailKey=key;
  const header=el('div','detail-header');header.append(el('h3','',title(c)),el('span','detail-badge'+(failed(c)?' failed':''),failed(c)?'Widerspruch gefunden':verified(c)?'Zielbereich erfüllt':'Originalbelege offen'));
  const nodes=[header,el('p','detail-summary',candidateSummary(c)),el('div','change-explanation',changeDescription(c))];
  const firstFailure=c.checks?.find(x=>x.status==='FAIL'||x.status==='ERROR');
  if(firstFailure)nodes.push(el('div','failure-explanation',(firstFailure.label||firstFailure.id)+': '+(firstFailure.detail||c.first_failure||'Der aufgezeichnete Prüfer hat diesen Kandidaten zurückgewiesen.')));
  const diff=c.tested_difference||c.trial?.tested_difference;
  if(diff && Object.keys(diff).length) {
    const table=el('table','difference-table');table.append(el('caption','','Gegenüber dem gespeicherten Vergleichsprogramm'));
    const head=el('tr');['Antwort','Vorher','Jetzt','Änderung'].forEach(x=>head.append(el('th','',x)));const th=el('thead');th.append(head);table.append(th);
    const body=el('tbody');['Z0','ZU','ZY','ZUY'].forEach(k=>{if(!diff[k])return;const row=el('tr');[k,pretty(diff[k].before),pretty(diff[k].after),pretty(diff[k].delta)].forEach(x=>row.append(el('td','',x)));body.append(row);});table.append(body);
    nodes.push(table,el('p','result-explainer','Eine Änderung von 0 bedeutet: Dieser Eingriff verändert diese Antwort im gewählten Modell nicht. Auch das wird als Ergebnis gespeichert.'));
  }
  const path=el('ol','path-list');(c.human_steps||[]).forEach(text=>path.append(el('li','',text)));nodes.push(path);
  if(calculated(c)) {
    const table=el('table','values-table');table.append(el('caption','','Die vier Antworten und ihr Zusammenhang'));
    const body=el('tbody');
    [['Quelle ohne Markierung','Z0'],['Mit Markierung U','ZU'],['Mit Markierung Y','ZY'],['Mit beiden Markierungen','ZUY'],['Verbundenes Signal','connected']].forEach(([name,k])=>{const row=el('tr');row.append(el('td','',name),el('td','',pretty(c.computed[k])));body.append(row);});
    table.append(body);nodes.push(table,el('p','result-explainer','U und Y sind die vorgeschlagenen Probeeinfügungen in dieselbe Quelle. Ihre ursprüngliche TFPT Bedeutung ist für dieses Versuchsprogramm noch zu belegen.'));
  }
  const proofs=el('ul','proof-list');
  const calculation=c.checks?.find(x=>x.id==='constructed_common_values');
  proofs.append(proofItem(calculation?.status==='PASS'?'pass':calculation?.status==='FAIL'?'failed':'open',calculation?.status==='PASS'?'Die vier Antworten wurden unabhängig nachgerechnet.':calculation?.status==='FAIL'?'Die unabhängige Nachrechnung widerspricht dem Ergebnis.':'Die vollständige Nachrechnung ist noch offen.'));
  const identity=c.checks?.find(x=>x.id==='primitive_original_source_identity');
  proofs.append(proofItem(identity?.status==='PASS'?'pass':identity?.status==='FAIL'?'failed':'open',identity?.status==='PASS'?'Die ursprüngliche Quellenidentität ist in diesem Zielbereich belegt.':'Die Identität mit der ursprünglichen Quelle ist noch nicht belegt.'));
  const time=c.checks?.find(x=>x.id==='original_marked_time_word_match');
  proofs.append(proofItem(time?.status==='PASS'?'pass':'open',time?.status==='PASS'?'Ursprüngliche Zeit und Einfügungen passen.':'Ursprüngliche Zeit und Einfügungen müssen noch zugeordnet werden.'));
  nodes.push(proofs);
  const assumptions=(c.hypotheses||[]).map(x=>el('p','',x));
  if(c.proposed_origin_argument)assumptions.push(el('h4','','Vorgeschlagene Erklärung des Ursprungs'),el('p','',c.proposed_origin_argument),el('p','','Diese Modellargumentation ist ein Vorschlag. Sie wird durch ihre Anzeige nicht zum Originalbeweis.'));
  if(assumptions.length)nodes.push(disclosure('Annahmen und Ursprungsargument ansehen',assumptions,'assumptions'));
  const checks=el('table','check-table');const cb=el('tbody');
  (c.checks||[]).forEach(x=>{const tr=el('tr');const status=x.status==='PASS'?'Bestanden':x.status==='FAIL'?'Widerspruch':x.status==='ERROR'?'Prüffehler':'Noch offen';const checkCopy=el('td');checkCopy.append(el('strong','',x.label||x.id));if(x.detail)checkCopy.append(el('p','check-detail',x.detail));tr.append(checkCopy,el('td','check-status '+(x.status==='PASS'?'pass':x.status==='FAIL'?'fail':''),status));cb.append(tr);});checks.append(cb);
  nodes.push(disclosure('Alle Prüfungen und das Programm ansehen',[el('p','','Vorhandene Teilbeweise bleiben erhalten. Hier wird geprüft, ob sie genau diesen Kandidaten tragen.'),checks,el('h4','','Tatsächlich ausgeführtes Programm'),el('pre','technical-code',JSON.stringify(c.ast,null,2)),el('p','','Zielbereich: '+(c.scope||'nicht aufgezeichnet')),el('p','','Damals aufgezeichnetes Modell: '+(ui.data?.run?.model||ui.data?.run?.requested_model||'nicht aufgezeichnet'))],'technical'));
  target.replaceChildren(...nodes);
  target.querySelectorAll('details').forEach(d=>{if(opened.has(d.dataset.disclosure))d.open=true;});
}
function renderEvents() {
  const events=[...(ui.data?.events||[])].reverse().slice(0,16);
  const key=JSON.stringify([events.map(e=>[e.id,e.title,e.description]),isResearch()?ui.data?.active_call:null]);if(key===ui.eventKey)return;ui.eventKey=key;
  if(!events.length){$('eventList').replaceChildren(el('li','empty-state','Noch keine Ereignisse gespeichert.'));return;}
  const fragment=document.createDocumentFragment();
  events.forEach(e=>{
    const li=el('li'), time=el('time','event-time',clock(e.time));if(date(e.time))time.dateTime=date(e.time).toISOString();
    const body=el('div','event-body');const severity=['invalid','error','failed'].includes(e.kind)?'error':['repeat','repeated','partial','backend_change'].includes(e.kind)?'warning':'';
    const researchActive=isResearch()&&e.kind==='research_activity';
    const eventTitle=researchActive?'Sol Aufruf aktiv':e.title||'Suchschritt gespeichert';
    const description=researchActive?researchActivityCopy(ui.data):isResearch()&&e.kind==='research_proposal'?plainResearchQuestion(e.description):isResearch()&&e.kind==='started'?'Der Ablauf verbindet einen Herleitungsaufruf mit einer unabhängigen Gegenprüfung.':e.description||'';
    body.append(el('strong','event-kind '+severity,eventTitle),el('p','',description));
    if(e.candidate_id && (candidateById(e.candidate_id)||researchRound(e.candidate_id))) {const b=el('button','event-link',isResearch()?'Forschungsrunde ansehen':'Versuchsprogramm ansehen');b.addEventListener('click',()=>{selectCandidate(e.candidate_id,true);$('detailTitle').scrollIntoView({behavior:'smooth',block:'start'});});body.append(b);}
    li.append(time,body);fragment.append(li);
  });$('eventList').replaceChildren(fragment);
}
function render(data) {
  ui.data=data;
  const candidates=isResearch()?data.research.rounds||[]:data.candidates||[];
  if(ui.follow && candidates.length)ui.selected=candidates[candidates.length-1].id;
  else if(!(isResearch()?researchRound(ui.selected):candidateById(ui.selected))&&candidates.length)ui.selected=candidates[0].id;
  renderHeader(data);if(isResearch())renderResearchHeader(data);renderStep();renderCandidates();renderDetail();renderEvents();
}
async function loadRuns() {
  const response=await fetch('/api/runs',{cache:'no-store'});if(!response.ok)throw Error('Laufliste nicht verfügbar');
  const body=await response.json();const select=$('runSelect');
  select.replaceChildren(...body.runs.map(r=>{const o=el('option','',(r.is_default?'Aktuell · ':'Historisch · ')+(r.label||r.id));o.value=r.id;return o;}));
  ui.defaultRun=body.default_run;
  ui.run=body.runs.some(r=>r.id===ui.run)?ui.run:body.default_run||body.runs[0]?.id;
  if(ui.run)select.value=ui.run;
}
async function poll() {
  if(ui.fetching)return;
  window.clearTimeout(ui.timer);
  ui.fetching=true;
  try {
    if(!ui.run)await loadRuns();
    const requestedRun=ui.run;
    const response=await fetch('/api/state'+(requestedRun?'?run='+encodeURIComponent(requestedRun):''),{cache:'no-store',signal:AbortSignal.timeout(10000)});
    if(!response.ok)throw Error('Suchdaten nicht verfügbar');
    const data=await response.json();
    if(requestedRun!==ui.run)return;
    $('errorBanner').hidden=true;
    $('connectionDot').className='connection-dot'+(ui.frozen?'':' connected');
    $('connectionText').textContent=ui.frozen?'Ansicht eingefroren':'Datenverbindung aktiv';
    if(!ui.frozen)render(data);
  } catch(error) {
    console.error('Suchmonitor:', error.message);
    $('errorBanner').hidden=false;$('connectionDot').className='connection-dot error';$('connectionText').textContent='Verbindung unterbrochen';
  } finally {ui.fetching=false;ui.timer=window.setTimeout(poll,2000);}
}
document.querySelectorAll('[data-step]').forEach(b=>b.addEventListener('click',()=>{ui.step=Number(b.dataset.step);renderStep();}));
document.querySelectorAll('[data-filter]').forEach(b=>b.addEventListener('click',()=>{ui.filter=b.dataset.filter;document.querySelectorAll('[data-filter]').forEach(x=>x.setAttribute('aria-pressed',String(x===b)));renderCandidates();}));
$('helpButton').addEventListener('click',()=>{const open=$('helpButton').getAttribute('aria-expanded')!=='true';$('helpButton').setAttribute('aria-expanded',String(open));$('flowHelp').hidden=!open;});
$('freezeButton').addEventListener('click',()=>{ui.frozen=!ui.frozen;$('freezeButton').setAttribute('aria-pressed',String(ui.frozen));$('freezeButton').querySelector('span').textContent=ui.frozen?'Ansicht fortsetzen':'Ansicht einfrieren';$('connectionText').textContent=ui.frozen?'Ansicht eingefroren':'Daten werden aktualisiert';if(!ui.frozen)poll();});
$('followButton').addEventListener('click',()=>{ui.follow=!ui.follow;$('followButton').setAttribute('aria-pressed',String(ui.follow));if(ui.follow&&ui.data)render(ui.data);});
function switchRun(id){ui.run=id;$('runSelect').value=id;ui.selected=null;ui.follow=true;ui.filter='all';ui.listKey=ui.detailKey=ui.eventKey='';ui.lastCount=null;document.querySelectorAll('[data-filter]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.filter==='all')));$('followButton').setAttribute('aria-pressed','true');poll();}
$('runSelect').addEventListener('change',()=>switchRun($('runSelect').value));
$('liveRunButton').addEventListener('click',()=>{if(ui.defaultRun)switchRun(ui.defaultRun);});
$('retryButton').addEventListener('click',poll);
$('downloadButton').addEventListener('click',()=>{if(!ui.data)return;const link=document.createElement('a');link.href=URL.createObjectURL(new Blob([JSON.stringify(ui.data,null,2)],{type:'application/json'}));link.download='TFPT-Suchstand-'+ui.run+'.json';link.click();setTimeout(()=>URL.revokeObjectURL(link.href),1000);});


const initialResearchQuestion='Resolve one actual original regulator block: derive nu_Vn and D_Vn in the same field chart and its unnormalized 1 and u=H^dagger H whole marked fibres with all fermion/rest/reference factors; alternatively evaluate the identical marked original IR transfer caller and prove its R113 source/time dictionary. Read the owning definitions before calculating. If a primitive input cannot yet be derived, obtain a new source-specific equation, exclusion or countertest that advances that exact derivation, not another finite kernel.';
function plainResearchQuestion(value){
  if(typeof value!=='string')return value;
  if(value.trim()!==initialResearchQuestion)return value;
  return 'Gesucht sind die ursprünglichen Feldgewichte νVn und DVn in derselben Feldkarte sowie die unnormierten Markierungen 1 und u = H†H mit allen Faktoren für Fermionen, übrige Felder und Referenz. Alternativ wird derselbe markierte ursprüngliche Transfer im Infrarotbereich berechnet und seine Zuordnung von Quelle und Zeit gemäß R113 belegt. Zuerst werden die maßgeblichen Definitionen gelesen. Fehlt eine primitive Eingabe, muss eine neue quellenbezogene Gleichung, ein Ausschluss oder eine Gegenprobe genau diese Herleitung voranbringen.';
}

function normalizedResearchCopy(value){return String(value||'').trim().replace(/\s+/g,' ').replace(/[.!?]+$/,'').toLowerCase();}
function briefResearchQuestion(value){
  if(typeof value!=='string')return value;
  if(value.trim()===initialResearchQuestion)return 'Wie entstehen die ursprünglichen Feldgewichte und die gemeinsamen Markierungen derselben Quelle?';
  if(value==='Für einen Originalblock p_mn identifizieren und zwischen Restriktion, transportierter Mittelung und effektivem zusammengesetztem Vertex entscheiden. Anschließend Delta_u im identischen vollständigen unnormierten Caller auswerten; bei Restriktion zuerst den lokalen Pullback prüfen.')return 'Welche ursprüngliche Abbildung fasst die Felder zusammen und erhält dabei die Higgsmarkierung?';
  if(value==='Die hinter Metrology005:166–168 beanspruchte D_Vn Konstruktion für einen benannten Originalblock verfolgen. Dort zuerst den Higgs und Metrikpullback prüfen; anschließend sämtliche Quellenableitungen und einen ganzen unnormierten Operatoranker pro zusammenhängender Komponente auswerten.')return 'Welche ursprüngliche Dirac Konstruktion bestimmt den Higgsreadout, die Metrik und den ganzen Quellenanker?';
  return 'Sol prüft jetzt die nächste offene Frage aus der Gegenprüfung an den Originalquellen. Die genaue Aufgabe kannst du darunter öffnen.';
}
function plainResearchStatus(value){return ({new_constraint:'Zusätzliche notwendige Bedingung',scoped_result_survives:'Bedingtes Ergebnis besteht die Gegenprüfung'})[value]||value;}
function researchRoundTitle(r){return r.proposal?.question||briefResearchQuestion(r.question)||'Forschungsrunde '+r.ordinal;}
function isKnownFirstRound(r){return r.id==='research-000001'&&r.proposal?.new_information?.startsWith('Die konkrete Markierungsbedingung lautet Δu(x)=')&&r.review?.new_information?.startsWith('Die bedingte Varianzidentität besteht. Eine rohe Mittelung');}
function researchOutcome(r){
  if(isKnownFirstRound(r))return 'Die bedingte Rechnung besteht; die ursprüngliche Feldabbildung bleibt laut Gegenprüfung offen.';
  const info=r.review?.new_information||r.proposal?.new_information;
  if(typeof info==='string')return info.match(/^.*?[.!?](?:\s|$)/)?.[0]?.trim()||info;
  return plainResearchStatus(r.critic_verdict)||'Ergebnis und Kritik ansehen.';
}

function plainResearchRole(value){return typeof value==='string'?value.replace(/\bproposer\b/gi,'Herleitung').replace(/\bcritic\b/gi,'Gegenprüfung'):value;}
function researchActivityCopy(data){
  const call=data.active_call||{};
  if(call.error)return plainResearchRole(String(call.error));
  const count=call.completed_commands;
  if(typeof count==='number'&&count>0)return 'Sol arbeitet an diesem Aufruf. '+number(count)+' Befehle sind beendet'+(typeof call.failed_commands==='number'?', davon '+number(call.failed_commands)+' mit Fehler':'')+'.';
  if(typeof call.tool_events==='number'&&call.tool_events>0)return 'Sol arbeitet an diesem Aufruf. '+number(call.tool_events)+' Werkzeugereignisse sind protokolliert.';
  return 'Sol arbeitet an diesem Aufruf.';
}

function isResearch(){return ui.data?.run?.kind==='original_source_research';}
function researchRound(id){return ui.data?.research?.rounds?.find(r=>r.id===id);}
const recordedNumber=v=>typeof v==='number'?number(v):'Nicht erfasst';
const researchSteps=[
  ['Ansatz aus Originalquellen','Sol liest die ursprünglichen Prinzipien und den bisherigen Stand. Der Ansatz muss die Herkunft der gemeinsamen Integrationsvorschrift erklären und eine entscheidende Frage nennen.'],
  ['Quellen lesen und rechnen','Der tatsächliche Aufruf kann Originalstellen lesen und Rechnungen in einem getrennten Arbeitsbereich ausführen. Werkzeugereignisse und beobachtete Ergebnisse bleiben einzeln nachvollziehbar.'],
  ['Unabhängige Kritik','Ein eigener Sol Aufruf prüft die tragende Identifikation, zusätzliche Annahmen und den Geltungsbereich. Seine Bewertung bleibt eine Modellinterpretation.'],
  ['Nächste Forschungsfrage','Aus dem Ergebnis und der Kritik entsteht die nächste präzise Aufgabe. Ein akzeptierter Vorschlag ist noch kein anerkannter physikalischer Beweis.']
];
function renderRunMode(data){
  const research=data.run?.kind==='original_source_research', historical=data.run?.historical===true;
  $('historyNotice').hidden=!historical;
  $('historyDescription').textContent='Dieser gespeicherte Lauf zeigt frühere Ergebnisse. Die aktuelle Originalquellenforschung verwendet OpenAI Sol 6.1.';
  $('liveRunButton').hidden=!ui.defaultRun||ui.defaultRun===data.run?.id;
  $('researchQuestion').hidden=!research;
  $('actualPhase').hidden=!research;$('pipeline').classList.toggle('research-mode',research);
  $('pageDescription').textContent=research?'Originalquellen lesen, entscheidend rechnen, unabhängig kritisieren und die nächste Forschungsfrage verfolgen.':'Gespeicherte Rechenversuche, ihre Ergebnisse und offene Originalbelege nachvollziehen.';
  $('candidateTitle').textContent=research?'Forschungsrunden':'Versuche im Vergleich';
  $('detailTitle').textContent=research?'Was diese Runde ergeben hat':'Was dieser Versuch verändert';
  const labels=research?['Forschungspaare','Modellaufrufe','Werkzeugereignisse','beendete Befehle']:['ganze Programme','Modellrunden','verschiedene Quellenrechnungen','gleiche Quellenrechnungen'];
  ['countProgramsLabel','countRoundsLabel','countNovelLabel','countRepeatsLabel'].forEach((id,i)=>$(id).textContent=labels[i]);
  const filters=research?['Alle','Mit Werkzeugen','Kritik mit Einwand','Beweis offen']:['Alle','Berechnet','Widerlegt','Noch offen'];
  document.querySelectorAll('[data-filter]').forEach((b,i)=>b.textContent=filters[i]);
  const flow=research?[['Ansatz entwickeln','Aus Originalprinzipien'],['Lesen und rechnen','Echte Quellen und Werkzeuge'],['Unabhängig kritisieren','Eigener Sol Aufruf'],['Aufgabe fortführen','Präziser nächster Schritt']]:[['Programme erzeugen','Rechenregeln kombinieren'],['Quelle ausrechnen','Vier Antworten gemeinsam'],['Gegenprüfen','Widersprüche erkennen'],['Ursprung belegen','Mit der Originalquelle verbinden']];
  document.querySelectorAll('[data-step]').forEach((b,i)=>{b.querySelector('strong').textContent=flow[i][0];b.querySelector('span:not(.node-circle):not(.node-stat)').textContent=flow[i][1];});
  $('flowHelpMethod').textContent=research?'Jede Runde beginnt bei einer konkreten Herkunftsfrage. Vorschlag, tatsächlich ausgeführte Quellenprüfung oder Rechnung und unabhängige Kritik werden getrennt gespeichert.':'Ein Versuchsprogramm ist eine ausführbare Rechenanweisung. Danach werden die vier Antworten derselben Quelle erneut berechnet.';
  $('flowHelpReview').textContent=research?'Die nächste Aufgabe folgt aus den belegten Ergebnissen und der Kritik. Ein Modellurteil ersetzt keinen Originalbeweis.':'Die Gegenprüfung sucht Fehler. Danach muss die Identität mit der ursprünglichen TFPT Quelle belegt werden.';
  if(historical){$('runStatus').dataset.historical='true';}
}
function renderResearchStep(){
  document.querySelectorAll('[data-step]').forEach(b=>b.setAttribute('aria-pressed',String(Number(b.dataset.step)===ui.step)));
  const [head,copy]=researchSteps[ui.step];$('stepExplanation').replaceChildren(el('strong','','Schritt '+(ui.step+1)+' erklärt: '+head),el('p','',copy));
}
function renderResearchHeader(data){
  const counts=data.summary?.counts||{}, proc=data.process||{}, rounds=data.research?.rounds||[];
  const question=data.research?.question,nextTask=data.research?.latest_next_task;
  $('questionText').textContent=briefResearchQuestion(question)||'Die konkrete Frage ist noch nicht protokolliert.';
  const sameNext=normalizedResearchCopy(question)===normalizedResearchCopy(nextTask);
  $('questionNext').textContent=nextTask&&!sameNext?'Nächste Aufgabe aus der Gegenprüfung: '+briefResearchQuestion(nextTask):'';
  $('questionOrigin').textContent=nextTask&&sameNext?'Aus der letzten Gegenprüfung':'';
  const exact=$('exactResearchQuestion'),wasExactOpen=exact.querySelector('details')?.open;
  const exactDisclosure=disclosure('Genaue Forschungsfrage',[el('p','',question||'Noch nicht protokolliert.')],'exact-question');exactDisclosure.open=!!wasExactOpen;exact.replaceChildren(exactDisclosure);
  const phase=String(data.activity?.id||data.active_call?.role||'').toLowerCase();
  const phaseIndex=/critic|review/.test(phase)?2:/read|calcul|tool/.test(phase)?1:/propos/.test(phase)?0:/next|task/.test(phase)?3:null;
  const phaseLabels=['Herleitung','Lesen und rechnen','Unabhängige Gegenprüfung','Nächste Forschungsaufgabe'];
  const isActive=proc.status==='running'&&phaseIndex!==null;
  $('actualPhase').textContent=isActive?'Jetzt: Schritt '+(phaseIndex+1)+' · '+phaseLabels[phaseIndex]:'Aktueller Lauf: '+(plainResearchRole(data.activity?.label)||data.process?.message||'Status wird gelesen');
  document.querySelectorAll('[data-step]').forEach((b,i)=>{
    b.dataset.active=String(isActive&&i===phaseIndex);if(isActive&&i===phaseIndex)b.setAttribute('aria-current','step');else b.removeAttribute('aria-current');
    let marker=b.querySelector('.node-live');if(!marker){marker=el('span','node-live','Jetzt aktiv');b.append(marker);}marker.hidden=!(isActive&&i===phaseIndex);
  });
  const labels=['rounds','model_calls','tool_events','completed_commands'];
  ['countPrograms','countRounds','countNovel','countRepeats'].forEach((id,i)=>$(id).textContent=recordedNumber(counts[labels[i]]));
  $('flowProposals').textContent=number(rounds.length)+' gespeicherte Ansätze';
  $('flowCalculated').textContent=recordedNumber(counts.tool_events)+' Werkzeugereignisse';
  $('flowChecks').textContent=number(counts.reviewed_rounds)+' Kritiken';
  $('flowProofs').textContent=number(rounds.filter(r=>r.next_task).length)+' Folgeaufgaben';
  $('modelLabel').textContent='Vorschlag und Kritik';
  const activeEvidence=$('activeResearchEvidence'), evidenceOpen=activeEvidence.querySelector('details')?.open;
  const activeDisclosure=disclosure('Aktueller Aufruf und Wissensstand',researchFields({original_question:data.research?.question,active_call:data.active_call,backend_receipt:data.run?.backend_receipt,knowledge_context:data.run?.knowledge_context}),'active-evidence');activeDisclosure.open=!!evidenceOpen;activeEvidence.replaceChildren(activeDisclosure);
  const actual=data.run?.model, requested=data.run?.requested_model;
  $('modelName').textContent=actual==='gpt-6.1-sol'?'OpenAI Sol 6.1 · Antwort protokolliert':requested==='gpt-6.1-sol'?'OpenAI Sol 6.1 · Aufruf angefordert':data.run?.model_label||'Modellbeleg fehlt';
  $('modelName').title='gpt-6.1-sol. Eine abgeschlossene Antwort bestätigt den Aufruf mit diesem eingestellten Modell. Sie ist keine unabhängige Serverattestation.';
  $('roundLabel').textContent='Abgeschlossene Forschungspaare';$('roundValue').textContent=recordedNumber(data.summary?.finished_rounds??counts.rounds);
  $('liveLabel').textContent=proc.status==='running'?'Aktiv in der Originalquelle':'Letzter Stand der Forschung';
  $('nowTitle').textContent=plainResearchRole(data.activity?.label||data.activity?.id)||'Sol Aufruf aktiv';
  $('nowDescription').textContent=proc.status==='running'?researchActivityCopy(data):data.summary?.explanation||'Der nächste tatsächliche Schritt wird hier protokolliert.';
  $('scaleTitle').textContent='Was der Ablauf tatsächlich prüft';
  $('scaleExplanation').textContent='Die Forschung setzt bei den Gewichten der Feldkarten νVn und DVn sowie beim markierten Transfer im Infrarotbereich an. Sol liest die ursprünglichen Regeln und bisherigen Ausschlüsse, prüft einen Ansatz und übergibt ihn an einen unabhängigen Kritiker.';
  $('scaleDetail').textContent=recordedNumber(counts.model_calls)+' Modellaufrufe · '+recordedNumber(counts.completed_commands)+' beendete Befehle · '+recordedNumber(counts.failed_commands)+' Befehle mit Fehler. Diese Zahlen messen beobachtete Arbeitsschritte und keinen Fortschrittsanteil der Theorie.';
  $('boundaryTitle').textContent='Vorschlag, Kritik und physikalischer Beweis haben ihren eigenen Status.';
  $('boundaryDescription').textContent='Die Runden zeigen Modellvorschläge, beobachtete Werkzeugergebnisse und unabhängige Kritik. Der Ablauf erkennt daraus keine fundamentale Lösung automatisch an. Offene Herkunftsfragen bleiben in jeder Runde sichtbar.';
  $('proofStatus').textContent='Beweis offen';
  $('knowledgeCoverage').textContent=data.run?.knowledge_context?'Der Aufruf erhält den gespeicherten Kontext mit Originalprinzipien, bisherigen Ergebnissen, Voraussetzungen und Ausschlüssen. Das konkrete Manifest ist in den Rundendetails einsehbar.':'Das Wissensmanifest wird noch eingelesen.';
  $('coverageText').textContent=data.coverage?.detail||'Alle gespeicherten Forschungsrunden bleiben sichtbar.';
}
function critiqueObjection(r){return /reject|fail|invalid|unsupported|objection|revise|einwand|widerspruch/i.test(r.critic_verdict||'');}
function hasTools(r){return Number(r.proposal_call?.tool_events||0)+Number(r.review_call?.tool_events||0)>0;}
function renderResearchRounds(){
  const rounds=ui.data?.research?.rounds||[];
  const visible=rounds.filter(r=>ui.filter==='all'||ui.filter==='calculated'&&hasTools(r)||ui.filter==='failed'&&critiqueObjection(r)||ui.filter==='open');
  $('candidateCount').textContent=number(visible.length)+' von '+number(rounds.length)+' Runden';
  const key=JSON.stringify([visible,ui.selected,ui.filter]);if(key===ui.listKey)return;ui.listKey=key;
  if(!visible.length){$('candidateList').replaceChildren(el('div','empty-state',rounds.length?'Für diesen Filter ist noch keine Runde gespeichert.':'Die erste Forschung mit unabhängiger Kritik läuft. Sobald beide Antworten gespeichert sind, erscheint die Runde hier.'));return;}
  $('candidateList').replaceChildren(...[...visible].reverse().map(r=>{
    const b=el('button','candidate-row'+(r.id===ui.selected?' selected':''));b.type='button';b.setAttribute('aria-pressed',String(r.id===ui.selected));
    const copy=el('div'),name=el('div','candidate-name');name.append(el('span','candidate-number',number(r.ordinal)),el('span','',researchRoundTitle(r)));copy.append(name,el('div','candidate-description',researchOutcome(r)));
    const status=el('div','candidate-row-status');status.append(el('span','status-label '+(r.status==='CRITIC_REVIEWED'?'calculated':'open'),r.status==='CRITIC_REVIEWED'?'Kritik gespeichert':'Vorgeschlagen'),el('span','status-label open','Beweis offen'));b.append(copy,status);b.addEventListener('click',()=>selectCandidate(r.id,true));return b;
  }));
}
const researchLabels={new_information:'Neue Erkenntnis',open_obligations:'Offene Nachweise',proposed_status:'Vorgeschlagener Status',critic_verdict:'Kritikerurteil',original_question:'Originalfrage',question:'Forschungsfrage',summary:'Ergebnis',title:'Ansatz',hypothesis:'Hypothese',hypotheses:'Hypothesen',result:'Ergebnis',findings:'Befunde',what_was_learned:'Was gelernt wurde',what_remains:'Was offen bleibt',remaining:'Offene Schritte',next_task:'Nächste Forschungsaufgabe',sources:'Originalquellen',citations:'Quellenstellen',assumptions:'Zusätzliche Annahmen',verdict:'Bewertung',scope:'Geltungsbereich',proof_status:'Beweisstatus',status:'Status',critique:'Kritik',commands:'Befehle',artifacts:'Artefakte',exit_code:'Rückgabewert',output:'Beobachtete Ausgabe',reason:'Begründung',evidence:'Belege',objections:'Einwände',derivation:'Herleitung',missing_derivation:'Fehlende Herleitung'};
function researchFields(value,depth=0){
  if(value===null||value===undefined)return [el('p','','Nicht aufgezeichnet.')];
  if(typeof value!=='object')return [el('p','',String(value))];
  if(depth>5)return [el('pre','technical-code',JSON.stringify(value,null,2))];
  const nodes=[];
  Object.entries(value).forEach(([key,v])=>{
    if(v===null||v===undefined||Array.isArray(v)&&!v.length)return;
    const group=el('div','research-field');group.append(el('h4','',Array.isArray(value)?'Eintrag '+(Number(key)+1):researchLabels[key]||key.replaceAll('_',' ')));
    group.append(...researchFields(v,depth+1));nodes.push(group);
  });return nodes.length?nodes:[el('p','','Keine Angaben gespeichert.')];
}
function renderResearchDetail(){
  const target=$('candidateDetail'),r=researchRound(ui.selected);
  if(!r){target.replaceChildren(el('div','empty-state','Noch keine vollständige Forschungspaarung gespeichert. Die aktuelle Frage und die tatsächliche Tätigkeit stehen oben.'));return;}
  const key=JSON.stringify(r);if(key===ui.detailKey)return;ui.detailKey=key;
  const opened=new Set([...target.querySelectorAll('details[open]')].map(d=>d.dataset.disclosure));
  const head=el('div','detail-header');head.append(el('h3','',researchRoundTitle(r)),el('span','detail-badge','Beweis offen'));
  const nodes=[head,el('p','detail-summary','Runde '+r.ordinal+' · '+(r.status==='CRITIC_REVIEWED'?'Vorschlag und unabhängige Kritik gespeichert.':'Vorschlag gespeichert.'))];
  if(r.proposed_status)nodes.push(el('p','research-status','Vorschlag: '+plainResearchStatus(r.proposed_status)));
  if(r.critic_verdict)nodes.push(el('div',critiqueObjection(r)?'failure-explanation':'change-explanation','Gegenprüfung: '+plainResearchStatus(r.critic_verdict)));
  const knownFirst=isKnownFirstRound(r);
  if(knownFirst)nodes.push(el('p','detail-summary','Laut Gegenprüfung zeigt die Rechnung eine zusätzliche Bedingung beim Zusammenfassen der Felder. Die bedingte Rechnung besteht; einfache Mittelung hat ein Problem mit der lokalen Eichsymmetrie. Welche Abbildung die Originaltheorie vorgibt, bleibt offen.'));
  const summary=knownFirst?'Beim Zusammenfassen der Felder muss auch die markierte Quelle passen. Die bereits bekannte unmarkierte Bedingung legt diese zusätzliche Zuordnung noch nicht fest.':r.proposal?.new_information||r.proposal?.summary||r.proposal?.what_was_learned||r.proposal?.result;
  if(summary)nodes.push(el('h4','','Ergebnis des Vorschlags'),...researchFields(summary));
  const reviewed=knownFirst?'Die Gegenprüfung bestätigt die bedingte Varianzrechnung. Eine einfache Mittelung verschiedener Higgsorte erhält die lokale Eichsymmetrie ohne zusätzlich begründeten gemeinsamen Rahmen nicht.':r.review?.new_information;
  if(reviewed)nodes.push(el('h4','','Neue Erkenntnis laut Gegenprüfung'),...researchFields(reviewed));
  const remains=r.review?.open_obligations||r.proposal?.open_obligations||r.review?.what_remains||r.review?.remaining||r.review?.missing_derivation||r.proposal?.what_remains;
  if(remains)nodes.push(el('h4','','Was noch offen ist'),...researchFields(remains));
  if(r.next_task)nodes.push(el('h4','','Nächste konkrete Forschungsaufgabe'),el('p','detail-summary',briefResearchQuestion(r.next_task)));
  nodes.push(disclosure('Ansatz, Hypothesen und Originalstellen',[el('h4','','Originalfrage'),el('p','',r.question),el('h4','','Originalstatus'),el('p','',r.proposed_status),...researchFields(r.proposal)],'proposal'));
  nodes.push(disclosure('Unabhängige Kritik und offene Beweisschritte',[el('p','',r.critic_verdict),...researchFields(r.review)],'review'));
  nodes.push(disclosure('Tatsächliche Werkzeugergebnisse und Artefakte',researchFields(r.actual_result),'actual'));
  nodes.push(disclosure('Aufrufbelege und verwendeter Wissensstand',[el('p','','Die Pfade benennen lokale Belege. Diese Ansicht lädt keine Originaldateien oder externen Inhalte.'),...researchFields({proposal_call:r.proposal_call,review_call:r.review_call,knowledge_context:ui.data.run.knowledge_context,physics_admission:r.physics_admission})],'receipts'));
  target.replaceChildren(...nodes);target.querySelectorAll('details').forEach(d=>{if(opened.has(d.dataset.disclosure))d.open=true;});
}

renderStep();poll();
