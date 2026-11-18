(() => {
  'use strict';

  const $ = (selector, root = document) => root.querySelector(selector);
  const view = $('#view');
  const side = $('#sidebar-content');
  const sideTitle = $('#sidebar-title');
  const sideCount = $('#sidebar-count');
  const notice = $('#notice');
  const searchInput = $('#search-input');
  const state = { topics: [], topicsLoaded: false, current: null, results: [], searchQuery: '', sourceReturn: null, serial: 0, graphSerial: 0, statusSerial: 0, statusTimer: null, activeThemeTab: 'solved', activeThemeTopicId: '', activeThemeTabUserSelected: false, statusFingerprint: '', language: localStorage.getItem('tfpt-atlas-language')==='en'?'en':'de' };
  const tx=(text)=>window.TFPT_I18N?.t(text)||text;
  const english=()=>state.language==='en';
  function localized(object,key,fallback=''){if(!object)return fallback;const suffix=english()?'_en':'_de',preferred=object[`${key}${suffix}`];return preferred!=null?preferred:object[key]??fallback;}
  const localField=(object,key)=>localized(object,key,object?.translations?.en?.[key]??'');

  const topicCatalog = [
    ['topic:origin', 'Ursprung', '◉'], ['topic:e8', 'E₈ & Codes', '⌘'], ['topic:clocks', 'Verzahnte Clocks', '◷'],
    ['topic:flavor', 'Standardmodell & Flavor', '✳'], ['topic:alpha', 'Feinstrukturkonstante', 'α'], ['topic:quantum', 'Quanten- & Randstruktur', '◌'],
    ['topic:gravity', 'Gravitation', '⌁'], ['topic:cosmology', 'Kosmologie', '◍'], ['topic:rh', 'Riemann-Hypothese', 'ζ'], ['topic:validation', 'Prüfungen & Vorhersagen', '✓'],
    ['topic:universalraum', 'Universalraum', '◇'], ['topic:metrology', 'Metrologie', '∿'], ['topic:spacetime', 'Raumzeit', '⊞'], ['topic:process', 'Prozess & Kompilation', '⇢']
  ];
  const topicNamesEn = {
    'topic:origin':'Origin','topic:e8':'E₈ & codes','topic:clocks':'Interlocking clocks','topic:flavor':'Standard Model & flavor',
    'topic:alpha':'Fine-structure constant','topic:quantum':'Quantum & boundary structure','topic:gravity':'Gravity','topic:cosmology':'Cosmology',
    'topic:rh':'Riemann hypothesis','topic:validation':'Tests & predictions','topic:universalraum':'Universal space','topic:metrology':'Metrology','topic:spacetime':'Spacetime','topic:process':'Process & compilation'
  };
  const topicCopy = {
    'topic:origin': 'Die Ausgangsprinzipien und die Frage, was die Theorie eindeutig festlegt.',
    'topic:e8': 'Endliche Codestrukturen, E₈ und ihre möglichen Beziehungen im Theoriegraphen.',
    'topic:clocks': 'Wie ineinander verschränkte Uhren und ihre Relationen beschrieben werden.',
    'topic:flavor': 'Ladungen, Teilchenfamilien und die Einordnung des Standardmodells.',
    'topic:alpha': 'Herleitungen, Werte und Prüfungen rund um die Feinstrukturkonstante.',
    'topic:quantum': 'Quantenstruktur, Randbedingungen und die Auswahl physikalischer Zustände.',
    'topic:gravity': 'Aussagen und offene Herleitungen zu Raumzeit und Gravitation.',
    'topic:cosmology': 'Kosmologische Zusammenhänge, Konsequenzen und beobachtbare Fragen.',
    'topic:rh': 'Verbindungen zur Riemann-Hypothese und der zugehörigen Zahlentheorie.',
    'topic:validation': 'Prüfprogramme, Vorhersagen, Skripte und formale Lean-Verknüpfungen.'
  };
  const icons = { Claim: 'A', Result: 'R', Derivation: 'D', Script: '⌘', Lean: 'L', Source: '§', Topic: '◉', Prediction: '↗', Observation: '◌' };

  const formatControl = /[\u0000-\u0008\u000B\u000C\u000E-\u001F\u007F]/g;
  function cleanText(value) { return String(value ?? '').replace(formatControl, '[Zeichenfehler]').replace(/\\u007f/gi,'[Zeichenfehler]'); }
  function containsFormatControl(value) {
    if(typeof value==='string')return /[\u0000-\u0008\u000B\u000C\u000E-\u001F]/.test(value);
    if(Array.isArray(value))return value.some(containsFormatControl);
    if(value&&typeof value==='object')return Object.values(value).some(containsFormatControl);
    return false;
  }
  function esc(value) { return cleanText(value).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c])); }
  function asArray(value) { return Array.isArray(value) ? value : value == null ? [] : [value]; }
  function attrs(node) { return node?.attrs || node?.attributes || {}; }
  function semantic(node) { return attrs(node).semantic || {}; }
  function isClaimNode(node) { return /claim/i.test(String(node?.type || node?.kind || attrs(node).type || '')); }
  function freshEnglish(node,key) {
    const a=attrs(node),s=semantic(node),tr=a.translations?.en||{};
    // Prefer English text authored alongside the current semantic record. The
    // service only includes its translation overlay when its fingerprint
    // matches, so it remains a safe fallback when no authored field is present.
    return s[`${key}_en`]||a[`${key}_en`]||node?.[`${key}_en`]||tr[key]||'';
  }
  function assessment(node) { return attrs(node).assessment || node?.assessment || {}; }
  function ledgerStatus(node) { const a=attrs(node),l=a.ledger||{};return a.status||a.verdict||a.ledger_status||l.status||l.verdict||l.canonical_status||l.claim_status||node?.ledger_status||node?.status||node?.verdict||''; }
  function typeName(node) { return node?.type || node?.kind || attrs(node).type || 'Eintrag'; }
  function typeLabel(node) {
    const raw=String(typeof node==='string'?node:typeName(node)),key=raw.toLowerCase().replace(/[^a-z]/g,'');
    if(key.includes('claim'))return tx('Aussage');
    if(key.includes('result'))return tx('Ergebnis');
    if(key.includes('derivation'))return tx('Herleitung');
    if(key.includes('script'))return tx('Prüfprogramm');
    if(key.includes('theorypattern')||key==='pattern')return tx('Forschungsmuster');
    if(key.includes('lean')||key.includes('theorem'))return tx('Lean-Satz');
    if(key.includes('theoryaspect')||key.includes('aspect'))return tx('Themenaspekt');
    if(key.includes('object'))return tx('Objekt');
    if(key.includes('check'))return tx('Prüfung');
    if(key.includes('premise')||key.includes('assumption'))return tx('Voraussetzung');
    if(key.includes('topic'))return tx('Thema');
    if(key.includes('structure'))return tx('Struktur');
    if(key.includes('goal'))return tx('Ziel');
    if(key.includes('attempt')||key.includes('contract')||attrs(node).original_contract)return tx('Versuch / Befund');
    if(key.includes('source'))return tx('Quelle');
    if(key.includes('input'))return tx('Eingabe');
    return raw;
  }
  function label(node) { const s=semantic(node),a=attrs(node),tr=a.translations?.en||{},englishTitle=isClaimNode(node)?freshEnglish(node,'title'):(node?.label_en||a.label_en||s.title_en||a.title_en||tr.title),candidate=(english()&&(englishTitle||tr.title))||(node?.label||node?.title||a.label||a.title||s.title_de||a.statement||node?.id||'Ohne Titel');const id=String(node?.id||'').replace(/^claim:/i,'');return String(candidate)===id&&(s.title_en||s.title_de||tr.title)?(english()?englishTitle||tr.title||s.title_de:s.title_de||s.title_en||tr.title):String(candidate); }
  function summary(node) { const a=attrs(node),s=semantic(node),tr=a.translations?.en||{},englishOverview=isClaimNode(node)?freshEnglish(node,'overview'):(node?.summary_en||node?.overview_en||a.summary_en||s.overview_en||a.overview_en||a.interpretation_en||tr.overview);return english()?(englishOverview||node?.summary||node?.overview||a.summary||a.description||s.overview_de||a.interpretation_de||''):(node?.summary||node?.overview||a.summary||a.description||s.overview_de||a.interpretation_de||''); }
  function excerpt(value,limit=280){const text=cleanText(value).replace(/\s+/g,' ').trim();return text.length>limit?`${text.slice(0,limit-1).trimEnd()}…`:text;}
  function wrapGraphText(value,max=32,limit=3){const words=String(value??'').split(/\s+/).filter(Boolean),out=[];let line='',cut=false;for(const word of words){if(word.length>max){if(line){out.push(line);line='';}for(let i=0;i<word.length;i+=max){if(out.length>=limit){cut=true;break;}out.push(word.slice(i,i+max));}if(cut)break;continue;}const next=line?`${line} ${word}`:word;if(next.length>max){out.push(line);line=word;if(out.length>=limit){cut=true;break;}}else line=next;}if(line&&out.length<limit)out.push(line);if(cut&&out.length)out[out.length-1]=`${out[out.length-1].slice(0,max-1)}…`;return out;}
  function categoryLabel(value){const c=String(value||'').toLowerCase();if(/methodological[_ -]?rule|prüfregel/.test(c))return tx('Geltende Prüfregel');if(/established|solved|proved|recorded_established/.test(c))return tx('belegt im Rahmen');if(/conditional|bedingt/.test(c))return tx('bedingt');if(/open|offen/.test(c))return tx('offen');if(/mixed|gemischt/.test(c))return tx('gemischt');if(/excluded|ausgeschlossen/.test(c))return tx('ausgeschlossen');if(/assumed|axiom|angenommen/.test(c))return tx('angenommen');if(/historical|historisch/.test(c))return tx('historisch');if(/unknown|unbekannt/.test(c))return tx('Status nicht automatisch zugeordnet');return '';}
  function statusClassFor(category,authority){const c=String(category||'').toLowerCase();let status='status-unknown';if(/methodological[_ -]?rule|prüfregel/.test(c))status='status-methodological';else if(/established|solved|proved/.test(c))status='status-established';else if(/conditional/.test(c))status='status-conditional';else if(/open/.test(c))status='status-open';else if(/mixed/.test(c))status='status-mixed';else if(/excluded/.test(c))status='status-excluded';else if(/assumed/.test(c))status='status-assumed';return authority==='ai_unreviewed'?`${status} status-unreviewed`:status;}
  function isUnreviewed(node){return assessment(node).authority==='ai_unreviewed'||attrs(node).provenance?.reviewed===false||node?.provenance?.reviewed===false;}
  function edgeLabel(value){return ({requires:'benötigt als Voraussetzung',required_together:'gemeinsam erforderlich',derives:'leitet her',has_result:'hat Ergebnis',has_derivation:'hat Herleitung',implements:'implementiert',checks:'prüft',formalization_candidate:'Formalisierungskandidat',uses_claim_under_scope:'verwendet im Geltungsbereich',sets_input:'setzt Eingabe',interprets_as:'ordnet als Eingabe ein'})[String(value||'').toLowerCase()]||String(value||'').replaceAll('_',' ');}
  function enrichment(node){
    const a=attrs(node),s=semantic(node),state=String(a.enrichment_state||s.enrichment_state||a.semantic_state||node?.enrichment_state||'').toLowerCase();
    if(state==='error')return {label:tx('Auswertung fehlgeschlagen'),className:'enrichment-error'};
    if(state==='pending')return {label:tx('Auswertung ausstehend'),className:'enrichment-pending'};
    if(state==='enriched'||s.overview_de||node?.overview||node?.overview_en)return {label:tx('KI-Auswertung · nicht fachlich abgenommen'),className:'ai'};
    if(a.semantic)return {label:tx('Auswertung ausstehend'),className:'enrichment-pending'};
    return null;
  }
  function sourceId(source) { return typeof source === 'string' ? source : source?.id || source?.source_id || source?.sourceId || ''; }
  function readableValue(value){
    if(value==null)return '';
    if(typeof value!=='object')return String(value);
    if(Array.isArray(value))return value.map(readableValue).filter(Boolean).join(' · ');
    const preferred=['text','label','statement','description','summary','value','condition','step','requirement','name'];
    if(english())for(const key of preferred)if(value[`${key}_en`]!=null)return readableValue(value[`${key}_en`]);
    for(const key of preferred)if(value[key]!=null)return readableValue(value[key]);
    return Object.entries(value).map(([key,val])=>`${key.replaceAll('_',' ')}: ${readableValue(val)}`).filter(Boolean).join(' · ');
  }
  function readableText(value){const text=readableValue(value).trim();return text||'';}
  function textBlock(value) {
    if (value == null || value === '') return '';
    if (Array.isArray(value)) return `<ul>${value.map(item => `<li data-original-content>${esc(readableValue(item))}</li>`).join('')}</ul>`;
    if (typeof value === 'object') return `<p data-original-content>${esc(readableValue(value))}</p>`;
    return `<p data-original-content>${esc(value)}</p>`;
  }
  function objectCards(items, kind, allowEnglish=true) {
    const entries=asArray(items); if(!entries.length)return '';
    const local=(obj,key,fallback='')=>english()&&allowEnglish?localized(obj,key,fallback):(obj?.[key]??fallback);
    const field=(obj,key)=>asArray(english()&&allowEnglish?(obj?.[`${key}_en`]??obj?.[key]):obj?.[key]).map(x=>{if(typeof x==='string')return x;if(!x||typeof x!=='object')return String(x??'');const preferred=english()&&allowEnglish?(x.label_en||x.text_en):'';return preferred||x.label||x.text||x.id||JSON.stringify(x);}).filter(Boolean);
    const renderLines=(title,values)=>values.length?`<div class="entry-detail"><b>${title}</b><span>${values.map(esc).join(' · ')}</span></div>`:'';
    const refs=(obj)=>{const refs=asArray(english()?(obj?.source_refs_en??obj?.source_refs):obj?.source_refs).filter(Boolean);return refs.length?`<div class="entry-sources"><b>Quellen</b>${refs.map(ref=>{const id=sourceId(ref);return id?`<button type="button" class="inline-source" data-source="${esc(id)}" data-line-start="${esc(ref?.line_start??'')}" data-line-end="${esc(ref?.line_end??'')}">${esc(typeof ref==='string'?ref:sourceLabel(ref))}</button>`:''}).join('')}</div>`:'';};
    return `<div class="entry-cards">${entries.map(item=>{
      if(typeof item==='string')return `<article class="entry-card"><p>${esc(item)}</p></article>`;
      if(kind==='result')return `<article class="entry-card"><h4>${esc(local(item,'label',item.id||tx('Ergebnis')))}</h4>${local(item,'statement')?`<p class="entry-statement">${esc(local(item,'statement'))}</p>`:''}${local(item,'interpretation')?`<p>${esc(local(item,'interpretation'))}</p>`:''}${renderLines('Einordnung',field(item,'epistemic'))}${renderLines('Geltungsbereich',field(item,'scope'))}${renderLines('Begrenzungen',field(item,'limitations'))}${refs(item)}</article>`;
      if(kind==='premise')return `<article class="entry-card compact"><h4>${esc(local(item,'kind',tx('Voraussetzung')))}</h4><p>${esc(local(item,'text',local(item,'label',item.id||'')))}</p>${renderLines('Herkunft',field(item,'origin_claim_id'))}${refs(item)}</article>`;
      if(kind==='derivation')return `<article class="entry-card"><h4>${esc(local(item,'label',item.id||tx('Herleitung')))}</h4>${local(item,'argument')?`<p>${esc(local(item,'argument'))}</p>`:''}${renderLines('Voraussetzungen',field(item,'premise_ids'))}${renderLines('Vorgelagerte Aussagen',field(item,'upstream_claim_ids'))}${renderLines('Ergebnisse',field(item,'result_ids'))}${renderLines('Geltungsbereich',field(item,'scope'))}${refs(item)}</article>`;
      if(kind==='script')return `<article class="entry-card"><h4>${esc(local(item,'purpose',item.purpose_de||item.script_id||tx('Prüfprogramm')))}</h4>${renderLines('Eingaben',field(item,'inputs'))}${renderLines('Schritte',field(item,'operations'))}${Array.isArray(item.checks)&&item.checks.length?`<div class="script-checks"><b>Prüfungen</b>${item.checks.map(check=>`<div><strong>${esc(local(check,'name','Prüfung'))}</strong><span>${esc(check.relation||'')}${check.method?` · ${esc(local(check,'method'))}`:''}${check.scope?` · ${esc(local(check,'scope'))}`:''}</span></div>`).join('')}</div>`:''}${renderLines('Ausgaben',field(item,'outputs'))}${renderLines('Begrenzungen',field(item,'limitations'))}${refs(item)}</article>`;
      if(kind==='lean')return `<article class="entry-card compact"><h4>${esc(item.theorem_id||item.id||'Lean-Verknüpfung')} <span class="relation-tag">${esc(item.relation||'')}</span></h4>${item.explanation?`<p>${esc(item.explanation)}</p>`:''}${renderLines('Voraussetzungen',field(item,'premises'))}${renderLines('Bezug zu Ergebnissen',field(item,'result_ids'))}${refs(item)}</article>`;
      return `<article class="entry-card"><h4>${esc(item.label||item.id||'Eintrag')}</h4><p>${esc(item.text||item.statement||item.argument||item.explanation||JSON.stringify(item))}</p></article>`;
    }).join('')}</div>`;
  }
  function structuredSection(title,icon,items,kind,allowEnglish=true){return asArray(items).length?`<section class="dossier-section"><h3><span class="section-mark">${icon}</span>${title}</h3>${objectCards(items,kind,allowEnglish)}</section>`:'';}
  function attemptDossierSection(node){
    if(!/attempt|contract/i.test(typeName(node))&&!attrs(node).original_contract)return '';
    const a=attrs(node),c=a.original_contract||semantic(node).original_contract||{},pick=(key,...fallback)=>{if(english()&&c[`${key}_en`]!=null)return c[`${key}_en`];if(english()&&a[`${key}_en`]!=null)return a[`${key}_en`];return c[key]??fallback.find(value=>value!=null);},t8=pick('no_t1_t8_closure'),fields=[['Fragestellung',pick('question',a.question)],['Erfolgskriterium',pick('success_criterion',a.success_criterion)],['Ergebnis',pick('result',c.outcome,a.result,a.outcome)],['Status',pick('verdict',c.status,c.verdict_raw,a.verdict,a.verdict_raw,a.status)],['Verbleibende Lücke',pick('remaining_gap',c.open_gap,a.remaining_gap)],['Bedingte Ergebnisse',pick('analytic_conditional_results',a.analytic_conditional_results)],['Ausschlusskriterien / Kill-Kriterien',pick('kills',c.kill_criteria,a.kills)],['T1–T8-Ziele',t8===true?tx('T1–T8 bleiben offen; kein Gesamtabschluss.'):t8===false?tx('T1–T8-Gesamtabschluss im Vertrag nicht als offen markiert.'):t8]];
    const rows=fields.filter(([,v])=>v!=null&&v!=='').map(([name,value])=>`<div class="attempt-field"><b>${esc(name)}</b><span>${esc(readableValue(value))}</span></div>`).join('');
    const shown=new Set(['question','success_criterion','result','outcome','verdict','status','remaining_gap','open_gap','analytic_conditional_results','kills','kill_criteria','no_t1_t8_closure']);
    const rest=Object.entries(c).filter(([key])=>!shown.has(key)&&!['sources','source_refs','id','claim_id'].includes(key));
    return `<section class="dossier-section attempt-dossier"><h3><span class="section-mark">⌕</span>Versuch / Befund</h3>${rows?`<div class="attempt-fields">${rows}</div>`:'<p>Für diesen Versuch sind noch keine Vertragsfelder eingetragen.</p>'}${rest.length?`<details><summary class="tiny-label">Weitere Vertragsangaben</summary><div class="attempt-fields">${rest.map(([key,value])=>`<div class="attempt-field"><b>${esc(key.replaceAll('_',' '))}</b><span>${esc(readableValue(value))}</span></div>`).join('')}</div></details>`:''}</section>`;
  }
  function claimAttemptsSection(response,node){
    if(!response)return '';
    const statusLabel=value=>{const raw=String(value||'').trim(),normalized=raw.toUpperCase();if(normalized==='PARTIAL')return tx('Teilweise');if(normalized==='PASS')return tx('Bestanden');if(normalized==='FAIL')return tx('Fehlgeschlagen');return raw;};
    const entries=asArray(response.attempts),total=Number(response.total||entries.length),offset=Number(response.offset||0),cards=entries.map(attempt=>{
      const a=attrs(attempt),c=a.original_contract||{},question=localField(c,'question')||localField(a,'question')||summary(attempt)||label(attempt),criterion=localField(c,'success_criterion')||localField(a,'success_criterion'),method=localField(c,'method')||localField(c,'approach')||localField(c,'procedure')||localField(a,'method'),outcome=localField(c,'result')||localField(c,'outcome')||localField(a,'result')||localField(a,'outcome'),rawStatus=a.verdict_raw||a.verdict||c.verdict||c.verdict_raw||c.status||a.status||a.assessment?.label_en||a.assessment?.label,status=statusLabel(rawStatus)||tx('Befund nicht klassifiziert'),scope=localField(c,'scope')||localField(c,'evidence_class')||localField(a,'scope'),findings=c.resolved||c.findings||a.findings,open=c.open||c.remaining_gap||a.remaining_gap,gap=localField(c,'remaining_gap')||localField(a,'remaining_gap'),questionText=readableText(question)||label(attempt),criterionText=readableText(criterion),methodText=readableText(method),outcomeText=readableText(outcome),scopeText=readableText(scope),findingsText=readableText(findings),openText=readableText(open),gapText=readableText(gap);
      const refs=[...asArray(a.source_refs),...asArray(c.source_refs),...asArray(a.sources)].filter(Boolean).filter((ref,index,all)=>all.findIndex(other=>sourceId(other)===sourceId(ref))===index);
      return `<article class="entry-card attempt-card"><div class="attempt-card-head"><h4>${esc(label(attempt))}</h4><span class="badge">${esc(status)}</span></div><p><b>${esc(tx('Fragestellung'))}:</b> ${esc(questionText)}</p>${criterionText?`<p><b>${esc(tx('Erfolgskriterium'))}:</b> ${esc(criterionText)}</p>`:''}${methodText?`<p><b>${esc(tx('Vorgehen'))}:</b> ${esc(methodText)}</p>`:''}${outcomeText?`<p><b>${esc(tx('Ergebnis'))}:</b> ${esc(outcomeText)}</p>`:''}${scopeText?`<p><b>${esc(tx('Geltungsbereich'))}:</b> ${esc(excerpt(scopeText,480))}</p>`:''}${findingsText?`<p><b>${esc(tx('Befunde'))}:</b> ${esc(excerpt(findingsText,520))}</p>`:''}${openText?`<p><b>${esc(tx('Offene Punkte'))}:</b> ${esc(excerpt(openText,420))}</p>`:''}${gapText&&gapText!==openText?`<p><b>${esc(tx('Offene Lücke'))}:</b> ${esc(gapText)}</p>`:''}<button type="button" class="inline-node" data-node="${esc(attempt.id)}">${esc(tx('Versuchsdossier öffnen'))}</button>${refs.length?`<div class="entry-sources"><b>${esc(tx('Originalquellen'))}</b>${sourceButtons(refs)}</div>`:''}</article>`;
    }).join('');
    return `<section class="dossier-section claim-attempts"><h3><span class="section-mark">⌕</span>${esc(tx('Versuche und Befunde'))} <span class="muted">${total}</span></h3><p class="gap-intro">${esc(tx('Dokumentierte Versuche zur Aussage; Einträge zählen nicht automatisch als unabhängige Experimente.'))}</p>${cards?`<div class="entry-cards">${cards}</div>`:`<p>${esc(tx('Noch keine explizit verknüpften Versuche oder Befunde.'))}</p>`}${offset+entries.length<total?`<button type="button" class="theme-more" data-more-attempts>${esc(tx('Weitere Versuche laden'))}</button>`:''}</section>`;
  }
  const inputOrigin = {
    upstream_result: ['vorgelagerte Quelle','upstream','Graph-Verknüpfung zu einer vorgelagerten Quelle. Das Label bestätigt für sich allein weder Beweisstatus noch fachliche Abnahme.'],
    declared_axiom: ['eingesetzt · Axiom','set','Als Axiom oder deklarierte Voraussetzung eingesetzt.'],
    external_data: ['extern','external','Wert stammt laut Audit aus externen Daten.'],
    hardcoded: ['eingesetzt','set','Wert ist im betrachteten Aufrufpfad fest eingetragen.'],
    free_parameter: ['eingesetzt · freier Parameter','set','Wert wird als freier Parameter vorgegeben.'],
    derived_in_script: ['im Skript abgeleitet','derived','Der Audit ordnet den Wert einer Ableitung innerhalb eines Skripts zu; das allein zeigt nicht, dass die Ableitung aus den Theorieprinzipien folgt.'],
    unknown: ['Herkunft offen','unknown','Die Herkunft ist in den verfügbaren Angaben nicht geklärt.']
  };
  const inputImpact = {
    benign: ['kein markierter Ableitungsmangel','benign'],
    conditional: ['bedingt','conditional'],
    circularity_risk: ['Zirkularitätsrisiko','risk'],
    missing_derivation: ['Herleitung fehlt','missing']
  };
  const selectionStatusLabels={forced_by_supplied_derivation:'durch angegebene Herleitung festgelegt',assumed_or_selected:'angenommen oder ausgewählt',not_applicable:'nicht anwendbar',unknown:'unklar'};
  function sourceButtons(refs){return asArray(refs).filter(Boolean).map((ref,index)=>{const id=sourceId(ref),text=typeof ref==='string'?`${tx('Quelle')} ${index+1}`:sourceLabel(ref);return id?`<button type="button" class="inline-source" data-source="${esc(id)}" data-line-start="${esc(ref?.line_start??'')}" data-line-end="${esc(ref?.line_end??'')}" title="${esc(typeof ref==='string'?id:sourceLabel(ref))}">${esc(text)}</button>`:''}).join('');}
  function researchContextPanel(id){
    const details=document.createElement('details');details.className='research-context';
    const summaryEl=document.createElement('summary');summaryEl.textContent=tx('Forschungsstand · vorhandene Evidenz und nächste Frage');
    const body=document.createElement('div');body.className='research-context-body';
    details.append(summaryEl,body);
    let loaded=false;
    const routeCurrent=()=>details.isConnected&&(isCurrentTopicRoute(id)||isCurrentNodeRoute(id));
    details.addEventListener('toggle',async()=>{
      if(!details.open||loaded)return;loaded=true;body.innerHTML=`<p class="loading-inline">${esc(tx('Forschungsstand wird zusammengestellt …'))}</p>`;
      try{
        const packet=await api(`/api/context?id=${encodeURIComponent(id)}&max_chars=12000`);
        if(!routeCurrent())return;
        const section=(title,items,render=value=>`<p>${esc(readableValue(value))}</p>`)=>{
          const rows=asArray(items);return rows.length?`<section><h4>${esc(title)}</h4><div class="context-items">${rows.map(render).join('')}</div></section>`:'';
        };
        const attempts=asArray(packet.prior_attempts?.items),seeds=asArray(packet.recorded_shared_inputs?.seeds),objects=asArray(packet.recorded_shared_inputs?.objects),connections=asArray(packet.recorded_connections),negatives=asArray(packet.scoped_negatives),candidates=asArray(packet.candidate_connections),discriminators=asArray(packet.next_discriminators),assumptions=asArray(packet.assumption_records);
        const attemptCard=item=>{const conclusions=item.recorded_conclusions||{},question=item.question_or_scope||item.label||item.id,status=item.verdict,conclusionRows=[['Befunde',conclusions.findings],['Bedingte Ergebnisse',conclusions.conditional_results],['Ausgeschlossen',conclusions.rules_out],['Offen geblieben',conclusions.missing_or_open]].map(([title,values])=>{const clean=asArray(values).map(readableValue).filter(value=>value&&value!=='{}');return clean.length?`<small><b>${esc(tx(title))}:</b> ${esc(clean.join(' · '))}</small>`:''}).join('');return `<article class="context-card"><strong>${esc(question||tx('Versuch / Befund'))}</strong>${status?`<span class="context-status">${esc(status)}</span>`:''}${conclusionRows}${item.notice?`<small>${esc(tx(item.notice))}</small>`:''}${sourceButtons(item.source_ids)}</article>`;};
        const relationCard=item=>`<article class="context-card"><strong>${esc(tx(edgeLabel(item.edge_type||item.relation||tx('Verknüpfung'))))}</strong><p>${esc(item.src||'')} → ${esc(item.dst||'')}</p><small>${esc(item.src_type||'')} → ${esc(item.dst_type||'')} · ${esc(tx('gespeicherte Beziehung, kein Beweis'))}</small></article>`;
        const negativeCard=item=>{const rules=asArray(item.rules_out).map(readableValue).filter(Boolean);return `<article class="context-card"><strong>${esc(item.id||tx('Negativer Befund'))} · ${esc(item.verdict||'')}</strong>${item.scope?`<p><b>${esc(tx('Geltungsbereich'))}:</b> ${esc(item.scope)}</p>`:''}${rules.length?`<small>${esc(tx('Ausgeschlossen'))}: ${esc(rules.join(' · '))}</small>`:''}${item.notice?`<small>${esc(tx(item.notice))}</small>`:''}${sourceButtons(item.source_ids)}</article>`;};
        const assumptionCard=item=>{const origin=inputOrigin[item.origin]||inputOrigin.unknown,impact=inputImpact[item.impact];return `<article class="context-card flagged"><strong>${esc(item.name||tx('Eingabe'))} · ${esc(item.claim_id||'')}</strong><span class="context-status">${esc(tx(origin[0]))}${impact?` · ${esc(tx(impact[0]))}`:''}</span>${item.used_for?`<small>${esc(tx('Verwendet für'))}: ${esc(item.used_for)}</small>`:''}${sourceButtons(item.source_ids)}</article>`;};
        const candidateCard=item=>{const candidateLabel=english()?(item.label_en||item.label):(item.label||item.label_en),candidateSummary=english()?(item.summary_en||item.summary):(item.summary||item.summary_en);return `<article class="context-card ${item.reviewed===false?'unreviewed':''}"><strong>${esc(candidateLabel||item.id||tx('Kandidat'))}</strong><span class="context-status">${esc(item.reviewed===false?tx('Kandidat · ungeprüft'):tx('Mögliche Verbindung · kein Beweis'))}</span>${candidateSummary?`<p>${esc(candidateSummary)}</p>`:''}${asArray(item.conditions).length?`<small>${esc(tx('Bedingungen'))}: ${esc(asArray(item.conditions).map(readableValue).join(' · '))}</small>`:''}${asArray(item.counterevidence).length?`<small>${esc(tx('Gegenbefunde'))}: ${esc(asArray(item.counterevidence).map(readableValue).join(' · '))}</small>`:''}${sourceButtons(item.source_ids)}</article>`;};
        const sharedCard=(item,key)=>`<article class="context-card"><strong>${esc(key==='seeds'?item.seed_id:key==='objects'?item.object_id:key==='structures'?item.structure_id:item.id||item.label||tx('Objekt'))}${item.label?` · ${esc(item.label)}`:''}</strong><small>${esc((item.claim_ids||[]).join(' · '))}</small>${item.relation?`<small>${esc(item.relation.replaceAll('_',' '))} · ${esc(tx('keine unabhängige Bestätigung'))}</small>`:''}${sourceButtons(item.source_ids)}</article>`;
        const sourceIds=asArray(packet.source_ids);
        const contextQuestion=english()&&topicNamesEn[id]?topicNamesEn[id]:packet.question;
        body.innerHTML=`<p class="context-notice">${esc(tx('Diese Übersicht zeigt gespeicherte Zusammenhänge und frühere Befunde; sie fügt keinen Beweis hinzu.'))}</p>${english()?`<p class="language-fallback">${esc(tx('Detailtexte bleiben im erfassten Originalwortlaut, wenn keine englische Übersetzung vorliegt.'))}</p>`:''}${contextQuestion?`<section class="context-question"><h4>${esc(tx('Leitfrage'))}</h4><p>${esc(contextQuestion)}</p></section>`:''}${packet.scope?`<p class="context-scope"><b>${esc(tx('Geltungsbereich'))}:</b> ${esc(packet.scope)}</p>`:`<p class="context-scope"><b>${esc(tx('Geltungsbereich'))}:</b> ${esc(tx('Geltungsbereich nicht angegeben'))}</p>`}${packet.coverage?.claims_total!=null?`<p class="context-notice">${esc(tx('Verknüpfte Aussagen im Kontext'))}: ${Number(packet.coverage.claims_in_page||0)}/${Number(packet.coverage.claims_total||0)}${packet.coverage.claims_total>packet.coverage.claims_in_page?` · ${esc(tx('Gemeinsame Eingaben werden nur im gezeigten Ausschnitt verglichen.'))}`:''}</p>`:''}${packet.what_we_have?section(tx('Bisher festgehalten'),packet.what_we_have):''}${section(tx('Nächste trennscharfe Fragen'),discriminators)}${section(tx('Versuche und Befunde'),attempts,attemptCard)}${packet.prior_attempts?.total?`<p class="context-notice">${attempts.length}/${packet.prior_attempts.total} ${tx('Versuche in der Auswahl')}</p>`:''}${section(tx('Eingaben mit offenen Herkunfts- oder Wirkungsfragen'),assumptions,assumptionCard)}${section(tx('Gemeinsam genutzte Eingaben · gleiche gespeicherte IDs'),[...seeds.map(item=>({...item,_kind:'seeds'})),...objects.map(item=>({...item,_kind:'objects'})),...asArray(packet.recorded_shared_inputs?.structures).map(item=>({...item,_kind:'structures'}))],item=>sharedCard(item,item._kind))}${packet.recorded_shared_inputs?.identity_notice?`<p class="context-notice">${esc(tx(packet.recorded_shared_inputs.identity_notice))}</p>`:''}${section(tx('Gespeicherte Beziehungen · keine Beweiskette'),connections,relationCard)}${section(tx('Begrenzte negative Befunde'),negatives,negativeCard)}${section(tx('Ungeprüfte Anschlusskandidaten'),candidates,candidateCard)}${sourceIds.length?`<section><h4>${esc(tx('Quellen für den Einstieg'))}</h4><div class="context-items">${sourceButtons(sourceIds)}</div></section>`:''}${packet.truncated?`<p class="context-truncated">${esc(tx('Der Forschungsrahmen ist begrenzt; weitere Einträge können über den vollständigen Wissensgraphen geöffnet werden.'))}</p>`:''}`;
        if(packet.prior_attempts?.truncated){const attemptsSection=[...body.querySelectorAll('section')].find(node=>node.querySelector('h4')?.textContent===tx('Versuche und Befunde'));if(attemptsSection){const more=document.createElement('button');more.type='button';more.className='context-more';more.textContent=`${tx('Weitere Versuche laden')} · ${packet.prior_attempts.total}`;more.addEventListener('click',async()=>{more.disabled=true;more.textContent=tx('Wird geladen …');try{const offset=Number(packet.prior_attempts.next_offset||0),extra=await api(`/api/context?id=${encodeURIComponent(id)}&max_chars=12000&offset=${offset}`);if(!routeCurrent())return;const grid=attemptsSection.querySelector('.context-items');for(const item of asArray(extra.prior_attempts?.items)){grid?.insertAdjacentHTML('beforeend',attemptCard(item));}attemptsSection.querySelectorAll('[data-source]').forEach(button=>{if(button.dataset.contextBound)return;button.dataset.contextBound='true';button.addEventListener('click',()=>loadSource(button.dataset.source,button,{line_start:button.dataset.lineStart?Number(button.dataset.lineStart):undefined,line_end:button.dataset.lineEnd?Number(button.dataset.lineEnd):undefined}));});if(extra.prior_attempts?.truncated){packet.prior_attempts=extra.prior_attempts;more.disabled=false;more.textContent=`${tx('Weitere Versuche laden')} · ${extra.prior_attempts.total}`;}else more.remove();}catch{more.disabled=false;more.textContent=tx('Laden fehlgeschlagen · erneut versuchen');}});attemptsSection.append(more);}}
        body.querySelectorAll('[data-source]').forEach(button=>{if(button.dataset.contextBound)return;button.dataset.contextBound='true';button.addEventListener('click',()=>loadSource(button.dataset.source,button,{line_start:button.dataset.lineStart?Number(button.dataset.lineStart):undefined,line_end:button.dataset.lineEnd?Number(button.dataset.lineEnd):undefined}));});
      }catch(error){if(!routeCurrent())return;loaded=false;body.innerHTML=`<p class="error-inline">${esc(error.message||tx('Forschungsstand konnte nicht geladen werden.'))} <button type="button" class="mini-button" data-context-retry>${esc(tx('Erneut laden'))}</button></p>`;body.querySelector('[data-context-retry]')?.addEventListener('click',()=>{details.open=false;setTimeout(()=>{details.open=true;},0);});}
    });
    return details;
  }
  function inputAuditSection(items,contextNode){
    const entries=asArray(items);if(!entries.length)return '';
    const attention=entries.filter(item=>item.origin==='unknown'||['circularity_risk','missing_derivation'].includes(item.impact)).length;
    return `<section class="dossier-section input-audit ${attention?'has-attention':''}"><div class="audit-intro"><div><h3><span class="section-mark">!</span>Woher kommen die Eingaben?</h3><p>Herkunft und Einsatz jeder erfassten Eingabe. Diese Einordnung beschreibt den Datenpfad; sie ist kein automatischer Beweis.</p></div>${attention?`<span class="audit-count">${attention} Prüfhinweis${attention===1?'':'e'}</span>`:''}</div><div class="entry-cards">${entries.map(item=>{
      const origin=inputOrigin[item.origin]||inputOrigin.unknown, impact=inputImpact[item.impact];
      const value=item.value==null?tx('(Wert nicht angegeben)'):typeof item.value==='string'?item.value:JSON.stringify(item.value);
      const claimId=normalizeGraphId(item.upstream_claim_id,'claim',contextNode);
      const selection=item.selection_status_en||(english()?selectionStatusLabels[item.selection_status]&&tx(selectionStatusLabels[item.selection_status]):selectionStatusLabels[item.selection_status])||item.selection_status, dependencies=item.dependency_names_en||item.dependency_names;
      return `<article class="entry-card input-card ${impact?.[1]==='risk'||impact?.[1]==='missing'?'flagged':''}"><div class="input-card-title"><div><h4>${esc(item.name||tx('Eingabe'))}</h4><code>${esc(value)}</code></div><span class="origin-badge ${origin[1]}">${esc(tx(origin[0]))}</span></div><p class="origin-explanation">${esc(english()?(item.explanation_en||item.explanation||origin[2]):(item.explanation||origin[2]))}</p>${selection?`<div class="entry-detail"><b>${esc(tx('Auswahlstatus'))}</b><span>${esc(selection)}</span></div>`:''}${asArray(dependencies).length?`<div class="entry-detail"><b>${esc(tx('Abhängige Eingaben'))}</b><span>${asArray(dependencies).map(esc).join(', ')}</span></div>`:''}${item.upstream_claim_id?`<div class="entry-detail"><b>${esc(tx('Vorgelagerte Aussage'))}</b><span><button type="button" class="inline-node" data-node="${esc(claimId)}">${esc(claimId)}</button></span></div>`:''}${item.derivation?`<div class="entry-detail"><b>${esc(tx('Ableitung'))}</b><span>${esc(english()?(item.derivation_en||item.derivation):item.derivation)}</span></div>`:''}${item.used_for?`<div class="entry-detail"><b>${esc(tx('Verwendet für'))}</b><span>${esc(english()?(item.used_for_en||item.used_for):item.used_for)}</span></div>`:''}${impact?`<div class="impact-line"><b>${esc(tx('Prüfstatus'))}</b><span class="impact-badge ${esc(impact[1])}">${esc(tx(impact[0]))}</span></div>`:''}${sourceButtons(item.source_refs)?`<div class="entry-sources"><b>${esc(tx('Quellen'))}</b>${sourceButtons(item.source_refs)}</div>`:''}</article>`;
    }).join('')}</div><p class="audit-caveat">Ein gesetzter oder externer Wert kann eine Rechnung reproduzierbar machen, ohne seine Herleitung aus TFPT zu begründen. Ein vorgelagerter Verweis übernimmt nicht automatisch den Beweisstatus der verknüpften Aussage.</p></section>`;
  }
  function proofGapsSection(items,contextNode){
    const gaps=asArray(items);if(!gaps.length)return '';
    return `<section class="dossier-section proof-gaps"><h3><span class="section-mark">?</span>Fehlende Beweisschritte</h3><p class="gap-intro">Diese Lücken markieren, wo ein Prüfprogramm die theoretische Herleitung nicht ersetzt.</p><div class="entry-cards">${gaps.map(gap=>`<article class="entry-card flagged"><h4>${esc(gap.label||'Offener Beweisschritt')}</h4>${gap.missing_step?`<div class="entry-detail"><b>Was fehlt</b><span>${esc(gap.missing_step)}</span></div>`:''}${gap.why_checks_insufficient?`<div class="entry-detail"><b>Warum Checks nicht genügen</b><span>${esc(gap.why_checks_insufficient)}</span></div>`:''}${asArray(gap.affected_result_ids).length?`<div class="entry-detail"><b>Betroffene Ergebnisse</b><span>${asArray(gap.affected_result_ids).map(id=>{const resultId=normalizeGraphId(id,'result',contextNode);return `<button type="button" class="inline-node" data-node="${esc(resultId)}">${esc(resultId)}</button>`}).join(' ')}</span></div>`:''}${sourceButtons(gap.source_refs)?`<div class="entry-sources"><b>Quellen</b>${sourceButtons(gap.source_refs)}</div>`:''}</article>`).join('')}</div></section>`;
  }
  function staticInputSection(data,note=''){
    const raw=data?.syntactic_input_candidates||[];const candidates=asArray(raw);if(!candidates.length&&!note)return '';
    const noteSection=note?`<section class="dossier-section script-audit-note"><h3><span class="section-mark">i</span>Hinweis zur Skriptauswertung</h3><p>${esc(note)}</p></section>`:'';
    const inputsSection=candidates.length?`<section class="dossier-section static-inputs"><h3><span class="section-mark">⌘</span>Statische Eingabekandidaten</h3><p class="gap-intro">Syntaxbeobachtungen aus dem Skriptbestand – noch keine semantische Bewertung von Herkunft, Verwendung oder Beweiswirkung.</p><div class="entry-cards">${candidates.map(item=>`<article class="entry-card compact"><h4>${esc(item.name||'Eingabekandidat')}</h4>${item.value!=null?`<code>${esc(item.value)}</code>`:''}<div class="static-input-meta">${item.kind?`<span>${esc(item.kind)}</span>`:''}${item.line!=null?`<span>Zeile ${esc(item.line)}</span>`:''}</div></article>`).join('')}</div></section>`:'';
    return `${noteSection}${inputsSection}`;
  }
  function staticCandidateOriginsSection(node){
    const a=attrs(node),s=semantic(node),normalizations=asArray(a.annotation_normalizations||s.annotation_normalizations);
    const rows=normalizations.filter(item=>item?.kind==='static_candidate_source_origins').flatMap(item=>asArray(item.cases));
    const cases=[...new Map(rows.filter(item=>item&&typeof item==='object').map(item=>[item.id||item.candidate_id||`${item.script_id}:${item.candidate_name}:${item.source?.line}`,item])).values()];
    if(!cases.length)return '';
    const roleNames={
      author_chosen_diagnostic_threshold:'diagnostischer Schwellenwert aus dem Skript',declared_scan_tolerance:'deklarierte Scantoleranz',
      finite_matrix_test_size:'Testgröße einer endlichen Matrix',finite_cylinder_instance:'endliche Zylindermodell-Instanz',
      conditional_integer_split_anchor:'Anker unter ganzzahliger Aufspaltung',external_verification_height:'extern verifizierte Höhe',
      unused_cited_constant_declaration:'zitierte Konstante, hier ungenutzt',declared_zero_position_error_budget:'deklarierte Positionsfehlerschranke',
      declared_tail_margin:'deklarierter Randsicherheitsfaktor',numerical_discretization_setting:'numerische Diskretisierung',
      quadrature_order:'Quadraturordnung',one_finite_model_instance:'endliche Modell-Instanz',finite_structure_entry:'Komponente einer endlichen Struktur',
      function_parameter_default:'Standardwert eines Funktionsparameters',local_initializer_or_formula:'lokale Initialisierung oder Formel',
      counter_or_diagnostic:'Zähler oder Diagnose',computed_check:'berechnete Prüfung',source_integrity_readback:'Quelltextintegritätsprüfung'
    };
    const basisNames={
      author_chosen_test_threshold:'im Skript gewählter Testschwellenwert',conditional_on_integer_split_and_typing:'bedingt durch ganzzahlige Aufspaltung und Typisierung',
      external_cited_verification_input:'zitierter externer Verifikationswert',classical_source_cited_but_unused_here:'klassische Quelle zitiert, hier ungenutzt',
      declared_diagnostic_budget:'deklarierte Diagnosetoleranz',declared_numerical_budget_without_per_root_certificate:'numerisches Budget ohne Zertifikat je Nullstelle',
      declared_numerical_margin:'deklarierter numerischer Spielraum',numerical_discretization_setting:'numerische Diskretisierung',
      finite_numerical_test_instance:'endliche numerische Testinstanz',one_finite_model_instance:'eine endliche Modell-Instanz',
      algorithmic_or_structural_definition:'algorithmische oder strukturelle Definition',computed_or_diagnostic:'berechneter oder diagnostischer Wert',explicit_in_source:'explizite Quellcode-Angabe'
    };
    const sourceOriginNames={literal_in_source:'Literal im Quellcode',function_default_in_source:'Standardwert im Quellcode',source_integrity_readback:'eingelesener Quelltext',computed_from_source_expression:'aus Quellcode-Ausdruck berechnet'};
    const friendly=value=>String(value||'').replaceAll('_',' ');
    const authored=(item,key)=>{
      if(english())return {text:item[`${key}_en`]||item[key]||'',original:!item[`${key}_en`]&&!!item[key]};
      return {text:item[`${key}_de`]||item[key]||'',original:!item[`${key}_de`]&&!!item[key]};
    };
    const sourceRefs=item=>{
      const refs=[...asArray(item.source_refs),...asArray(item.sources)];
      if(item.source?.source_ref)refs.push({id:item.source.source_ref,path:item.source.path,line_start:item.source.line,line_end:item.source.line});
      if(item.source?.context_ref){const range=asArray(item.source.context_lines);refs.push({id:item.source.context_ref,path:item.source.path,line_start:range[0],line_end:range[1]});}
      for(const ref of asArray(item.claim_source_contexts))if(ref?.id)refs.push({id:ref.id,path:ref.path,line_start:ref.line_start,line_end:ref.line_end});
      return refs.filter((ref,index,all)=>sourceId(ref)&&all.findIndex(other=>sourceId(other)===sourceId(ref))===index);
    };
    const cards=cases.map(item=>{
      const name=item.candidate_name||item.candidate_id||item.id||tx('Wert'),role=item.semantic_role||item.role||'',candidateKind=item.candidate_kind||'',basis=item.selection_basis||'',purpose=authored(item,'purpose'),rationale=authored(item,'selection_rationale'),originNote=authored(item,'source_origin_note'),limit=authored(item,'interpretation_limit');
      const source=item.source||{},location=[source.path,source.line!=null?`${tx('Zeile')} ${source.line}`:''].filter(Boolean).join(' · '),refs=sourceRefs(item);
      const badgeText='Originalwortlaut · nicht übersetzt';
      const field=(title,value,original=false)=>value?`<div class="entry-detail code-origin-detail"><b>${esc(tx(title))}</b><span>${esc(value)}${original?` <em class="origin-language-note">${esc(tx(badgeText))}</em>`:''}</span></div>`:'';
      return `<article class="entry-card code-origin-card"><div class="code-origin-heading"><h4>${esc(name)}</h4>${item.expression!=null?`<code>${esc(item.expression)}</code>`:''}</div><div class="static-input-meta">${item.script_id?`<span>${esc(item.script_id)}</span>`:''}${location?`<span>${esc(location)}</span>`:''}</div>${role||basis||candidateKind?`<div class="code-origin-tags">${role?`<span>${esc(tx(roleNames[role]||friendly(role)))}</span>`:candidateKind?`<span>${esc(tx('Syntaxklassifikation'))}: ${esc(tx(friendly(candidateKind)))}</span>`:''}${basis?`<span>${esc(tx(basisNames[basis]||friendly(basis)))}</span>`:''}</div>`:''}${field('Rolle im Skript',purpose.text,purpose.original)}${field('Auswahlgrund',rationale.text,rationale.original)}${field('Ursprungshinweis',originNote.text,originNote.original)}${field('Grenze dieser Einordnung',limit.text,limit.original)}${item.source_origin?`<small class="origin-source-kind">${esc(tx(sourceOriginNames[item.source_origin]||friendly(item.source_origin)))}</small>`:''}${refs.length?`<div class="entry-sources"><b>${esc(tx('Quellen'))}</b>${sourceButtons(refs)}</div>`:''}</article>`;
    }).join('');
    return `<details class="dossier-section code-origins"><summary><span class="section-mark">⌘</span><strong>${esc(tx('Geprüfte Code-Herkunft'))}</strong><span class="code-origin-count">${cases.length} ${esc(tx('Einträge'))}</span></summary><p class="code-origin-intro">${esc(tx('Die Einträge erklären, woher ein Skriptwert stammt und welche Rolle er dort hat. Sie belegen keine physikalische Herleitung oder allgemeine Gültigkeit.'))}</p>${english()&&cases.some(item=>!item.purpose_en&&!item.selection_rationale_en&&!item.purpose_de&&!item.selection_rationale_de)?`<p class="language-fallback">${esc(tx('Originalwortlaut · nicht übersetzt'))}</p>`:''}<div class="entry-cards">${cards}</div></details>`;
  }
  function renderAuditRelatedClaims(claims,contextNode){
    const values=asArray(claims);if(!values.length)return '';
    return `<section class="dossier-section related-audits"><h3><span class="section-mark">↗</span>Verknüpfte Aussagen im Audit</h3><p class="gap-intro">Zu diesen Aussagen enthält der Skript-Audit eigene Herkunftshinweise und Beweislücken.</p><div class="entry-cards">${values.map(item=>{
      const id=normalizeGraphId(item.id||item.claim_id,'claim',contextNode),inputs=asArray(item.input_audit),gaps=asArray(item.proof_gaps);
      return `<article class="entry-card related-audit-card"><div class="related-audit-heading"><h4><button type="button" class="inline-node" data-node="${esc(id)}">${esc(item.label||id||'Aussage')}</button></h4><span class="audit-state">${esc(item.enrichment_state||'')}</span></div>${inputs.length?`<div class="related-audit-inputs"><b>Eingaben</b>${inputs.map(input=>{const origin=inputOrigin[input.origin]||inputOrigin.unknown,impact=inputImpact[input.impact];return `<div class="related-audit-row"><span><strong>${esc(input.name||'Eingabe')}</strong><code>${esc(input.value==null?'(Wert nicht angegeben)':typeof input.value==='string'?input.value:JSON.stringify(input.value))}</code></span><span class="origin-badge ${origin[1]}">${esc(origin[0])}</span>${impact?`<span class="impact-badge ${impact[1]}">${esc(impact[0])}</span>`:''}${input.explanation?`<p>${esc(input.explanation)}</p>`:`<p>${esc(origin[2])}</p>`}${sourceButtons(input.source_refs)?`<div class="entry-sources"><b>Quellen</b>${sourceButtons(input.source_refs)}</div>`:''}</div>`;}).join('')}</div>`:''}${gaps.length?`<div class="related-audit-gaps"><b>Offene Beweisschritte</b>${gaps.map(gap=>`<div><strong>${esc(gap.label||'Offener Schritt')}</strong>${gap.missing_step?`<span>${esc(gap.missing_step)}</span>`:''}${gap.why_checks_insufficient?`<p>${esc(gap.why_checks_insufficient)}</p>`:''}${sourceButtons(gap.source_refs)}</div>`).join('')}</div>`:''}${!inputs.length&&!gaps.length?'<p>Für diese verknüpfte Aussage sind hier keine Eingaben oder Beweislücken verzeichnet.</p>':''}</article>`;
    }).join('')}</div></section>`;
  }
  function listSection(title, icon, value, cls = '') {
    if (value == null || (Array.isArray(value) && !value.length) || value === '') return '';
    const content = Array.isArray(value) ? `<ul>${textBlock(value)}</ul>` : textBlock(value);
    return `<section class="dossier-section ${cls}"><h3><span class="section-mark">${icon}</span>${title}</h3>${content}</section>`;
  }
  function statusSection(value,status=''){
    const entries=[];if(status)entries.push(['Status',status]);
    if(typeof value==='string'&&value)entries.push(['Einordnung',value]);
    else if(Array.isArray(value))value.forEach(item=>entries.push([item?.label||item?.name||'Einordnung',readableValue(item?.text||item?.statement||item?.value||item)]));
    else if(value&&typeof value==='object')Object.entries(value).forEach(([key,val])=>entries.push([key.replaceAll('_',' '),readableValue(val)]));
    if(!entries.length)return '';
    return `<section class="dossier-section"><h3><span class="section-mark">◌</span>Erkenntnisstatus</h3><div class="status-facts">${entries.map(([key,val])=>`<div><b>${esc(key)}</b><span>${esc(val)}</span></div>`).join('')}</div></section>`;
  }
  function setNotice(message = '', isError = false) {
    notice.hidden = !message;
    notice.textContent = message;
    notice.classList.toggle('error', isError);
  }
  async function api(path, options = {}) {
    const response = await fetch(path, { headers: { 'Accept': 'application/json', ...(options.body ? { 'Content-Type': 'application/json' } : {}) }, ...options });
    if (!response.ok) throw new Error(`${tx('Die Anfrage ist fehlgeschlagen')} (${response.status}).`);
    return response.json();
  }
  function plural(count, one, many = `${one}e`) { return `${count} ${count === 1 ? one : many}`; }
  function normalizedTopics(data) {
    const incoming = Array.isArray(data) ? data : data?.topics || data?.results || [];
    const fallbackById=new Map(topicCatalog.map(([id,name,glyph])=>[id,{name,glyph}]));
    const glyphFor=(node,index)=>attrs(node).glyph||attrs(node).symbol||fallbackById.get(node.id)?.glyph||['◉','⌘','◷','✳','α','◌','⌁','◍','ζ','✓','◇','∿','⊞','◈'][index%14];
    const topicFromNode=(node,index)=>({id:node.id,label:english()?(node.label_en||semantic(node).title_en||attrs(node).title_en||attrs(node).translations?.en?.title||topicNamesEn[node.id]||label(node)):label(node),glyph:glyphFor(node,index),node,count:attrs(node).count||attrs(node).node_count||attrs(node).topic_count||''});
    const listed=incoming.filter(node=>node?.id).map(topicFromNode);
    const listedIds=new Set(listed.map(topic=>topic.id));
    const legacy=topicCatalog.filter(([id])=>!listedIds.has(id)).map(([id,name,glyph])=>({id,label:name,glyph,node:null,count:''}));
    return [...listed,...legacy];
  }
  function nodeHref(node) { return `#node/${encodeURIComponent(node?.id || '')}`; }
  function parseSourceRoute(value){const raw=String(value||'').replace(/^#?source\//,''),split=raw.indexOf('?'),encodedId=split<0?raw:raw.slice(0,split),params=new URLSearchParams(split<0?'':raw.slice(split+1)),start=params.get('line_start'),end=params.get('line_end');return {id:decodeURIComponent(encodedId),line_start:start==null?undefined:Number(start),line_end:end==null?undefined:Number(end)};}
  function currentRoute(){try{return decodeURIComponent(location.hash.slice(1));}catch{return location.hash.slice(1);}}
  function isCurrentTopicRoute(id){return currentRoute()===`topic/${id}`;}
  function isCurrentNodeRoute(id){const route=currentRoute();return route===`node/${id}`||route===`chain/${id}`||(route.startsWith('source/')&&state.current?.node?.id===id);}
  function sourceRouteHref(id,range={}){const params=new URLSearchParams();if(range.line_start!=null)params.set('line_start',range.line_start);if(range.line_end!=null)params.set('line_end',range.line_end);const query=params.toString();return `#source/${encodeURIComponent(id)}${query?`?${query}`:''}`;}
  function topicFor(node) {
    const ids = asArray(semantic(node).topic_ids || attrs(node).topic_ids);
    return topicCatalog.find(t => ids.includes(t[0]))?.[0] || '';
  }
  function goNode(id, replace = false) {
    if (!id) return;
    if(String(id).startsWith('topic:')){const next=`#topic/${encodeURIComponent(id)}`;if(replace)history.replaceState(null,'',next);else if(location.hash!==next)location.hash=next;else loadTopic(id);return;}
    const next = `#node/${encodeURIComponent(id)}`;
    const current=decodeURIComponent(location.hash.slice(6));
    if (replace) { history.replaceState(null, '', next); loadNode(id); }
    else if (location.hash.startsWith('#node/') && current===id) loadNode(id);
    else location.hash = next;
  }
  function openChain(id){if(!id)return;const next=`#chain/${encodeURIComponent(id)}`;if(location.hash===next){loadNode(id,true);return;}location.hash=next;}
  function goTopic(topic) {
    const match = state.topics.find(t => t.id === topic);
    if(match){const next=`#topic/${encodeURIComponent(topic)}`;if(location.hash===next)loadTopic(topic);else location.hash=next;return;}
    search(topicCatalog.find(t => t[0] === topic)?.[1] || topic);
  }
  function topicName(id) {
    if (english() && topicNamesEn[id]) return topicNamesEn[id];
    // API labels may be English even when the UI is German. Prefer the curated
    // local catalog for known topics so breadcrumbs, sidebars and related links
    // all follow the selected UI language consistently.
    return topicCatalog.find(t => t[0] === id)?.[1] || state.topics.find(t => t.id === id)?.label || id;
  }
  function normalizeGraphId(id,kind,contextNode){
    const raw=String(id||'');if(!raw||raw.includes(':'))return raw;
    if(kind==='claim')return `claim:${raw}`;
    if(kind==='result'){
      const claimId=String(contextNode?.id||'').replace(/^claim:/,'');
      return claimId?`result:${claimId}:${raw}`:raw;
    }
    return raw;
  }
  function updateSidebarTopics() {
    sideTitle.textContent = tx('Themen entdecken');
    sideCount.textContent = `${state.topics.length} ${tx('Themen')}`;
    side.innerHTML = `<div class="topic-group-label">${esc(tx('Themenfelder'))}</div>`;
    for (const topic of state.topics) {
      const button = document.createElement('button');
      button.type = 'button'; button.className = 'topic-button'; button.dataset.topic = topic.id;
      button.innerHTML = `<span class="topic-glyph">${esc(topic.glyph)}</span><span class="topic-name"></span>${topic.count ? `<span class="topic-count">${esc(topic.count)}</span>` : ''}`;
      $('.topic-name', button).textContent = topicName(topic.id);
      button.addEventListener('click', () => goTopic(topic.id)); side.append(button);
    }
  }
  function setSidebarResults(query, results, count) {
    state.searchQuery=query;
    sideTitle.textContent = tx('Suchtreffer'); sideCount.textContent = String(count ?? results.length);
    side.replaceChildren();
    const meta = document.createElement('div'); meta.className = 'sidebar-search-meta'; meta.textContent = `${tx('Ergebnisse für')} „${query}“`; side.append(meta);
    for (const node of results.slice(0, 30)) {
      const button = document.createElement('button'); button.type = 'button'; button.className = 'result-button';
      const kind = document.createElement('span'); kind.className = 'result-type'; kind.textContent = typeLabel(node);
      const copy = document.createElement('span'); copy.className = 'result-copy';
      const title = document.createElement('strong'); title.textContent = excerpt(label(node), 110);
      const desc = document.createElement('small'); desc.textContent = excerpt(summary(node), 90);
      copy.append(title, desc); button.append(kind, copy); button.addEventListener('click', () => goNode(node.id)); side.append(button);
    }
    const back = document.createElement('button'); back.type='button'; back.className='topic-button'; back.innerHTML=`<span class="topic-glyph">←</span><span class="topic-name">${esc(tx('Alle Themen'))}</span>`; back.addEventListener('click', updateSidebarTopics); side.append(back);
  }
  function statusKind(node) {
    const assess=assessment(node),status=String(ledgerStatus(node)||'').toLowerCase(),category=String(assess.category||'').toLowerCase();
    if (/open|unresolved|offen|ungeklärt/.test(status)||/open|unresolved|offen|ungeklärt/.test(category)) return 'open';
    if (/axiom|assum|declared input|gesetzte voraussetzung/.test(status)||/assum/.test(category)) return 'assumed';
    if (/proof|proved|verified|abgenommen|accepted/.test(status)) return 'verified';
    return '';
  }
  function badge(node, kind) {
    if (kind === 'ai') {const value=enrichment(node);return value?`<span class="badge ${value.className}">${esc(value.label)}</span>`:'';}
    const status = ledgerStatus(node);
    if (kind === 'ledger') {if(!status&&!/claim/i.test(typeName(node)))return '';return `<span class="badge ${statusKind(node)}">${esc(tx('Originalstatus'))}: ${esc(status || tx('siehe Ledger'))}</span>`;}
    return '';
  }
  function maybeJson(value) { return Array.isArray(value) ? value : value == null ? [] : [value]; }
  function sourceLabel(source) {
    if(typeof source==='string')return source;
    const page=source?.page,base=source?.title||source?.label||source?.path||source?.id||tx('Quelle');
    return page!=null?`${tx('PDF-Seite')} ${page} · ${base}`:base;
  }
  function sourceRegistry(sources) {
    const records = new Map();
    for (const source of asArray(sources)) {
      const id = sourceId(source);
      if (id && typeof source === 'object') records.set(id, source);
    }
    return records;
  }
  function resolveSource(source, registry) {
    const id = sourceId(source);
    const record = registry?.get(id);
    if (!record) return source;
    const ref = typeof source === 'object' && source ? source : {};
    return {
      ...record,
      ...ref,
      id: id || record.id,
      path: ref.path || record.path,
      title: ref.title || record.title,
      label: ref.label || record.label,
      page: ref.page ?? record.page,
      line_start: ref.line_start ?? record.line_start,
      line_end: ref.line_end ?? record.line_end,
    };
  }
  function renderSources(sources, registry) {
    const list = asArray(sources); if (!list.length) return '<p class="aside-empty">Keine direkt verknüpften Quellen gefunden.</p>';
    return list.map(ref => {
      const id = sourceId(ref), source = resolveSource(ref, registry);
      const path = source && typeof source === 'object' ? (source.path || source.title || source.label || '') : '';
      const start = source?.line_start, end = source?.line_end;
      const range = start != null ? `${tx('Zeile')} ${start}${end != null && end !== start ? `–${end}` : ''}` : '';
      const label = path || (source?.page != null ? `${tx('PDF-Seite')} ${source.page}` : id || tx('Quelle'));
      const detail = [path && range, !path && source?.page != null && range].filter(Boolean).join(' · ') || (path ? '' : id);
      return `<button type="button" class="source-row" data-source="${esc(id)}" data-line-start="${esc(start??'')}" data-line-end="${esc(end??'')}" title="${esc(id)}"><strong>${esc(label)}</strong><small>${esc(detail)}</small></button>`;
    }).join('');
  }
  function renderRelated(nodes) {
    const list = asArray(nodes); if (!list.length) return '<p class="aside-empty">Noch keine verknüpften Einträge.</p>';
    return list.slice(0, 18).map(node => `<button type="button" class="related-row" data-node="${esc(node.id)}"><strong>${esc(String(node.id||'').startsWith('topic:')?topicName(node.id):label(node))}</strong><small>${esc(typeLabel(node))}</small></button>`).join('');
  }
  function goalReadingMarkup(data,node){
    const ui=(en,de)=>english()?en:de;
    const a=attrs(node),steps=asArray(english()?a.required_steps_en:(a.required_steps_de||a.required_steps)),related=asArray(data.related_nodes),edges=asArray(data.edges),claims=related.filter(item=>/claim/i.test(typeName(item))),topics=related.filter(item=>String(item.id||'').startsWith('topic:'));
    const mapped=claims.filter(claim=>edges.some(edge=>(edge.src===claim.id||edge.source===claim.id)&&(edge.dst===node.id||edge.target===node.id)&&/addresses_open_goal/.test(edge.type||'')));
    const linked=mapped.length?mapped:claims;
    const claimCard=claim=>{const ca=attrs(claim),assessment=ca.assessment||{},englishOverview=english()&&freshEnglish(claim,'overview'),overview=english()?(englishOverview||summary(claim)):summary(claim),category=String(assessment.category||'').toLowerCase(),statusLabels={established:['belegt im Rahmen','Recorded within scope'],recorded_established:['belegt im Rahmen','Recorded within scope'],conditional:['bedingt','Conditional'],mixed:['gemischt','Mixed'],open:['offen','Open'],excluded:['ausgeschlossen','Excluded'],assumed:['angenommen','Assumed'],historical:['historisch','Historical'],unknown:['Status nicht automatisch zugeordnet','Status not automatically classified'],methodological_rule:['Geltende Prüfregel','Active methodological rule']},rawStatus=ledgerStatus(claim),status=(english()&&assessment.label_en)||statusLabels[category]?.[english()?1:0]||(english()?(rawStatus?`Original status: ${rawStatus}`:ui('Status not automatically classified','Status nicht automatisch zugeordnet')):(assessment.label||rawStatus||'Status nicht automatisch zugeordnet')),relation=mapped.includes(claim)?ui('Addresses this goal','Bezieht sich auf dieses Ziel'):ui('Related in the dossier','Im Dossier verknüpft'),summaryText=overview?String(overview):'',summaryBlock=summaryText.length>360?`<p>${esc(excerpt(summaryText,340))}</p><details class="reading-claim-full"><summary>${esc(ui('Full claim summary','Vollständige Aussage'))}</summary><p>${esc(summaryText)}</p></details>`:summaryText?`<p>${esc(summaryText)}</p>`:'',refs=sourceButtons(asArray(ca.source_refs||semantic(claim).source_refs).slice(0,3));return `<article class="reading-linked-claim"><div class="reading-card-head"><strong>${esc(label(claim))}</strong><span class="reading-status">${esc(status)}</span></div>${summaryBlock}${english()&&!englishOverview?`<small class="language-fallback">${esc(ui('No English summary is available; the authored text is shown.',''))}</small>`:''}<small>${esc(relation)}${edges.some(edge=>(edge.src===claim.id||edge.source===claim.id)&&(edge.dst===node.id||edge.target===node.id)&&edge.attrs?.method==='curated_source_mapping')?` · ${esc(ui('curated source mapping','kuratierte Quellenzuordnung'))}`:''}</small>${refs?`<div class="reading-sources">${refs}</div>`:''}<button type="button" class="reading-open-link" data-node="${esc(claim.id)}">${esc(ui('Open claim','Aussage öffnen'))} ↗</button></article>`;};
    const stepList=steps.length?`<ol class="reading-steps">${steps.map(step=>`<li>${esc(readableValue(step))}</li>`).join('')}</ol>`:`<p>${esc(ui('No required steps are recorded for this goal.','Für dieses Ziel sind noch keine erforderlichen Schritte erfasst.'))}</p>`;
    const topicsMarkup=topics.length?`<div class="reading-topic-links">${topics.map(topic=>`<button type="button" data-node="${esc(topic.id)}">${esc(topicName(topic.id))} ↗</button>`).join('')}</div>`:'';
    const nextStep=(english()&&a.next_step_en)||a.next_step||steps[0]||'';
    return `<section class="goal-reading-story" aria-label="${esc(ui('Goal overview','Zielübersicht'))}"><p class="reading-eyebrow">${esc(ui('WHERE IT SITS','EINORDNUNG'))}</p>${topicsMarkup}<section class="reading-story-part"><h3>${esc(ui('What this goal requires','Was dieses Ziel verlangt'))}</h3>${stepList}</section><section class="reading-story-part"><h3>${esc(ui('Related work and its status','Verknüpfte Arbeit und ihr Status'))}</h3>${linked.length?`<div class="reading-claim-list">${linked.slice(0,3).map(claimCard).join('')}</div>${linked.length>3?`<details class="reading-more"><summary>${esc(english()?`Show all ${linked.length} linked claims`:`Alle ${linked.length} verknüpften Aussagen anzeigen`)}</summary><div class="reading-claim-list">${linked.slice(3).map(claimCard).join('')}</div></details>`:''}`:`<p>${esc(ui('No claims are linked to this goal yet.','Noch keine Aussagen sind diesem Ziel zugeordnet.'))}</p>`}<p class="reading-caveat">${esc(ui('A topic or goal link records relevance; it does not mean a requirement is satisfied.','Eine Themen- oder Zielverknüpfung dokumentiert Relevanz, nicht die Erfüllung eines Schritts.'))}</p></section>${nextStep?`<section class="reading-story-part reading-next"><h3>${esc(ui('Next stated requirement','Nächster angegebener Schritt'))}</h3><p>${esc(readableValue(nextStep))}</p></section>`:''}</section>`;
  }
  function renderDossier(data) {
    const node = data.node || data;
    const s = semantic(node), a = attrs(node), ledger = a.ledger || {};
    const audit=data.audit||{};
    const topicId = topicFor(node);
    const overview = (english()&&(isClaimNode(node)?freshEnglish(node,'overview'):(s.overview_en||a.overview_en||a.translations?.en?.overview)))||s.overview_de||summary(node);
    const field=(key,fallback='')=>(english()&&(s[`${key}_en`]||a[`${key}_en`]))||s[key]||a[key]||fallback;
    const whatWeHave=field('what_we_have');
    const whatWeSeek=field('what_we_seek');
    const requiredSteps=field('required_steps',[]);
    const conditions=field('conditions',[]);
    const counterevidence=field('counterevidence',[]);
    const epistemic=field('epistemic');
    const epistemicStatus=s.status||a.epistemic_status||a.status||a.verdict||'';
    const hasTextFormatError=containsFormatControl(s);
    const ownSourceList=asArray(s.source_refs||a.source_refs||a.original_contract?.source_refs);
    const relatedSourceList=asArray(data.sources||a.sources).filter(source=>!ownSourceList.some(own=>sourceId(own)&&sourceId(own)===sourceId(source)));
    const sourceList=ownSourceList.length?ownSourceList:asArray(a.sources||data.sources);
    const sourceRecords=sourceRegistry([...asArray(data.sources),...asArray(a.sources)]);
    const related = data.related_nodes || [];
    const statement = (english()&&(!isClaimNode(node)&&(a.statement_en||a.claim_en||a.formula_en)))||a.statement || a.claim || a.formula || ledger.statement || ledger.claim || ledger.formula || '';
    const meanings = (english()&&(isClaimNode(node)?freshEnglish(node,'interpretation'):(s.interpretation_en||a.interpretation_en||a.translations?.en?.interpretation||a.meaning_en)))||s.interpretation_de || a.interpretation || a.meaning || '';
    const extractionQuality=a.extraction_quality||{};
    const bases = s.premises || s.objects || a.premises || a.assumptions || [];
    const results = s.results || a.results || a.consequences || [];
    const derivations = s.derivations || a.derivations || [];
    const unknowns = s.unknowns || a.unknowns || a.open_questions || [];
    const scripts = s.script_analysis || a.scripts || [];
    const lean = s.lean_links || a.lean_links || a.lean || [];
    const relationships = s.relationships || [];
    const inputAudits=s.input_audit||a.input_audit||audit.input_audit||[];
    const proofGaps=s.proof_gaps||a.proof_gaps||audit.proof_gaps||[];
    const staticInputs=a.static_inputs||audit;
    const relatedAuditClaims=asArray(audit.related_claims);
    const relatedInputs=relatedAuditClaims.flatMap(claim=>asArray(claim.input_audit));
    const relatedGaps=relatedAuditClaims.flatMap(claim=>asArray(claim.proof_gaps));
    const auditAttention=[...asArray(inputAudits),...relatedInputs].filter(item=>item.origin==='unknown'||['circularity_risk','missing_derivation'].includes(item.impact)).length;
    const proofGapCount=asArray(proofGaps).length+relatedGaps.length;
    const rawLedger = Object.keys(ledger).length ? ledger : Object.fromEntries(Object.entries(a).filter(([k]) => !['semantic', 'statement', 'summary', 'sources'].includes(k)));
    const goalReading=/goal/i.test(typeName(node));
    const patternInference=/theory_pattern/i.test(typeName(node))&&(a.reviewed===false||a.kind==='cross_source_inference');
    view.innerHTML = `<article class="dossier">
      <header class="dossier-top"><div class="breadcrumb"><button type="button" data-home>Wissensatlas</button><span>›</span><span>${esc(topicId ? topicName(topicId) : 'Eintrag')}</span><span>›</span><span>${esc(typeLabel(node))}</span></div>
        <div class="dossier-type">${esc(typeLabel(node))}</div><h2>${esc(label(node))}</h2><div class="dossier-id">${esc(node.id || '')}</div><div class="dossier-badges">${badge(node,'ledger')}${patternInference?`<span class="badge audit-warning">${esc(tx('Arbeitsinferenz · ungeprüft'))}</span>`:''}${badge(node,'ai')}${statusKind(node)==='open' ? `<span class="badge open">${esc(tx('Offene Frage / ungeklärt'))}</span>` : ''}${auditAttention?`<span class="badge audit-warning">${auditAttention} ${esc(tx('Herkunfts-/Herleitungshinweis'))}</span>`:''}${proofGapCount?`<span class="badge open">${proofGapCount} ${esc(tx('Beweisschritte offen'))}</span>`:''}${extractionQuality.state==='needs_review'?`<span class="badge audit-warning">${esc(tx('Extraktion braucht Prüfung'))}</span>`:''}</div>
        ${hasTextFormatError?'<p class="format-warning">Die KI-Textdarstellung enthält Formatierungsfehler. Für Formeln bitte Originalquellen prüfen.</p>':''}
        ${extractionQuality.state==='needs_review'?`<section class="quality-warning"><strong>${esc(tx('Extraktion braucht Prüfung'))}</strong><p>${esc(tx('Ein Qualitätsaudit hat mögliche Probleme in der KI-Extraktion markiert. Diese Hinweise sind keine fachliche Bewertung des ursprünglichen Anspruchs.'))}</p>${extractionQuality.scope?`<p><b>${esc(tx('Prüfumfang'))}:</b> ${esc(extractionQuality.scope)}</p>`:''}${asArray(extractionQuality.issues).length?`<ul>${extractionQuality.issues.map(issue=>{const refs=[...asArray(issue.source_refs),...asArray(issue.evidence?.source_refs),...asArray(issue.evidence)].filter(ref=>typeof ref==='string'?ref.startsWith('source:'):!!sourceId(ref)).slice(0,4);return `<li>${esc((english()&&issue.message_en)||issue.message||issue.code||readableValue(issue))}${refs.length?`<div class="quality-issue-sources">${sourceButtons(refs)}</div>`:''}</li>`}).join('')}</ul>`:''}</section>`:''}
        ${english()&&!(isClaimNode(node)?freshEnglish(node,'overview'):(s.overview_en||a.overview_en||a.summary_en||a.translations?.en?.overview))?`<p class="language-fallback">${esc(tx('Für diesen Eintrag ist keine englische Übersetzung vorhanden; der Originaltext bleibt sichtbar.'))}</p>`:''}
        ${english()&&isClaimNode(node)&&[results,derivations,bases,scripts,lean,unknowns,relationships,inputAudits,proofGaps].some(value=>asArray(value).length>0)?`<p class="language-fallback">${esc(tx('Detailed annotations are shown in their original language where a translation is unavailable.'))}</p>`:''}
      </header>
      <nav class="view-tabs" aria-label="Eintragsansichten"><button class="view-tab active" data-tab="detail">Einordnung</button><button class="view-tab" data-tab="graph">Lösungskette ansehen</button></nav>
      <div class="tab-content" id="tab-content"><div class="dossier-layout"><div class="dossier-body">
        ${overview ? `<section class="dossier-section semantic-note"><h3><span class="section-mark">i</span>Einordnung</h3>${textBlock(overview)}${enrichment(node)?.className==='ai' ? `<span class="badge ai">${esc(enrichment(node).label)}</span>` : ''}</section>` : ''}
        ${attemptDossierSection(node)}
        ${listSection('Was wir haben','✓',whatWeHave)}
        ${listSection('Was wir suchen','⌕',whatWeSeek)}
        ${listSection('Erforderliche Schritte','→',requiredSteps)}
        ${listSection('Bedingungen','◇',conditions)}
        ${listSection('Gegenbefunde','⊣',counterevidence)}
        ${statusSection(epistemic,epistemicStatus)}
        ${statement ? listSection('Aussage im Original', '§', statement, 'ledger-box') : ''}
        ${structuredSection('Was wurde gezeigt?', '✓', results, 'result',!isClaimNode(node))}
        ${listSection('Was bedeutet das?', '↗', meanings)}
        ${inputAuditSection(inputAudits,node)}
        ${proofGapsSection(proofGaps,node)}
        ${structuredSection('Worauf baut es auf?', '↑', bases, 'premise',!isClaimNode(node))}
        ${structuredSection('Herleitungen', '∴', derivations, 'derivation',!isClaimNode(node))}
        ${listSection('Zusammenhänge', '⌁', relationships)}
        ${listSection('Was bleibt offen?', '?', unknowns, 'open-section')}
        ${claimAttemptsSection(data.claim_attempts,node)}
        ${structuredSection('Prüfprogramme & Skripte', '⌘', scripts, 'script',!isClaimNode(node))}
        ${staticInputSection(staticInputs,audit.interpretation_note)}
        ${staticCandidateOriginsSection(node)}
        ${renderAuditRelatedClaims(audit.related_claims,node)}
        ${structuredSection('Lean-Verknüpfungen', 'L', lean, 'lean',!isClaimNode(node))}
        ${/claim/i.test(typeName(node))?`<section class="dossier-section"><h3><span class="section-mark">§</span>Original-Ledger</h3><div class="ledger-box"><p>${esc(rawLedger.status || rawLedger.verdict || rawLedger.claim_status || 'Originalstatus im Ledger nachsehen.')}</p>${rawLedger.rationale ? `<p>${esc(rawLedger.rationale)}</p>` : ''}${rawLedger.evidence ? `<p>${esc(rawLedger.evidence)}</p>` : ''}<details><summary class="tiny-label">Ledger-Felder einblenden</summary><dl class="ledger-fields">${Object.entries(rawLedger).map(([k,v]) => `<div><dt>${esc(k.replaceAll('_',' '))}</dt><dd>${esc(typeof v === 'string' ? v : JSON.stringify(v))}</dd></div>`).join('')}</dl></details></div></section>`:''}
        <div id="source-viewer"></div>
      </div><aside class="dossier-aside"><div class="aside-block"><h3>Verknüpfte Einträge</h3>${renderRelated(related)}</div><div class="aside-block"><h3>${ownSourceList.length?'Direkte Originalbelege':'Originalquellen'}</h3>${renderSources(sourceList,sourceRecords)}${relatedSourceList.length?`<details class="related-source-details"><summary class="tiny-label">Weitere verknüpfte Quellen · ${relatedSourceList.length}</summary>${renderSources(relatedSourceList,sourceRecords)}</details>`:''}</div></aside></div></div>
    </article>`;
    if(goalReading){
      const article=view.querySelector('.dossier'),body=view.querySelector('.dossier-body'),aside=view.querySelector('.dossier-aside');
      article?.classList.add('goal-reading-dossier');
      if(body){
        const details=document.createElement('details');details.className='goal-specialist-details';
        const summaryEl=document.createElement('summary');summaryEl.textContent=english()?'Full dossier · sources, checks, and specialist details':'Vollständiges Dossier · Quellen, Prüfungen und Fachdetails';
        const full=document.createElement('div');full.className='goal-specialist-content';
        while(body.firstChild)full.append(body.firstChild);
        const viewTabs=view.querySelector('.view-tabs');if(viewTabs)full.prepend(viewTabs);
        if(aside)full.append(aside);
        details.append(summaryEl,full);
        const story=document.createElement('div');story.className='goal-reading';story.innerHTML=goalReadingMarkup(data,node);
        body.append(story,details);
      }
    }
    const frontier=researchContextPanel(node.id),overviewSection=view.querySelector('.dossier-body .semantic-note');
    if(goalReading)view.querySelector('.goal-specialist-content')?.prepend(frontier);
    else if(overviewSection)overviewSection.insertAdjacentElement('afterend',frontier);else view.querySelector('.dossier-body')?.prepend(frontier);
    view.querySelector('[data-home]').addEventListener('click', home);
    view.querySelectorAll('[data-node]').forEach(el => el.addEventListener('click', () => goNode(el.dataset.node)));
    view.querySelectorAll('[data-source]').forEach(el => el.addEventListener('click', () => loadSource(el.dataset.source, el, {line_start: el.dataset.lineStart ? Number(el.dataset.lineStart) : undefined, line_end: el.dataset.lineEnd ? Number(el.dataset.lineEnd) : undefined})));
    view.querySelector('[data-more-attempts]')?.addEventListener('click',async event=>{const button=event.currentTarget;button.disabled=true;button.textContent='Wird geladen …';try{const response=data.claim_attempts||{},offset=Number(response.offset||0)+asArray(response.attempts).length,next=await api(`/api/attempts?id=${encodeURIComponent(node.id)}&limit=100&offset=${offset}`);data.claim_attempts={...next,attempts:[...asArray(response.attempts),...asArray(next.attempts)]};state.current=data;renderDossier(data);}catch{button.disabled=false;button.textContent='Laden fehlgeschlagen · erneut versuchen';}});
    view.querySelectorAll('[data-tab]').forEach(el => el.addEventListener('click', () => {
      view.querySelectorAll('[data-tab]').forEach(tab => tab.classList.toggle('active', tab === el));
      if (el.dataset.tab === 'graph') {
        const claimInputs=relatedAuditClaims.flatMap(claim=>{const claimId=normalizeGraphId(claim.id||claim.claim_id,'claim',node);return asArray(claim.input_audit).map((item,index)=>({...item,upstream_claim_id:claimId,graph_id:`audit-input:${claimId||'claim'}:${item.id||item.name||index}`}));});
        const graphAudits=[...asArray(inputAudits),...claimInputs];
        const graphNode=graphAudits.length?{...node,attrs:{...a,semantic:{...s,input_audit:graphAudits}}}:node;
        const claimNodes=relatedAuditClaims.filter(claim=>claim.id||claim.claim_id).map(claim=>{const id=normalizeGraphId(claim.id||claim.claim_id,'claim',node);return {id,label:claim.label||id,type:'Claim',attrs:{}};});
        loadTrace(graphNode,[...related,...claimNodes],data.edges||[]);
      }
      else renderDossier(data);
    }));
  }
  async function loadTrace(root, related, fallbackEdges, depth=4) {
    const holder=$('#tab-content'); if(!holder)return;
    const request=++state.graphSerial,rootId=root?.id;
    const isCurrent=()=>request===state.graphSerial&&holder.isConnected&&$('#tab-content')===holder&&state.current?.node?.id===rootId;
    holder.innerHTML='<div class="loading-block"><span class="spinner"></span><span>Verknüpfungen werden verfolgt …</span></div>';
    try {
      let trace;
      if(/claim/i.test(typeName(root))){try{trace=await api(`/api/chain?id=${encodeURIComponent(root.id)}&depth=${depth}`);}catch{/* Older servers still provide trace. */}}
      if(!trace)trace=await api(`/api/trace?id=${encodeURIComponent(root.id)}&depth=${depth}`);
      if(!isCurrent())return;
      renderGraph({root:trace.root||root,nodes:[root,...asArray(trace.nodes)],edges:[...asArray(trace.edges),...asArray(trace.references)],truncated:trace.truncated,node_limit:trace.node_limit,depth:trace.depth,counts:trace.counts,reason:trace.reason,truncation_reasons:trace.truncation_reasons,affected_derivations:trace.affected_derivations,traversal_complete:trace.traversal_complete,all_premises_complete:trace.all_premises_complete});
    } catch {if(!isCurrent())return;const connected=new Set(asArray(fallbackEdges).flatMap(e=>[e.source||e.from||e.from_id||e.src,e.target||e.to||e.to_id||e.dst]).filter(Boolean));renderGraph({root,nodes:[root,...related.filter(n=>connected.has(n.id))],edges:fallbackEdges}); }
  }
  function renderGraph(data) {
    const holder = $('#tab-content'); if (!holder) return;
    const rawNodes=asArray(data.nodes).filter(n=>n?.id),root=data.root||rawNodes[0];
    const rootNode=typeof root==='object'?root:rawNodes.find(n=>n.id===root)||rawNodes[0];
    const audits=asArray(semantic(rootNode).input_audit||attrs(rootNode).input_audit);
    const auditParentIds=[...new Set(audits.map(item=>normalizeGraphId(item.upstream_claim_id,'claim',rootNode)).filter(Boolean))];
    const parentNodes=auditParentIds.map(id=>rawNodes.find(n=>n.id===id)||{id,type:'Claim',label:`Aussage ${id}`,attrs:{}});
    const auditNodeId=(item,index)=>item.graph_id||item.id||`audit-input:${item.upstream_claim_id||rootNode?.id||'node'}:${item.name||index}`;
    const auditTargets=audits.map((item,index)=>{const id=auditNodeId(item,index),name=String(item.name||'').trim().toLowerCase();const existing=rawNodes.find(n=>n.id===id)||rawNodes.find(n=>/assumption|input|premise/i.test(typeName(n))&&name&&String(label(n)).toLowerCase().includes(name));return existing?.id||id;});
    const auditNodes=audits.map((item,index)=>{const id=auditTargets[index],existing=rawNodes.find(n=>n.id===id),auditAttrs={origin:item.origin,impact:item.impact,upstream_claim_id:item.upstream_claim_id,audit_target:item.upstream_claim_id||rootNode?.id};return existing?{...existing,attrs:{...attrs(existing),...auditAttrs}}:{id,type:'Eingabe',label:`${item.name||'Eingabe'}${item.value!=null?`: ${typeof item.value==='string'?item.value:JSON.stringify(item.value)}`:''}`,attrs:auditAttrs};});
    const allNodes=[rootNode,...parentNodes,...auditNodes,...rawNodes].filter(Boolean).filter((node,index,list)=>list.findIndex(other=>other.id===node.id)===index);
    const nodeMap=new Map(allNodes.map(n=>[n.id,n]));
    const coreTypes=new Set(['requires','required_together','derives','has_result','has_derivation','implements','checks','formalization_candidate','uses_claim_under_scope','sets_input','interprets_as']);
    const relationForOrigin={upstream_result:'requires',declared_axiom:'sets_input',external_data:'sets_input',hardcoded:'sets_input',free_parameter:'sets_input',derived_in_script:'interprets_as'};
    const auditEdges=audits.map((item,index)=>({source:normalizeGraphId(item.upstream_claim_id,'claim',rootNode)||rootNode?.id,target:auditTargets[index],type:relationForOrigin[item.origin]})).filter(e=>e.type);
    const allEdges=[...auditEdges,...asArray(data.edges)].map(e=>({source:e.source||e.from||e.from_id||e.src,target:e.target||e.to||e.to_id||e.dst,type:String(e.type||e.relation||e.label||'')})).filter(e=>nodeMap.has(e.source)&&nodeMap.has(e.target)&&!/topic/i.test(e.type)&&!/(topic|subject|theme)/i.test(typeName(nodeMap.get(e.source)))&&!/(topic|subject|theme)/i.test(typeName(nodeMap.get(e.target)))).filter((edge,index,all)=>all.findIndex(other=>other.source===edge.source&&other.target===edge.target&&other.type===edge.type)===index);
    const isCross=e=>!coreTypes.has(e.type.toLowerCase());
    const coreEdges=allEdges.filter(e=>!isCross(e)),proofTypes=new Set(['required_together','derives']),proofEdges=coreEdges.filter(e=>proofTypes.has(e.type.toLowerCase()));
    const levels=new Map(allNodes.map(n=>[n.id,0]));
    const rank=n=>{const t=String(typeName(n)).toLowerCase();if(/claim/.test(t))return 0;if(/premise|assumption|input|source/.test(t))return 1;if(/derivation|script/.test(t))return 2;if(/result|lean/.test(t))return 3;return 1;};
    for(const n of allNodes)levels.set(n.id,rank(n));
    const indegree=new Map(allNodes.map(n=>[n.id,0]));for(const e of proofEdges)if(e.source!==e.target)indegree.set(e.target,(indegree.get(e.target)||0)+1);
    const queue=allNodes.filter(n=>(indegree.get(n.id)||0)===0).map(n=>n.id),ordered=new Set();
    while(queue.length){const id=queue.shift();ordered.add(id);for(const e of proofEdges.filter(edge=>edge.source===id&&edge.target!==id)){levels.set(e.target,Math.min(3,Math.max(levels.get(e.target)||0,(levels.get(id)||0)+1)));indegree.set(e.target,indegree.get(e.target)-1);if(indegree.get(e.target)===0)queue.push(e.target);}}
    // Cyclic leftovers use their semantic type layer; cycles cannot stretch the canvas.
    for(const n of allNodes)if(!ordered.has(n.id))levels.set(n.id,rank(n));
    const columns=[...levels.values()].length?Math.max(...levels.values())+1:1,groups=Array.from({length:columns},()=>[]);
    for(const n of allNodes)groups[Math.min(levels.get(n.id)||0,columns-1)].push(n);
    const cardW=258,cardH=112,colGap=92,rowGap=34,pad=38;
    const width=Math.max(760,pad*2+columns*cardW+(columns-1)*colGap),height=Math.max(250,...groups.map(g=>pad*2+g.length*cardH+Math.max(0,g.length-1)*rowGap));
    const places=new Map();groups.forEach((group,col)=>group.forEach((n,row)=>places.set(n.id,[pad+cardW/2+col*(cardW+colGap),pad+cardH/2+row*(cardH+rowGap)])));
    let svg=`<svg class="graph-canvas" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" role="img" aria-label="${esc(tx('Gerichtete Lösungskette'))}">`;
    svg+='<defs><marker id="chain-arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="#74938a"/></marker></defs>';
    const structuralTypes=new Set(['requires','has_result','has_derivation']);
    for(const e of allEdges){const [x1,y1]=places.get(e.source),[x2,y2]=places.get(e.target);const core=!isCross(e),structural=structuralTypes.has(e.type.toLowerCase()),d=`M ${x1+cardW/2-3} ${y1} C ${x1+cardW/2+colGap/2} ${y1}, ${x2-cardW/2-colGap/2} ${y2}, ${x2-cardW/2+3} ${y2}`,labelText=tx(edgeLabel(e.type));svg+=`<path class="chain-edge ${structural?'structural-edge':''} ${core?'':'cross-edge'}" d="${d}" marker-end="url(#chain-arrow)"><title>${esc(labelText)}</title></path><text class="graph-edge-label ${structural?'structural-label':''} ${core?'':'cross-label'}" x="${(x1+x2)/2}" y="${(y1+y2)/2-7}">${esc(labelText)}</text>`;}
    for(let col=0;col<groups.length;col++)for(const n of groups[col]){const [x,y]=places.get(n.id),a=attrs(n),assessment=a.assessment||{},category=assessment.category||a.ledger_status||'',statusClass=statusClassFor(category,isUnreviewed(n)?'ai_unreviewed':assessment.authority),isInput=typeLabel(n)==='Eingabe',isPremise=/premise|assumption/i.test(typeName(n)),rawTitle=String(label(n)),localId=String(a.id||''),detailTitle=(english()&&(a.statement_en||a.text_en||a.description_en||a.translations?.en?.statement||a.translations?.en?.text))||a.statement||a.text||a.description||summary(n)||'',genericLabel=/^(?:P\d+|A\d+|premise|assumption)$/i.test(rawTitle)||rawTitle===localId||rawTitle===n.id,owner=a.priority_owner||a.owner||a.origin_claim_id||rootNode?.id||'',title=String(isPremise&&genericLabel?(detailTitle||`${tx('Voraussetzung')} · ${owner||n.id}`):rawTitle),lines=wrapGraphText(title,32,3),targetId=a.audit_target||a.upstream_claim_id||n.id,origin=a.origin?` · ${tx(inputOrigin[a.origin]?.[0]||a.origin)}`:'',badgeText=isUnreviewed(n)?tx('KI-Auswertung · ungeprüft'):categoryLabel(category)||String(category).replaceAll('_',' '),scope=(english()&&(a.scope_en||assessment.scope_en))||a.scope||assessment.scope||'',ownerNote=owner?` · ${tx('Zugehörige Aussage')}: ${owner}`:'',localNote=isPremise&&localId?` · ${localId}`:'';svg+=`<g class="chain-node ${statusClass} ${n.id===rootNode?.id?'main':''} ${isInput?'input-node':''}" transform="translate(${x-cardW/2},${y-cardH/2})" data-node="${esc(targetId)}" tabindex="0" role="button" aria-label="${esc(`${title}${localNote}${scope?` · ${scope}`:''}`)}"><title>${esc(`${title}${localNote}${scope?` · ${scope}`:assessment.notice?` · ${assessment.notice}`:''}${ownerNote}`)}</title><rect width="${cardW}" height="${cardH}" rx="12"/><text class="chain-type" x="14" y="20">${esc(typeLabel(n))}${localNote?esc(localNote):''}${origin?esc(origin):''}</text><text class="chain-assessment" x="${cardW-12}" y="20" text-anchor="end">${esc(badgeText)}</text>${lines.map((line,i)=>`<text class="chain-title" x="14" y="43" dominant-baseline="hanging" transform="translate(0 ${i*17})">${esc(line)}</text>`).join('')}${scope?`<text class="chain-scope" x="14" y="${cardH-12}">${esc(excerpt(scope,46))}</text>`:''}${/derivation/i.test(typeName(n))?`<text class="conjunction" x="${cardW-12}" y="${cardH-12}" text-anchor="end">${esc(tx('Alle Voraussetzungen'))}</text>`:''}</g>`;}
    svg+='</svg>';
    const crossCount=allEdges.filter(isCross).length,reasonText=[data.reason,...asArray(data.truncation_reasons)].filter(Boolean).join(' · '),affected=asArray(data.affected_derivations).filter(item=>item.incomplete||item.truncated).map(item=>item.label||item.id).filter(Boolean),truncated=data.truncated?`<div class="graph-limit"><span>${esc(tx('Begrenzter Ausschnitt'))}: ${allNodes.length} ${esc(tx('Einträge'))}, ${esc(tx('Tiefe'))} ${esc(data.depth||4)}. ${esc(reasonText||tx('Die Anzeige ist begrenzt; der vollständige Umfang ist hier nicht bestimmt.'))}${data.traversal_complete===false?` ${esc(tx('Einige Voraussetzungen und Herleitungsschritte sind hier nicht vollständig angezeigt. Daraus folgt keine Aussage über ihre Beweiskraft.'))}`:''}${affected.length?` ${esc(tx('Betroffen'))}: ${esc(affected.slice(0,5).join(', '))}${affected.length>5?' …':''}.`:''}</span> <button type="button" class="mini-button" data-expand-chain>${esc(tx('Weitere Stufe laden'))}</button></div>`:'';
    const rootScope=(english()&&(attrs(rootNode).assessment?.scope_en||attrs(rootNode).scope_en))||attrs(rootNode).assessment?.scope||attrs(rootNode).scope||'';
    const chainNeedsLanguageNote=english()&&allNodes.some(n=>{const a=attrs(n),s=semantic(n),t=a.translations?.en||s.translations?.en||{};return !(isClaimNode(n)?freshEnglish(n,'title'):(n.label_en||a.label_en||a.title_en||s.title_en||t.title))&&!(isClaimNode(n)?freshEnglish(n,'overview'):(n.summary_en||a.summary_en||a.overview_en||s.overview_en||t.overview));});
    const chainLanguageNote=chainNeedsLanguageNote&&!document.querySelector('.dossier-top .language-fallback');
    holder.innerHTML=`<header class="chain-root"><span>${esc(tx('Ausgewählte Aussage'))}</span><h3>${esc(label(rootNode))}</h3>${rootScope?`<p><b>${esc(tx('Geltungsbereich'))}:</b> ${esc(rootScope)}</p>`:`<p class="scope-missing">${esc(tx('Geltungsbereich nicht angegeben'))}</p>`}</header><div class="chain-toolbar"><strong>${allNodes.length} ${esc(tx('Einträge'))} · ${proofEdges.length} ${esc(tx('Herleitungsschritte'))}</strong><span><button type="button" data-zoom="-1" aria-label="${esc(tx('Verkleinern'))}">−</button><button type="button" data-zoom="1" aria-label="${esc(tx('Vergrößern'))}">+</button><button type="button" data-fit>${esc(tx('Einpassen'))}</button><label><input type="checkbox" data-cross> ${esc(tx('Querverweise (keine Herleitung)'))} · ${crossCount}</label></span></div><div class="chain-legend"><span class="legend-established">${esc(tx('Belegt im angegebenen Rahmen'))}</span><span class="legend-conditional">${esc(tx('Bedingt'))}</span><span class="legend-open">${esc(tx('Offen'))}</span><span class="legend-unreviewed">${esc(tx('KI-Auswertung · ungeprüft'))}</span></div>${chainLanguageNote?`<p class="language-fallback">${esc(tx('Detaillierte Erläuterungen bleiben mangels Übersetzung in ihrer Originalsprache sichtbar.'))}</p>`:''}${truncated}<div class="graph-viewport"><div class="graph-world">${svg}</div></div><p class="graph-caption">${esc(tx('Pfeile zeigen die tatsächliche Richtung. „Benötigt als Voraussetzung“, „hat Ergebnis“ und „hat Herleitung“ beschreiben Strukturbezüge; nur „gemeinsam erforderlich“ und „leitet her“ zählen hier als Herleitungsschritte. Knoten öffnen ihr Dossier.'))}</p>`;
    const viewport=$('.graph-viewport',holder),world=$('.graph-world',holder);let zoom=1;
    const applyZoom=()=>{world.style.transform=`scale(${zoom})`;world.style.transformOrigin='top left';};
    applyZoom();holder.querySelectorAll('[data-zoom]').forEach(btn=>btn.addEventListener('click',()=>{zoom=Math.max(.55,Math.min(1.7,zoom+Number(btn.dataset.zoom)*.15));applyZoom();}));
    $('[data-fit]',holder)?.addEventListener('click',()=>{zoom=Math.min(1,(viewport.clientWidth-16)/width);applyZoom();viewport.scrollTo({left:0,top:0,behavior:'smooth'});});
    $('[data-cross]',holder)?.addEventListener('change',e=>holder.classList.toggle('show-cross',e.target.checked));
    $('[data-expand-chain]',holder)?.addEventListener('click',()=>loadTrace(rootNode,[],[],Math.min(16,(data.depth||4)+2)));
    holder.querySelectorAll('.chain-node[data-node]').forEach(el=>{el.addEventListener('click',()=>goNode(el.dataset.node));el.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();goNode(el.dataset.node)}})});
  }
  async function loadNode(id,openGraph=false) {
    if(String(id).startsWith('topic:')){loadTopic(id);return;}
    if(!isCurrentNodeRoute(id))return;
    const n=++state.serial; setNotice(''); view.innerHTML='<div class="loading-block"><span class="spinner"></span><span>Eintrag und Verknüpfungen werden geladen …</span></div>';
    void refreshStatus();
    try {
      const data=await api(`/api/node?id=${encodeURIComponent(id)}`); if(n!==state.serial)return;
      const node=data.node||data;
      if(String(node.type||'').startsWith('research_')){
        if(!isCurrentNodeRoute(id))return;
        const questions=asArray(data.related_nodes).filter(item=>item.type==='research_problem');
        const target=node.type==='research_problem'?id:questions.length===1?questions[0].id:id;
        history.replaceState(null,'',`#research/${encodeURIComponent(target)}`);handleHash();return;
      }
      if(/claim|script/i.test(typeName(node))){try{data.audit=await api(`/api/audit?id=${encodeURIComponent(id)}`);}catch{/* Audit annotations are optional when absent for this node. */}}
      if(/claim/i.test(typeName(node))){try{data.claim_attempts=await api(`/api/attempts?id=${encodeURIComponent(id)}&limit=100&offset=0`);}catch{/* Attempt inventory may be unavailable for this claim. */}}
      if(n!==state.serial||!isCurrentNodeRoute(id))return; state.current=data; renderDossier(data);if(openGraph)setTimeout(()=>{if(isCurrentNodeRoute(id))view.querySelector('[data-tab="graph"]')?.click()},0);
    }
    catch(error) {if(n!==state.serial||!isCurrentNodeRoute(id))return; view.innerHTML=`<div class="error-state">${esc(error.message)} <button class="mini-button" id="retry-node">Erneut laden</button></div>`; $('#retry-node')?.addEventListener('click',()=>loadNode(id));}
  }
  async function loadSource(id, trigger, range={}) {
    if(!id)return;
    const container=$('#source-viewer'); if(!container)return;
    container.innerHTML=`<div class="loading-block"><span class="spinner"></span><span>${esc(tx('Quelle wird geladen …'))}</span></div>`;
    try {
      const params=new URLSearchParams({id});
      if(range.line_start!=null)params.set('line_start',range.line_start);
      if(range.line_end!=null)params.set('line_end',range.line_end);
      const source=await api(`/api/source?${params}`);
      container.innerHTML=`<section class="source-viewer"><div class="source-viewer-head"><div><h3>${esc(sourceLabel(source))}</h3><div class="dossier-id">${esc(source.path||'')}${source.page!=null?` · ${esc(tx('PDF-Seite'))} ${esc(source.page)}`:''}${source.line_start!=null?` · ${esc(tx('Zeilen'))} ${esc(source.line_start)}–${esc(source.line_end??source.line_start)}`:''}</div></div><div class="source-actions"><button class="mini-button" id="copy-source-link">${esc(tx('Quellenlink kopieren'))}</button></div></div><pre></pre></section>`;
      $('pre',container).textContent=cleanText(source.text||'Für diese Quelle wurde kein Text zurückgegeben.');
      $('#copy-source-link',container).addEventListener('click',async()=>{const url=`${location.origin}${location.pathname}${sourceRouteHref(id,{line_start:range.line_start??source.line_start,line_end:range.line_end??source.line_end})}`;try{await navigator.clipboard.writeText(url);setNotice(tx('Quellenlink wurde kopiert.'));}catch{setNotice(url);}});
      container.scrollIntoView({behavior:'smooth',block:'start'});
    } catch(error) {container.textContent=error.message;}
  }
  function renderResults(query, data) {
    const results=data.results||[]; state.results=results;
    setSidebarResults(query,results,data.count);
    const head=document.createElement('div'); head.className='topic-result-head';
    const text=document.createElement('div'); text.innerHTML=`<h2>${esc(tx('Suchergebnisse'))}</h2><p></p>`; $('p',text).textContent=`${data.count ?? results.length} ${tx('Treffer für')} „${query}“ · ${tx('Nach Relevanz sortiert')}`;
    const count=document.createElement('span'); count.className='tiny-label'; count.textContent=tx('Originalstatus bleibt maßgeblich'); head.append(text,count);
    view.replaceChildren(head);
    if(!results.length){const empty=document.createElement('div');empty.className='empty-state';empty.textContent=tx('Keine passenden Einträge gefunden. Probiere einen anderen Begriff oder ein Themenfeld.');view.append(empty);return;}
    const list=document.createElement('div');list.className='result-list';
    for(const node of results){
      const item=document.createElement('div');item.className='search-result-item';
      const button=document.createElement('button');button.type='button';button.className='search-result';
      const line=document.createElement('div');line.className='result-topline';
      const type=document.createElement('span');type.className='result-type';type.textContent=typeLabel(node);
      const meta=document.createElement('span');meta.className='result-meta';
      if(ledgerStatus(node)) {const b=document.createElement('span');b.className=`badge ${statusKind(node)}`;b.textContent=`${tx('Ledger')}: ${ledgerStatus(node)}`;meta.append(b);}
      const enrichmentStatus=enrichment(node);if(enrichmentStatus){const b=document.createElement('span');b.className=`badge ${enrichmentStatus.className}`;b.textContent=enrichmentStatus.label;meta.append(b);}
      line.append(type,meta);
      const title=document.createElement('h3');title.className='result-title';title.textContent=excerpt(label(node),220);
      const desc=document.createElement('p');desc.className='result-summary';desc.dataset.originalContent='';desc.textContent=excerpt(summary(node),280);
      const id=document.createElement('div');id.className='result-id';id.textContent=node.id;
      button.append(line,title,desc,id);button.addEventListener('click',()=>goNode(node.id));item.append(button);
      const sourceIds=Array.isArray(node.source_ids)?node.source_ids:Array.isArray(node.source_refs)?node.source_refs:[];
      if(sourceIds.length){const sources=document.createElement('div');sources.className='search-source-links';sourceIds.slice(0,3).forEach((ref,index)=>{const sourceId=typeof ref==='string'?ref:ref?.source_id;if(!sourceId)return;const sourceButton=document.createElement('button');sourceButton.type='button';sourceButton.className='source-chip';sourceButton.textContent=`${tx('Quelle')} ${index+1}`;sourceButton.title=typeof ref==='object'?(ref.path||sourceId):sourceId;sourceButton.addEventListener('click',event=>{event.stopPropagation();location.hash=sourceRouteHref(sourceId,typeof ref==='object'?ref:{})});sources.append(sourceButton)});if(sources.childElementCount)item.append(sources)}
      list.append(item);
    }
    view.append(list);
  }
  async function search(query, type='') {
    query=query.trim();if(!query)return;
    state.searchQuery=query;document.body.classList.remove('research-mode');
    const mode=document.querySelector('[name="search-mode"][value="entries"]');if(mode)mode.checked=true;
    searchInput.value=query; const n=++state.serial;setNotice('');view.innerHTML=`<div class="loading-block"><span class="spinner"></span><span>${esc(tx('Der Wissensgraph wird durchsucht …'))}</span></div>`;
    const url=new URL('/api/search',location.origin);url.searchParams.set('q',query);url.searchParams.set('limit','40');if(type)url.searchParams.set('type',type);
    try {const data=await api(url.pathname+url.search);if(n!==state.serial)return;renderResults(query,data);}
    catch(error){if(n!==state.serial)return;view.innerHTML=`<div class="error-state">${esc(error.message)} <button class="mini-button" id="retry-search">Erneut suchen</button></div>`;$('#retry-search')?.addEventListener('click',()=>search(query,type));}
  }
  function showHomeData() {
    const v=document.createElement('div');
    const head=document.createElement('div');head.className='welcome-head';head.innerHTML='<div><h2>Wo möchtest du anfangen?</h2><p>Wähle ein Themenfeld oder suche direkt nach einer Aussage.</p></div><span class="tiny-label">THEMENLANDKARTE</span>';v.append(head);
    const grid=document.createElement('div');grid.className='topic-grid';
    for(const t of state.topics){const card=document.createElement('button');card.type='button';card.className='topic-card';
      const glyph=document.createElement('span');glyph.className='topic-card-glyph';glyph.textContent=t.glyph;
      const name=document.createElement('strong');name.textContent=topicName(t.id);
      const desc=document.createElement('p');const topicDescription=topicCopy[t.id]||summary(t.node)||'Aussagen, Herleitungen und offene Fragen zu diesem Feld.';desc.textContent=topicDescription;if(!topicCopy[t.id])desc.dataset.originalContent='';
      const count=document.createElement('span');count.className='card-count';count.textContent=t.count?`${t.count} Einträge`:'Themenfeld öffnen';
      const arrow=document.createElement('span');arrow.className='arrow';arrow.textContent='↗';card.append(glyph,name,desc,count,arrow);card.addEventListener('click',()=>goTopic(t.id));grid.append(card);}
    v.append(grid);
    const orientation=document.createElement('section');orientation.className='orientation-panel';
    orientation.innerHTML=`<div class="orientation-copy"><p class="eyebrow">${esc(tx('Orientierung'))}</p><h3>${esc(tx('Vom Compiler zur Physik'))}</h3><p>${esc(tx('Der Compiler bündelt endliche Relationen, Codes und Prüfungen. Eine vollständige physikalische Theorie erfordert zusätzlich die Herleitung von Raumzeit, Dynamik und beobachtbaren Größen aus den Ausgangsprinzipien.'))}</p><p>${esc(tx('Resultate, die gemeinsame Seeds oder dieselben Grunddaten verwenden, sind korreliert und zählen nicht automatisch als unabhängige Bestätigungen.'))}</p><a class="orientation-link" href="#orientation">${esc(tx('Aufbau und Ziele ansehen'))} <span aria-hidden="true">↗</span></a><a class="map-home-link" href="#map">${esc(tx('Gesamtkarte ansehen'))} <span aria-hidden="true">↗</span></a></div><button type="button" class="universal-card"><span class="universal-symbol">∞</span><span><strong>${esc(topicName('topic:universalraum'))}</strong><small>${esc(tx('Zum Themenfeld'))}</small></span><span class="universal-arrow">↗</span></button>`;
    orientation.querySelector('.universal-card').addEventListener('click',()=>goTopic('topic:universalraum'));
    v.append(orientation);
    const recent=state.results.slice(0,4);if(recent.length){const section=document.createElement('section');section.className='recent-section';section.innerHTML='<div class="section-heading"><h3>Zuletzt angesehen</h3><span>aus dieser Sitzung</span></div>';const list=document.createElement('div');list.className='recent-list';for(const node of recent){const row=document.createElement('button');row.type='button';row.className='recent-row';row.innerHTML=`<span class="recent-kind">${esc(typeLabel(node))}</span><span class="recent-copy"><strong></strong><small></small></span><span class="recent-arrow">↗</span>`;$('.recent-copy strong',row).textContent=excerpt(label(node),110);$('.recent-copy small',row).textContent=excerpt(summary(node),95);row.addEventListener('click',()=>goNode(node.id));list.append(row);}section.append(list);v.append(section);}
    view.replaceChildren(v);
  }
  async function showOrientation(){
    const request=++state.serial;setNotice('');updateSidebarTopics();
    view.innerHTML='<div class="loading-block"><span class="spinner"></span><span>Theorieaufbau und Ziele werden geladen …</span></div>';
    const ids=['topic:origin',...Array.from({length:8},(_,i)=>`goal:T${i+1}`)];
    const loaded=await Promise.all(ids.map(async id=>{try{return await api(`/api/node?id=${encodeURIComponent(id)}`)}catch{return null}}));
    if(request!==state.serial)return;
    const originData=loaded[0],origin=originData?.node||originData;
    const results=document.createElement('div');results.className='orientation-page';
    const intro=document.createElement('header');intro.className='orientation-page-head';
    intro.innerHTML='<p class="eyebrow">THEORIEAUFBAU</p><h2>Aufbau verstehen</h2><p>Beginne bei den Ausgangsprinzipien und folge dann den acht offenen Zielrichtungen. Der Graph hält vorhandene Resultate und noch fehlende Herleitungen getrennt.</p>';
    const originCard=document.createElement('button');originCard.type='button';originCard.className='origin-start-card';
    const glyph=document.createElement('span');glyph.className='topic-card-glyph';glyph.textContent='◉';const copy=document.createElement('span');copy.className='origin-start-copy';const title=document.createElement('strong');title.textContent=origin?label(origin):'Ursprung';const desc=document.createElement('small');desc.textContent=origin?summary(origin):'Ausgangsprinzipien und die Frage, was die Theorie eindeutig festlegt.';copy.append(title,desc);const arrow=document.createElement('span');arrow.textContent='↗';originCard.append(glyph,copy,arrow);originCard.addEventListener('click',()=>goTopic('topic:origin'));
    const goalsHeading=document.createElement('div');goalsHeading.className='section-heading goals-heading';goalsHeading.innerHTML='<h3>Zielrichtungen</h3><span>Goal-Knoten aus dem Theoriegraphen</span>';
    const goals=document.createElement('div');goals.className='goal-grid';
    for(let i=0;i<8;i++){
      const data=loaded[i+1],node=data?.node||data,id=`goal:T${i+1}`;
      const card=document.createElement('button');card.type='button';card.className='goal-card';
      const number=document.createElement('span');number.className='goal-number';number.textContent=`T${i+1}`;
      const body=document.createElement('span');body.className='goal-copy';
      const name=document.createElement('strong');name.textContent=node?label(node):id;
      const desc=document.createElement('small');desc.textContent=node?summary(node):'Zieldossier wird im Graphen geöffnet.';
      body.append(name,desc);const go=document.createElement('span');go.className='goal-arrow';go.textContent='↗';card.append(number,body,go);card.addEventListener('click',()=>goNode(node?.id||id));goals.append(card);
    }
    results.append(intro,originCard,goalsHeading,goals);view.replaceChildren(results);
  }
  function renderTopicDashboard(data,node){
    const root=data.topic||node,topicId=root?.id||node?.id,summaryData=data.summary||{};
    const orientation=attrs(root).topic_orientation||semantic(root).topic_orientation||{};
    const ui=(en,de)=>english()?en:de;
    const allTabs=[['solved','Belegte Resultate','Recorded results'],['mixed','Teilresultate / gemischter Status','Partial / mixed status'],['conditional','Bedingte Anschlüsse','Conditional connections'],['open','Offene Fragen','Open questions'],['excluded','Ausgeschlossen','Excluded'],['assumed','Angenommen','Assumed'],['historical','Historisch','Historical'],['methodological_rule','Geltende Prüfregel','Active methodological rule'],['unknown','Status nicht automatisch zugeordnet','Status not automatically classified']];
    const categories={solved:data.solved||[],mixed:data.mixed||[],conditional:data.conditional||[],open:data.open||[],excluded:data.excluded||[],assumed:data.assumed||[],historical:data.historical||[],methodological_rule:data.methodological_rule||[],unknown:data.unknown||[],attempts:data.attempts||[],structures:data.structures||[],goals:data.goals||[]};
    categories.solved=asArray(categories.solved).filter(entry=>{const a=attrs(entry),assessment=a.assessment||{},state=String(a.enrichment_state||'').toLowerCase();return assessment.authority==='ledger'||(assessment.authority!=='ai_unreviewed'&&!['pending','error'].includes(state));});
    const categoryCount=key=>summaryData[key]??summaryData[key==='solved'?'recorded_established':key]??data.totals?.[key]??asArray(categories[key]).length;
    const orientationValue=key=>english()?(orientation[`${key}_en`]||''):(orientation[`${key}_de`]||orientation[key]||'');
    const host=document.createElement('article');host.className='topic-dashboard reading-dashboard';
    const head=document.createElement('header');head.className='topic-dashboard-head reading-topic-head';
    const title=document.createElement('h2');title.textContent=label(root);
    const purpose=document.createElement('p');purpose.className='reading-purpose';purpose.textContent=orientationValue('role')||(english()?(attrs(root).summary_en||semantic(root).summary_en||''):(summary(root)||summaryData.description||''));
    head.append(title,purpose);
    const researchLink=document.createElement('a');researchLink.className='topic-research-link';researchLink.href=`#research/${encodeURIComponent(topicId)}`;researchLink.textContent=ui('Open reviewed findings and next tasks →','Eingeordnete Befunde und nächste Aufgaben öffnen →');head.append(researchLink);
    const relatedToggle=document.createElement('label');relatedToggle.className='topic-related-toggle';const relatedInput=document.createElement('input');relatedInput.type='checkbox';relatedInput.checked=localStorage.getItem(`tfpt-topic-related:${topicId}`)==='true';const relatedText=document.createElement('span');relatedText.textContent=`${ui('Include related topics','Auch angrenzende Themen zeigen')}${data.related_claims_total?` · ${data.related_claims_total}`:''}`;relatedToggle.append(relatedInput,relatedText);head.append(relatedToggle);
    if(data.membership_notice&&!relatedInput.checked){const membership=document.createElement('small');membership.className='topic-membership-note';membership.textContent=english()?'Showing the core of this topic. Related findings remain available through the switch; their evidence status stays unchanged.':'Hier siehst du den Themenkern. Angrenzende Einträge bleiben über den Schalter erreichbar; ihr Belegstatus bleibt unverändert.';head.append(membership);}
    const story=document.createElement('div');story.className='topic-reading-story';
    const englishNotice=english()&&(!orientation.role_en||(!orientation.have_en&&!!orientation.have_de)||(!orientation.connection_en&&!!orientation.connection_de)||(!orientation.next_bridge_en&&!orientation.seek_en&&(!!orientation.next_bridge_de||!!orientation.seek_de)));
    if(englishNotice){const note=document.createElement('p');note.className='language-fallback';note.textContent='Some topic annotations have no English version; those items remain in the specialist details in their authored language.';story.append(note);}
    const statusText=(entry,key)=>{const a=attrs(entry),assessment=a.assessment||{},category=String(assessment.category||key||'').toLowerCase(),labels={established:['belegt im Rahmen','recorded within scope'],recorded_established:['belegt im Rahmen','recorded within scope'],solved:['belegt im Rahmen','recorded within scope'],mixed:['gemischt','mixed'],conditional:['bedingt','conditional'],open:['offen','open'],excluded:['ausgeschlossen','excluded'],assumed:['angenommen','assumed'],historical:['historisch','historical'],unknown:['Status nicht automatisch zugeordnet','Status not automatically classified'],methodological_rule:['geltende Prüfregel','active methodological rule']};return (english()&&assessment.label_en)||labels[category]?.[english()?1:0]||(english()?assessment.label_en:assessment.label)||ui('Status not automatically classified','Status nicht automatisch zugeordnet');};
    const shortClaim=(entry,key)=>{
      const a=attrs(entry),assessment=a.assessment||{},card=document.createElement('article');card.className=`reading-claim ${statusClassFor(assessment.category||key,isUnreviewed(entry)?'ai_unreviewed':assessment.authority)}`;
      const top=document.createElement('div');top.className='reading-card-head';const name=document.createElement('strong');name.textContent=label(entry);const badge=document.createElement('span');badge.className='reading-status';badge.textContent=statusText(entry,key);top.append(name,badge);card.append(top);
      const englishText=english()&&(isClaimNode(entry)?freshEnglish(entry,'overview'):(entry.summary_en||a.summary_en||a.overview_en||a.statement_en||semantic(entry).overview_en));
      const content=english()?(englishText||summary(entry)||a.statement):(a.statement||semantic(entry).overview_de||a.overview_de||summary(entry));
      if(content){const full=String(content),p=document.createElement('p');p.textContent=full.length>360?excerpt(full,340):full;card.append(p);if(full.length>360){const more=document.createElement('details');more.className='reading-claim-full';const moreSummary=document.createElement('summary');moreSummary.textContent=ui('Full explanation','Vollständige Erläuterung');const fullText=document.createElement('p');fullText.textContent=full;more.append(moreSummary,fullText);card.append(more);}if(english()&&!englishText){const note=document.createElement('small');note.className='language-fallback';note.textContent='No English claim summary is available; the authored text is shown.';card.append(note);}}
      const scope=(english()&&assessment.scope_en)||assessment.scope;if(scope){const p=document.createElement('p');p.className='reading-scope';p.textContent=`${tx('Geltungsbereich')}: ${scope}`;card.append(p);}
      const refs=asArray(a.source_refs||a.sources||assessment.source_refs||semantic(entry).source_refs).slice(0,3);if(refs.length){const sources=document.createElement('div');sources.className='reading-sources';refs.forEach((ref,index)=>{const id=sourceId(ref);if(!id)return;const button=document.createElement('button');button.type='button';button.textContent=`${tx('Quelle')} ${index+1}`;button.title=sourceLabel(ref);button.addEventListener('click',()=>loadSource(id,button,{line_start:ref?.line_start,line_end:ref?.line_end}));sources.append(button);});card.append(sources);}
      const open=document.createElement('button');open.type='button';open.className='reading-open-link';open.textContent=/claim/i.test(typeName(entry))?ui('Open claim','Aussage öffnen'):ui('Open details','Details öffnen');open.addEventListener('click',()=>goNode(entry.id));card.append(open);return card;
    };
    const addStoryPart=(titleText,copy,entries,key)=>{
      const part=document.createElement('section');part.className='reading-story-part';const heading=document.createElement('h3');heading.textContent=titleText;part.append(heading);
      if(copy){const paragraph=document.createElement('p');paragraph.className='reading-story-copy';paragraph.textContent=readableValue(copy);part.append(paragraph);}
      const items=asArray(entries).slice(0,3);if(items.length){const list=document.createElement('div');list.className='reading-claim-list';items.forEach(item=>list.append(shortClaim(item,key)));part.append(list);}
      if(!copy&&!items.length){const empty=document.createElement('p');empty.className='reading-empty';empty.textContent=ui('No entry is assigned here yet.','Hier ist noch kein Eintrag zugeordnet.');part.append(empty);}
      story.append(part);
    };
    const have=orientationValue('have')||(english()?'':semantic(root).what_we_have||attrs(root).what_we_have||'');
    const conditionalEntries=[...asArray(categories.conditional),...asArray(categories.mixed)].filter((item,index,all)=>all.findIndex(other=>other.id===item.id)===index);
    const nextEntries=asArray(categories.open);
    const nextQuestion=orientationValue('next_bridge')||orientationValue('seek')||(english()?'':semantic(root).what_we_seek||attrs(root).what_we_seek||'');
    addStoryPart(ui('What is recorded?','Was ist festgehalten?'),have,categories.solved,'solved');
    addStoryPart(ui('What remains conditional or mixed?','Was gilt nur bedingt oder bleibt gemischt?'),orientationValue('connection'),conditionalEntries,'conditional');
    addStoryPart(ui('Next specific question','Nächste konkrete Frage'),nextQuestion,nextEntries,'open');
    if(data.membership_notice){const note=document.createElement('p');note.className='reading-caveat';note.textContent=ui('Topic links organize material; they are not evidence by themselves.','Die Zuordnung zu einem Thema ordnet Material ein; sie ist selbst kein Beweis.');story.append(note);}

    const attemptSection=document.createElement('section');attemptSection.className='reading-attempts';
    const attempts=asArray(categories.attempts),attemptTotal=Number(data.totals?.attempts??summaryData.attempts??attempts.length);let shown=0;
    const attemptHead=document.createElement('div');attemptHead.className='reading-attempts-head';const attemptTitle=document.createElement('h3');attemptTitle.textContent=ui('Attempts and findings','Versuche und Befunde');const attemptCount=document.createElement('span');attemptCount.textContent=`${attemptTotal}`;attemptHead.append(attemptTitle,attemptCount);attemptSection.append(attemptHead);
    const attemptNote=document.createElement('p');attemptNote.textContent=ui('Indexed attempts and findings are not necessarily independent experiments.','Verzeichnete Versuche und Befunde sind nicht automatisch unabhängige Experimente.');attemptNote.className='reading-caveat';attemptSection.append(attemptNote);
    const attemptList=document.createElement('div');attemptList.className='reading-attempt-list';attemptSection.append(attemptList);
    const attemptCard=entry=>{const a=attrs(entry),contract=a.original_contract||{},card=document.createElement('article');card.className='reading-attempt-card';const title=document.createElement('strong');title.textContent=label(entry);const raw=String(a.verdict_raw||a.verdict||contract.verdict||contract.status||'').trim();const status=document.createElement('span');status.className='reading-status';status.textContent=raw||(ui('Finding not classified','Befund nicht klassifiziert'));const kind=document.createElement('small');kind.className='reading-type';kind.textContent=ui('Attempt / finding','Versuch / Befund');const header=document.createElement('div');header.className='reading-card-head';header.append(title,status);card.append(header,kind);const summaryText=summary(entry);if(summaryText){const p=document.createElement('p');p.textContent=summaryText.length>360?excerpt(summaryText,340):summaryText;card.append(p);if(summaryText.length>360){const more=document.createElement('details');const moreSummary=document.createElement('summary');moreSummary.textContent=ui('Full attempt details','Vollständige Versuchsdaten');const full=document.createElement('p');full.textContent=summaryText;more.append(moreSummary,full);card.append(more);}if(english()&&!(entry.summary_en||entry.overview_en||a.summary_en||a.overview_en||semantic(entry).overview_en)){const note=document.createElement('small');note.className='language-fallback';note.textContent='No English summary is available; the authored text is shown.';card.append(note);}}
      const findings=asArray(contract.findings||contract.conditional_results||a.findings||a.positive_findings);const negative=asArray(contract.rules_out||contract.missing_or_open||a.negative_findings);const finding=findings[0]||negative[0];if(finding&&(!summaryText||['PARTIAL','PASS','FAIL','OPEN'].includes(summaryText.trim().toUpperCase()))){const detail=document.createElement('p');detail.className='reading-finding';detail.textContent=`${english()?(findings.length?'Recorded finding':'Recorded limit'):(findings.length?'Erfasster Befund':'Erfasste Grenze')}: ${readableValue(finding)}`;card.append(detail);}
      const basis=String(a.topic_link_basis||'').toLowerCase();if(basis){const link=document.createElement('small');link.textContent=/explicit|direct|topic_link/.test(basis)?ui('Direct topic link','Direkt mit dem Thema verknüpft'):/inherit|linked_claim|derived/.test(basis)?ui('Inherited through a linked claim','Über verknüpfte Aussage zugeordnet'):'';if(link.textContent)card.append(link);}
      const open=document.createElement('button');open.type='button';open.className='reading-open-link';open.textContent=ui('Open attempt dossier','Versuchsdossier öffnen');open.addEventListener('click',()=>goNode(entry.id));card.append(open);return card;};
    const showAttempts=()=>{const until=Math.min(attempts.length,shown+3);for(;shown<until;shown++)attemptList.append(attemptCard(attempts[shown]));attemptControls();};
    let loadingAttempts=false;
    const attemptControls=()=>{attemptSection.querySelector('.reading-attempt-controls')?.remove();if(shown>=attemptTotal)return;const controls=document.createElement('div');controls.className='reading-attempt-controls';if(shown<attempts.length){const more=document.createElement('button');more.type='button';more.textContent=`${ui('Show more attempts','Weitere Versuche anzeigen')} (${Math.min(3,attempts.length-shown)})`;more.addEventListener('click',showAttempts);controls.append(more);}else if(attempts.length<attemptTotal){const more=document.createElement('button');more.type='button';more.disabled=loadingAttempts;more.textContent=loadingAttempts?ui('Loading…','Wird geladen …'):ui('Load more attempts','Weitere Versuche laden');more.addEventListener('click',async()=>{if(loadingAttempts)return;loadingAttempts=true;attemptControls();try{const extra=await api(`/api/theme?id=${encodeURIComponent(topicId)}&offset=${attempts.length}&limit=100${relatedInput.checked?'&include_related=true':''}`);const known=new Set(attempts.map(item=>item.id));attempts.push(...asArray(extra.attempts).filter(item=>!known.has(item.id)));categories.attempts=attempts;loadingAttempts=false;showAttempts();}catch{loadingAttempts=false;more.disabled=false;more.textContent=ui('Loading failed · retry','Laden fehlgeschlagen · erneut versuchen');controls.append(more);}});controls.append(more);}if(controls.childElementCount)attemptSection.append(controls);};
    showAttempts();

    const info=document.createElement('div');info.className='topic-dashboard-info';
    for(const [key,de,en] of [['have','Was wir haben','What we have'],['seek','Was wir suchen','What we seek'],['required_steps','Erforderliche Schritte','Required steps'],['conditions','Bedingungen','Conditions']]){const value=english()?(orientation[`${key}_en`]||semantic(root)[`${key}_en`]||attrs(root)[`${key}_en`]||''):(orientation[`${key}_de`]||orientation[key]||semantic(root)[`${key}_de`]||semantic(root)[key]||attrs(root)[`${key}_de`]||attrs(root)[key]);if(value){const section=document.createElement('section');const heading=document.createElement('h3');heading.textContent=ui(en,de);const paragraph=document.createElement('div');paragraph.textContent=readableValue(value);section.append(heading,paragraph);info.append(section);}}
    const specialist=document.createElement('details');specialist.className='topic-specialist-details';const specialistSummary=document.createElement('summary');specialistSummary.textContent=ui('All classifications, topic links, and specialist details','Alle Einträge, Themenverknüpfungen und Fachdetails');specialist.append(specialistSummary);
    const tabsEl=document.createElement('div');tabsEl.className='theme-tabs';const panel=document.createElement('div');panel.className='theme-panel';
    const activeTopicChanged=state.activeThemeTopicId!==topicId,hasTab=allTabs.some(([key])=>key===state.activeThemeTab);if(activeTopicChanged||!state.activeThemeTabUserSelected||!hasTab){state.activeThemeTab=allTabs.find(([key])=>categoryCount(key)>0)?.[0]||'solved';state.activeThemeTabUserSelected=false;}state.activeThemeTopicId=topicId;
    const draw=(key)=>{panel.replaceChildren();const entries=asArray(categories[key]),count=categoryCount(key),heading=document.createElement('div');heading.className='theme-panel-heading';const h=document.createElement('h3');h.textContent=english()?allTabs.find(tab=>tab[0]===key)?.[2]||key:allTabs.find(tab=>tab[0]===key)?.[1]||key;const number=document.createElement('span');number.textContent=`${count} ${ui('entries','Einträge')}`;heading.append(h,number);panel.append(heading);if(!entries.length){const empty=document.createElement('p');empty.className='reading-empty';empty.textContent=ui('No entries are recorded in this group yet.','Für diese Gruppe sind derzeit keine Einträge erfasst.');panel.append(empty);}else{const list=document.createElement('div');list.className='theme-entry-list';entries.forEach(entry=>list.append(shortClaim(entry,key)));panel.append(list);}if(entries.length<count){const more=document.createElement('button');more.type='button';more.className='theme-more';more.textContent=ui('Load more entries','Weitere Einträge laden');more.addEventListener('click',async()=>{more.disabled=true;more.textContent=ui('Loading…','Wird geladen …');try{const offset=entries.length,extra=await api(`/api/theme?id=${encodeURIComponent(topicId)}&offset=${offset}&limit=100${relatedInput.checked?'&include_related=true':''}`),known=new Set(entries.map(item=>item.id));categories[key].push(...asArray(extra[key]).filter(item=>!known.has(item.id)));draw(key);}catch{more.disabled=false;more.textContent=ui('Loading failed · retry','Laden fehlgeschlagen · erneut versuchen');}});panel.append(more);}};
    for(const [key,de,en] of allTabs){const button=document.createElement('button');button.type='button';button.className='theme-tab';button.textContent=`${ui(en,de)} · ${categoryCount(key)}`;button.addEventListener('click',()=>{state.activeThemeTab=key;state.activeThemeTabUserSelected=true;tabsEl.querySelectorAll('.theme-tab').forEach(item=>item.classList.toggle('active',item===button));draw(key);});tabsEl.append(button);if(key===state.activeThemeTab)button.classList.add('active');}
    draw(state.activeThemeTab);specialist.append(tabsEl,panel,info,researchContextPanel(topicId));
    const sources=document.createElement('section');sources.className='topic-source-list';const sourceRefs=asArray(semantic(root).source_refs||attrs(root).source_refs);if(sourceRefs.length){const h=document.createElement('h3');h.textContent=ui('Original sources for this topic','Originalquellen des Themas');sources.append(h);sourceRefs.slice(0,12).forEach((ref,index)=>{const id=sourceId(ref);if(!id)return;const button=document.createElement('button');button.type='button';button.textContent=sourceLabel(ref)||`${tx('Quelle')} ${index+1}`;button.addEventListener('click',()=>loadSource(id,button,{line_start:ref?.line_start,line_end:ref?.line_end}));sources.append(button);});specialist.append(sources);}
    const extras=document.createElement('div');extras.className='theme-extra';for(const [key,de,en] of [['structures','Strukturkarte · keine Beweiskette','Structure map · not a proof chain'],['goals','Zielrichtungen','Goals']]){const entries=asArray(categories[key]);if(!entries.length)continue;const section=document.createElement('section');const h=document.createElement('h3');h.textContent=ui(en,de);section.append(h);const list=document.createElement('div');list.className='reading-claim-list';entries.forEach(item=>list.append(shortClaim(item,key)));section.append(list);extras.append(section);}if(extras.childElementCount)specialist.append(extras);
    const sourceViewer=document.createElement('div');sourceViewer.id='source-viewer';
    host.append(head,story,attemptSection,specialist,sourceViewer);view.replaceChildren(host);
    relatedInput.addEventListener('change',()=>{localStorage.setItem(`tfpt-topic-related:${topicId}`,String(relatedInput.checked));loadTopic(topicId);});
  }

  async function loadTopic(topicId) {
    if(!isCurrentTopicRoute(topicId))return;
    const topic=state.topics.find(t=>t.id===topicId);if(!topic){if(!state.topicsLoaded){view.innerHTML=`<div class="loading-block"><span class="spinner"></span><span>${english()?'Loading topics before opening this page…':'Themen werden geladen; die gewünschte Seite bleibt geöffnet …'}</span></div>`;return;}view.innerHTML=`<div class="error-state">${english()?'This topic is not in the current topic list.':'Dieses Themenfeld ist in der aktuellen Themenliste nicht vorhanden.'} <button class="mini-button" id="retry-topic-list">${english()?'Retry':'Erneut laden'}</button></div>`;$('#retry-topic-list')?.addEventListener('click',()=>location.reload());return;}
    const request=++state.serial;
    view.innerHTML=`<div class="loading-block"><span class="spinner"></span><span>${english()?'Opening topic…':'Themenfeld wird geöffnet …'}</span></div>`;
    void refreshStatus();
    try {const related=localStorage.getItem(`tfpt-topic-related:${topicId}`)==='true',themePath=`/api/theme?id=${encodeURIComponent(topicId)}${related?'&include_related=true':''}`;const [data,theme]=await Promise.all([api(`/api/node?id=${encodeURIComponent(topicId)}`),api(themePath)]);if(request!==state.serial||!isCurrentTopicRoute(topicId))return;const node=data.node||data;renderTopicDashboard(theme,node);}
    catch(error) {if(request!==state.serial||!isCurrentTopicRoute(topicId))return;view.innerHTML=`<div class="error-state">${esc(error.message|| (english()?'Could not load the topic overview.':'Themenübersicht konnte nicht geladen werden.'))} <button class="mini-button" id="retry-topic">${english()?'Retry':'Erneut laden'}</button></div>`;$('#retry-topic')?.addEventListener('click',()=>loadTopic(topicId));}
  }
  async function showResearch(query='') {
    const request=++state.serial,route=location.hash,language=state.language;state.current=null;state.searchQuery='';setNotice('');
    document.body.classList.add('research-mode');document.body.classList.remove('detail-mode','map-mode','research-detail-mode');
    if(state.researchRoute!==route)window.scrollTo({top:0,behavior:'instant'});state.researchRoute=route;
    const mode=document.querySelector('[name="search-mode"][value="research"]');if(mode)mode.checked=true;
    searchInput.value=/^(problem:|topic:|goal:)/.test(query)?'':query;
    await window.TFPTResearch.open({host:view,query,api,isCurrent:()=>request===state.serial&&route===location.hash&&language===state.language});
  }
  function home() {if(location.hash)history.replaceState(null,'',location.pathname+location.search);state.searchQuery='';updateSidebarTopics();void showResearch();setNotice('');}

  function handleHash() {
    const rawHash=location.hash.slice(1),hash=decodeURIComponent(rawHash);
    document.body.classList.toggle('detail-mode',/^(node\/|topic\/|chain\/|source\/)/.test(hash)||hash==='map');
    document.body.classList.toggle('map-mode',hash==='map');
    document.body.classList.toggle('research-mode',!hash||hash.startsWith('research/'));
    if(hash.startsWith('research/')){void showResearch(hash.slice(9));return;}
    if(hash.startsWith('search/')){search(hash.slice(7));return;}
    if(hash==='discover'){state.serial++;state.current=null;state.searchQuery='';updateSidebarTopics();showHomeData();return;}
    if(hash==='map'){state.current=null;window.TFPTMap?.open();return;}
    if(hash==='orientation'){showOrientation();return;}
    if(hash.startsWith('chain/')){loadNode(hash.slice(6),true);return;}
    if(hash.startsWith('node/')){loadNode(hash.slice(5));return;}
    if(hash.startsWith('topic/')){loadTopic(hash.slice(6));return;}
    if(hash.startsWith('source/')){
      const sourceRoute=parseSourceRoute(rawHash);
      if(state.current)renderDossier(state.current);
      else view.innerHTML=`<button class="back-button" id="source-home">← ${esc(tx('Zur Themenübersicht'))}</button><div id="source-viewer"></div>`;
      $('#source-home')?.addEventListener('click',home);
      setTimeout(()=>loadSource(sourceRoute.id,undefined,sourceRoute),0);return;
    }
    home();
  }
  function applyStatus(status){
    state.lastStatus=status;
    const counts=status.counts||{},meta=status.meta||{},coverage=meta.coverage||{};
    const fingerprint=JSON.stringify([counts,status.source_count,coverage]);if(state.statusFingerprint&&state.statusFingerprint!==fingerprint){const hint=$('#refresh-hint');if(hint){hint.hidden=false;hint.textContent=tx('Neue Daten verfügbar');}}state.statusFingerprint=fingerprint;
    const nodeCount=counts.nodes??counts.node_count??meta.node_count??counts.total??'—';
    const claimCount=counts.claims??counts.claim_count??coverage.claims_total??'—';
    const sourceCount=status.source_count??counts.sources??counts.source_count??'—';
    $('#stat-nodes').textContent=nodeCount;$('#stat-claims').textContent=claimCount;$('#stat-sources').textContent=sourceCount;
    const coverageEl=$('#coverage-detail');
    if(coverage.claims_total!=null){
      const pieces=[`${coverage.claims_total} ${tx('Aussagen insgesamt')}`];
      if(coverage.claims_enriched!=null)pieces.push(`${coverage.claims_enriched} ${tx('KI-ausgewertet')}`);
      if(coverage.claims_pending!=null)pieces.push(`${coverage.claims_pending} ${tx('ausstehend')}`);
      if(coverage.claims_error!=null)pieces.push(`${coverage.claims_error} ${tx('mit Extraktionsfehler')}`);
      if(coverage.scripts_total!=null)pieces.push(`${coverage.scripts_total} ${tx('Skripte')}`);
      if(coverage.lean_theorems_total!=null)pieces.push(`${coverage.lean_theorems_total} ${tx('statisch erfasste Lean-Sätze')}`);
      coverageEl.textContent=pieces.join(' · ');coverageEl.hidden=false;
    }
  }
  async function refreshStatus(){
    const request=++state.statusSerial;
    try{
      const status=await api('/api/status');if(request!==state.statusSerial)return;
      applyStatus(status);
      $('.live-dot')?.classList.remove('error');$('.live-dot')?.classList.add('ready');
      const connection=$('#connection-label');if(connection)connection.textContent=tx('Wissensgraph verbunden');
    }catch{
      // Keep the last successful coverage numbers visible while the graph is temporarily unavailable.
    }
  }
  async function refreshCurrent(){
    const button=$('#refresh-view'),hint=$('#refresh-hint'),scroll=window.scrollY,hash=location.hash,graphActive=!!view.querySelector('[data-tab="graph"].active');
    if(button){button.disabled=true;button.textContent=tx('Wird aktualisiert …');}
    try{
      const [status,topics]=await Promise.all([api('/api/status'),api('/api/topics')]);state.topics=normalizedTopics(topics);state.topicsLoaded=true;applyStatus(status);updateSidebarTopics();
      if(hash.startsWith('#topic/'))await loadTopic(decodeURIComponent(hash.slice(7)));
      else if(hash.startsWith('#chain/'))await loadNode(decodeURIComponent(hash.slice(7)),true);
      else if(hash.startsWith('#node/'))await loadNode(decodeURIComponent(hash.slice(6)));
      else if(hash.startsWith('#source/')){if(state.current?.node?.id)await loadNode(state.current.node.id);else if(state.current)renderDossier(state.current);const sourceRoute=parseSourceRoute(hash.slice(1));await loadSource(sourceRoute.id,undefined,sourceRoute);}
      else if(hash.startsWith('#research/'))await showResearch(decodeURIComponent(hash.slice(10)));
      else if(hash.startsWith('#search/'))await search(decodeURIComponent(hash.slice(8)));
      else if(hash==='#discover')showHomeData();
      else if(hash==='#orientation')await showOrientation();
      else if(hash==='#map')await window.TFPTMap?.open();
      else if(state.searchQuery)await search(state.searchQuery);
      else {await showResearch();}
      if(graphActive&&hash.startsWith('#node/'))view.querySelector('[data-tab="graph"]')?.click();
      if(hint){hint.hidden=false;hint.textContent=tx('Aktualisiert');}setTimeout(()=>{if(hint)hint.hidden=true;},2200);
    }catch{if(hint){hint.hidden=false;hint.textContent=tx('Aktualisierung fehlgeschlagen');}}
    finally{if(button){button.disabled=false;button.textContent=tx('Aktualisieren');}setTimeout(()=>window.scrollTo({top:scroll,behavior:'instant'}),60);}
  }
  async function rerenderForLanguage(){
    const scroll=window.scrollY,hash=location.hash,graphActive=!!view.querySelector('[data-tab="graph"].active');
    if(state.lastStatus)applyStatus(state.lastStatus);
    if(!state.searchQuery)updateSidebarTopics();
    if(hash.startsWith('#topic/'))await loadTopic(decodeURIComponent(hash.slice(7)));
    else if(hash.startsWith('#chain/'))await loadNode(decodeURIComponent(hash.slice(7)),true);
    else if(hash.startsWith('#node/'))await loadNode(decodeURIComponent(hash.slice(6)));
    else if(hash.startsWith('#research/'))await showResearch(decodeURIComponent(hash.slice(10)));
    else if(hash.startsWith('#search/'))await search(decodeURIComponent(hash.slice(8)));
    else if(hash==='#discover')showHomeData();
    else if(hash==='#orientation')await showOrientation();
    else if(hash==='#map')await window.TFPTMap?.open();
    else if(hash.startsWith('#source/')){if(state.current)renderDossier(state.current);const sourceRoute=parseSourceRoute(hash.slice(1));await loadSource(sourceRoute.id,undefined,sourceRoute);}
    else if(state.searchQuery)await search(state.searchQuery);
    else await showResearch();
    if(graphActive&&hash.startsWith('#node/'))view.querySelector('[data-tab="graph"]')?.click();
    setTimeout(()=>window.scrollTo({top:scroll,behavior:'instant'}),100);
  }
  async function boot() {
    $('#search-form').addEventListener('submit',e=>{e.preventDefault();const query=searchInput.value.trim();if(!query)return;const mode=document.querySelector('[name="search-mode"]:checked')?.value;const hash=`#${mode==='entries'?'search':'research'}/${encodeURIComponent(query)}`;if(location.hash===hash)handleHash();else location.hash=hash;});
    searchInput.addEventListener('keydown',e=>{if(e.key==='Escape'){searchInput.value='';home();}});
    $('#refresh-view')?.addEventListener('click',refreshCurrent);
    document.addEventListener('tfpt-language-change',event=>{state.language=event.detail.lang;rerenderForLanguage();});
    const dot=$('.live-dot'), connection=$('#connection-label');
    try {
      const [status,topics]=await Promise.all([api('/api/status'),api('/api/topics')]);
      state.topics=normalizedTopics(topics);state.topicsLoaded=true;updateSidebarTopics();
      applyStatus(status);
      dot.classList.add('ready');connection.textContent=tx('Wissensgraph verbunden');
      state.statusTimer=setInterval(()=>{void refreshStatus();},30000);
      handleHash();
    } catch(error) {
      dot.classList.add('error');connection.textContent=tx('Verbindung zum Wissensgraphen fehlgeschlagen');
      side.innerHTML=`<div class="empty-state">${esc(tx('Themen konnten nicht geladen werden.'))}</div>`;
      view.innerHTML=`<div class="error-state">${esc(tx('Die Atlas-Oberfläche ist bereit, aber der lokale Wissensdienst antwortet nicht.'))} ${esc(error.message)}<button class="mini-button" id="retry-boot">${esc(tx('Erneut verbinden'))}</button></div>`;
      $('#retry-boot')?.addEventListener('click',()=>location.reload());
    }
  }
  window.addEventListener('hashchange',handleHash);
  boot();
})();
