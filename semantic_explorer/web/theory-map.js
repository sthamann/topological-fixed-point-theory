(() => {
  'use strict';
  const esc = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const en = () => document.documentElement.lang === 'en';
  const text = (de, english) => en() ? english : de;
  const field = (obj, key) => obj?.[`${key}_${en() ? 'en' : 'de'}`] || obj?.[key] || '';
  const state = {data:null, stage:'closure', selected:'map:e8-glue', serial:0};
  const stages = [
    {id:'inputs', title:['Ausgangspunkte','Starting points'], subtitle:['Randkern · Normierung · Träger','Kernel · normalization · carrier'], description:['P1 ist der strukturierte Randkern. Die c₃ Normierung und die P2 Reduktionen besitzen eigene Herleitungen unter benannten Voraussetzungen.','P1 is the structured boundary kernel. The c₃ normalization and P2 reductions have their own derivations under stated premises.'], status:['Vorausgesetzt','Assumed'], nodes:['map:p1','map:c3-normalization','map:p2','map:p2-reductions'], focus:'map:p1'},
    {id:'construction',title:['Strukturen aufbauen','Build the structures'],subtitle:['Naht · Double-Cover · D₅ / A₃','Seam · double cover · D₅ / A₃'],description:['Dasselbe Viermarkenobjekt liefert Familienrang drei und Funktionsraumdimension fünf. Die konkreten physikalischen Zuordnungen behalten ihre Voraussetzungen.','The same four-mark object gives family rank three and function-space dimension five. The physical identifications retain their stated premises.'],status:['Mit festgelegten Zuordnungen','Under specified dictionaries'],nodes:['map:seam','map:shared-fourmark-origin','map:double-cover','map:carrier-d5','map:family-a3'],focus:'map:double-cover'},
    {id:'closure',title:['Gemeinsame Algebra','Shared algebra'],subtitle:['D₅ ⊕ A₃ → E₈','D₅ ⊕ A₃ → E₈'],description:['Die festgelegte Verklebung schließt zum E₈-Gitter. Der Anker bündelt interne Zahlbeziehungen.','The specified glue closes to the E₈ lattice. The anchor organizes internal numerical relations.'],status:['Im gewählten Rahmen gezeigt','Shown within the chosen framework'],nodes:['map:e8-glue','map:primitive-integral-completion','map:anchor','map:backbootstrap'],focus:'map:e8-glue'},
    {id:'shared',title:['Gemeinsame Ordnung','Shared organization'],subtitle:['Clocks · Anker · φ₀','Clocks · anchor · φ₀'],description:['Clockwirkungen und gemeinsame Formeleingaben verbinden mehrere Auslesungen.','Clock actions and shared formula inputs connect several readouts.'],status:['Interne Beziehungen','Internal relations'],nodes:['map:clock-network','map:shared-seed'],focus:'map:shared-seed'},
    {id:'readouts',title:['Was daraus gelesen wird','What is read out'],subtitle:['α · Flavor · Quanten- und Gravitationsmodelle','α · flavor · quantum and gravity models'],description:['Die Zweige liefern Formeln und Modellresultate mit jeweils eigenen Voraussetzungen.','These branches supply formulas and model results, each with its own assumptions.'],status:['Bedingte Ergebnisse','Conditional results'],nodes:['map:alpha','map:flavor','map:quantum','map:gravity-readouts'],focus:'map:flavor'},
    {id:'physics',title:['Physisch vervollständigen','Complete the physical theory'],subtitle:['Gemeinsame Quelle · Raumzeit · Dynamik','Common source · spacetime · dynamics'],description:['Die operative Saat enthält Zustand und Zeit bereits. Die vollständige gemeinsame physische Quellenauswertung und ihre Feldübergänge bleiben zu leisten.','The operative seed already contains state and time. Full evaluation of the common physical source and its field transitions remains to be completed.'],status:['Gemeinsame Herleitung offen','Joint derivation remains open'],nodes:['map:state-selection','map:local-net','map:geometry-selection','map:continuum','map:local-boost','map:gravity-variation','map:area-normalization','map:process-law','map:toe-goal'],focus:'map:state-selection'}
  ];
  const status = n => ({assumed:text('Gesetzte Eingabe','Declared input'),established_scoped:text('Im angegebenen Rahmen gezeigt','Shown within stated scope'),conditional:text('Gilt unter Voraussetzungen','Conditional on assumptions'),open:text('Herleitung offen','Derivation open'),goal:text('Forschungsziel','Research goal')}[n?.state] || text('Einordnung offen','Status unresolved'));
  const topicNames={origin:['Ursprung','Origin'],e8:['E₈ und Codes','E₈ and codes'],clocks:['Verzahnte Clocks','Interlocking clocks'],flavor:['Standardmodell und Flavor','Standard Model and flavor'],alpha:['Feinstrukturkonstante','Fine-structure constant'],quantum:['Quanten- und Randstruktur','Quantum and boundary structure'],gravity:['Gravitation','Gravity'],cosmology:['Kosmologie','Cosmology'],universalraum:['Universalraum','Universal space'],process:['Prozess und Compilation','Process and compilation'],spacetime:['Raumzeit','Spacetime'],metrology:['Metrologie','Metrology'],rh:['Riemann-Hypothese','Riemann hypothesis'],predictions:['Tests und Vorhersagen','Tests and predictions']};
  const nodeById = id => state.data?.nodes?.find(n => n.id === id);
  const select = (id, stage) => { state.selected=id; state.stage=stage || stages.find(s=>s.nodes.includes(id))?.id || 'context'; render(); document.querySelector('.story-focus')?.scrollIntoView({block:'start',behavior:'instant'}); };
  const labelId = id => id.startsWith('claim:') ? `${text('Aussage','Claim')} ${id.slice(6)}` : id.startsWith('topic:') ? `${text('Thema','Topic')}: ${id.slice(6)}` : `${text('Baustein','Structure')}: ${id.split(':').slice(1).join(':')}`;
  function sourceLinks(refs) {
    return (refs || []).map(s => {
      const id=s.source_id || s.id; const path=s.path || text('Originalquelle','Original source');
      if (!id) return `<span>${esc(path)}</span>`;
      const q=new URLSearchParams(); if(s.line_start != null) q.set('line_start',s.line_start); if(s.line_end != null)q.set('line_end',s.line_end);
      const where=s.line_start != null ? ` · ${text('Zeile','line')} ${s.line_start}${s.line_end && s.line_end!==s.line_start ? '–'+s.line_end : ''}` : '';
      return `<a href="#source/${encodeURIComponent(id)}${q.size?'?'+q:''}">${esc(path+where)}</a>`;
    }).join('');
  }
  function nodeButton(n, extra='') {
    if(!n) return '';
    return `<button class="story-node-button ${extra}" data-story-node="${esc(n.id)}"><span class="story-node-state is-${esc(n.state)}">${esc(status(n))}</span><strong>${esc(field(n,'label'))}</strong></button>`;
  }
  function relationEntry(edge, inbound) {
    const other=nodeById(inbound ? edge.from : edge.to); if(!other)return '';
    const consistency=edge.from==='map:backbootstrap' || edge.to==='map:backbootstrap';
    const open=edge.kind==='open_requirement';
    const conditional=edge.kind==='conditional_bridge';
    const shared=edge.kind==='shared_seed';
    const kind=consistency?text('Konsistenzprüfung','Consistency check'):open?text('Offene Verbindung','Open connection'):conditional?text('Bedingter Anschluss','Conditional connection'):shared?text('Gemeinsame Eingabe','Shared input'):text('Konstruktion im Modell','Construction within the model');
    const relation=consistency ? text('Interne Konsistenzprüfung','Internal consistency check') : field(edge,'label');
    return `<div class="story-relation ${open?'is-open':''} ${consistency?'is-consistency':''} ${conditional?'is-conditional':''} ${shared?'is-shared':''}"><span class="story-relation-kind">${esc(kind)}</span>${inbound?nodeButton(other):''}<p class="story-connector"><span aria-hidden="true">${consistency?'↔':open?'⋯ →':'→'}</span>${esc(relation)}</p>${!inbound?nodeButton(other):''}${edge.conditions?.length?`<details class="story-edge-conditions"><summary>${text('Voraussetzungen ansehen','View conditions')}</summary><p>${esc(edge.conditions.map(c=>typeof c==='string'?c:field(c,'label')).join(' '))}</p></details>`:''}</div>`;
  }
  function focusMarkup() {
    const n=nodeById(state.selected); if(!n)return '';
    const incoming=state.data.edges.filter(e=>e.to===n.id && e.kind!=='scope_boundary');
    const outgoing=state.data.edges.filter(e=>e.from===n.id && e.kind!=='scope_boundary');
    const group=stages.find(s=>s.id===state.stage);
    const contextNodes=n.id==='map:rh-branch'?['map:rh-branch']:['map:universalraum','map:ur-to-tfpt','map:tfpt-to-ur'];
    const groupNodes=(group?.nodes || contextNodes).map(nodeById).filter(Boolean);
    const scope=field(n,'scope_note');
    return `<section class="story-focus" aria-labelledby="story-focus-title">
      <header><button class="story-back" data-story-back>${text('↑ Zur Gesamtübersicht','↑ Back to overview')}</button><p class="story-kicker">${text('Den Zusammenhang im Detail verstehen','Understand the connection')}</p><h2 id="story-focus-title">${esc(field(n,'label'))}</h2></header>
      <nav class="story-node-tabs" aria-label="${text('Bausteine dieser Stufe','Components in this stage')}">${groupNodes.map(x=>`<button data-story-node="${esc(x.id)}" aria-pressed="${x.id===n.id}">${esc(field(x,'label'))}</button>`).join('')}</nav>
      <div class="story-neighborhood">
        <section class="story-upstream"><h3>${text('Baut auf …','Builds on …')}</h3>${incoming.length ? incoming.map(e=>relationEntry(e,true)).join('') : `<p class="story-empty">${n.state==='assumed'?text('Dieser Baustein wird als Ausgangspunkt gesetzt. Seine Herkunft ist eine eigene offene Frage.','This component is declared as a starting point. Its origin is a separate open question.'):text('Hier ist keine direkte vorgelagerte Herleitung erfasst. Der Themenbezug allein liefert keinen Beweis.','No direct upstream derivation is recorded here. A topic link alone supplies no proof.')}</p>`}</section>
        <section class="story-selected is-${esc(n.state)}"><span class="story-state">${esc(status(n))}</span><h3>${esc(field(n,'label'))}</h3><p>${esc(field(n,'summary'))}</p>${(n.related_ids||[]).filter(id=>id.startsWith('topic:')).length?`<div class="story-topic-links">${(n.related_ids||[]).filter(id=>id.startsWith('topic:')).map(id=>`<a href="#topic/${encodeURIComponent(id)}">${text('Thema','Topic')}: ${esc(topicNames[id.slice(6)]?.[en()?1:0]||id.slice(6))} ↗</a>`).join('')}</div>`:''}${scope?`<details><summary>${text('Genauer Geltungsbereich','Exact scope')}</summary><p>${esc(scope)}</p></details>`:''}<details class="story-evidence"><summary>${text('Aussagen und Originalquellen öffnen','Open claims and original sources')}</summary><div class="story-proof-links">${(n.related_ids||[]).filter(x=>x.startsWith('claim:')).map(id=>`<a href="#node/${encodeURIComponent(id)}">${esc(labelId(id))}</a>`).join('')}${sourceLinks(n.source_refs)}</div></details></section>
        <section class="story-downstream"><h3>${text('Führt weiter zu …','Leads on to …')}</h3>${outgoing.length ? outgoing.map(e=>relationEntry(e,false)).join('') : `<p class="story-empty">${text('Ein Endpunkt dieses Ausschnitts. Voraussetzungen und Belege stehen im zugehörigen Dossier.','An endpoint in this view. Assumptions and evidence are available in its dossier.')}</p>`}</section>
      </div>
      <p class="story-edge-key">${text('Durchgezogene Verbindungen: benannte Konstruktion oder Modellbeziehung. Gepunktete Verbindungen: noch zu begründender Übergang. Rückbezüge: Konsistenzprüfung.','Solid links: named construction or model relation. Dotted links: a transition still requiring justification. Back-links: consistency checks.')}</p>
    </section>`;
  }
  function render() {
    const host=document.getElementById('view'); if(!host || !state.data)return;
    document.body.classList.add('story-map-mode');
    host.innerHTML=`<article class="story-map">
      <header class="story-header"><p class="story-kicker">${text('TFPT · Das Gesamtbild','TFPT · The big picture')}</p><h1>${text('Von den Eingaben zur physikalischen Theorie','From the inputs to a physical theory')}</h1><p>${text('TFPT verbindet wenige Ausgangsannahmen mit einer gemeinsamen algebraischen Struktur. Daraus werden verschiedene Größen ausgelesen. Die gemeinsame physikalische Herkunft dieser Konstruktionen bleibt der entscheidende offene Teil.','TFPT connects a small set of starting assumptions to a shared algebraic structure. Different quantities are read from that structure. A common physical origin for these constructions remains the central open part.')}</p></header>
      <nav class="story-spine" aria-label="${text('Hauptlinie der Theorie in sechs Schritten','Main theory line in six steps')}">${stages.map((s,i)=>`<button data-story-stage="${s.id}" aria-pressed="${state.stage===s.id}" class="story-stage ${s.id==='physics'?'story-stage-open':''}"><span class="story-step">${i+1}</span><h2>${esc(s.title[en()?1:0])}</h2><strong>${esc(s.subtitle[en()?1:0])}</strong><span class="story-stage-status">${esc(s.status[en()?1:0])}</span></button>`).join('')}</nav>
      <p class="story-reading-order"><strong>${esc(stages.find(s=>s.id===state.stage)?.description[en()?1:0] || '')}</strong> ${text('Die sechs Schritte ordnen den Aufbau. Wähle einen Schritt: Darunter siehst du die tatsächlich erfassten Verbindungen und ihre Voraussetzungen.','The six steps organize the architecture. Choose a step to see the recorded connections and their assumptions below.')}</p>
      <div class="story-crosslinks"><button data-story-node="map:universalraum"><strong>${text('Wo steht der Universalraum?','Where does Universalraum fit?')}</strong><span>${text('Er fragt nach der gemeinsamen Quelle und dem vollständigen Prozess — über den gesamten Aufbau hinweg.','It asks about the common source and complete process, across the entire architecture.')}</span></button><button data-story-node="map:local-boost"><strong>${text('Wie schließt Gravitation an?','How does gravity connect?')}</strong><span>${text('Lokale Randstruktur → geometrischer Fluss → Flächenantwort. Die Übergänge sind eigene offene Aufgaben.','Local boundary structure → geometric flow → area response. The transitions are separate open tasks.')}</span></button></div>
      ${focusMarkup()}
      <details class="story-other"><summary>${text('Weitere Verbindungen und Forschungszweige','Further connections and research branches')}</summary><div class="story-secondary">${nodeButton(nodeById('map:rh-branch'))}${nodeButton(nodeById('map:backbootstrap'))}${nodeButton(nodeById('map:toe-goal'))}</div><p>${text('Rückbootstrap prüft die Konsistenz innerhalb der gewählten Architektur. Der RH-Zweig verfolgt ein eigenes mathematisches Ziel. Beide werden getrennt von der physikalischen Gesamtherleitung eingeordnet.','Back-bootstrap checks consistency within the chosen architecture. The RH branch pursues a separate mathematical goal. Both are distinguished from the complete physical derivation.')}</p></details>
    </article>`;
    host.querySelector('[data-story-back]')?.addEventListener('click',()=>host.querySelector('.story-header').scrollIntoView({block:'start',behavior:'instant'}));
    host.querySelectorAll('[data-story-stage]').forEach(b=>b.addEventListener('click',()=>{const s=stages.find(x=>x.id===b.dataset.storyStage);select(s.focus,s.id);}));
    host.querySelectorAll('[data-story-node]').forEach(b=>b.addEventListener('click',()=>select(b.dataset.storyNode)));
  }
  async function open() {
    document.body.classList.add('detail-mode','story-map-mode');
    const host=document.getElementById('view');const serial=++state.serial;
    host.innerHTML=`<div class="loading-block">${text('Der Aufbau wird geladen …','Loading the architecture …')}</div>`;
    try {
      const r=await fetch('/api/theory-map');if(!r.ok)throw new Error(`HTTP ${r.status}`);const payload=await r.json();
      if(serial!==state.serial || location.hash!=='#map')return;
      state.data=payload.theory_map || payload.map || payload;render();
    } catch(e) {if(location.hash==='#map')host.innerHTML=`<p>${text('Die Karte konnte nicht geladen werden. Bitte erneut aktualisieren.','The map could not be loaded. Please refresh.')}</p>`;}
  }
  window.addEventListener('hashchange',()=>{if(location.hash!=='#map')document.body.classList.remove('story-map-mode');});
  window.TFPTMap={open,focus:id=>{state.selected=id;state.stage=stages.find(s=>s.nodes.includes(id))?.id || 'context';if(location.hash==='#map')render();else location.hash='#map';}};
})();
