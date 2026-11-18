(() => {
  'use strict';
  const esc = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const arr = value => Array.isArray(value) ? value : [];
  const en = () => document.documentElement.lang === 'en';
  const t = (de, english) => en() ? english : de;
  const href = (kind, id) => `#${kind}/${encodeURIComponent(id)}`;
  const researchHref = id => href('research', id);
  const topicTitle = id => {
    const names = {'topic:origin':['Ursprung','Origin'],'topic:e8':['E₈ & Codes','E₈ & codes'],'topic:clocks':['Verzahnte Clocks','Interlocking clocks'],'topic:flavor':['Standardmodell & Flavor','Standard Model & flavor'],'topic:alpha':['Feinstrukturkonstante','Fine-structure constant'],'topic:quantum':['Quanten- & Randstruktur','Quantum & boundary structure'],'topic:gravity':['Gravitation','Gravity'],'topic:cosmology':['Kosmologie','Cosmology'],'topic:rh':['Riemann-Hypothese','Riemann hypothesis'],'topic:validation':['Prüfungen & Vorhersagen','Tests & predictions'],'topic:universalraum':['Universalraum','Universal space'],'topic:metrology':['Metrologie','Metrology'],'topic:spacetime':['Raumzeit','Spacetime'],'topic:process':['Prozess & Kompilation','Process & compilation']};
    return names[id]?.[en()?1:0] || id;
  };
  const title = item => item?.title || item?.label || item?.id || t('Ohne Titel', 'Untitled');
  const kind = finding => {
    const status = String(finding.status || '').toLowerCase();
    if (['excluded_scoped','excluded','rejected','negative','refuted','counterexample'].includes(status)) return 'excluded';
    if (['resolved','established','proved','established_scoped','recorded_established','solved','exact'].includes(status)) return 'established';
    if (['conditional','conditional_result'].includes(status)) return 'conditional';
    return 'recorded';
  };
  const findingLabel = finding => finding.source_state === 'stale' ? t('Quellenprüfung nötig', 'Source review needed') : {
    established: t('Belegt im angegebenen Rahmen', 'Established within scope'),
    excluded: t('Im angegebenen Rahmen ausgeschlossen', 'Excluded within scope'),
    conditional: t('Bedingter Befund', 'Conditional finding'),
    recorded: t('Erfasster Befund', 'Recorded finding')
  }[kind(finding)];
  const strategyLabel = outcome => ({
    reuse_result: t('Vorhandenes Ergebnis weiterverwenden', 'Reuse the existing result'),
    excluded_scoped: t('Unter diesen Voraussetzungen ausgeschlossen', 'Excluded under these premises'),
    insufficient_for_goal: t('Für das Gesamtziel nicht hinreichend', 'Insufficient for the overall goal')
  }[outcome] || t('Einordnung im Datensatz prüfen', 'Check the recorded assessment'));
  const connectionLabel = kind => ({
    requires_result: t('Benötigt ein Ergebnis', 'Requires a result'),
    shares_unresolved_input: t('Teilt eine noch ungeklärte Eingabe', 'Shares an unresolved input'),
    constrains_under: t('Schränkt unter Voraussetzungen ein', 'Constrains under premises')
  }[kind] || t('Gespeicherte Beziehung', 'Recorded relationship'));
  const taskStatus = task => task.effective_status || task.status || 'needs_review';
  const taskLabel = task => ({
    ready: t('Jetzt bearbeitbar', 'Ready to investigate'), blocked: t('Voraussetzung fehlt', 'Prerequisite missing'),
    resolved: t('Teilaufgabe gelöst', 'Subtask resolved'), excluded: t('Weg ausgeschlossen', 'Approach excluded'),
    needs_review: t('Prüfung nötig', 'Review needed'), open: t('Offen', 'Open')
  }[taskStatus(task)] || t('Status prüfen', 'Check status'));
  const nextTask = problem => arr(problem.tasks).find(task => taskStatus(task) === 'ready') || arr(problem.tasks).find(task => !['resolved','excluded'].includes(taskStatus(task)));
  function refs(items, index, css = '') {
    return arr(items).map(id => `<a class="research-ref ${css}" href="${esc(href('node', id))}">${esc(index.get(id) ? title(index.get(id)) : id)}</a>`).join('');
  }
  function sources(items) {
    return arr(items).map((source, i) => {
      const id = typeof source === 'string' ? source : source.source_id || source.id;
      if (!id) return '';
      const range = new URLSearchParams();
      if (source.line_start != null) range.set('line_start', source.line_start);
      if (source.line_end != null) range.set('line_end', source.line_end);
      const path = source.path || id;
      const rangeLabel = source.line_start != null ? ` · ${source.line_start}${source.line_end != null ? `–${source.line_end}` : ''}` : '';
      return `<a class="research-source" href="${esc(href('source', id))}${range.size ? `?${range}` : ''}" title="${esc(path + rangeLabel)}"><span>§ ${i + 1}</span> ${esc(path.split('/').slice(-2).join('/') + rangeLabel)}</a>`;
    }).join('');
  }
  function findingCard(finding, index) {
    return `<article class="research-finding is-${kind(finding)} ${finding.source_state === 'stale' ? 'is-stale' : ''}">
      <span class="research-pill">${esc(findingLabel(finding))}</span><h4>${esc(title(finding))}</h4>
      <p>${esc(finding.statement || '')}</p>
      <p class="research-scope"><strong>${t('Gilt für', 'Scope')}:</strong> ${esc(finding.scope || t('Geltungsbereich nicht erfasst; Originalquelle prüfen.', 'Scope not recorded; inspect the original source.'))}</p>
      ${arr(finding.premise_ids).length ? `<div class="research-requires"><span>${t('Unter diesen Voraussetzungen gemeinsam', 'With all these premises together')}</span>${refs(finding.premise_ids, index)}</div>` : ''}
      <div class="research-sources">${sources(finding.sources)}</div>
      <details class="research-record-links"><summary>${t('Befund und ursprüngliche Versuche', 'Finding and original attempts')}</summary><div>${refs([finding.id, ...arr(finding.attempt_ids)], index)}</div></details>
    </article>`;
  }
  function taskCard(task, index, number) {
    const status = taskStatus(task);
    return `<li class="research-task is-${esc(status)}" id="research-task-${esc(task.id)}">
      <span class="research-step" aria-hidden="true">${status === 'resolved' ? '✓' : status === 'excluded' ? '×' : number + 1}</span>
      <div><span class="research-pill">${esc(taskLabel(task))}</span><h4><a href="${esc(href('node', task.id))}">${esc(title(task))}</a></h4>
      ${arr(task.depends_on).length ? `<div class="research-requires"><span>${t('Benötigt gemeinsam', 'Requires together')}</span>${refs(task.depends_on, index)}</div>` : `<p class="research-task-note">${t('Keine vorgelagerte Teilaufgabe eingetragen.', 'No prerequisite subtask recorded.')}</p>`}
      ${task.success_criterion ? `<p><strong>${t('Erfolgskriterium', 'Success criterion')}:</strong> ${esc(task.success_criterion)}</p>` : ''}
      ${!['resolved','excluded'].includes(status) && task.next_action ? `<p><strong>${t('Nächste Handlung', 'Next action')}:</strong> ${esc(task.next_action)}</p>` : ''}
      ${arr(task.required_premise_ids).length ? `<div class="research-requires"><span>${t('Zusätzliche Voraussetzungen', 'Additional premises')}</span>${refs(task.required_premise_ids, index)}</div>` : ''}
      ${arr(task.resolution_ids).length ? `<div class="research-requires"><span>${t('Abschluss durch', 'Resolved by')}</span>${refs(task.resolution_ids, index)}</div>` : ''}</div>
    </li>`;
  }
  function problemCard(problem) {
    const findings = arr(problem.findings), next = nextTask(problem);
    const established = findings.filter(item => kind(item) === 'established' && item.source_state !== 'stale');
    const conditional = findings.filter(item => kind(item) === 'conditional' && item.source_state !== 'stale');
    const excluded = findings.filter(item => kind(item) === 'excluded' && item.source_state !== 'stale');
    const stale = findings.filter(item => item.source_state === 'stale');
    return `<article class="research-problem-card">
      <div class="research-card-meta"><span>${established.length} ${t('gesicherte Teilresultate', 'established subresults')}</span>${conditional.length ? `<span>${conditional.length} ${t('bedingte Teilresultate', 'conditional subresults')}</span>` : ''}<span>${excluded.length} ${t('begrenzte Ausschlüsse', 'scoped exclusions')}</span>${stale.length ? `<span class="research-attention">${stale.length} ${t('zu prüfen', 'need review')}</span>` : ''}</div>
      <h3><a href="${esc(researchHref(problem.id))}">${esc(title(problem))}</a></h3><p class="research-question">${esc(problem.question || '')}</p>
      <div class="research-card-result"><span>${conditional.length && !established.length ? t('Bedingtes Teilresultat', 'Conditional subresult') : t('Darauf können wir aufbauen', 'What we can build on')}</span><p>${esc(established[0]?.statement || established[0]?.title || conditional[0]?.statement || conditional[0]?.title || t('Noch kein aktuelles gesichertes Teilresultat in dieser Auswahl.', 'No current established subresult in this selection yet.'))}</p>${established.length > 1 ? `<small>+ ${established.length - 1} ${t('weitere Teilresultate', 'more subresults')}</small>` : ''}</div>
      <div class="research-card-next"><span>${next ? esc(taskLabel(next)) : t('Kein nächster Schritt eingetragen', 'No next step recorded')}</span><p>${esc(next ? title(next) : t('Aus abgeschlossenen Teilaufgaben folgt kein automatischer Gesamtabschluss.', 'Completed subtasks do not automatically establish the overall goal.'))}</p></div>
      <a class="research-open" href="${esc(researchHref(problem.id))}">${t('Ergebnisse und nächsten Schritt öffnen', 'Open findings and next step')} <span aria-hidden="true">→</span></a>
    </article>`;
  }

  function indexCard(problem) {
    const next = arr(problem.ready_tasks)[0];
    const count = typeof problem.findings === 'number' ? problem.findings : arr(problem.findings).length;
    return `<article class="research-problem-card is-index"><div class="research-card-meta"><span>${count} ${t('eingeordnete Befunde', 'reviewed findings')}</span><span>${problem.resolved_tasks ?? '—'} ${t('gelöste Teilaufgaben', 'resolved subtasks')}</span></div><h3><a href="${esc(researchHref(problem.id))}">${esc(title(problem))}</a></h3><p class="research-question">${esc(problem.question || t('Ergebnisse, Voraussetzungen und Ausschlüsse dieser Frage gemeinsam öffnen.', 'Open the findings, premises and exclusions for this question together.'))}</p>${next ? `<div class="research-card-next"><span>${t('Jetzt bearbeitbar', 'Ready to investigate')}</span><p>${esc(title(next))}</p></div>` : ''}<a class="research-open" href="${esc(researchHref(problem.id))}">${t('Ergebnisse und nächsten Schritt öffnen', 'Open findings and next step')} <span aria-hidden="true">→</span></a></article>`;
  }
  function problemDetail(problem, packet) {
    const findings = arr(problem.findings), tasks = arr(problem.tasks), premises = arr(problem.premises || problem.accepted_premises), strategies = arr(problem.strategies);
    const activeTaskHandles = arr(problem.active_task_handles);
    const accepted = new Set([...arr(problem.accepted_premise_ids), ...arr(problem.accepted_premises).map(item => typeof item === 'string' ? item : item.id)]);
    const index = new Map([problem, ...findings, ...tasks, ...activeTaskHandles, ...premises, ...strategies].map(item => [item.id, item]));
    const next = nextTask(problem) || activeTaskHandles[0], pending = problem.pending_attempts || {};
    const partial = problem.packet_state === 'partial', detail = problem.partial_coverage || {}, continuation = packet?.budget?.continuation || {};
    const established = findings.filter(item => kind(item) === 'established' || kind(item) === 'conditional');
    const exclusions = findings.filter(item => kind(item) === 'excluded');
    const other = findings.filter(item => kind(item) === 'recorded');
    const group = (heading, items, empty) => `<section class="research-evidence-section"><h3>${heading} <span>${items.length}</span></h3>${items.length ? `<div class="research-findings">${items.map(item => findingCard(item, index)).join('')}</div>` : `<p class="research-empty-inline">${empty}</p>`}</section>`;
    return `<article class="research-problem-detail">
      <header class="research-detail-head"><a class="research-back" href="#">← ${t('Alle Forschungsfragen', 'All research questions')}</a><p class="research-eyebrow">${t('Forschungsfrage', 'Research question')}</p><h2>${esc(title(problem))}</h2><p>${esc(problem.question || '')}</p><div class="research-detail-links">${refs([problem.id], index)}${arr(problem.topics).map(id => `<a href="${esc(href('topic', typeof id === 'string' ? id : id.id))}">${esc(typeof id === 'string' ? topicTitle(id) : title(id))} ↗</a>`).join('')}</div></header>
      ${partial ? `<section class="research-notice"><strong>${t('Teilansicht eines großen Forschungsstands', 'Partial view of a large research record')}</strong><p>${esc(problem.notice || '')}</p><p>${t('Datensätze', 'Records')} ${esc(Number(detail.detail_offset || 0) + 1)}–${esc(Number(detail.detail_offset || 0) + Number(detail.detail_returned || 0))} ${t('von', 'of')} ${esc(detail.detail_total ?? '—')}.</p>${continuation.detail_offset != null ? `<button type="button" class="mini-button" data-research-detail-next="${esc(continuation.detail_offset)}">${t('Nächsten vollständigen Abschnitt laden', 'Load the next complete section')}</button>` : ''}</section>` : ''}
      ${next ? `<section class="research-next-step"><div class="research-next-label"><span class="research-step-mark">→</span><span>${t('Der nächste entscheidende Schritt', 'The next decisive step')}<small>${esc(taskLabel(next))}</small></span></div><div><h3>${esc(title(next))}</h3><p>${esc(next.next_action || next.success_criterion || t('Die konkrete Handlung ist noch zu präzisieren.', 'The concrete action still needs to be specified.'))}</p>${next.success_criterion ? `<p class="research-next-criterion"><strong>${t('Entschieden, wenn', 'Decided when')}:</strong> ${esc(next.success_criterion)}</p>` : ''}${arr(next.depends_on).length ? `<div class="research-requires"><span>${t('Voraussetzungen', 'Prerequisites')}</span>${refs(next.depends_on, index)}</div>` : ''}</div></section>` : `<p class="research-notice">${t('Für diese Frage ist keine offene Folgeaufgabe eingetragen. Das allein belegt keinen vollständigen physikalischen Abschluss.', 'No open follow-up task is recorded for this question. That alone does not establish complete physical closure.')}</p>`}
      <details class="research-assumptions"><summary>${t('Gesetzte Voraussetzungen', 'Accepted premises')} <span>${premises.length}</span><small>${t('Der Geltungsbereich aller folgenden Befunde', 'The scope of the findings below')}</small></summary><div class="research-premise-list">${premises.length ? premises.map(item => `<article><h4><a href="${esc(href('node', item.id))}">${esc(title(item))}</a></h4><p>${esc(item.statement || '')}</p><small>${accepted.has(item.id) ? t('Für diese Frage angenommen', 'Accepted for this question') : t('Zusätzliche Voraussetzung einzelner Ansätze', 'Additional premise of specific approaches')}</small></article>`).join('') : `<p>${t('Keine Voraussetzungen in dieser Antwort erfasst. Daraus folgt keine voraussetzungslose Herleitung.', 'No premises recorded in this response. This does not imply an assumption-free derivation.')}</p>`}</div></details>
      <div class="research-detail-grid"><div>${group(t('Was bereits geklärt ist', 'What is already established'), established, t('Noch kein bestätigtes Teilresultat in dieser Auswahl.', 'No established subresult in this selection yet.'))}${group(t('Welche Wege begrenzt ausgeschlossen sind', 'Which approaches are excluded within scope'), exclusions, t('Kein geprüfter Ausschluss in dieser Auswahl.', 'No reviewed exclusion in this selection.'))}${other.length ? group(t('Weitere erfasste Befunde', 'Other recorded findings'), other, '') : ''}</div>
      <aside class="research-chain"><h3>${t('Von den Voraussetzungen zum Ziel', 'From prerequisites to the goal')}</h3><p>${t('Jede Aufgabe nennt, was vorher gemeinsam vorliegen muss. Gelöste Teilaufgaben bleiben gelöst.', 'Each task names what must hold together first. Resolved subtasks remain resolved.')}</p><ol>${tasks.map((item, i) => taskCard(item, index, i)).join('')}</ol>${!tasks.length ? `<p>${t('Noch keine Aufgabenfolge erfasst.', 'No task sequence recorded yet.')}</p>` : ''}</aside></div>
      ${arr(problem.connections).length ? `<section class="research-connections"><h3>${t('Verbindung zu anderen Fragen', 'Connections to other questions')}</h3><div>${arr(problem.connections).map(connection => `<article><span class="research-pill">${esc(connectionLabel(connection.kind))}</span><h4><a href="${esc(researchHref(connection.other_problem_id || connection.target_problem_id))}">${esc(connection.other_title || connection.target_title || connection.other_problem_id || connection.target_problem_id)} →</a></h4><p>${esc(connection.explanation || '')}</p>${connection.source_state === 'stale' ? `<p class="research-connection-warning">${t('Eine zugrunde liegende Quelle wurde verändert; die Verbindung muss erneut geprüft werden.', 'An underlying source has changed; this connection needs review.')}</p>` : ''}<details><summary>${t('Tragende Voraussetzungen und Befunde', 'Supporting premises and findings')}</summary><div class="research-requires">${refs(connection.premise_ids, index)}${refs(connection.finding_ids, index)}</div><p>${t('Diese Beziehung ordnet die Forschung; sie liefert keine zusätzliche Herleitung.', 'This relationship organizes the research; it supplies no additional derivation.')}</p></details></article>`).join('')}</div></section>` : ''}
      ${strategies.length ? `<section class="research-strategies"><h3>${t('Vorhandene Ansätze vor einem neuen Versuch prüfen', 'Check existing approaches before a new attempt')}</h3><div>${strategies.map(strategy => `<details><summary>${esc(title(strategy))}<span>${esc(strategyLabel(strategy.outcome))}</span></summary><p>${esc(strategy.explanation || '')}</p>${arr(strategy.premise_ids).length ? `<div class="research-requires"><span>${t('Voraussetzungen', 'Premises')}</span>${refs(strategy.premise_ids, index)}</div>` : ''}<div class="research-requires">${refs(strategy.finding_ids, index)}</div>${strategy.changes_required ? `<p><strong>${t('Für einen neuen Versuch muss sich ändern', 'A new attempt requires a change to')}:</strong> ${esc(Array.isArray(strategy.changes_required) ? strategy.changes_required.join(' · ') : strategy.changes_required)}</p>` : ''}</details>`).join('')}</div></section>` : ''}
      <details class="research-pending"><summary>${t('Historische Versuche noch nicht in den Forschungsstand eingeordnet', 'Historical attempts not yet reconciled into the research state')} <span>${pending.total ?? '—'}</span></summary><p>${t('Diese Einträge sind noch keine bestätigten Ergebnisse oder neuen offenen Aufgaben. Ihre Quellen und Voraussetzungen müssen zuerst eingeordnet werden.', 'These records are not yet established findings or new open tasks. Their sources and assumptions need to be reconciled first.')}</p><div>${refs(arr(pending.items).map(item => item.id), new Map(arr(pending.items).map(item => [item.id,item])))}</div>${Number(pending.total) > arr(pending.items).length ? `<p>${pending.total - arr(pending.items).length} ${t('weitere Versuche hier nicht angezeigt.', 'more attempts not displayed here.')}</p>` : ''}</details>
    </article>`;
  }
  function coverageMarkup(packet) {
    const coverage = packet.coverage || {}, budget = packet.budget || {};
    const rows = [
      [coverage.reviewed_problems, t('eingeordnete Fragen', 'reviewed questions')],
      [coverage.reviewed_findings, t('geprüfte Befunde', 'reviewed findings')],
      [coverage.pending_attempts, t('historische Versuche noch einzuordnen', 'historical attempts awaiting review')],
      [coverage.stale_findings, t('Befunde mit veränderten Quellen', 'findings with changed sources')]
    ];
    const parts = rows.filter(([value]) => value != null).map(([value,label]) => `<span><strong>${esc(value)}</strong> ${label}</span>`).join('');
    const missing = rows.slice(0,3).some(([value]) => value == null);
    const omitted = budget.omitted_problems ?? budget.problems_omitted ?? coverage.omitted_problems ?? budget.index_omitted;
    const truncated = budget.truncated || packet.truncated || Number(omitted) > 0;
    return `<section class="research-coverage" aria-label="${t('Abdeckung des Forschungsstands', 'Research state coverage')}"><div>${parts || `<span>${t('Abdeckung nicht gemeldet', 'Coverage not reported')}</span>`}</div><p>${t('Gezielt eingeordneter Forschungsstand; der gesamte Atlas bleibt unter Themen und Quellen zugänglich.', 'A curated research state; the full atlas remains accessible through topics and sources.')}${coverage.corpus_complete === false ? ` ${t('Der gesamte Korpus ist noch nicht in diese Entscheidungsebene eingeordnet.', 'The entire corpus has not yet been reconciled into this decision layer.')}` : ''}${missing ? ` ${t('Einzelne Abdeckungszahlen fehlen.', 'Some coverage counts are missing.')}` : ''}${truncated ? ` <strong>${t('Begrenzter Ausschnitt', 'Limited selection')}${Number(omitted) > 0 ? ` · ${esc(omitted)} ${t('weitere Fragen', 'more questions')}` : t(' der Detaildaten; alle Fragen bleiben über den Index erreichbar', ' of detail records; all questions remain accessible through the index')}.</strong>` : ''}</p></section>`;
  }
  async function open({host, query = '', detailOffset = 0, expectedRevision = '', api, isCurrent = () => true}) {
    host.innerHTML = `<div class="loading-block"><span class="spinner"></span><span>${t('Forschungsstand wird geladen …', 'Loading research state …')}</span></div>`;
    try {
      const params = new URLSearchParams({query, lang: en() ? 'en' : 'de', max_chars: '24000'});
      if (detailOffset) params.set('detail_offset', detailOffset);
      if (expectedRevision) params.set('expected_revision', expectedRevision);
      let packet = await api(`/api/research?${params}`);
      if (query && arr(packet.problem_index).some(item => item.id === query) && !arr(packet.problems).some(item => item.id === query) && Number(packet.budget?.required_chars_for_next_problem) <= 60000) {
        params.set('max_chars', '60000');packet = await api(`/api/research?${params}`);
      }
      if (!isCurrent()) return;
      const problems = arr(packet.problems);
      const exact = problems.find(problem => problem.id === query) || problems.find(problem => problem.packet_state === 'partial');
      document.body?.classList.toggle('research-detail-mode', Boolean(exact));
      const fullById = new Map(problems.map(problem => [problem.id,problem]));
      const problemIndex = arr(packet.problem_index).length ? packet.problem_index : problems;
      const topicCoverage = arr(packet.topic_coverage || packet.coverage?.topic_coverage).filter(topic => topic.reviewed === false || topic.reviewed_problems === 0);
      host.innerHTML = `<div class="research-navigator" data-original-content>
        ${!exact ? `<header class="research-heading"><div><p class="research-eyebrow">${t('Forschungsnavigator', 'Research navigator')}</p><h2>${query ? `${t('Forschungsstand zu', 'Research state for')} „${esc(query)}“` : t('Was wissen wir — und was folgt daraus?', 'What do we know — and what follows?')}</h2><p>${t('Von gesicherten Teilresultaten über ihre Voraussetzungen zur nächsten entscheidenden Frage.', 'From established subresults and their premises to the next decisive question.')}</p></div><a class="research-discover" href="#discover">${t('Alle Themen entdecken', 'Explore all topics')} ↗</a></header>` : ''}
        ${!exact && !query ? `<p class="research-full-search"><a href="/website-audit.html">${t('Website, Atlas und Kernessenz: geprüfte Verbindungen und konkrete Korrekturen', 'Website, Atlas and Kernessenz: reviewed connections and concrete corrections')} ↗</a></p>` : ''}
        ${coverageMarkup(packet)}
        ${exact ? problemDetail(exact, packet) : problemIndex.length ? `<div class="research-problem-grid">${problemIndex.map(item => fullById.has(item.id) ? problemCard(fullById.get(item.id)) : indexCard(item)).join('')}</div>` : `<div class="research-empty"><h3>${t('Für diese Suche liegt noch kein eingeordneter Forschungsstand vor.', 'No curated research state is available for this search yet.')}</h3><p>${t('Das bedeutet nicht, dass es keine Ergebnisse gibt. Suche im vollständigen Atlas nach Aussagen, Versuchen und Originalquellen.', 'This does not mean there are no results. Search the full atlas for claims, attempts and original sources.')}</p><a href="${esc(href('search', query))}">${t('Im vollständigen Atlas suchen', 'Search the full atlas')} →</a></div>`}
        ${!exact && topicCoverage.length ? `<section class="research-unreviewed-topics"><h3>${t('Weitere Themen im Atlas', 'More topics in the atlas')}</h3><p>${t('Für diese Themen ist der geprüfte Entscheidungsstand noch nicht zusammengestellt. Ihre Aussagen, Ergebnisse und Originalquellen sind bereits zugänglich.', 'The reviewed decision state for these topics has not yet been compiled. Their claims, results and original sources are already accessible.')}</p><div>${topicCoverage.map(topic => `<a href="${esc(href('topic',topic.id))}">${esc(title(topic))}<small>${t('Einordnung ausstehend', 'Review pending')} ↗</small></a>`).join('')}</div></section>` : ''}
        ${query && !exact && problemIndex.length ? `<a class="research-full-search" href="${esc(href('search', query))}">${t('Auch einzelne Aussagen und Quellen durchsuchen', 'Also search individual claims and sources')} ↗</a>` : ''}
      </div>`;
      host.querySelector('[data-research-detail-next]')?.addEventListener('click', event => open({host, query, detailOffset: Number(event.currentTarget.dataset.researchDetailNext || 0), expectedRevision: packet.revision, api, isCurrent}));
    } catch (error) {
      if (!isCurrent()) return;
      host.innerHTML = `<div class="research-empty"><h2>${t('Forschungsstand derzeit nicht verfügbar', 'Research state currently unavailable')}</h2><p>${t('Der Wissensgraph bleibt über die Themen und die Suche erreichbar.', 'The knowledge graph remains accessible through topics and search.')}</p><p>${esc(error.message)}</p><button type="button" class="mini-button" data-research-retry>${t('Erneut laden', 'Retry')}</button> <a href="#discover">${t('Themen öffnen', 'Open topics')} ↗</a></div>`;
      host.querySelector('[data-research-retry]').addEventListener('click', () => open({host,query,detailOffset,expectedRevision,api,isCurrent}));
    }
  }
  window.TFPTResearch = {open};
})();
