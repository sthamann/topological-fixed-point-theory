(() => {
  'use strict';
  const dictionaries = {
    en: {
      'Zum Inhalt springen':'Skip to content','TFPT Wissensatlas':'TFPT Knowledge Atlas','TFPT Wissensatlas Startseite':'TFPT Knowledge Atlas home','Wissensatlas':'Knowledge Atlas',
      'Aufbau verstehen':'Understand the structure','Ein Weg durch die Theorie':'A path through the theory','Finde den Zusammenhang.':'Find the connection.',
      'Frage nach einer Aussage, folge ihrer Herleitung und sieh direkt, was belegt ist – und was noch offen bleibt.':'Ask about a claim, follow its reasoning, and see what is supported and what remains open.',
      'Den Wissensgraph durchsuchen':'Search the knowledge graph','Suchen':'Search','Enter zum Suchen':'Press Enter to search',
      'Suche in Aussagen, Ergebnissen und Quellen':'Search claims, results, and sources','Zum Beispiel: Wie entsteht die Feinstrukturkonstante?':'For example: How does the fine-structure constant arise?',
      'Umfang des Wissensgraphen':'Knowledge graph coverage','Der aktuelle Wissensstand':'Current graph contents','Extrahierte Inhalte aus dem Theoriegraphen':'Extracted content from the theory graph',
      'Aussagen insgesamt':'claims in total','structure-link':'Understand the structure','full-map-link':'Full map','coverage-link':'How to read this','reading-title':'A guide, not an endorsement',
      'reading-p1':'The atlas organizes claims, arguments, and test material from the existing theory graph. AI analysis helps locate and summarize material; it is',
      'reading-bold':'not expert review or proof','reading-p2':'. The original ledger status remains authoritative. Open questions are explicitly marked as open.',
      'Einträge':'Entries','Aussagen':'Claims','Quellen':'Sources','Wie ist das zu lesen?':'How should I read this?',
      'Forschungsmuster':'Research pattern','Arbeitsinferenz · ungeprüft':'Working inference · unreviewed','Gegenbefunde':'Counterevidence',
      'Suchergebnisse':'Search results','Treffer für':'matches for','Nach Relevanz sortiert':'sorted by relevance','Originalstatus bleibt maßgeblich':'The original status remains authoritative','Ledger':'Ledger','Quelle':'Source','Themen':'Topics','Der Wissensgraph wird durchsucht …':'Searching the knowledge graph …',
      'Keine passenden Einträge gefunden. Probiere einen anderen Begriff oder ein Themenfeld.':'No matching entries found. Try another term or topic.','Themenaspekt':'Theory aspect','Objekt':'Object','Prüfung':'Check','Orientierung':'Orientation',
      'Teilresultate / gemischter Status':'Partial / mixed results','Originaler Hinweis':'Original note',
      'Diese Übersicht zeigt gespeicherte Zusammenhänge und frühere Befunde; sie fügt keinen Beweis hinzu.':'This view surfaces recorded connections and earlier findings; it adds no proof.','Only repeated explicit seed IDs or non-claim-local object IDs are grouped. Equal labels, symbols, or numeric values do not establish a shared object.':'Only repeated explicit seed IDs or non-claim-local object IDs are grouped. Matching labels, symbols, or numeric values do not establish shared identity.',
      'Weitere Versuche laden':'Load more attempts','Befunde':'Findings','Bedingte Ergebnisse':'Conditional results','Ausgeschlossen':'Ruled out','Offen geblieben':'Left open','keine unabhängige Bestätigung':'not an independent confirmation','Stored attempt record; its verdict is not a proof of a linked claim.':'Stored attempt record; its verdict is not proof of a linked claim.','Recorded within this attempt\'s stated scope; not a general exclusion or proof.':'Recorded within this attempt\'s stated scope; not a general exclusion or proof.',
      'Eingaben mit offenen Herkunfts- oder Wirkungsfragen':'Inputs with unresolved origin or impact',
      'Mögliche Verbindung · kein Beweis':'Possible connection · not proof',
      'Versuche in der Auswahl':'attempts in this list',
      'Detailtexte bleiben im erfassten Originalwortlaut, wenn keine englische Übersetzung vorliegt.':'Detailed context text stays in its recorded language when an English translation is unavailable.',
      'Verknüpfte Aussagen im Kontext':'Linked claims in this context','Gemeinsame Eingaben werden nur im gezeigten Ausschnitt verglichen.':'Shared inputs are compared only within this displayed sample.',
      'Forschungsstand · vorhandene Evidenz und nächste Frage':'Research frontier · recorded evidence and next question','Forschungsstand wird zusammengestellt …':'Preparing the research context …','Forschungsstand konnte nicht geladen werden.':'Could not load research context.','Orientierung für die Recherche; gespeicherte Verknüpfungen und Kandidaten sind kein Beweis.':'Research orientation only. Stored links and candidates are not proof.','Research orientation only. Stored edges, common inputs, annotations, and authored patterns are surfaced with their recorded scope; no new proof or physical interpretation is inferred.':'Research orientation only. Stored edges, common inputs, annotations, and authored patterns are surfaced with their recorded scope; no new proof or physical interpretation is inferred.','Leitfrage':'Guiding question','Bisher festgehalten':'Recorded so far','Nächste trennscharfe Fragen':'Next discriminating questions','Gemeinsam genutzte Eingaben · gleiche gespeicherte IDs':'Shared inputs · matching recorded IDs','Gespeicherte Beziehungen · keine Beweiskette':'Recorded relationships · not a proof chain','gespeicherte Beziehung, kein Beweis':'recorded relationship, not proof','Versuche und Befunde':'Attempts and findings','Kandidat · ungeprüft':'Candidate · unreviewed','Ungeprüfte Anschlusskandidaten':'Unreviewed connection candidates','Begrenzte negative Befunde':'Scoped negative findings','Quellen für den Einstieg':'Starting sources','Der Forschungsrahmen ist begrenzt; weitere Einträge können über den vollständigen Wissensgraphen geöffnet werden.':'This research context is bounded; open the full knowledge graph for more entries.',
      'Themen und Suchergebnisse':'Topics and search results','Themen entdecken':'Explore topics','Themen werden geladen …':'Loading topics …',
      'Wissensgraph wird geladen':'Loading knowledge graph','Ein Kompass, kein Gütesiegel':'A guide, not a seal of approval',
      'Wo möchtest du anfangen?':'Where would you like to start?','Wähle ein Themenfeld oder suche direkt nach einer Aussage.':'Choose a topic or search for a claim.',
      'THEMENLANDKARTE':'TOPIC MAP','Einträge insgesamt':'entries in total','KI-ausgewertet':'AI analyzed','ausstehend':'pending','mit Extraktionsfehler':'with extraction errors','Skripte':'scripts',
      'Wissensgraph verbunden':'Knowledge graph connected','Verbindung zum Wissensgraphen fehlgeschlagen':'Could not connect to the knowledge graph','Themen konnten nicht geladen werden.':'Could not load topics.',
      'Der TFPT-Wissensatlas macht Aussagen, Herleitungen, Prüfprogramme und offene Fragen durchsuchbar.':'The TFPT Knowledge Atlas makes claims, derivations, test programs, and open questions searchable.',
      'Die Atlas-Oberfläche ist bereit, aber der lokale Wissensdienst antwortet nicht.':'The atlas is ready, but the local knowledge service is not responding.',
      'Erneut verbinden':'Reconnect','Suchtreffer':'Search results','Alle Themen':'All topics','Themenfelder':'Topics','Ergebnisse für':'Results for',
      'Theorieaufbau und Ziele werden geladen …':'Loading theory structure and goals …','Zielrichtungen':'Goals','Goal-Knoten aus dem Theoriegraphen':'Goal nodes from the theory graph',
      'Beginne bei den Ausgangsprinzipien und folge dann den acht offenen Zielrichtungen. Der Graph hält vorhandene Resultate und noch fehlende Herleitungen getrennt.':'Start with the initial principles, then follow the eight open goals. The graph separates recorded results from missing derivations.',
      'Zieldossier wird im Graphen geöffnet.':'Opening goal details in the graph.','Ausgangsprinzipien und die Frage, was die Theorie eindeutig festlegt.':'Initial principles and the question of what the theory determines uniquely.',
      'Die Ausgangsprinzipien und die Frage, was die Theorie eindeutig festlegt.':'Initial principles and the question of what the theory determines uniquely.',
      'Vom Compiler zur Physik':'From compiler to physics','Der Compiler bündelt endliche Relationen, Codes und Prüfungen. Eine vollständige physikalische Theorie erfordert zusätzlich die Herleitung von Raumzeit, Dynamik und beobachtbaren Größen aus den Ausgangsprinzipien.':'The compiler collects finite relations, codes, and checks. A complete physical theory also requires deriving spacetime, dynamics, and observables from the initial principles.',
      'Resultate, die gemeinsame Seeds oder dieselben Grunddaten verwenden, sind korreliert und zählen nicht automatisch als unabhängige Bestätigungen.':'Results that share seeds or underlying data are correlated and do not automatically count as independent confirmations.',
      'Aufbau und Ziele ansehen':'Explore structure and goals','Zum Themenfeld':'Open topic','Themenfeld öffnen':'Open topic','Zuletzt angesehen':'Recently viewed','aus dieser Sitzung':'from this session',
      'Der Atlas ordnet Aussagen, Argumente und Prüfmaterial aus dem vorhandenen Theoriegraphen. Eine KI-Auswertung hilft beim Auffinden und Zusammenfassen; sie ist':'The atlas organizes claims, arguments, and test material from the existing theory graph. AI analysis helps locate and summarize material; it is',
      'keine fachliche Abnahme und kein Beweis':'not expert review or proof','Der Status im ursprünglichen Ledger bleibt maßgeblich. Offene Fragen sind ausdrücklich als offen markiert.':'The original ledger status remains authoritative. Open questions are marked as open.',
      'im Original verzeichnet':'Recorded in the original','KI-Auswertung':'AI analysis','offen / ungeklärt':'Open / unresolved',
      'Orientierung im Theoriegraphen':'Navigate the theory graph','Eintrag':'Entry','Einordnung':'Overview','Lösungskette ansehen':'View reasoning chain',
      'Versuch / Befund':'Attempt / finding','Versuche und Befunde':'Attempts and findings','Belegte Resultate':'Recorded results','Teilresultate · offen':'Partial results · open',
      'Bedingte Anschlüsse':'Conditional links','Offene Fragen':'Open questions','Ausgeschlossen':'Excluded','Angenommen':'Assumed','Historisch':'Historical','Status unklar':'Status unclear',
      'Strukturkarte · keine Beweiskette':'Structure map · not a proof chain','Zielrichtungen':'Goals','Originalquellen dieses Themas':'Original sources for this topic',
      'Was wir haben':'What we have','Was wir suchen':'What we seek','Erforderliche Schritte':'Required steps','Bedingungen':'Conditions',
      'Was wurde gezeigt?':'What was shown?','Was bedeutet das?':'What does it mean?','Worauf baut es auf?':'What does it build on?',
      'Was bleibt offen?':'What remains open?','Prüfprogramme & Skripte':'Test programs & scripts','Lean-Verknüpfungen':'Lean links',
      'Original-Ledger':'Original ledger','Ledger-Felder einblenden':'Show ledger fields','Originalstatus im Ledger nachsehen.':'Check the original ledger for status.',
      'Versuche und Befunde sind dokumentierte Arbeitsschritte; sie zählen nicht automatisch als unabhängige Experimente.':'Attempts and findings are documented work; they do not automatically count as independent experiments.',
      'Dokumentierte Versuche zur Aussage; Einträge zählen nicht automatisch als unabhängige Experimente.':'Documented attempts related to this claim; entries do not automatically count as independent experiments.',
      'Noch keine explizit verknüpften Versuche oder Befunde.':'No explicitly linked attempts or findings yet.',
      'Weitere Versuche laden':'Load more attempts','Versuchsdossier öffnen':'Open attempt details','Fragestellung':'Question','Erfolgskriterium':'Success criterion',
      'Ergebnis':'Result','Status':'Status','Verbleibende Lücke':'Remaining gap','Bedingte Ergebnisse':'Conditional results','Geltungsbereich':'Scope','Befunde':'Findings','Offene Punkte':'Open points','Teilweise':'Partial','Bestanden':'Passed','Fehlgeschlagen':'Failed',
      'Ausschlusskriterien / Kill-Kriterien':'Exclusion / kill criteria','T1–T8-Ziele':'T1–T8 goals','T1–T8 bleiben offen; kein Gesamtabschluss.':'T1–T8 remain open; there is no overall completion.',
      'T1–T8-Gesamtabschluss im Vertrag nicht als offen markiert.':'The contract does not mark overall T1–T8 completion as open.',
      'Weitere Vertragsangaben':'More contract details','Vorgehen':'Method','Offene Lücke':'Open gap','Originalquellen':'Original sources',
      'Woher kommen die Eingaben?':'Where do the inputs come from?','Fehlende Beweisschritte':'Missing proof steps','Die Herkunft ist in den verfügbaren Angaben nicht geklärt.':'The origin is not clear from the available information.',
      'eingesetzt':'assigned','extern':'external','vorgelagerte Quelle':'upstream reference','vorgelagert belegt':'upstream reference','Herkunft offen':'origin unknown','im Skript abgeleitet':'derived in script',
      'eingesetzt · Axiom':'assigned · axiom','eingesetzt · freier Parameter':'assigned · free parameter','kein markierter Ableitungsmangel':'no derivation gap flagged',
      'Zirkularitätsrisiko':'circularity risk','Herleitung fehlt':'missing derivation','Prüfstatus':'Review status','Vorgelagerte Aussage':'Upstream claim','Ableitung':'Derivation',
      'Verwendet für':'Used for','Quellen':'Sources','Prüfhinweis':'review note','Fehlender Beweisschritt':'Missing proof step','Was fehlt':'What is missing',
      'Warum Checks nicht genügen':'Why checks are insufficient','Betroffene Ergebnisse':'Affected results','Eingabekandidaten':'Input candidates','Statische Eingabekandidaten':'Static input candidates',
      'Geprüfte Code-Herkunft':'Reviewed code origins','Die Einträge erklären, woher ein Skriptwert stammt und welche Rolle er dort hat. Sie belegen keine physikalische Herleitung oder allgemeine Gültigkeit.':'These entries describe where a script value comes from and the role it plays there. They do not establish a physical derivation or general validity.','Rolle im Skript':'Role in script','Auswahlgrund':'Selection rationale','Ursprungshinweis':'Origin note','Grenze dieser Einordnung':'Limits of this interpretation','Originalwortlaut · nicht übersetzt':'Original wording · not translated','Syntaxklassifikation':'Syntax classification',
      'diagnostischer Schwellenwert aus dem Skript':'diagnostic threshold chosen in the script','deklarierte Scantoleranz':'declared scan tolerance','Testgröße einer endlichen Matrix':'finite matrix test size','endliche Zylindermodell-Instanz':'finite cylinder-model instance','Anker unter ganzzahliger Aufspaltung':'anchor conditional on an integer split','extern verifizierte Höhe':'externally verified height','zitierte Konstante, hier ungenutzt':'cited constant, unused here','deklarierte Positionsfehlerschranke':'declared position-error budget','deklarierter Randsicherheitsfaktor':'declared tail margin','numerische Diskretisierung':'numerical discretization','Quadraturordnung':'quadrature order','endliche Modell-Instanz':'finite model instance','Komponente einer endlichen Struktur':'component of a finite structure','Standardwert eines Funktionsparameters':'function-parameter default','lokale Initialisierung oder Formel':'local initializer or formula','Zähler oder Diagnose':'counter or diagnostic','berechnete Prüfung':'computed check','Quelltextintegritätsprüfung':'source-integrity readback',
      'im Skript gewählter Testschwellenwert':'test threshold chosen in the script','bedingt durch ganzzahlige Aufspaltung und Typisierung':'conditional on integer-split and typing assumptions','zitierter externer Verifikationswert':'cited external verification input','klassische Quelle zitiert, hier ungenutzt':'classical source cited, unused here','deklarierte Diagnosetoleranz':'declared diagnostic budget','numerisches Budget ohne Zertifikat je Nullstelle':'numerical budget without a per-zero certificate','deklarierter numerischer Spielraum':'declared numerical margin','eine endliche Modell-Instanz':'one finite model instance','algorithmische oder strukturelle Definition':'algorithmic or structural definition','berechneter oder diagnostischer Wert':'computed or diagnostic value','explizite Quellcode-Angabe':'explicit code declaration','numerische Zuweisung':'numeric assignment','Literal im Quellcode':'literal in source code','Standardwert im Quellcode':'default value in source code','eingelesener Quelltext':'source-text read-back','aus Quellcode-Ausdruck berechnet':'computed from source expression',
      'Hinweis zur Skriptauswertung':'Script analysis note','Syntaxbeobachtungen aus dem Skriptbestand – noch keine semantische Bewertung von Herkunft, Verwendung oder Beweiswirkung.':'Syntax observations from the script corpus; not a semantic assessment of origin, use, or proof value.',
      'Ein gesetzter oder externer Wert kann eine Rechnung reproduzierbar machen, ohne seine Herleitung aus TFPT zu begründen. Ein vorgelagerter Verweis übernimmt nicht automatisch den Beweisstatus der verknüpften Aussage.':'An assigned or external value can make a calculation reproducible without deriving it from TFPT. An upstream reference does not automatically inherit the linked claim’s proof status.',
      'Diese Lücken markieren, wo ein Prüfprogramm die theoretische Herleitung nicht ersetzt.':'These gaps mark where a test program does not replace a theoretical derivation.',
      'Herkunft und Einsatz jeder erfassten Eingabe. Diese Einordnung beschreibt den Datenpfad; sie ist kein automatischer Beweis.':'Origin and use of each recorded input. This describes the data path; it is not proof by itself.',
      'belegt im Rahmen':'supported within scope','bedingt':'conditional','offen':'open','gemischt':'mixed','ausgeschlossen':'excluded','angenommen':'assumed','historisch':'historical','unklar':'unclear',
      'KI-Auswertung · nicht fachlich abgenommen':'AI analysis · not expert-reviewed','Auswertung fehlgeschlagen':'Analysis failed','Auswertung ausstehend':'Analysis pending',
      'Originalstatus aus Ledger':'Status from original ledger','Status aus Original-Ledger':'Status from original ledger','Originalstatus:':'Original status:',
      'Ergebnis im angegebenen Rahmen':'Result within the stated scope','Aufbau':'Structure','Eintrag wird geladen …':'Loading entry …','Eintragsansichten':'Entry views',
      'Einordnung':'Overview','Graphansicht':'Graph view','Eintrag':'Entry','Aussage im Original':'Original claim','Was wurde gezeigt?':'What was shown?',
      'Auf Originalquellen prüfen':'Check original sources','Der Eintrag enthält keine englische Zusammenfassung.':'No English summary is available for this entry.',
      'Aussage':'Claim','Herleitung':'Derivation','Prüfprogramm':'Test program','Lean-Satz':'Lean theorem','Voraussetzung':'Premise','Thema':'Topic','Struktur':'Structure','Ziel':'Goal','Quelle':'Source','Eingabe':'Input',
      'Prüfstatus':'Review status','Abhängige Eingaben':'Dependent inputs','Auswahlstatus':'Selection status','(Wert nicht angegeben)':'(value not specified)','Eingabe':'Input',
      'durch angegebene Herleitung festgelegt':'fixed by the stated derivation','angenommen oder ausgewählt':'assumed or selected','nicht anwendbar':'not applicable','unklar':'unclear',
      'Die Ausgangsprinzipien und die Frage, was die Theorie eindeutig festlegt.':'Initial principles and what the theory determines uniquely.',
      'Endliche Codestrukturen, E₈ und ihre möglichen Beziehungen im Theoriegraphen.':'Finite code structures, E₈, and their possible relationships in the theory graph.',
      'Wie ineinander verschränkte Uhren und ihre Relationen beschrieben werden.':'How interlocking clocks and their relations are described.',
      'Ladungen, Teilchenfamilien und die Einordnung des Standardmodells.':'Charges, particle families, and the place of the Standard Model.',
      'Herleitungen, Werte und Prüfungen rund um die Feinstrukturkonstante.':'Derivations, values, and checks concerning the fine-structure constant.',
      'Quantenstruktur, Randbedingungen und die Auswahl physikalischer Zustände.':'Quantum structure, boundary conditions, and selection of physical states.',
      'Aussagen und offene Herleitungen zu Raumzeit und Gravitation.':'Claims and open derivations concerning spacetime and gravity.',
      'Kosmologische Zusammenhänge, Konsequenzen und beobachtbare Fragen.':'Cosmological connections, consequences, and observable questions.',
      'Verbindungen zur Riemann-Hypothese und der zugehörigen Zahlentheorie.':'Connections to the Riemann hypothesis and related number theory.',
      'Prüfprogramme, Vorhersagen, Skripte und formale Lean-Verknüpfungen.':'Test programs, predictions, scripts, and formal Lean links.',
      'benötigt als Voraussetzung':'requires as a premise','gemeinsam erforderlich':'jointly required','leitet her':'derives','hat Ergebnis':'has result','hat Herleitung':'has derivation',
      'implementiert':'implements','prüft':'checks','Formalisierungskandidat':'formalization candidate','verwendet im Geltungsbereich':'uses within scope','setzt Eingabe':'sets input','ordnet als Eingabe ein':'interprets as input',
      'Belegt im angegebenen Rahmen':'Supported within the stated scope','Geltende Prüfregel':'Active methodological rule','Bedingt':'Conditional','Offen':'Open','KI-Auswertung · ungeprüft':'AI analysis · unreviewed',
      'benötigt':'requires','Verknüpfte Einträge':'Related entries','Direkte Originalbelege':'Direct original sources','Weitere verknüpfte Quellen':'More related sources',
      'Eintrag und Verknüpfungen werden geladen …':'Loading entry and links …','Themenfeld wird geöffnet …':'Opening topic …','Verknüpfungen werden verfolgt …':'Following links …',
      'Quelle wird geladen …':'Loading source …','Erneut laden':'Reload','Erneut suchen':'Search again','Details öffnen':'Open details','Weitere Einträge laden':'Load more entries',
      'Für diese Gruppe sind derzeit keine Einträge erfasst.':'There are currently no entries recorded for this group.',
      'Daten aktualisieren':'Refresh data','Aktualisieren':'Refresh','Neue Daten verfügbar':'New data available','Keine Änderungen':'No changes','Wird aktualisiert …':'Refreshing …',
      'Aktualisierung fehlgeschlagen':'Refresh failed','Suchtreffer':'Search results','Originalstatus':'Original status','siehe Ledger':'see ledger',
      'Für dieses Themenfeld ist keine englische Übersetzung vorhanden; der Originaltext bleibt sichtbar.':'No English translation is available for this topic; the original wording is shown.',
      'Wird geladen …':'Loading …','Laden fehlgeschlagen · erneut versuchen':'Loading failed · try again','Befund nicht klassifiziert':'Finding not classified','Vorgehen':'Method','Offene Lücke':'Open gap',
      'Offene Frage / ungeklärt':'Open question / unresolved','Herkunfts-/Herleitungshinweis':'origin / derivation note','Beweisschritte offen':'proof steps open',
      'Extraktion braucht Prüfung':'Extraction needs review','Ein Qualitätsaudit hat mögliche Probleme in der KI-Extraktion markiert. Diese Hinweise sind keine fachliche Bewertung des ursprünglichen Anspruchs.':'A quality audit flagged possible issues in the AI extraction. These notes are not an expert assessment of the original claim.',
      'Für diesen Eintrag ist keine englische Übersetzung vorhanden; der Originaltext bleibt sichtbar.':'No English translation is available for this entry; the original wording is shown.',
      'Quelle':'Source','PDF-Seite':'PDF page','Zeilen':'Lines','Quellenlink kopieren':'Copy source link','Lösungskette':'Reasoning chain',
      'Herkunft und Einsatz jeder erfassten Eingabe. Diese Einordnung beschreibt den Datenpfad; sie ist kein automatischer Beweis.':'Origin and use of each recorded input. This describes the data path; it is not proof by itself.',
      'Diese Lücken markieren, wo ein Prüfprogramm die theoretische Herleitung nicht ersetzt.':'These gaps mark where a test program does not replace a theoretical derivation.',
      'Syntaxbeobachtungen aus dem Skriptbestand – noch keine semantische Bewertung von Herkunft, Verwendung oder Beweiswirkung.':'Syntax observations from the script corpus; not a semantic assessment of origin, use, or proof value.',
      'Aussage im Original':'Original claim','Zusammenhänge':'Connections','Erkenntnisstatus':'Epistemic status',
      'Was wurde geprüft?':'What was checked?','Geltungsbereich':'Scope','Warum belegt? / Was fehlt?':'Why supported? / What is missing?',
      'Weitere Stufe laden':'Load another level','Einpassen':'Fit view','Vergrößern':'Zoom in','Verkleinern':'Zoom out',
      'Querverweise (keine Herleitung)':'Cross-references (not derivations)','Begrenzter Ausschnitt':'Limited view','Alle Voraussetzungen':'all premises',
      'Geltungsbereich nicht angegeben':'Scope not specified',
      'Ausgewählte Aussage':'Selected claim','Herleitungsschritte':'reasoning steps','Tiefe':'depth','Einige Voraussetzungen und Herleitungsschritte sind hier nicht vollständig angezeigt. Daraus folgt keine Aussage über ihre Beweiskraft.':'Some premises and reasoning steps are not fully displayed here. This does not assess whether they are proved.',
      'Betroffen':'Affected','Die Anzeige ist begrenzt; der vollständige Umfang ist hier nicht bestimmt.':'This view is limited; the full extent is not determined here.',
      'Pfeile zeigen die tatsächliche Richtung. „Benötigt als Voraussetzung“, „hat Ergebnis“ und „hat Herleitung“ beschreiben Strukturbezüge; nur „gemeinsam erforderlich“ und „leitet her“ zählen hier als Herleitungsschritte. Knoten öffnen ihr Dossier.':'Arrows show the stored direction. “Requires as a premise,” “has result,” and “has derivation” describe structural links; only “jointly required” and “derives” count as reasoning steps here. Open a node to view its details.',
      'Die Anfrage ist fehlgeschlagen':'The request failed','Prüfumfang':'Review scope','statisch erfasste Lean-Sätze':'statically cataloged Lean theorems',
      'Gesamtkarte wird geladen …':'Loading the full map …','Gesamtkarte ansehen':'View the full map','Gesamtkarte':'Full map','Alle Stufen':'All stages','Fokus auf Stufe':'Focus on stage',
      'Knoten anklicken: Erklärung und Quellen anzeigen':'Select a node to see its explanation and sources','Knoten':'Nodes','Verknüpfungen':'Links','Gegenbefunde':'Counterevidence',
      'Wähle einen Knoten, um Details zu sehen.':'Select a node to see details.','Die Gesamtkarte konnte nicht geladen werden.':'Could not load the full map.',
      'Für diese Karte ist keine englische Übersetzung vorhanden; der Originaltext bleibt sichtbar.':'No English translation is available for this map entry; the original wording is shown.',
      'Begrenzt':'Limited','Die Anzeige ist begrenzt; der vollständige Umfang ist hier nicht bestimmt.':'This view is limited; the full extent is not determined here.',
      'Ein Kompass, kein Gütesiegel':'A guide, not an endorsement','Aktualisiert':'Updated',
      'Details unten':'Details below','Geltungsbereichsnotiz mangels englischer Übersetzung im Originalwortlaut.':'Scope note shown in its original wording because no English translation is available.','Geltungsbereichsnotiz ist nur im Originalwortlaut verfügbar.':'Scope note is available only in its original wording.',
      'Auch angrenzende Themen zeigen':'Include related topics','Verbindung zur Theorie':'Connection to the theory','Nächste fehlende Brücke':'Next missing bridge','Status nicht automatisch zugeordnet':'Status not automatically classified',
      'Gerichtete Lösungskette':'Directed reasoning chain','Detaillierte Erläuterungen bleiben mangels Übersetzung in ihrer Originalsprache sichtbar.':'Detailed annotations are shown in their original language where a translation is unavailable.','Zugehörige Aussage':'Owning claim',
      'angezeigt':'shown','öffnen':'open','Berichteter Ledgerstatus im genannten Geltungsbereich; keine neue fachliche Abnahme oder Ausführung.':'Ledger-reported status within the stated scope; no new expert review or execution.',
      'Zur Themenübersicht':'Back to topic overview','Quellenlink kopieren':'Copy source link','Quellenlink wurde kopiert.':'Source link copied.','Zeile':'Line','Zeilen':'Lines','PDF-Seite':'PDF page'
    }
  };
  const de = Object.fromEntries(Object.entries(dictionaries.en).map(([g,e])=>[e,g]));
  const directLabels={
    'research-link':{de:'Forschungsstand',en:'Research state'},
    'research-hero-eyebrow':{de:'Vom Befund zur nächsten Frage',en:'From findings to the next question'},
    'research-hero-title':{de:'Baue auf dem auf, was wir wissen.',en:'Build on what we already know.'},
    'research-hero-intro':{de:'Sieh, was bereits geklärt ist, welche Wege unter welchen Voraussetzungen ausscheiden und welcher Schritt jetzt weiterführt.',en:'See what is established, which approaches fail under which assumptions, and what the next step can resolve.'},
    'research-search-hint':{de:'Eine Frage oder ein Thema eingeben',en:'Enter a question or topic'},
    'research-search-mode':{de:'Forschungsstand',en:'Research state'},
    'research-entries-mode':{de:'Alle Aussagen & Quellen',en:'All claims & sources'},
    'structure-link':{de:'Aufbau verstehen',en:'Understand the structure'},'full-map-link':{de:'Gesamtkarte',en:'Full map'},'coverage-link':{de:'Wie ist das zu lesen?',en:'How to read this'},
    'reading-title':{de:'Ein Kompass, kein Gütesiegel',en:'A guide, not an endorsement'},
    'reading-p1':{de:'Der Atlas ordnet Aussagen, Argumente und Prüfmaterial aus dem vorhandenen Theoriegraphen. Eine KI-Auswertung hilft beim Auffinden und Zusammenfassen; sie ist',en:'The atlas organizes claims, arguments, and test material from the existing theory graph. AI analysis helps locate and summarize material; it is'},
    'reading-bold':{de:'keine fachliche Abnahme und kein Beweis',en:'not expert review or proof'},
    'reading-p2':{de:'. Der Status im ursprünglichen Ledger bleibt maßgeblich. Offene Fragen sind ausdrücklich als offen markiert.',en:'. The original ledger status remains authoritative. Open questions are explicitly marked as open.'}
  };
  const attrs = {'placeholder':{'Zum Beispiel: Wie entsteht die Feinstrukturkonstante?':'For example: How does the fine-structure constant arise?'},'aria-label':{'Sprache / Language':'Language'},'title':{'Daten aktualisieren':'Refresh data'},'content':{'Der TFPT-Wissensatlas macht Aussagen, Herleitungen, Prüfprogramme und offene Fragen durchsuchbar.':'The TFPT Knowledge Atlas makes claims, derivations, test programs, and open questions searchable.'}};
  let lang = localStorage.getItem('tfpt-atlas-language') === 'en' ? 'en' : 'de';
  document.documentElement.lang=lang;
  const sourceSelector = '.source-viewer, .source-text, pre, code.source-content, [data-original-content], .result-summary, .theme-entry-content, .recent-copy small, .entry-statement';
  function translateText(value) {
    const text=String(value??''); if(lang==='de') return de[text]||text;
    if(dictionaries.en[text]) return dictionaries.en[text];
    let out=text;
    out=out.replace(/^(\d+) Einträge$/, '$1 entries').replace(/^(\d+) Quellen$/, '$1 sources').replace(/^(\d+) Aussagen$/, '$1 claims').replace(/^(\d+) Themen$/, '$1 topics').replace(/^Quelle (\d+)$/, 'Source $1');
    const counted=out.match(/^(.*) · (\d+)$/);if(counted&&dictionaries.en[counted[1]])out=`${dictionaries.en[counted[1]]} · ${counted[2]}`;
    out=out.replace(/^Ergebnisse für /, 'Results for ').replace(/^Originalquellen dieses Themas · /, 'Original sources for this topic · ').replace(/^Weitere (\d+) Quellen anzeigen$/, 'Show $1 more sources').replace(/^Weitere (\d+) Einträge anzeigen$/, 'Show $1 more entries').replace(/^Belegte Resultate · /, 'Recorded results · ').replace(/^Teilresultate · offen · /, 'Partial results · open · ').replace(/^Bedingte Anschlüsse · /, 'Conditional links · ').replace(/^Offene Fragen · /, 'Open questions · ').replace(/^Versuche und Befunde · /, 'Attempts and findings · ');
    return out;
  }
  function translateTree(root=document) {
    if(!root)return;
    const direct=root.matches?.('[data-i18n]')?[root]:[];for(const el of [...direct,...(root.querySelectorAll?.('[data-i18n]')||[])]){const key=el.getAttribute('data-i18n'),entry=directLabels[key];if(entry&&el.textContent!==entry[lang])el.textContent=entry[lang];}
    const walker=document.createTreeWalker(root,NodeFilter.SHOW_TEXT);
    let n; while((n=walker.nextNode())) { if(n.parentElement?.closest(sourceSelector)||n.parentElement?.closest('[data-i18n]'))continue;const v=n.nodeValue,t=translateText(v);if(t!==v)n.nodeValue=t; }
    const elements=root.querySelectorAll?.('[placeholder],[title],[aria-label],meta[name="description"]')||[];
    for(const el of elements) for(const attr of ['placeholder','title','aria-label',...(el.matches('meta[name="description"]')?['content']:[])]) if(el.hasAttribute(attr)){const v=el.getAttribute(attr),t=translateText(v);if(t!==v)el.setAttribute(attr,t);}
  }
  function setLanguage(value) {
    lang=value==='en'?'en':'de';localStorage.setItem('tfpt-atlas-language',lang);document.documentElement.lang=lang;translateTree();document.dispatchEvent(new CustomEvent('tfpt-language-change',{detail:{lang}}));
  }
  document.addEventListener('DOMContentLoaded',()=>{
    const select=document.getElementById('language-select'); if(select){select.value=lang;select.addEventListener('change',()=>setLanguage(select.value));}
    translateTree();
    const observer=new MutationObserver(records=>{for(const record of records){if(record.type==='childList')for(const node of record.addedNodes)if(node.nodeType===Node.ELEMENT_NODE||node.nodeType===Node.TEXT_NODE){const target=node.nodeType===Node.TEXT_NODE?node.parentElement:node;if(target?.isConnected)translateTree(target);}}});
    observer.observe(document.body,{subtree:true,childList:true});
  });
  window.TFPT_I18N={get language(){return lang;},t:translateText,setLanguage};
})();
