/* TFPT Landkarte · kuratierte Inhalte, Teil 2: Muster, Ursprung von P1/P2, Universalraum,
   Origin Idee, OPH und andere Theorien, noch nicht ausgeführte Tests.
   Lesehilfe, kein Beweis. Maßgeblich bleiben Ledger, Originalquellen und der lebende Atlas. */
(function () {
  const N = (o) => window.TFPT.nodes.push(o);
  const FP = "docs/research/fundamental-plan/";
  const RAD = "experiments/theory-contracts/radical-global-hypotheses-20261118/README.md";
  const F1S = "research/originlab/F1-ORIGIN-SCENARIOS-20261116.txt";

  /* ───────────── c) Muster und Verbindungen ───────────── */
  const P = (o) => N(Object.assign({ k: "pat" }, o));

  P({ id: "pat-code", grp: "Information", t: "E₈ ist ein Fehlerkorrekturcode", st: "belegt",
    s: "Der Hamming Code [8,4,4] baut über den gaußschen Zahlen genau das E₈ Gitter.",
    b: ["Construction A über Z[i]: Der Quotient L/(1+i)L ist F₂⁴, also vier Bits. Die 240 Wurzeln fallen in 15 Fasern zu je 16.",
      "Idee dahinter: Die Struktur, die die Teilchenladungen ordnet, ist dieselbe, die Fehler erkennt und behebt. Ob die Natur diese Eigenschaft benutzt, ist nicht gezeigt."],
    n: ["8", "4", "16", "15", "240"], g: ["Code", "Hamming", "Gitter", "Information"],
    l: [["qcode", "karte"], ["pat-qubit", "karte"]], src: ["note_e8_gaussian_code.tex"] });

  P({ id: "pat-qubit", grp: "Information", t: "Zwei Qubits im Gitter", st: "belegt",
    s: "60 gaußsche Linien sind die 60 Stabilisatorzustände zweier Qubits, geordnet in 15 Messkontexten.",
    b: ["Die 15 Kontexte bilden das „Doily“ (verallgemeinertes Viereck GQ(2,2)). Die Spiegelungsgruppe G₃₁ bildet auf die Cliffordgruppe ab.",
      "Grenze: Stabilisatorzustände allein sind klassisch simulierbar. Ein universeller Quantenrechner folgt daraus nicht.",
      "Geordnete E₈ Ereignisse auf dem Doily geben Erfolgsraten 1/3 gegen 1, je nach Reihenfolge."],
    n: ["60", "15", "2", "1/3", "46080"], g: ["Qubit", "Clifford", "Kontext", "Information"],
    l: [["qcode", "karte"], ["clocks", "gleich"], ["pat-galois", "gleich"]], src: ["note_e8_gaussian_code.tex"] });

  P({ id: "pat-loss", grp: "Information", t: "Auslesen verliert Information", st: "belegt",
    s: "Warnmuster: Zwei verschiedene Ursprünge können denselben Schatten werfen.",
    b: ["Mittelung über Geschichten, Projektion auf einen Quotienten und Spurbildung zerstören genau die Unterschiede, die man für den Ursprung braucht.",
      "Folge für die Suche: Aus einer fertigen Auslese lässt sich die Quelle nicht eindeutig zurückrechnen. Man braucht Zukunftsantworten oder die ursprüngliche Reihenfolge."],
    n: [], g: ["Warnung", "Schatten", "Information", "Gedächtnis"],
    l: [["ur-shadow", "karte"], ["dyn-wordsource", "karte"]], src: [FP + "TFPT-Gesamtkonsolidierung-2026-11-12.md"] });

  P({ id: "pat-prime", grp: "Information", t: "Primzahlen und π als Ereignisprotokoll", st: "aus",
    s: "Idee: Die Ziffern von π oder die Primzahlen sind das Logbuch der Ausführung. Der Test vom 8. September 2026 fand kein Signal.",
    b: ["π lieferte 14 Treffer, der Nullvergleich im Mittel 24,5 (95 % Bereich 16 bis 34). e gab 24, √2 gab 27. π liegt also eher unter dem Zufall.",
      "Der Test prüfte nur kleine Faktoren. Die Riemann Arbeiten laufen als eigenes Projekt (Prime Front) und zählen nicht als Beleg für TFPT.",
      "Was als Bild bleibt: Ein Ereignisprotokoll braucht eine Regel, die Ereignisse erzeugt. Genau die fehlt."],
    n: ["14", "24,5"], g: ["Primzahlen", "Ereignis", "Protokoll", "Zufall"],
    l: [["ext-prime", "karte"], ["dyn-wordsource", "gleich"]], src: ["experiments/pi-prime-event-log-2026-09-08/run.log"] });

  P({ id: "pat-hyper", grp: "Geometrie", t: "Ikosaeder, McKay und (2,3,5)", st: "belegt",
    s: "Dieselben drei Zahlen 2, 3, 5 erscheinen als Drehachsen des Ikosaeders, als Arme des E₈ Diagramms und als Primteiler von 30.",
    b: ["Das Ikosaeder als 3 uniformer Hypergraph hat 120 Symmetrien. Die McKay Marken (1,2,3,4,5,6,4,2,3) sind der Attraktor der Regel m ← (A+2I)/4·m.",
      "Die Reihe T(2,3,r) hat Determinante 3, 2, 1, 0, −1: Genau bei r = 5 wird sie 1, also unimodular. Davor zu locker, danach unendlich."],
    n: ["2", "3", "5", "30", "120"], g: ["Ikosaeder", "McKay", "Hypergraph", "Attraktor"],
    l: [["dyn-hypergraph", "karte"], ["ext-oph", "karte"], ["pat-growth", "karte"]], src: ["hypergraph.md"] });

  P({ id: "pat-growth", grp: "Geometrie", t: "Wachstum stoppt genau bei E₈", st: "belegt",
    s: "Lässt man ein Diagramm mit drei Armen wachsen und verlangt Endlichkeit, ist bei E₈ Schluss.",
    b: ["Spektralradius ρ ≤ 2 trennt endliche von unendlichen Diagrammen. E₈ ist das letzte endliche in der Reihe.",
      "Gleichwertig gelesen: „kein freier abelscher Sektor“ ⟺ det K = 1."],
    n: ["8", "2", "1"], g: ["Wachstum", "Endlichkeit", "Netzwerk"],
    l: [["e8", "gleich"], ["org-p2-det", "karte"]], src: ["hypergraph.md", "tfpt_3_e8_audit_bootstrap.tex"] });

  P({ id: "pat-picard", grp: "Geometrie", t: "Picardkurve y³ = x⁴ − 1", st: "belegt",
    s: "Eine Kurve vom Geschlecht 3 mit 48 Symmetrien, von denen die markierte Auslese nur 12 sieht.",
    b: ["Sie verbindet die Dreierdrehung (Familie) mit der Viererdrehung (Marken). Auch Kummerfläche und K3 treten in derselben Nachbarschaft auf."],
    n: ["3", "4", "48", "12"], g: ["Picard", "Kurve", "Auslese"],
    l: [["dyn-reader", "karte"], ["viermarken", "gleich"]], src: [F1S] });

  P({ id: "pat-galois", grp: "Zahlen", t: "Getriebe: Galois und Invarianten", st: "belegt",
    s: "φ(30) = 8 Exponenten, Invariantengrade 8, 12, 20, 24, und die Vier als Schaltwerk des Fünfertakts.",
    b: ["Die acht zu 30 teilerfremden Zahlen 1, 7, 11, 13, 17, 19, 23, 29 sind die Exponenten von E₈. G₃₁ teilt mit E₈ die Grade 8, 12, 20, 24."],
    n: ["30", "8", "12", "20", "24"], g: ["Galois", "Coxeter", "Clock"],
    l: [["clocks", "karte"], ["pat-qubit", "gleich"]], src: ["origin_theory.tex"] });

  P({ id: "pat-sweep", grp: "Zahlen", t: "Reihenprüfer: 141 Hypothesen auf einmal", st: "belegt",
    s: "Ein Programm bewertet ganze Hypothesenfamilien automatisch. Die Zahlen der Theorie sind im erklärten Suchraum einzig.",
    b: ["Gebaut am 19.11.: sweep.py prüft Familien zu Schnittgeometrie, Zustand, Exponent, Singularität, Markenlage, Zahlen und Kubik in zwei Sekunden gegen dieselben festen Zielgrößen.",
      "Ergebnis: Unter 95000 Kombinationen liegt nur φ_tree = (4/3)c₃ mit δ_top = 48c₃⁴ innerhalb 3 σ. Von 25 Gitterklassen trägt nur die Spinorklasse von D₅. In der Kubik tragen nur 2 und 8·41/10. Von 46 Singularitätsschnitten gibt nur Z₅ die Fünferuhr.",
      "Grenzen: „trägt“ heißt widerspruchsfrei, nicht bewiesen. Die Einzigkeit gilt nur im erklärten Suchraum. Eigene Hypothesen lassen sich über eine Datei ergänzen."],
    n: ["141", "95000", "4/3", "48"], g: ["Suche", "Einzigkeit", "Alpha"],
    l: [["pat-a4cut", "gleich"], ["alpha", "dockt"]], src: ["experiments/theory-contracts/hypothesis-sweep-20261119/README.md"] });

  P({ id: "pat-a4cut", grp: "Zahlen", t: "Die Fünferuhr sitzt am A₄ Schnitt", st: "belegt",
    s: "Der Phasenterm der originalen P1 Wirkung gibt nur am Z₅ Schnitt eine Uhr. E₈ selbst gibt keine.",
    b: ["Gerechnet am 19.11.: Die primitive Wirkung enthält die Phase (iπ/2)·Δη. Für den Schnitt S³/Z_p ist η = (p²−1)/(6p). Die Phase hat Ordnung 5 bei Z₅, Ordnung 32 bei Z₄, Ordnung 16 bei RP³.",
      "Alle anderen Singularitäten: E₆, E₇, E₈ geben Ordnung 576, 1152, 2880. Die D Reihe gibt keine Uhr. Nur A₄ trifft, passend zu „Ordnung 5 braucht Rang 4“.",
      "Muster am Z₄ Schnitt: Die zwei Spinstrukturen geben 5/8 und 3/8. Das sind die Gewichte der Klebefelder von D₅ und A₃. Ihre Summe ist nur bei p = 4 genau 1. Am Z₅ Schnitt ist η = 4/5, die Norm der A₄ Klebeklasse. Eine gemeinsame Regel für beide fehlt, deshalb nur Muster.",
      "Grenzen: Nicht gezeigt ist, dass P1 an einem solchen Schnitt sitzt. Der gewichtete Originalschnitt bleibt unberechnet."],
    n: ["5", "4/5", "5/8", "3/8", "2880"], g: ["Uhr", "Singularität", "Phase"],
    l: [["pat-clockrank", "gleich"], ["dyn-reset", "dockt"]], src: ["experiments/theory-contracts/p1-kernel-ground-state-20261119/README.md"] });

  P({ id: "pat-clockrank", grp: "Zahlen", t: "Jede Uhr braucht ihren Rang", st: "belegt",
    s: "In der Ebene gehen nur die Takte 2, 3, 4, 6. Der Fünfertakt braucht Rang 4, und in E₈ verschmelzen je zwei Uhren, alle drei nie.",
    b: ["Gerechnet am 19.11.: Die Takte 15, 20, 24 und 30 sind in E₈ jeweils eine einzige primitive Uhr vom Rang 8. Ein Takt 60 (drei, vier und fünf zugleich) existiert in E₈ nicht; als eine Uhr bräuchte er Rang 16.",
      "Deshalb greift die Vier als Getriebe in die Fünf ein: Sie ist deren Galoisgruppe, und die Getrieberelation ist die Frobeniusgruppe F₂₀ in S₅.",
      "Die Vier wählt das Quadrat, die Drei wählte das Sechseck. Am Sechseck wäre c₃ um 2/√3 größer und α⁻¹ = 104,4. Kohärent verdoppelte Keime gäben 37 oder 512."],
    n: ["4", "5", "30", "8", "16", "60", "20"], g: ["Clock", "Galois", "Auswahl", "Gitter"],
    l: [["clocks", "karte"], ["taui", "karte"], ["pat-galois", "gleich"], ["pat-picard", "gleich"]], src: ["experiments/theory-contracts/landkarte-origin-scenarios-20261119/README.md"] });

  P({ id: "pat-41", grp: "Zahlen", t: "Die 41", st: "belegt",
    s: "41 = aᵀKa = 10·b₁ = Summe der Leptonen und Higgsbeiträge = 5² + 4².",
    b: ["Dieselbe Zahl als Operatorwert, als Laufkoeffizient der Hyperladung und als gaußsche Norm von 5 + 4i. Drei Lesarten, eine Abbildung zwischen ihnen ist nur teilweise gezeigt.",
      "Dazu die Nachbarn 7, 13, 55 als Normen."],
    n: ["41", "5", "4", "10"], g: ["Alpha", "Norm", "Hyperladung"],
    l: [["alpha", "karte"], ["famop", "karte"]], src: ["tfpt_2_standard_model.tex", RAD] });

  P({ id: "pat-54", grp: "Zahlen", t: "Der Exponent 5/4", st: "bedingt",
    s: "5/4 = Norm des D₅ Klebevektors = 2h mit h = 5/8 = ein Viertel je Cartanrichtung.",
    b: ["Als Gitterzahl exakt. Der Term erscheint im Archiv erstmals in Fassung 2.5 (27.01.2026). Das ist eine Reihenfolge, kein Fit: Der angepasste Wert wäre 1,2404 ± 0,0051, und 5/4 liegt 1,9 σ daneben.",
      "Offene Frage dahinter: Warum liest der Operator die Zeit am ganzen Klebevektor, die Oszillatorgrade aber nur an den fünf Trägerrichtungen?"],
    n: ["5/4", "5/8", "3/4", "5"], g: ["Alpha", "Verklebung", "Chronologie"],
    l: [["alpha", "karte"], ["e8", "gleich"], ["gap-operator", "offen"]], src: [RAD, "experiments/theory-contracts/origin-ideas-alpha-recheck-20261119/README.md"] });

  P({ id: "pat-48", grp: "Zahlen", t: "Die 48 und die Acht", st: "belegt",
    s: "48 = 16·3 Zustände, 48c₃⁴ = 1/(2⁸·Vol S⁷), und die Acht ist dreifach erzwungen.",
    b: ["Acht als Rang von E₈, als 1/c₃ = 8π und als p₄ − p₃. Dazu 48 Symmetrien der Picardkurve. Gleiche Zahl heißt hier noch nicht gleiche Sache."],
    n: ["48", "16", "3", "8"], g: ["Seed", "Topologie", "Zählung"],
    l: [["phi0", "karte"], ["sm", "gleich"], ["pat-picard", "gleich"]], src: ["tfpt_2_standard_model.tex", "origin_theory.tex"] });

  P({ id: "pat-23", grp: "Zahlen", t: "Zwei Drittel überall", st: "bedingt",
    s: "2/3 als Rückgewinnungsrate, als Nariai Verhältnis und in der Koide Formel.",
    b: ["Drei verschiedene Messgrößen mit demselben Wert. Eine gemeinsame Ursache ist eine Vermutung, keine Folge."],
    n: ["2/3", "1/3"], g: ["Attraktor", "Koide", "Horizont"],
    l: [["attraktor", "karte"], ["ext-koide", "gleich"]], src: ["tfpt_horizon_readouts.tex"] });

  P({ id: "pat-1920", grp: "Zahlen", t: "1920 = 3·4·8·20", st: "belegt",
    s: "Das Produkt der vier Familiendeterminanten ist die Ordnung der Weylgruppe von D₅.",
    b: ["Flavor und Träger zählen dieselbe Symmetrie. Als Zählgleichheit exakt, als Abbildung offen."],
    n: ["1920", "3", "4", "8", "20"], g: ["Flavor", "Weyl", "Zählung"],
    l: [["famop", "karte"], ["d5", "gleich"]], src: ["tfpt_2_standard_model.tex"] });

  P({ id: "pat-modular", grp: "Zahlen", t: "c₃ als Modulform", st: "belegt",
    s: "c₃ = E₂(i)/24 und 1/(6π) = E₂(i)/18. Die Normierung ist ein Wert am quadratischen Torus.",
    b: ["Das macht aus einer Kugelzahl eine Toruszahl. Die Brücke zwischen Kugel und Torus über dieselbe Viermarkenüberlagerung ist noch nicht gebaut."],
    n: ["1/(8π)", "24", "18", "1/(6π)"], g: ["Modul", "Torus", "Normierung"],
    l: [["org-c3-modular", "karte"], ["taui", "karte"], ["gap-torus", "offen"]], src: [RAD] });

  P({ id: "pat-seed", grp: "Warnmuster", t: "Gemeinsamer Seed, korrelierte Treffer", st: "belegt",
    s: "Warnmuster: Zehn Treffer aus einem Seed sind kein zehnfacher Beleg.",
    b: ["Fast alle Flavor und Kosmologiezahlen lesen φ₀. Stimmt φ₀ ungefähr, stimmen sie gemeinsam. Unabhängig sind nur Tests, die verschiedene Eingaben prüfen."],
    n: [], g: ["Warnung", "Seed", "Korrelation"], l: [["phi0", "karte"]], src: ["semantic_explorer/skills/tfpt-knowledge/SKILL.md"] });

  P({ id: "pat-closure", grp: "Warnmuster", t: "Abschluss ist nicht Ursprung", st: "belegt",
    s: "Warnmuster: Dass alles zusammenpasst, erklärt nicht, woher es kommt.",
    b: ["Der Rückbootstrap zeigt, dass g = 5 die einzige stimmige Wahl ist. Er zeigt nicht, welcher Vorgang sie trifft. Diagnose vom 16. November: P1, P2 und E₈ legen die vier Marken nicht fest. Der Produzent liegt außerhalb der erklärten Eingaben."],
    n: ["5", "4"], g: ["Warnung", "Bootstrap", "Ursprung"],
    l: [["bootstrap", "karte"], ["org-diagnose", "karte"]], src: [FP + "TFPT-Ursprungshypothesen-2026-11-12.md"] });

  /* ───────────── d) Ursprung von P1 und P2 ───────────── */
  const O = (o) => N(Object.assign({ k: "org" }, o));

  O({ id: "org-c3-geo", lane: "p1", t: "c₃ aus reiner Geometrie", st: "bedingt",
    s: "8π = |Z₂| · 2π · χ(S²): zwei Ufer, ein Umlauf, Eulerzahl der Kugel.",
    b: ["Einfachste Lesart. Jeder Faktor ist benennbar. Offen ist, warum genau diese drei Faktoren multipliziert werden."],
    n: ["8", "2", "1/(8π)"], g: ["Normierung", "Kugel", "Naht"], l: [["p1", "testet"]], src: ["origin_theory.tex"] });

  O({ id: "org-c3-thermo", lane: "p1", t: "c₃ aus Horizontthermodynamik", st: "bedingt",
    s: "Jacobson: Verlangt man Entropie proportional zur Fläche, steht 8π in der Einstein Gleichung.",
    b: ["Setzt Raumzeit und lokale Horizonte voraus. Verwandt: Replika Rechnung nach Fursaev und Solodukhin."],
    n: ["8", "1/4"], g: ["Horizont", "Entropie", "Gravitation"],
    l: [["p1", "testet"], ["grav", "karte"], ["ext-jacobson", "karte"]], src: ["docs/OPEN_PROBLEMS.md"] });

  O({ id: "org-c3-modular", lane: "p1", t: "c₃ als Energie am selbstdualen Torus", st: "belegt",
    s: "c₃ = E₂(i)/24: Ladungsenergie je Cartanrichtung jedes selbstdualen Gitters bei τ = i.",
    b: ["Exakt als Energie. Als Erzeuger der Ortsfaktoren nicht gestützt: Am Torus tritt die lemniskatische Konstante auf, in den TFPT Ausgaben nur π."],
    n: ["1/(8π)", "24"], g: ["Modul", "Torus", "Normierung"],
    l: [["p1", "testet"], ["taui", "karte"], ["pat-modular", "karte"]], src: [RAD] });

  O({ id: "org-c3-sphere", lane: "p1", t: "c₃ als Randantwort der Kugel", st: "bedingt",
    s: "Der Dirichlet zu Neumann Kern gibt am Gegenpunkt R·C = 1/(8πκ).",
    b: ["Für κ = 1 ist das c₃, samt positiver radialer Zeitform. κ wird nicht ausgewählt. Verwandt: Grenzwert der Decktransformationsspur am Torus und das Lemma zum Inversionsgewicht 1/2."],
    n: ["1/(8π)", "1/2"], g: ["Holografie", "Kugel", "Rand"],
    l: [["p1", "testet"], ["dyn-holo", "karte"], ["gap-cnat", "offen"]], src: [F1S] });

  O({ id: "org-winding", lane: "p1", t: "Einheitswindung als Erzeuger markierter Geschichten", st: "bedingt",
    s: "Die Windung um die Naht als U(1) Charakter der geladenen Endlinien.",
    b: ["Topologisch bestanden. Aus dem Mittelstrom allein widerlegt: Der Strom in der Mitte trägt die Windung nicht. Sie muss an den Endlinien sitzen.",
      "Der entscheidende Test (Wurzelkontraktionen a und b) steht aus."],
    n: ["1"], g: ["Windung", "Endlinien", "Ladung"],
    l: [["p1", "testet"], ["dyn-openword", "karte"], ["gap-ab", "offen"]], src: [FP + "TFPT-Ursprungshypothesen-2026-11-12.md", "research/originlab/OPH-WINDING-MEMORY-20261119.md"] });

  O({ id: "org-current", lane: "p1", t: "Geometrie aus der Stromantwort zurückgewinnen", st: "bedingt",
    s: "Aus den gemessenen Stromkorrelationen die Fläche und ihre Marken rekonstruieren.",
    b: ["Gelingt für reine Zustände über die Pfaffsche Struktur. Eine Mischung verletzt sie. Damit ist die Rekonstruktion ein Test auf Reinheit, kein Ursprung."],
    n: [], g: ["Strom", "Rekonstruktion"], l: [["p1", "testet"], ["pat-loss", "karte"]], src: [FP + "TFPT-Gesamtkonsolidierung-2026-11-12.md"] });

  O({ id: "org-fourmark", lane: "both", t: "Ein Objekt für beide Eingaben", st: "belegt",
    s: "Die Kugel mit vier Marken liefert Rang 3 und Dimension 5 zugleich.",
    b: ["Das ist der bisher stärkste gemeinsame Ursprung: P1 (Naht mit Deckstruktur) und P2 (fünf Plätze) sind zwei Lesungen eines Objekts. Offen bleibt die Wahl des Divisors."],
    n: ["4", "3", "5"], g: ["Marken", "Riemann–Roch"],
    l: [["p1", "testet"], ["p2", "testet"], ["viermarken", "karte"]], src: ["experiments/theory-contracts/puncture-metric-joint-completion-20261102/ORIGIN-TOPOLOGY-RESULT.txt"] });

  O({ id: "org-anchor", lane: "both", t: "Ankerkompression (1, 1, 2)", st: "belegt",
    s: "Drei Zahlen kodieren beide Eingaben und die E₈ Zählung.",
    b: ["Kompression, kein Ursprung: Man erhält P1 und P2 zurück, wenn man den Anker schon hat."],
    n: ["1", "2", "4", "5"], g: ["Anker", "Kompression"], l: [["anker", "karte"], ["p1", "testet"], ["p2", "testet"]], src: ["docs/THEORY.md"] });

  O({ id: "org-fable", lane: "both", t: "Acht freie Fermionen", st: "bedingt",
    s: "Vorschlag: Acht freie Fermionen tragen zwei Operatoren D und T mit DT = TD, D² = T⁴.",
    b: ["Die Phasenräume zählen 60, 64, 60, 64, passend zu den Wurzelklassen. T hat Ordnung 8 und ist ein reguläres Element von G₃₁. Der Cartanraum trägt beide Eingaben.",
      "Bedingt: Die Wahl von D und T ist eine Setzung, bis eine Regel sie liefert."],
    n: ["8", "60", "64", "4"], g: ["Fermion", "Operator", "Cartan"],
    l: [["p1", "testet"], ["p2", "testet"], ["ext-oph", "gleich"]], src: [FP + "TFPT-Gesamtkonsolidierung-2026-11-12.md"] });

  O({ id: "org-split", lane: "both", t: "Algebraische 3+1 Aufteilung", st: "bedingt",
    s: "Vier kürzeste Logarithmen mit Normquadrat 3/4 teilen sich in drei plus eins.",
    b: ["Eine algebraische Aufteilung. Dass sie Raum und Zeit bedeutet, ist nicht gezeigt."],
    n: ["3", "1", "3/4", "4"], g: ["Raumzeit", "Aufteilung"], l: [["netz", "testet"]], src: [FP + "TFPT-Gesamtkonsolidierung-2026-11-12.md"] });

  O({ id: "org-p2-pascal", lane: "p2", t: "P2 aus Pascal und Rangfüllung", st: "belegt",
    s: "2^g = g² + g + 2 hat nur g = 5 als passende Lösung. In Lean geprüft.",
    b: ["Zusammen mit g + 3 = 8 und der Rückschleife über h = 30. Eindeutig in der Vergleichsklasse."],
    n: ["5", "8", "30"], g: ["Bootstrap", "Lean"], l: [["p2", "testet"], ["bootstrap", "karte"]], src: ["tfpt_3_e8_audit_bootstrap.tex"] });

  O({ id: "org-p2-r102", lane: "p2", t: "Primitive ganze Ladung erzwingt E₈", st: "bedingt",
    s: "R102: Verlangt man eine primitive ganzzahlige Ladung, bleibt nur das unimodulare Gitter.",
    b: ["Damit folgt die Verklebung aus einer Ladungsforderung. Die Forderung selbst ist die neue Annahme."],
    n: ["8", "1"], g: ["Ladung", "Gitter"], l: [["p2", "testet"], ["e8", "karte"]], src: [FP + "RESULTS.md"] });

  O({ id: "org-p2-det", lane: "p2", t: "Kein freier abelscher Sektor", st: "bedingt",
    s: "„Nichts bleibt ungebunden“ ist gleichwertig zu det K = 1.",
    b: ["Physikalisch lesbar: Es gibt keine freie Ladung, die nichts sieht. Als Prinzip schön, als Auswahlgrund noch nicht hergeleitet. Verwandt: Smith Normalform und McKay Wachstum."],
    n: ["1"], g: ["Ladung", "Endlichkeit"], l: [["p2", "testet"], ["pat-growth", "karte"]], src: ["hypergraph.md"] });

  O({ id: "org-klein", lane: "p2", t: "Kleinsche Quelle F₂²", st: "bedingt",
    s: "Fünf isotrope Punkte einer kleinen binären Geometrie als die fünf Plätze.",
    b: ["Braucht eine zusätzliche Fahne (S₅ → S₄), um 3+2 zu trennen. Diese Fahne ist die versteckte Wahl."],
    n: ["5", "4", "2"], g: ["Klein", "Binär"], l: [["p2", "testet"], ["pat-code", "gleich"]], src: [F1S] });

  O({ id: "org-fivebeat", lane: "p2", t: "Vier Raumrichtungen aus zwei Fünfertakten", st: "aus",
    s: "Zwei Fünfertakte sollten vier euklidische Richtungen aufspannen.",
    b: ["Ausgeschlossen: Sie erzeugen A₆, das in SO(3) × SO(3) nur trivial wirkt."],
    n: ["5", "4", "360"], g: ["Clock", "Raumzeit"], l: [["dyn-clockgear", "karte"]], src: [FP + "TFPT-Ursprungshypothesen-2026-11-12.md"] });

  O({ id: "org-tau", lane: "sel", t: "Auswahl des Quadrats τ = i", st: "bedingt",
    s: "Nur am Quadrat bleibt die Markenmenge unter der Vierteldrehung erhalten.",
    b: ["Die freie Reflexionspositivität unterscheidet t nicht. Die Vierteldrehung tut es. Offen: Woher kommt die Forderung nach der Vierteldrehung, wenn nicht aus μ₄ selbst?"],
    n: ["4", "i"], g: ["Auswahl", "Modul"], l: [["taui", "karte"], ["naht", "testet"]], src: [RAD] });

  O({ id: "org-diagnose", lane: "sel", t: "Diagnose vom 16. November", st: "belegt",
    s: "P1, P2 und E₈ zusammen legen die vier Marken nicht fest.",
    b: ["Der Produzent der Marken liegt außerhalb der erklärten Eingaben. Das ist die präziseste Angabe, wo der Ursprung fehlt: Gesucht ist ein Vorgang, kein weiteres Objekt."],
    n: ["4"], g: ["Diagnose", "Ursprung"], l: [["quelle", "karte"], ["pat-closure", "karte"]], src: [F1S] });

  O({ id: "org-sixroutes", lane: "sel", t: "OriginLab: sechs Routen", st: "bedingt",
    s: "Matrixnaht, Petz Rückgewinnung, Momentenabschluss, regionale Algebren, chiraler Index, Streuung.",
    b: ["Alle 57 nichttrivialen Kombinationen werden für 32 Kontrollvarianten geprüft. Eine Kombination notwendiger Prüfungen ist kein physisches Modell.",
      "Offen laut Labor: das gewichtete originale P1 Paket, Kontinuum, globale Anomalien, physische Spiegellücke, vollständige Streuung."],
    n: ["6", "57", "32"], g: ["OriginLab", "Suche"], l: [["quelle", "testet"], ["gap-r113", "karte"]], src: ["research/originlab/SIX-ROUTES.md"] });

  /* ───────────── d) Universalraum ───────────── */
  const U = (o) => N(Object.assign({ k: "ur" }, o));

  U({ id: "ur-two", t: "Zwei Richtungen", st: "offen",
    s: "TFPT → Universalraum: operative Vervollständigung. Universalraum → TFPT: Ursprung.",
    b: ["Richtung eins fragt: Welcher größte Raum führt alles aus, was TFPT beschreibt? Richtung zwei fragt: Aus welchem Raum fällt TFPT als Auslese heraus?",
      "Beide Richtungen dürfen einander nicht stillschweigend voraussetzen."],
    n: ["2"], g: ["Universalraum", "Ursprung"], l: [["quelle", "karte"], ["dyn-recursion", "karte"]], src: ["universal_room/README.md"] });

  U({ id: "ur-future", t: "Zukunftsäquivalenz", st: "belegt",
    s: "Zwei Zustände sind gleich, wenn keine künftige Frage sie unterscheidet.",
    b: ["Das ist die tragfähige Definition des Raums: der Quotient nach gleichen Zukunftsantworten. Er rekonstruiert, erzeugt aber nicht."],
    n: [], g: ["Universalraum", "Rekonstruktion", "Gedächtnis"], l: [["dyn-wordsource", "karte"], ["ur-shadow", "karte"]], src: ["universal_room/README.md"] });

  U({ id: "ur-shadow", t: "Gleicher Schatten, andere Zukunft", st: "belegt",
    s: "Zwei Räume mit derselben heutigen Auslese können sich morgen unterscheiden.",
    b: ["Deshalb reicht die Übereinstimmung von Zahlen, Spektren oder Gruppenordnungen nie als Beweis für denselben Ursprung."],
    n: [], g: ["Universalraum", "Schatten", "Warnung"], l: [["pat-loss", "karte"], ["dyn-lorentz", "karte"]], src: ["universal_room/README.md"] });

  U({ id: "ur-first", t: "Runden v1.4 bis v1.6.3 (September 2026)", st: "aus",
    s: "Kein Welttest wurde nativ bestanden.",
    b: ["Der vereinheitlichte Feldtyp lieferte nur ein antisymmetrisches Feld B_μν. Ein Spin 2 Feld (Gravitation) entstand nicht.",
      "Die Runden bleiben als sauber dokumentierte Grenze wertvoll."],
    n: ["2"], g: ["Universalraum", "Welttest", "Spin 2"], l: [["grav", "grenze"], ["dyn-unitary", "karte"]], src: ["universal_room/README.md"] });

  /* ───────────── e) f) Origin, OPH, andere Theorien ───────────── */
  const X = (o) => N(Object.assign({ k: "ext" }, o));

  X({ id: "ext-origin", sub: "Origin Idee", t: "Origin Theorie", st: "deutung",
    s: "Naht = Horizont. Ein α für Elektromagnetismus, kosmischen Horizont und Λ. Die Welt als sich selbst wiedergebende Lesung.",
    b: ["Die Origin Theorie ist die physische Erzählung zum TFPT Compiler. Ihr ganzzahliges Skelett folgt aus (g, N_fam) = (5, 3)."],
    same: ["Acht dreifach erzwungen", "eindeutiger Attraktor unter den gelückten Voraussetzungen", "Kinematik maschinengeprüft"],
    diff: ["Die zyklische, sich selbst reproduzierende Lesung ist ausdrücklich Deutung", "Familienfixpunkt ist kein exakter Quantenzyklus"],
    open: ["Übergangsdynamik am Kollaps", "warum die Naht ein physischer Horizont ist"],
    n: ["5", "3", "8", "137"], g: ["Origin", "Horizont", "Zyklus"],
    l: [["dyn-cyclic", "karte"], ["attraktor", "karte"], ["alpha", "karte"], ["org-sixroutes", "karte"]], src: ["origin_theory.tex", "documents/origin_theory/"] });

  X({ id: "ext-oph", sub: "OPH", t: "Observer Patch Holography", st: "bedingt",
    s: "Ein unabhängiges Programm, das von Beobachterflächen auf dem Ikosaeder ausgeht und an denselben Gittern ankommt.",
    b: ["Stand 19. November 2026. Alle Dynamikergebnisse gelten für erklärte endliche Modelle."],
    same: ["Affine McKay Konstruktion des OPH Spindubletts gibt E₈", "Löschen des vierdimensionalen Knotens gibt D₅ ⊕ A₃ mit Normen 5/4 und 3/4 und Z₄ Verklebung", "Wurzelklassen 52, 64, 60, 64", "Außenalgebra auf 3+2 gleich dem TFPT Materieträger samt Hyperladung"],
    diff: ["TFPT führt den neutralen Λ⁰ Zustand mit, OPH nicht", "Antipodaler Quotient hebt sich nicht zu primitiver C₄ Holonomie", "Keine harmonische Quadratmarkierung unter 45 Kombinationen", "Spinhebungen antikommutieren, Korrektur (7,3,3,3) → (6,4,4,2) nötig"],
    open: ["Gemeinsame positive Entwicklung auf derselben ursprünglichen Zeit und Referenz", "ursprüngliche Vakuummischung mit geladener Rückkehr"],
    n: ["5/4", "3/4", "60", "64", "52", "1/64"], g: ["OPH", "Ikosaeder", "McKay", "Holografie"],
    l: [["e8", "karte"], ["dyn-oph", "karte"], ["dyn-markov", "karte"], ["pat-hyper", "karte"], ["gap-oph", "offen"]],
    src: ["research/originlab/OPH-COHERENT-HISTORY-20261119.md", "research/originlab/OPH-WINDING-MEMORY-20261119.md"] });

  const T = (id, t, s, same, diff, open, n, g, l) =>
    X({ id, sub: "Andere Theorien", t, st: "bedingt", s, b: [], same, diff, open, n, g, l, src: ["docs/OPEN_PROBLEMS.md"] });

  T("ext-connes", "Connes–Chamseddine · Spektraltripel", "Das Standardmodell als Geometrie eines endlichen nichtkommutativen Raums.",
    ["Endliche Algebra, KO Dimension 6, 96 Dimensionen", "ein Higgs aus dem Materieoperator"], ["TFPT leitet die Algebra aus dem Träger her, statt sie zu setzen", "Faktor −12 gegen −2 in der Spurlesung"], ["Gewicht der Spektralwirkung aus der Quelle"],
    ["96", "6"], ["Spektral", "Higgs"], [["rand", "karte"], ["dyn-spectral", "karte"]]);
  T("ext-jacobson", "Jacobson · Gravitation als Zustandsgleichung", "Einstein Gleichung aus Thermodynamik lokaler Horizonte.",
    ["8π als Kopplung"], ["setzt Raumzeit voraus"], ["Herkunft der 4D Feldtheorie"], ["8"], ["Gravitation", "Horizont"], [["grav", "karte"], ["org-c3-thermo", "karte"]]);
  T("ext-bianconi", "Bianconi · entropische Gravitation", "Gravitation aus relativer Entropie zweier Metriken.",
    ["β′ = c₃/6 exakt"], ["R² Koeffizient verfehlt um 3·(8π)⁹"], [], ["1/6"], ["Gravitation", "Entropie"], [["dyn-entropic", "karte"]]);
  T("ext-kitaev", "Kitaev · E₈ Zustand und 16 fache Ordnung", "Ein gelückter 2+1D Zustand mit chiralem Rand c = 8.",
    ["c = 8, 248 Ströme", "Kondensation D₅⊕A₃ (16) → D₈ (4) → E₈ (1)"], ["Halbladungsfeld der Naht nicht identifiziert"], ["Kontinuumsgrenzwert"], ["8", "16", "248"], ["Anomalie", "Rand", "Kitaev"], [["dyn-anomaly", "karte"], ["rand", "karte"]]);
  T("ext-wolfram", "Wolfram · Hypergraph Physik", "Raumzeit als Grenzfall einer Ersetzungsregel auf Hypergraphen.",
    ["lokale Regel mit Attraktor", "Ikosaeder als 3 uniformer Hypergraph"], ["hier nur Skelettebene, keine Raumzeit", "Seed mit drei Zweigen bleibt Eingabe"], ["gewichtete Hyperkanten aus effektiven Vertices"], ["3", "120"], ["Hypergraph", "Netzwerk"], [["dyn-hypergraph", "karte"], ["gap-hyperedge", "offen"]]);
  T("ext-gut", "SO(10), Pati–Salam, SU(5)", "Große Vereinheitlichung der Eichgruppen.",
    ["16 von Spin(10)", "sin²θ_W = 3/8", "Hyperladungen"], ["E₈ ist hier Prüfobjekt, keine Eichgruppe", "kein Protonzerfall als Folge behauptet"], ["Zweischleifenlauf mit Schwellen"], ["16", "10", "3/8"], ["Standardmodell", "Ladung"], [["sm", "karte"], ["d5", "karte"], ["gap-twoloop", "offen"]]);
  T("ext-string", "Heterotische Gitter · Lisi E₈", "E₈ ⊕ E₈ und D₁₆⁺ als selbstduale Träger. Lisi: E₈ als Eichgruppe von allem.",
    ["gerade unimodulare Gitter", "II(9,1) = E₈ ⊕ U"], ["Torusdaten unterscheiden E₈⊕E₈ und D₁₆⁺ nicht", "gegen Lisi: E₈ ist hier nicht die Eichgruppe"], ["Geschlecht 5 Theta Differenz numerisch"], ["16", "8", "9"], ["Gitter", "Lorentz"], [["dyn-lorentz", "karte"], ["gap-theta", "offen"]]);
  T("ext-ccc", "Penrose CCC · Carroll Rekurrenz", "Zyklische und wiederkehrende Kosmologien.",
    ["Horizont als Übergang"], ["exakter Reset im selben Vakuum ausgeschlossen"], ["mehrere orthogonale Vakuumdarstellungen"], [], ["Zyklus", "Kosmos"], [["dyn-cyclic", "karte"], ["dyn-reset", "karte"]]);
  T("ext-koide", "Koide Formel", "Q = 2/3 für die geladenen Leptonen.",
    ["2/3 erscheint auch als Rückgewinnungsrate"], ["bewusst Schnittstelle, keine Compilerausgabe"], ["Transfer Quelle → Polmasse"], ["2/3"], ["Koide", "Flavor"], [["pat-23", "karte"], ["gap-lepton", "karte"]]);
  T("ext-prime", "Prime Front · Riemann", "Eigenes Projekt zu Primzahlen und Zetafunktion.",
    [], ["kein Beleg für TFPT, getrennt geführt"], [], [], ["Primzahlen"], [["pat-prime", "karte"]]);
  T("ext-starobinsky", "Starobinsky Inflation", "R + R² Inflation.",
    ["n_s = 1 − 2/N, r = 12/N²"], ["Skalaronmasse hier aus c₃^(7/2) statt angepasst"], ["Reheating, N"], ["7/2", "12"], ["Inflation", "Kosmologie"], [["kosmo", "karte"], ["grav", "karte"]]);

  /* ───────────── Noch nicht ausgeführt · nächste Entscheide ───────────── */
  const G = (o) => N(Object.assign({ k: "gap", st: "neu", b: [], n: [], src: [FP + "PLAN.md"] }, o));

  G({ id: "gap-f1", st: "offen", t: "F1: gemeinsame Quelle auswerten", s: "Die aktive Hauptaufgabe. Alles andere hängt daran.", g: ["F1", "Quelle"], l: [["quelle", "karte"], ["dyn-fieldint", "karte"]] });
  G({ id: "gap-r113", st: "offen", t: "R113 Blockade lösen", s: "Es fehlen das originale geblockte Faserpaket, die gewichtete Randform und die absolute Leerreferenz.", g: ["F1", "Quelle", "Referenz"], l: [["dyn-fieldint", "karte"]], src: [FP + "RESULTS.md"] });
  G({ id: "gap-ab", t: "Wurzelkontraktionen a und b am markierten P1 Streifen", s: "Entscheidet, ob die Windung an den Endlinien sitzt. Benannt, nicht ausgeführt.", g: ["Windung", "Endlinien", "Test offen"], l: [["org-winding", "testet"], ["dyn-openword", "testet"]], src: [FP + "TFPT-Ursprungshypothesen-2026-11-12.md"] });
  G({ id: "gap-cnat", t: "Native Operatorgleichheit C_nat", s: "Ist der Kugelrandkern derselbe Operator wie der originale P1 Kern? Priorität A.", b: ["Gerechnet am 19.11.: bedingt bestanden. Auf der runden S² trägt die positive Hälfte des originalen Dirac Randkerns je Kanal genau 1/(2π) = 4c₃ und gibt exp(−2α). Das Spektrum ist nur auf S² flach. Gezeigt sind Spektrum und Gewicht, keine volle Operatorgleichheit. Offen: runde S² und positive Hälfte herleiten."], g: ["Holografie", "Operator", "Test offen"], l: [["dyn-holo", "testet"], ["dyn-p1kern", "karte"]], src: [F1S, "experiments/theory-contracts/p1-kernel-ground-state-20261119/README.md"] });
  G({ id: "gap-torus", t: "Brücke Kugel ↔ Torus über dieselbe Überlagerung", s: "Verbindet c₃ als Kugelzahl mit c₃ als Toruszahl.", b: ["Gerechnet am 19.11.: Die Brücke ist exakt. E₂(q) = 1 − 12·Σ Z⁺(mβ), wobei Z⁺ die Wärmespur des originalen Randoperators auf S² ist und m die Umläufe zählt. Am quadratischen Torus tragen die Schwingungen −c₃ und die Nullmode +c₃.", "Folge: Kugelzahl, Toruszahl und KMS Zahl sind dieselbe Zahl T/4. Das ist ein Beleg, nicht drei. Das Quadrat ist ein Sattel und nur mit Spiegelung stabil."], g: ["Torus", "Kugel", "Modul", "Test offen"], l: [["pat-modular", "testet"], ["viermarken", "testet"], ["dyn-p1kern", "karte"]], src: [RAD, "experiments/theory-contracts/p1-kernel-ground-state-20261119/README.md"] });
  G({ id: "gap-operator", t: "Warum Zeit am ganzen Klebevektor, Grade nur an fünf Richtungen?", s: "Erklärte den Exponenten 5/4 und warum B−L mitzählt.", b: ["Gerechnet am 19.11.: (1−q)^(−5/4) ist die Zweipunktfunktion des D₅ Spinfelds (Gewicht 5/8, fünf Richtungen mal 1/4). D₄ verfehlt um −47 σ, D₆ um +51 σ. Das benennt denselben Exponenten wie q(D₅), ist also kein neuer Beleg, aber ein konkreter Operator. Ob dahinter ein volles Modenprodukt steht, kann α nicht entscheiden (0,03 σ)."], g: ["Alpha", "Verklebung", "Operator", "Test offen"], l: [["pat-54", "testet"], ["alpha", "testet"]], src: [RAD] });
  G({ id: "gap-alpha", t: "Vorausschauender α Test", s: "Eine Messung mit Unsicherheit ≤ 5·10⁻⁹, die um mehr als 5σ von 137,0359992168 abweicht, widerlegt die Gleichung.", g: ["Alpha", "Messung", "Test offen"], n: ["137"], l: [["alpha", "testet"]], src: [RAD] });
  G({ id: "gap-lepton", t: "Transfer Quellmasse → Polmasse der Leptonen", s: "+2,20 % und +0,21 % Abstand. Die schärfste empirische Stelle.", b: ["Gerechnet am 19.11.: Einschleifen QED Lauf auf eine gemeinsame Skala hat das falsche Vorzeichen (−0,98 % und −1,86 %). Nimmt man Koide als Schnittstelle, gibt das Quellverhältnis e/μ allein m_μ/m_τ auf 0,02 %. Dann bleibt nur die eine Lücke von 0,21 %; m_τ läge 4 σ zu hoch."], g: ["Flavor", "Koide", "Test offen"], l: [["flavor", "testet"]], src: [RAD, "experiments/theory-contracts/landkarte-origin-scenarios-20261119/README.md"] });
  G({ id: "gap-twoloop", t: "Zweischleifenlauf mit Schwellen (G4)", s: "Trägerskala gegen Skalaron: Verhältnis 2,97 bei einer Schleife. Zwei Schleifen fehlen.", g: ["Standardmodell", "Lauf", "Test offen"], n: ["2,97"], l: [["sm", "testet"], ["grav", "testet"]], src: [RAD] });
  G({ id: "gap-g2", t: "G2: Familie = Kanal im Leptonsektor", s: "Hypothese benannt, ohne eigene Probe.", g: ["Familie", "Flavor", "Test offen"], l: [["flavor", "testet"], ["a3", "testet"]], src: [RAD] });
  G({ id: "gap-g3", t: "G3: Horizont zuerst", s: "Thermofeldverdoppelter Zustand als Start. Am 19.11. geprüft.", b: ["Als Quantenzustand ausgeschlossen: Bei β = 2π trägt jede Quadratur mindestens 2c₃, nie c₃, und die Dämpfung wäre exp(−6,31α). Als klassische Nullmode bei T = 1/(2π) gibt der Horizont genau exp(−2α). Das ist dieselbe Zahl wie beim Grundzustand, kein zweiter Beleg."], g: ["Horizont", "Zustand", "Test offen"], l: [["dyn-selfdual", "testet"], ["grav", "testet"]], src: [RAD] });
  G({ id: "gap-ufield", t: "α als Feld mit Potential U(α)", s: "Am 19.11. geprüft: ein einziges Minimum.", b: ["Das Potential U mit U′ = F hat auf (0,1) genau ein Minimum, keine Barriere und kein zweites Vakuum. Krümmung 0,9885·3α². Ein kinetischer Term ist nicht bekannt, deshalb folgt daraus noch keine Dynamik."], g: ["Alpha", "Potential", "Test offen"], l: [["dyn-alphafield", "testet"]], src: [RAD] });
  G({ id: "gap-theta", t: "Geschlecht 5 Theta Differenz numerisch", s: "Unterscheidet E₈ ⊕ E₈ von D₁₆⁺. Erst dort sichtbar.", g: ["Gitter", "Modul", "Test offen"], n: ["5", "16"], l: [["ext-string", "testet"]], src: [F1S] });
  G({ id: "gap-picard", t: "Hilbertraum und Ladungsintertwiner für den Picard Wechsel", s: "Macht aus „P1 als wechselnder Leser“ ein prüfbares Modell. Priorität B.", g: ["Picard", "Auslese", "Test offen"], l: [["dyn-reader", "testet"]], src: [F1S] });
  G({ id: "gap-hyperedge", t: "Gewichtete Hyperkanten aus effektiven Vertices", s: "Die echten Vertices der Quelle rekursiv verketten statt ein Skelett wachsen zu lassen.", g: ["Hypergraph", "Rekursion", "Test offen"], l: [["dyn-hypergraph", "testet"], ["dyn-recursion", "testet"]], src: ["experiments/theory-contracts/recursive-mirrored-source-20261103/ANALYSIS.txt"] });
  G({ id: "gap-oph", t: "Ursprüngliche Vakuummischung mit geladener Rückkehr", s: "Der Anschluss von OPH an die TFPT Quelle auf derselben Zeit.", g: ["OPH", "Kohärenz", "Test offen"], l: [["dyn-oph", "testet"], ["ext-oph", "testet"]], src: ["research/originlab/OPH-COHERENT-HISTORY-20261119.md"] });

  /* ───────────── Beschriftungen ───────────── */
  window.TFPT.status = {
    gesetzt: ["Gesetzt", "Eingabe, die vorausgesetzt wird"],
    belegt: ["Belegt", "Im genannten Geltungsbereich gezeigt"],
    bedingt: ["Bedingt", "Gilt unter benannten Annahmen"],
    offen: ["Offen", "Benannte Lücke, in Arbeit"],
    aus: ["Ausgeschlossen", "Im geprüften Bereich gescheitert"],
    neu: ["Nicht probiert", "Benannt, aber noch nicht ausgeführt"],
    deutung: ["Deutung", "Erzählung, keine Folge"],
  };
  window.TFPT.kinds = { chain: "Herleitung", dyn: "Dynamik", pat: "Muster", org: "Ursprung", ur: "Universalraum", ext: "Verbindung", gap: "Offener Test" };
  window.TFPT.linkTypes = {
    leitet: "führt zu", bedingt: "führt bedingt zu", seed: "speist", offen: "fehlt für", rueck: "bestimmt rückwärts",
    gleich: "gleiche Zahl oder Form, keine Abbildung gezeigt", karte: "ausdrücklich verbunden", dockt: "dockt an",
    testet: "prüft", grenze: "begrenzt",
  };
})();
