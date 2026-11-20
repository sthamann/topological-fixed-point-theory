/* TFPT Landkarte · kuratierte Inhalte, Teil 1: Herleitung, Ausgaben, Dynamik.
   Stand der Lektüre: 19. November 2026. Dies ist eine Lesehilfe, kein Beweis.
   Maßgeblich bleiben Ledger, Originalquellen und der lebende Atlas. */
window.TFPT = window.TFPT || { nodes: [] };
(function () {
  const N = (o) => window.TFPT.nodes.push(o);
  const FP = "docs/research/fundamental-plan/";

  /* ───────────── a) Herleitungskette ───────────── */
  N({ id: "p1", k: "chain", col: 0, row: 1.2, t: "P1 · Randkern der Naht", st: "gesetzt",
    s: "Eine orientierte Naht mit reflexionspositivem Randkern, Einheitswindung und der Normierung c₃ = 1/(8π).",
    b: ["P1 ist mehr als eine Zahl. Die operative Saat im Original enthält ein lokales Netz, Zeit, Reflexion, Zustand, eine Nahtklasse und den elliptischen Kragen. c₃ ist die Normierung dieses Kerns.",
      "Bild: ein Rand, der zwei Seiten zusammennäht. Wie stark beide Seiten aneinander koppeln, sagt c₃.",
      "Status: gesetzte Eingabe. In mehreren Verträgen folgt c₃ unter Voraussetzungen (Geometrie, Thermodynamik, Modulformen), die Herkunft der ganzen Saat ist offen."],
    n: ["1/(8π)", "8", "2"], g: ["Naht", "Rand", "Normierung", "Reflexion"],
    l: [["naht", "leitet"], ["phi0", "seed"], ["grav", "bedingt"]],
    src: ["tfpt_1_architecture_e8.tex", "documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/"] });

  N({ id: "p2", k: "chain", col: 0, row: 4.2, t: "P2 · Träger mit fünf Plätzen", st: "gesetzt",
    s: "Eine Schnittstelle mit fünf Trägerplätzen, aufgeteilt in 3 + 2.",
    b: ["Drei Plätze werden zu Farbe, zwei zur schwachen Wechselwirkung. Aus dieser einen Aufteilung folgen später Hyperladung, Eichgruppe und die 16 Zustände einer Familie.",
      "Die Zahl 5 ist innerhalb ihrer Auswahlklasse mehrfach erzwungen (Pascal, Rangfüllung, Coxeter). Das ist Konsistenz. Warum die Natur diese Schnittstelle trägt, ist damit nicht hergeleitet."],
    n: ["5", "3", "2"], g: ["Träger", "3+2"],
    l: [["d5", "bedingt"]], src: ["tfpt_1_architecture_e8.tex", "docs/THEORY.md"] });

  N({ id: "naht", k: "chain", col: 1, row: 0, t: "Orientierte Naht · Doppeldeckung", st: "gesetzt",
    s: "Die Naht hat zwei Ufer (Z₂). Ihre Doppelüberlagerung und die μ₄ Deckstruktur ordnen alles Weitere.",
    b: ["Möbius und Doppeldeckung strukturieren die Naht. Vorsicht bei den Wörtern: Orientierungsüberlagerung, Spin(n)→SO(n), zwei Ufer und Zustandsverdopplung sind vier verschiedene Operationen. Ein Faktor zwei muss aus der tatsächlich benutzten Maßform kommen."],
    n: ["2", "4"], g: ["Naht", "Doppeldeckung", "Möbius"],
    l: [["viermarken", "bedingt"]], src: ["research/originlab/F1-ORIGIN-SCENARIOS-20261116.txt"] });

  N({ id: "taui", k: "chain", col: 1, row: 1.4, t: "Harmonischer Punkt τ = i", st: "bedingt",
    s: "Die vier Marken bilden ein Quadrat {1, i, −1, −i}. Nur dort bleibt der Divisor unter der Vierteldrehung erhalten.",
    b: ["In der Familie {1, −1, t, −t} bleiben 3+2, Hyperladung und der endliche Intertwiner für jedes t erhalten. Die Vierteldrehung z → iz erhält den Divisor genau für t = ±i.",
      "Neu am 18.11.: τ = i ist der einzige Modul, an dem die absolute Konstante der Stromzweipunktfunktion verschwindet (E₆(i) = 0). Das ist ein Kandidat für das fehlende Auswahlprinzip, als Produzent der Ortsfaktoren aber nicht gestützt."],
    n: ["4", "i"], g: ["Auswahl", "Modul", "Quadrat"],
    l: [["viermarken", "bedingt"], ["dyn-selfdual", "dockt"], ["org-tau", "karte"]],
    src: ["experiments/theory-contracts/radical-global-hypotheses-20261118/README.md", FP + "TFPT-Gesamtkonsolidierung-2026-11-12.md"] });

  N({ id: "viermarken", k: "chain", col: 1, row: 2.9, t: "Gemeinsames Viermarkenobjekt", st: "belegt",
    s: "Dieselbe Kugel mit vier Marken liefert zwei Zahlen: Rang 3 aus ihren Schleifen und Dimension 5 aus ihren Funktionen.",
    b: ["Schleifen um vier Löcher: H₁ ≅ Z³, also drei unabhängige Umläufe → Familienrang 3.",
      "Funktionen mit einfachen Polen an den vier Marken: Riemann–Roch gibt h⁰ = 5 → Trägerzahl 5. Die Halbdrehung wirkt darauf mit drei geraden und zwei ungeraden Richtungen → 3+2.",
      "Das ist bewiesen und wird nicht mehr gesucht. Offen bleibt, warum genau dieser Divisor gewählt wird."],
    n: ["4", "3", "5"], g: ["Marken", "Riemann–Roch", "Familie", "Träger"],
    l: [["a3", "leitet"], ["d5", "leitet"], ["org-fourmark", "karte"]],
    src: ["experiments/theory-contracts/puncture-metric-joint-completion-20261102/ORIGIN-TOPOLOGY-RESULT.txt"] });

  N({ id: "anker", k: "chain", col: 1, row: 5, t: "Anker a = (1, 1, 2)", st: "belegt",
    s: "Drei kleine Zahlen fassen beide Eingaben zusammen: e₁ = 4, e₂ = 5, e₃ = 2.",
    b: ["c₃ = 1/(2·e₁·π) und g = e₂. Die Potenzsummen pₙ = 2 + 2ⁿ geben 240 = p₁p₂p₃, 8 = p₄ − p₃ und 248 als Summe.",
      "Der Anker komprimiert die Architektur. Er erzeugt sie nicht aus dem Nichts."],
    n: ["1", "2", "4", "5", "240", "248", "8"], g: ["Anker", "Kompression"],
    l: [["p1", "rueck"], ["p2", "rueck"], ["famop", "leitet"], ["e8", "gleich"]],
    src: ["docs/THEORY.md"] });

  N({ id: "phi0", k: "chain", col: 2, row: 0.2, t: "Gemeinsamer Seed φ₀", st: "bedingt",
    s: "φ₀ = 1/(6π) + 48·c₃⁴ ≈ 0,0531720. Eine einzige kleine Zahl speist Flavor, α und Kosmologie.",
    b: ["Baumterm 1/(6π) = (4/3)·c₃, dazu der topologische Zusatz 48·c₃⁴ = 3/(256π⁴).",
      "Weil so viele Ausgaben denselben Seed lesen, sind ihre Treffer korreliert und keine unabhängigen Bestätigungen."],
    n: ["1/(6π)", "48", "4/3"], g: ["Seed", "Korrelation"],
    l: [["flavor", "seed"], ["alpha", "seed"], ["kosmo", "seed"], ["pat-seed", "gleich"], ["pat-48", "gleich"]],
    src: ["tfpt_2_standard_model.tex"] });

  N({ id: "a3", k: "chain", col: 2, row: 2.2, t: "Familienseite A₃", st: "bedingt",
    s: "Die drei Umläufe werden zum Gitter A₃ (su(4), Rang 3). Daher drei Familien.",
    b: ["Die ganzzahligen Meridiane führen über die Summenrelation und die deklarierte Zählmetrik zum A₃ Wurzel und Gewichtsgitter. Ein Meridian ist exakt das Familiengewicht des E₈ Verklebungsgenerators."],
    n: ["3", "4", "3/4"], g: ["Familie", "Gitter"],
    l: [["e8", "leitet"], ["famop", "leitet"]], src: [FP + "TFPT-Gesamtkonsolidierung-2026-11-12.md"] });

  N({ id: "d5", k: "chain", col: 2, row: 4.2, t: "Träger D₅ · Halbspinor 16", st: "belegt",
    s: "Aus fünf Plätzen entsteht S⁺ = Λ^gerade C⁵ mit 1 + 10 + 5 = 16 Zuständen: genau eine Familie.",
    b: ["Das ist der Halbspinor von Spin(10). Die Lesart ist an die gewählte Gruppenwirkung gebunden."],
    n: ["5", "16", "10", "5/4"], g: ["Träger", "Spinor"],
    l: [["e8", "leitet"], ["sm", "leitet"]], src: ["tfpt_1_architecture_e8.tex"] });

  N({ id: "clocks", k: "chain", col: 3, row: 0.9, t: "Ineinandergreifende Clocks", st: "belegt",
    s: "Fünfer, Sechser und Dreißigertakt, dazu μ₄ als Getriebe: G·C₅·G⁻¹ = C₅², G⁴ = 1.",
    b: ["h(E₈) = 30 = 2·3·5 und Z₃₀ ≅ Z₅ × Z₆. Die Vier ist kein Unterzyklus der 30 (4 teilt 30 nicht). Sie greift als Automorphismus in den Fünfertakt ein.",
      "Compilerclock c = Jσ mit c¹² = 1, σ = c⁴, J = c⁹ in G₃₁ (Ordnung 46080).",
      "Wichtig: Matrixordnung, Relaxationsrate, RG Skala und physische Dauer sind vier verschiedene Dinge."],
    n: ["30", "5", "6", "4", "12", "8"], g: ["Clock", "Galois", "Coxeter"],
    l: [["attraktor", "gleich"], ["dyn-clockgear", "dockt"], ["pat-galois", "karte"]],
    src: ["origin_theory.tex", "experiments/theory-contracts/clock-marking-audit/README.md"] });

  N({ id: "e8", k: "chain", col: 3, row: 3, t: "Verklebung D₅ ⊕ A₃ → E₈", st: "belegt",
    s: "Beide Gitter haben dieselbe Lücke Z₄. Klebt man sie dort zusammen (5/4 + 3/4 = 2), entsteht das Gitter E₈.",
    b: ["det = (4·4)/4² = 1: gerade, unimodular, Rang 8. 240 Wurzeln, Dimension 248.",
      "E₈ ist hier Abschluss und Prüfobjekt. Es ist nicht die physische Eichgruppe.",
      "Zerlegung: (45,1) ⊕ (1,15) ⊕ (16,4) ⊕ (16̄,4̄) ⊕ (10,6)."],
    n: ["8", "240", "248", "4", "5/4", "3/4", "2"], g: ["E8", "Gitter", "Verklebung"],
    l: [["clocks", "leitet"], ["qcode", "karte"], ["bootstrap", "leitet"], ["pat-growth", "gleich"], ["ext-oph", "karte"]],
    src: ["tfpt_3_e8_audit_bootstrap.tex", "verification/status_ledger.csv"] });

  N({ id: "qcode", k: "chain", col: 3, row: 4.4, t: "Endliche Quantenstruktur", st: "belegt",
    s: "Dasselbe Gitter ist ein Fehlerkorrekturcode und trägt die Zwei Qubit Geometrie.",
    b: ["E₈ über Z[i] entsteht per Construction A aus dem Hamming Code [8,4,4]. Der Quotient L/(1+i)L ist F₂⁴.",
      "60 gaußsche Linien entsprechen den Zwei Qubit Stabilisatorstrahlen, geordnet in 15 Viererkontexten."],
    n: ["8", "4", "16", "15", "60"], g: ["Code", "Qubit", "Clifford"],
    l: [["rand", "bedingt"], ["pat-code", "karte"], ["pat-qubit", "karte"]], src: ["note_e8_gaussian_code.tex"] });

  N({ id: "bootstrap", k: "chain", col: 3, row: 5.6, t: "Rückbootstrap", st: "bedingt",
    s: "Der Abschluss schaut zurück und lässt nur g = 5 und die Vier zu. Das ist Rückbestimmung als Konsistenz.",
    b: ["μ² − 5μ + 4 = 0 → μ = 4. Pascal: 2^g = g² + g + 2 → g = 5. Rangfüllung g + 3 = 8. Rückschleife g = 5 → E₈ → h = 30 = 2·3·5 → 5.",
      "„Eindeutig“ gilt immer innerhalb der Vergleichsklasse. Der Beweisgraph bleibt azyklisch: Eine Folge darf ihre eigene Annahme nicht begründen."],
    n: ["5", "4", "30", "8"], g: ["Bootstrap", "Rückwärts", "Fixpunkt"],
    l: [["p1", "rueck"], ["p2", "rueck"], ["pat-closure", "gleich"], ["org-p2-pascal", "karte"]],
    src: ["tfpt_3_e8_audit_bootstrap.tex"] });

  N({ id: "attraktor", k: "chain", col: 4, row: 0, t: "Transfer und Attraktor", st: "bedingt",
    s: "Ein Transferoperator mit Spektrum {1, (2/3)⁶, (1/3)⁶} läuft immer auf denselben Zustand zu.",
    b: ["Lücke Δ = 6·log(3/2). Unter den positiven, gelückten Voraussetzungen gibt es genau eine stationäre Richtung.",
      "Der Satz beschreibt das Verhalten eines gewählten T. Er wählt T nicht aus."],
    n: ["2/3", "1/3", "6"], g: ["Attraktor", "Transfer", "Lücke"],
    l: [["dyn-transfer", "dockt"], ["pat-23", "gleich"]], src: ["origin_theory.tex"] });

  N({ id: "alpha", k: "chain", col: 4, row: 1.2, t: "α als Fixpunkt", st: "bedingt",
    s: "Eine feste Gleichung hat genau eine positive Wurzel: α⁻¹ = 137,0359992168…",
    b: ["F(α) = α³ − 2c₃³α² − (4/5)·41·c₃⁶·log(1/φ(α)) = 0 mit φ = 1/(6π) + q(1−q)^(−5/4), q = 48c₃⁴·e^(−2α).",
      "In der Variable g = 32π²α lautet sie schlicht g²(g − c₃) = b₁·log(1/φ) mit b₁ = 41/10.",
      "Die Wurzel ist bestimmt und eindeutig. Offen sind der physische Ursprung des Funktionals und die Gleichsetzung mit dem gemessenen α. Der Exponent 5/4 ist kein angepasster Wert: Ein Fit gäbe 1,2404 ± 0,0051, und 5/4 = q(D₅) ist eine Gitterzahl. Unter den Exponenten n/4 trifft nur n = 5 (n = 4: −47 σ, n = 6: +51 σ)."],
    n: ["137", "41", "5/4", "48", "4/5"], g: ["Alpha", "Fixpunkt", "Ward"],
    l: [["kosmo", "bedingt"], ["pat-54", "karte"], ["pat-41", "karte"], ["gap-alpha", "offen"]],
    src: ["experiments/theory-contracts/radical-global-hypotheses-20261118/README.md", "verification/v3_*.py"] });

  N({ id: "flavor", k: "chain", col: 4, row: 2.4, t: "Flavor · Massen und Mischung", st: "bedingt",
    s: "Massen als Leiter in Potenzen von φ₀, Mischungswinkel aus demselben Seed.",
    b: ["λ_C = √(φ₀(1−φ₀)) ≈ 0,2244. Leptonen: m_μ/m_τ = (8/7)φ₀, m_e/m_μ = (12/7)φ₀². Quarks: 55/117, (34/47)/φ₀, (3/26)/φ₀².",
      "Das sind Quellmassen. Der Weg zu Polmassen (Renormierung, Schwellen) ist ein eigener, offener Transfer und derzeit die schärfste empirische Stelle."],
    n: ["8/7", "12/7", "55/117"], g: ["Flavor", "Yukawa", "CKM", "PMNS"],
    l: [["gap-lepton", "offen"]], src: ["tfpt_2_standard_model.tex"] });

  N({ id: "famop", k: "chain", col: 4, row: 3.6, t: "Familienoperatoren Q, K, R, L", st: "belegt",
    s: "Vier ganzzahlige Matrizen auf dem Familiengitter Z³ mit Determinanten 3, 4, 8, 20.",
    b: ["Produkt 1920 = |W(D₅)|. R hat die Zeilen (1,3,0), (1,5,2), (2,5,3). Derselbe K Block liefert 1ᵀK1 = 25 und aᵀKa = 41: Massenvolumen und elektromagnetisches Budget hängen an denselben Operatoren."],
    n: ["3", "4", "8", "20", "1920", "41", "25"], g: ["Flavor", "Operator"],
    l: [["flavor", "leitet"], ["alpha", "seed"], ["pat-1920", "gleich"]], src: ["tfpt_2_standard_model.tex"] });

  N({ id: "sm", k: "chain", col: 4, row: 4.8, t: "Standardmodell Gerüst", st: "belegt",
    s: "SU(3) × SU(2) × U(1), die richtigen Hyperladungen, 16 Zustände je Familie, drei Familien, ein Higgs.",
    b: ["Hyperladung −1/3 auf Farbplätzen, +1/2 auf schwachen Plätzen. k_Y = 5/3, sin²θ_W = 3/8 an der Trägerskala. Die Anomaliesummen passen.",
      "Das ist die Ladungs und Darstellungsstruktur. Eine konstruierte wechselwirkende chirale Quantenfeldtheorie in 3+1 Dimensionen ist es noch nicht."],
    n: ["3", "2", "16", "48", "96", "3/8", "5/3"], g: ["Standardmodell", "Ladung"],
    l: [["netz", "offen"], ["ext-gut", "gleich"]], src: ["tfpt_2_standard_model.tex"] });

  N({ id: "grav", k: "chain", col: 5, row: 0, t: "Gravitation als Randkanal", st: "bedingt",
    s: "Unter den Annahmen des Verschränkungsgleichgewichts folgt die Einstein Gleichung mit 1/c₃ = 8π.",
    b: ["Dazu ein R + R² Zweig mit Skalaron bei c₃^(7/2)·M̄ ≈ 3,06·10¹³ GeV.",
      "Die Route setzt die 4D Quantenfeldtheorie voraus, deren Verschränkung variiert wird. Eine quantisierte Metrik und ein universeller masseloser Spin 2 folgen daraus nicht."],
    n: ["8", "7/2"], g: ["Gravitation", "Einstein", "Horizont"],
    l: [["kosmo", "leitet"], ["dyn-entropic", "dockt"], ["ext-jacobson", "karte"]], src: ["docs/OPEN_PROBLEMS.md"] });

  N({ id: "kosmo", k: "chain", col: 5, row: 1.2, t: "Kosmologische Auslesungen", st: "bedingt",
    s: "Baryonenanteil, Doppelbrechungswinkel, Vakuumenergie und Inflation als Lesungen derselben Randdaten.",
    b: ["Ω_b = (1 − 1/4π)·φ₀, β = φ₀/(4π), ρ_Λ/M̄⁴ = 3·e^(−2/α)/(4π²), n_s = 1 − 2/N, r = 12/N².",
      "Jede Lesung trägt eine eigene Bedingung (Zweigwahl, Reheating, Kalibrierung)."],
    n: ["1/(4π)"], g: ["Kosmologie", "Vakuum", "Inflation"], l: [["vgeo", "bedingt"]],
    src: ["website/public/verification/predictions_frozen.json"] });

  N({ id: "rand", k: "chain", col: 5, row: 2.4, t: "Rand Quantenfeldtheorie", st: "bedingt",
    s: "Die Naht als ein Objekt: Algebra, Zustand, modulare Zeit und ein endliches Spektraltripel mit 96 Dimensionen.",
    b: ["Zielaussage SEAM.EQUIV: Die rohe Naht ist das holomorphe (E₈)₁ Netz bei τ = i. Geschlossen bis auf zitierte Kontinuumssätze (MMST, OS Rekonstruktion). Die mikroskopische geladene Naht ist damit nicht identifiziert."],
    n: ["96", "8", "16"], g: ["Rand", "KMS", "Spektraltripel"],
    l: [["netz", "offen"], ["dyn-anomaly", "dockt"], ["ext-connes", "karte"]], src: ["docs/THEORY.md", "docs/OPEN_PROBLEMS.md"] });

  N({ id: "netz", k: "chain", col: 5, row: 3.6, t: "Physisches lokales Netz in 3+1D", st: "offen",
    s: "Ein lokaler, unitärer, wechselwirkender, chiraler Parent in 3+1 Dimensionen fehlt.",
    b: ["Gebraucht werden: Lokalität, Lorentzsymmetrie, Spiegelentkopplung, kontrolliertes Infrarot. Kein Satz in 1+1 Dimensionen schließt eine Aussage in 3+1."],
    n: ["3+1"], g: ["Raumzeit", "Kontinuum", "Lokalität"], l: [["quelle", "offen"]], src: ["docs/OPEN_PROBLEMS.md"] });

  N({ id: "quelle", k: "chain", col: 5, row: 4.8, t: "Gemeinsame Quelle · Prozessregel", st: "offen",
    s: "Die eine Regel, die Naht, Träger, Zeit und Feldantwort zugleich erzeugt. Das ist die aktive Hauptaufgabe F1.",
    b: ["Gesucht: eine ursprünglich lokale Komposition, deren Reflexion einen positiven Zustand ergibt und deren zulässige Fortsetzungen denselben Zeittransport tragen.",
      "Konkret fehlt die vollständige S_adm/ω Darstellung mit allen Feldbeinen auf derselben Operator, Domain und Referenzkarte."],
    n: [], g: ["Quelle", "Prozess", "F1"],
    l: [["p1", "offen"], ["p2", "offen"], ["dyn-fieldint", "dockt"], ["gap-f1", "karte"]], src: [FP + "PLAN.md", FP + "GESAMTBILD.md"] });

  N({ id: "vgeo", k: "chain", col: 5, row: 6, t: "Eine Einheit v_geo", st: "gesetzt",
    s: "Der Compiler ist dimensionslos. Für GeV und Sekunden braucht es genau einen Maßstab von außen.",
    b: ["No Unit Satz: Der Compiler erzeugt keine eigene Skala. Dazu kommen typisierte Schnittstellen (Koide, η_B, Axion, m_p/m_e), die bewusst keine Compilerausgaben sind."],
    n: ["1"], g: ["Metrologie", "Anker"], l: [], src: ["docs/OPEN_PROBLEMS.md"] });

  /* ───────────── b) Wo steckt die Dynamik? ─────────────
     req: Z Zustand · T Zeit · P Positivität · L Lokalität/3+1D · W Wechselwirkung
          A Auswahl ohne Zusatzwahl · O an die Originalquelle gebunden
     j = liefert im Geltungsbereich · t = teilweise/bedingt · n = scheitert dort · o = nicht geprüft */
  const D = (o) => N(Object.assign({ k: "dyn" }, o));

  D({ id: "dyn-fieldint", grp: "Quelle auswerten", t: "Gemeinsame Feldintegration der Originalquelle (F1)", st: "offen",
    s: "Die ursprüngliche Kette S_adm → Z_rel → OS → Γ vollständig auswerten, mit allen Feldbeinen auf einer Karte.",
    b: ["Das ist der aktive Hauptweg. Vorhanden: bedingte Wirkungsteile, markierter Generator, tatsächliche Selektorkontakte (R112/R113), vollständige Clifford Feld und Zeitkarte für klassischen Hintergrund.",
      "Blockiert: Der R113 Lauf meldet BLOCKED. Es fehlen das original_blocked_fiber_packet, die vollständige gewichtete physische Randform und die gebundene absolute Leerreferenz."],
    req: { Z: "t", T: "t", P: "t", L: "o", W: "t", A: "o", O: "t" }, n: ["16", "96"], g: ["Quelle", "F1", "Feldintegration"],
    l: [["quelle", "dockt"], ["gap-r113", "offen"]], src: [FP + "PLAN.md", FP + "RESULTS.md", "research/originlab/OPH-WINDING-MEMORY-20261119.md"] });

  D({ id: "dyn-wordsource", grp: "Quelle auswerten", t: "Geordnete Wörter mit Gedächtnis", st: "belegt",
    s: "Die Quelle als Folge geordneter Operationen. Was man nicht sieht, wirkt als Erinnerung weiter.",
    b: ["Im Originalbaustein reichen exakt 68 innere Richtungen für alle Wortlängen: 5 sichtbar, 63 Gedächtnis. Zwei Quellenrückwirkungen geben 6144, dieselbe Folge mit J dazwischen gibt 0.",
      "Zukunftsquotient und positive GNS Rekonstruktion bauen aus vollständig gegebenen Antworten die minimale Prozessdarstellung. Die Regel, die diese Antworten erzeugt, fehlt."],
    req: { Z: "j", T: "t", P: "j", L: "t", W: "t", A: "n", O: "t" }, n: ["68", "5", "63", "6144"], g: ["Gedächtnis", "Wort", "Rekonstruktion"],
    l: [["quelle", "dockt"], ["pat-loss", "karte"], ["ur-future", "karte"]], src: [FP + "GESAMTBILD.md"] });

  D({ id: "dyn-openword", grp: "Quelle auswerten", t: "Regel für offene Übergänge", st: "neu",
    s: "Zwei Geschichten gelten als gleich, wenn sie dieselben geladenen Endlinien, Orte und Zukunftsantworten tragen.",
    b: ["Dann käme die Einheitswindung aus dem Endliniencharakter und die Ortsbewegung aus den tatsächlichen Inklusionen.",
      "Der entscheidende Test ist benannt und noch nicht ausgeführt: die beiden Wurzelkontraktionen a und b am vollständigen markierten P1 Streifen vergleichen."],
    req: { Z: "o", T: "o", P: "o", L: "o", W: "o", A: "o", O: "o" }, n: ["126", "72"], g: ["Windung", "Endlinien", "Test offen"],
    l: [["org-winding", "karte"], ["quelle", "dockt"]], src: [FP + "TFPT-Ursprungshypothesen-2026-11-12.md"] });

  D({ id: "dyn-transfer", grp: "Relaxation", t: "Transfer Attraktor", st: "bedingt",
    s: "Wiederholtes Anwenden eines positiven Transfers zieht alles auf einen Zustand zusammen.",
    b: ["Spektrum {1, (2/3)⁶, (1/3)⁶}, Lücke 6·log(3/2), eindeutige stationäre Richtung. Daraus der Satz „parameterfrei ist ein Theorem“ unter den gelückten Voraussetzungen.",
      "Grenze: Der Satz ist kein Selektor von T. Gleiche Kontraktionsklasse heißt nicht gleiche physische Uhr."],
    req: { Z: "j", T: "t", P: "j", L: "o", W: "o", A: "n", O: "t" }, n: ["2/3", "1/3", "6"], g: ["Attraktor", "Relaxation"],
    l: [["attraktor", "dockt"], ["dyn-relax", "gleich"]], src: ["origin_theory.tex", "docs/OPEN_PROBLEMS.md"] });

  D({ id: "dyn-relax", grp: "Relaxation", t: "Relaxations Clock · träge Paarregel", st: "bedingt",
    s: "Raten r(n) = −6·log(1 − n/3). Sechs Schritte erzeugen beide Lücken.",
    b: ["Träge Paarregel bleiben/springen/lecken = (1/2, 1/6, 1/3) innerhalb treuer Leiter, Paritäts und Positivitätsannahmen. Ketten zu Nariai, Born Quadrat und Entropie verbinden die Raten mit bedingten physikalischen Lesarten."],
    req: { Z: "t", T: "t", P: "j", L: "o", W: "o", A: "n", O: "o" }, n: ["2/3", "1/3", "1/2", "1/6", "6"], g: ["Relaxation", "Clock"],
    l: [["attraktor", "dockt"], ["pat-23", "gleich"]], src: ["tfpt_horizon_readouts.tex"] });

  D({ id: "dyn-markov", grp: "Offene Systeme", t: "Markov und Lindblad Zeit", st: "aus",
    s: "Den Transfer als kontinuierlichen offenen Prozess exp(t·log T) lesen.",
    b: ["Für die feste M₄ Karte ist die homogene GKSL Zeit ausgeschlossen. Ein GKLS q Intervall ist vorhanden und wiederverwendbar.",
      "Neu am 19.11.: In der kovarianten Klasse des Viererrings muss die Kohärenz mindestens so schnell zerfallen wie die Population mischt (γ ≥ m). Mischen ohne Kohärenzverlust ist dort unmöglich."],
    req: { Z: "t", T: "n", P: "j", L: "o", W: "o", A: "n", O: "n" }, n: ["4", "3/4"], g: ["Lindblad", "Kanal", "Kohärenz"],
    l: [["attraktor", "dockt"], ["ext-oph", "karte"]], src: ["research/originlab/OPH-WINDING-MEMORY-20261119.md"] });

  D({ id: "dyn-reader", grp: "Offene Systeme", t: "P1 als wechselnder Leser", st: "bedingt",
    s: "P1 ist die Wahl einer Auslese. Zwei Auslesen derselben Kurve vertragen sich nicht, und der Wechsel zwischen ihnen ist die Bewegung.",
    b: ["Die Picardkurve hat 48 Symmetrien, die markierte Auslese sieht nur 12. Zwei Deckprojektoren P, Q vom Rang 16 kommutieren nicht. Daraus entsteht ein vollständig positiver Kanal mit nachweisbarer Zustandsänderung.",
      "Grenzen: Der Prozess ist nicht primitiv. Die vier Quotienten sind als native Vierteluhr ausgeschlossen (ihre Permutationsgruppe ist A₄, ohne Vierzyklus). Hilbertraum, Messoperation und Rate sind Zusatzannahmen."],
    req: { Z: "t", T: "t", P: "j", L: "o", W: "o", A: "n", O: "o" }, n: ["48", "12", "16", "1/3"], g: ["Auslese", "Picard", "Kanal"],
    l: [["p1", "dockt"], ["pat-picard", "karte"]], src: ["research/originlab/F1-ORIGIN-SCENARIOS-20261116.txt"] });

  D({ id: "dyn-d4hist", grp: "Offene Systeme", t: "Mischung von D₄ Geschichten", st: "aus",
    s: "Die Welt als Mittel über mehrere mögliche Geschichten derselben Symmetrieklasse.",
    b: ["Geprüft und in der erklärten Klasse ausgeschlossen (History Gram und KMS für alle D₄). Mittelung über A₄ oder S₄ ergibt denselben vollständig depolarisierenden Kanal: Selbst eine vollständige reduzierte Zeitentwicklung wählt die markierte Geschichte nicht."],
    req: { Z: "n", T: "n", P: "j", L: "o", W: "o", A: "n", O: "t" }, n: ["4"], g: ["Geschichten", "Mittelung"],
    l: [["pat-loss", "karte"]], src: ["research/originlab/F1-HISTORY-GRAM-AND-ALL-D4-KMS-NO-GO-20261116.md", FP + "TFPT-Gesamtkonsolidierung-2026-11-12.md"] });

  D({ id: "dyn-modular", grp: "Zeit aus dem Zustand", t: "Modulare und thermische Zeit", st: "aus",
    s: "Der Zustand selbst erzeugt seine Zeit (Tomita–Takesaki, KMS mit β = 2π).",
    b: ["Als Baustein vorhanden und wiederverwendbar (Transfer/KMS Rekonstruktion). Als gemeinsame physische Uhr gescheitert: Die Viersektor Clock und Thermal Time Gleichsetzungen scheitern unter ihren Voraussetzungen. Einzelne Modularzeiten sind auf überlappenden Paaralgebren unverträglich, die gemeinsame absolute Modularzeit hat zweiseitiges Spektrum.",
      "Keine allgemeine Widerlegung relationaler Zeit."],
    req: { Z: "j", T: "n", P: "j", L: "n", W: "o", A: "t", O: "t" }, n: ["2π"], g: ["Zeit", "KMS", "Modular"],
    l: [["rand", "dockt"], ["clocks", "dockt"]], src: ["docs/OPEN_PROBLEMS.md", FP + "TFPT-Ursprungshypothesen-2026-11-12.md"] });

  D({ id: "dyn-selfdual", grp: "Zeit aus dem Zustand", t: "Selbstdualer Punkt als Auswahl (G1, G3)", st: "bedingt",
    s: "Die Naht steht am selbstdualen Modul τ = i. Dort braucht man keine absolute Referenz.",
    b: ["Exakt: c₃ = E₂(i)/24 ist die Ladungsenergie pro Cartanrichtung jedes selbstdualen Gitters. Das verklebte E₈ trifft c₃, das unverklebte D₅ + A₃ liefert 0,289·c₃.",
      "Gescheitert als Produzent der Ortsfaktoren der vier R113 Marken: Alle Punktwerte am Torus enthalten die lemniskatische Konstante ϖ, die TFPT Ausgaben nur π. Falls der Torus die Quelle ist, können die Auslesungen nur Energien sein.",
      "G3 (Horizont zuerst, thermofeldverdoppelter Zustand) hat noch keine eigene Probe."],
    req: { Z: "t", T: "t", P: "o", L: "n", W: "o", A: "j", O: "n" }, n: ["1/(8π)", "24", "12"], g: ["Auswahl", "Modul", "Torus"],
    l: [["taui", "dockt"], ["org-c3-modular", "karte"]], src: ["experiments/theory-contracts/radical-global-hypotheses-20261118/README.md"] });

  D({ id: "dyn-vacuum", grp: "Zeit aus dem Zustand", t: "Grundzustand der Vierteldrehung", st: "bedingt",
    s: "Die Vierteldrehung ist die Fouriertransformation. Ihr tiefster fester Zustand ist exakt gaußisch und gibt exp(−2α) ohne Rest.",
    b: ["Gerechnet am 19.11.: Der Grundzustand ψ₀ = 2^(1/4)·exp(−πx²) ist seine eigene Fouriertransformierte. Seine Nullpunktsenergie je Quadratur ist c₃, und der Erwartungswert des Verschiebeoperators entlang einer Wurzel ist genau exp(−2α).",
      "Gegenproben: Ein einzelnes E₈ Gitter verfehlt um +142 σ. Viele unabhängige Keime nähern sich mit 1/N, innerhalb 2 σ erst ab N = 1360. Der nächste feste Zustand (n = 4) verfehlt um mehr als 10⁵ σ.",
      "Grenzen: Die Konvention h = 1 und das Mitzählen beider Quadraturen wurden in Kenntnis des Zielwerts gewählt. Gaußische Momente enthalten nur π, das passt zu den Ausgaben. Warum die Quelle in diesem Zustand sitzt, ist nicht hergeleitet."],
    req: { Z: "j", T: "t", P: "j", L: "o", W: "o", A: "t", O: "o" }, n: ["1/(8π)", "4", "2", "1360"], g: ["Zustand", "Fourier", "Alpha", "Auswahl"],
    l: [["taui", "dockt"], ["alpha", "dockt"], ["dyn-selfdual", "gleich"], ["quelle", "dockt"]], src: ["experiments/theory-contracts/landkarte-origin-scenarios-20261119/README.md"] });

  D({ id: "dyn-p1kern", grp: "Zeit aus dem Zustand", t: "Der P1 Kern trägt den Grundzustand", st: "bedingt",
    s: "Auf einer runden Kugel als Schnitt liefert der originale Randoperator Vierteluhr, Periode 2π und c₃ = T/4 aus einem Stück.",
    b: ["Gerechnet am 19.11.: Der originale Dirac Randoperator hat auf S² ganzzahliges Spektrum. Daraus folgen eine Vierteluhr der Ordnung 4, die Eigenperiode 2π und ein flaches Kanalgewicht von 4c₃. Ein Kanal gibt exp(−2α).",
      "Kette: runde S² → Periode 2π → Vierteldrehung tauscht Zeitkreis und Raumkreis → quadratischer Torus → Nullmode trägt T/4 = c₃. Der Radius fällt heraus.",
      "Ehrliche Folge: KMS Weg, Grundzustand und Toruswert sind dieselbe Zahl. Das ist ein Beleg, nicht drei.",
      "Nachgeprüft am Original: Die positive Hälfte steht dort (Calderón Projektor, erster positiver Eigenwert). Die runde Kugel steht dort nicht und ist nicht nötig. Das Original hat eine flache Kugel mit vier Marken bei μ₄, und deren Doppelüberlagerung ist exakt der quadratische Torus.",
      "Warum liegen die Marken bei μ₄? Geprüft am 19.11.: Die Spiegelung allein reicht nicht. Auch vier Ecken auf der Naht geben nur ein Rechteck mit freiem Seitenverhältnis. Erst ein Element der Ordnung 4 erzwingt das Quadrat. Die Klebegruppe von D₅ ist zyklisch und könnte das leisten; dass sie auf den Marken wirkt, ist nicht gezeigt.",
      "Grenzen: Die Lage der Marken und genau eine chirale Nullmode sind Annahmen. δ_top und φ_tree erklärt das Bild nicht. Das Quadrat ist nur stabil, wenn die Spiegelung die Verdrehung einfriert."],
    req: { Z: "j", T: "j", P: "j", L: "o", W: "o", A: "t", O: "o" }, n: ["1/(8π)", "4", "2π"], g: ["Zustand", "Kern", "Torus", "Alpha"],
    l: [["dyn-vacuum", "gleich"], ["taui", "dockt"], ["alpha", "dockt"], ["quelle", "dockt"]], src: ["experiments/theory-contracts/p1-kernel-ground-state-20261119/README.md"] });

  D({ id: "dyn-unitary", grp: "Geschlossene Quantendynamik", t: "Unitäre Ausführung · Quantenautomat", st: "bedingt",
    s: "Die Welt als reversible Schrittfolge exp(−itH) oder als Kollisionsmodell.",
    b: ["Endliche unitäre, Kollisions und kontinuierliche Erweiterungen existieren im Korpus. Dieselbe Populationsregel lässt aber verschiedene kohärente Ausführungen zu. Welche gilt, ist genau die Ursprungsfrage.",
      "Reversible endliche Ausführung kehrt wieder (Rekurrenz). Ein zusätzlicher Oszillator oder Kopplungswert ist erst Fortschritt, wenn seine Auswahl folgt."],
    req: { Z: "t", T: "j", P: "j", L: "t", W: "t", A: "n", O: "n" }, n: [], g: ["Unitär", "Rekurrenz", "Automat"],
    l: [["quelle", "dockt"], ["ur-first", "karte"]], src: ["docs/TFPT_GEMEINSAME_AUSFUEHRUNG_2026-09-14_v1.0.md"] });

  D({ id: "dyn-hamilton", grp: "Geschlossene Quantendynamik", t: "Lokale Hamiltonquelle mit Gaußsektor", st: "bedingt",
    s: "Ein explizites lokales Gittermodell mit Zeitentwicklung, Grundzustand und Lücke.",
    b: ["Bewiesen: Für den entkoppelten DET16 Produktparent Spektrum {0…N}, eindeutiger Grundzustand, Lücke 1. Für eine gleichmäßig lokale Hamiltonklasse auf Z² existieren Zeitentwicklung, quasilokale Algebra, Grund und KMS Zustände.",
      "Ausgeschlossen: die direkte Gleichsetzung der alten Cliffordquelle mit einem separierten DET16 Transfer. Exaktes DET16 ist keine universelle TFPT Pflicht."],
    req: { Z: "j", T: "j", P: "j", L: "t", W: "t", A: "n", O: "n" }, n: ["16", "1"], g: ["Hamilton", "Spiegel", "Lücke"],
    l: [["netz", "dockt"], ["sm", "dockt"]], src: ["verification/status_ledger.csv"] });

  D({ id: "dyn-clockgear", grp: "Geschlossene Quantendynamik", t: "Das Clock Getriebe selbst als Dynamik", st: "aus",
    s: "Die diskreten Takte (5, 6, 30, 12, μ₄) direkt als Zeit oder Raum lesen.",
    b: ["Die Algebra ist exakt und bleibt. Als physische Zeit oder Raum scheitert die direkte Lesart: Eine native Clock im festen Hintergrund affin energiegleich zu machen ist ausgeschlossen. Aus zwei Fünfertakten vier euklidische Richtungen zu gewinnen scheitert (sie erzeugen A₆, das in SO(3) × SO(3) nur trivial wirkt).",
      "Periodische Anzeigen brauchen für mehrere Umläufe zusätzliche Zählinformation."],
    req: { Z: "o", T: "n", P: "o", L: "n", W: "o", A: "t", O: "j" }, n: ["5", "30", "12", "4", "360"], g: ["Clock", "Takt"],
    l: [["clocks", "dockt"], ["org-fivebeat", "karte"]], src: [FP + "TFPT-Ursprungshypothesen-2026-11-12.md"] });

  D({ id: "dyn-recursion", grp: "Rekursion und Netze", t: "Rekursive Bindung · Petz Rückkanal", st: "belegt",
    s: "Auslesen und Wiederherstellen als Schleife K = R₂E₂. Die Organisation bleibt stabil, während Zustände wechseln.",
    b: ["Der Rückkanal lebt auf End(C⁵) mit 25 Dimensionen. Die maximal gemischte Referenz I₅/5 bleibt fest. Vektorisiert ist das der verschränkte Paarvektor Ω₅.",
      "Kernidee: Ein Strukturfixpunkt kann Dynamik tragen. Auslesen, Paarbindung und physische Zeit bleiben aber verschiedene Anwendungen derselben Abbildung."],
    req: { Z: "j", T: "t", P: "j", L: "o", W: "o", A: "t", O: "t" }, n: ["25", "5"], g: ["Rekursion", "Petz", "Fixpunkt"],
    l: [["dyn-wordsource", "karte"], ["ur-two", "karte"]], src: [FP + "GESAMTBILD.md", "_newest2/TFPT_Gesamtdokumentation2_20260927.md"] });

  D({ id: "dyn-hypergraph", grp: "Rekursion und Netze", t: "Hypergraph und Netzwachstum", st: "bedingt",
    s: "Eine lokale Mittelungsregel plus Wachstum mit Balance Tor läuft von selbst auf E₈ zu und stoppt dort.",
    b: ["Belegt auf der Skelettebene: Die Marken (1,2,3,4,5,6,4,2,3) sind Attraktor der Regel m ← (A+2I)/4·m. Wachstum mit ρ ≤ 2 bleibt genau bis E₈ endlich.",
      "Grenzen: Es braucht einen Seed mit drei Zweigen, also (2,3,5), also wieder P2. Und es entsteht keine Raumzeit.",
      "Offener Anschluss: die tatsächlichen effektiven Vertices derselben Quelle als gewichtete Hyperkanten auswerten und rekursiv verketten."],
    req: { Z: "t", T: "t", P: "o", L: "n", W: "o", A: "t", O: "n" }, n: ["2", "3", "5", "120", "30"], g: ["Hypergraph", "Netzwerk", "Attraktor"],
    l: [["e8", "dockt"], ["pat-hyper", "karte"], ["ext-wolfram", "karte"]], src: ["hypergraph.md", "experiments/theory-contracts/recursive-mirrored-source-20261103/ANALYSIS.txt"] });

  D({ id: "dyn-spectral", grp: "Wirkung und Felder", t: "Spektralwirkung · Higgs aus der Materie", st: "bedingt",
    s: "Die Dynamik des Higgsfelds aus dem ganzen Materieoperator ablesen: zweite Spur → Bewegung, vierte → Quartik.",
    b: ["Im bedingten Raumzeitanschluss exakt ausgeführt. Dasselbe Higgs koppelt alle Familien, seine Integration erzeugt Vierfermionterme mit Produkten der bestehenden Yukawas.",
      "Benannte Lücke: Zwei Größen, die dieselbe Spur lesen müssten, stehen im Verhältnis −12 statt −2 (Faktor sechs). Die nackte Spektralwirkung unter Marginalisierung durch einen Skalar ist ausgeschlossen. Ein unbeschränktes Higgsvolumen lässt sich so nicht normieren: Das Gewicht muss die Quelle liefern."],
    req: { Z: "t", T: "o", P: "t", L: "t", W: "j", A: "n", O: "t" }, n: ["96", "6", "−12", "−2"], g: ["Higgs", "Spektral", "Wirkung"],
    l: [["flavor", "dockt"], ["rand", "dockt"], ["ext-connes", "karte"]], src: [FP + "GESAMTBILD.md"] });

  D({ id: "dyn-entropic", grp: "Wirkung und Felder", t: "Verschränkungsgleichgewicht · entropische Wirkung", st: "bedingt",
    s: "Gravitation als Thermodynamik der Verschränkung am Rand.",
    b: ["Liefert bedingt die Einstein Gleichung mit beiden Koeffizienten. Bianconis entropische Wirkung trifft die Konstante exakt (β′ = c₃/6), verfehlt aber den R² Koeffizienten um genau 3·(8π)⁹ ≈ 10¹³.",
      "Setzt die 4D Theorie voraus und kann deshalb nicht als Beleg für sie zählen."],
    req: { Z: "t", T: "t", P: "o", L: "t", W: "o", A: "n", O: "n" }, n: ["8", "1/6"], g: ["Gravitation", "Entropie"],
    l: [["grav", "dockt"], ["ext-bianconi", "karte"], ["ext-jacobson", "karte"]], src: ["docs/OPEN_PROBLEMS.md"] });

  D({ id: "dyn-master", grp: "Wirkung und Felder", t: "TFPT4D Masterroute · Gitterwirkung", st: "offen",
    s: "Eine explizite endliche 4D Gitterwirkung (Wilson + Ginsparg–Wilson + ein relevanter Sektor) als einziges Wachstumsprinzip.",
    b: ["Arbeitshypothese: Es fehlen nicht zwanzig Formeln, sondern eine Vervollständigung der Nahtalgebra, lokal, reflexionspositiv, näherungsweise quanten Markov, anomaliefrei.",
      "Dazu ein bedingter Dimensionswähler (d = 4 eindeutig unter vier Annahmen, μ₄ geht nirgends ein) und der störungstheoretische Epstein–Glaser Vertrag. Alles Konstruktive ist offen, mit registrierten Abbruchkriterien."],
    req: { Z: "o", T: "o", P: "o", L: "t", W: "t", A: "o", O: "o" }, n: ["4"], g: ["Gitter", "4D", "Masterroute"],
    l: [["netz", "dockt"]], src: ["docs/OPEN_PROBLEMS.md", "tfpt_research_contracts.tex"] });

  D({ id: "dyn-anomaly", grp: "Wirkung und Felder", t: "Anomalie Zufluss · Rand eines E₈ Zustands", st: "bedingt",
    s: "Die Naht als Rand eines gelückten 2+1D Zustands (Kitaevs E₈ Phase, c = 8).",
    b: ["Ein gelücktes p+ip Gittermodell trifft das Ziel auf jeder berechenbaren Ebene (c = 8, 248 Ströme, GSD 1, Reflexionspositivität). Extern bleibt der zitierte Kontinuumsgrenzwert.",
      "Die mikroskopische geladene Naht mit Halbladungsfeld ist damit nicht identifiziert."],
    req: { Z: "j", T: "t", P: "j", L: "t", W: "o", A: "t", O: "n" }, n: ["8", "16", "248"], g: ["Anomalie", "Rand", "Kitaev"],
    l: [["rand", "dockt"], ["ext-kitaev", "karte"]], src: ["docs/OPEN_PROBLEMS.md", "research/originlab/F1-E8-ANOMALY-INFLOW-ORIGIN-20261116.md"] });

  D({ id: "dyn-holo", grp: "Geometrische Ursprünge", t: "P1 als holografische Randantwort", st: "bedingt",
    s: "Aus einer einzigen freien Wirkung im Raum folgen Zahl und Zeitform am Rand der Kugel.",
    b: ["Der Dirichlet zu Neumann Kern gibt R·C(x, −x) = 1/(8πκ). Für κ = 1 ist das c₃, als Lesung eines ganzen Operators samt positiver radialer Zeitform.",
      "Grenzen: κ wird von Positivität und Geometrie nicht gewählt. Die antipodale Ablesung ist keine positive Reflexion. R³ ist keine hergeleitete Raumzeit.",
      "Nächster harter Entscheid (Priorität A): die native Operatorgleichheit prüfen."],
    req: { Z: "t", T: "j", P: "j", L: "t", W: "n", A: "n", O: "o" }, n: ["1/(8π)", "1/2"], g: ["Holografie", "Kugel", "Rand"],
    l: [["p1", "dockt"], ["org-c3-sphere", "karte"]], src: ["research/originlab/F1-ORIGIN-SCENARIOS-20261116.txt"] });

  D({ id: "dyn-lorentz", grp: "Geometrische Ursprünge", t: "Größerer lorentzscher Träger", st: "bedingt",
    s: "E₈ als Schatten eines Gitters mit einer Zeitrichtung: II(9,1) = E₈ ⊕ U.",
    b: ["Exakt: Für eine markierte Nullrichtung e ist e⊥/Ze = E₈. Aber jede Verschiebung im großen Gitter wirkt auf dem rekonstruierten E₈ als Identität. Bewegung im Größeren ist nach der Auslese unsichtbar.",
      "Eine negative Gitterrichtung ist keine hergeleitete Uhr."],
    req: { Z: "o", T: "n", P: "n", L: "o", W: "o", A: "n", O: "o" }, n: ["9", "1", "10", "8"], g: ["Lorentz", "Gitter", "Schatten"],
    l: [["e8", "dockt"], ["ext-string", "karte"], ["ur-shadow", "karte"]], src: ["research/originlab/F1-ORIGIN-SCENARIOS-20261116.txt"] });

  D({ id: "dyn-reset", grp: "Kosmische Lesarten", t: "P1 als Singularitätsnaht · Reset", st: "aus",
    s: "Die Naht als Neustart: alte Information bleibt erhalten, die sichtbare Welt beginnt neu.",
    b: ["Vier Varianten sind ausgeschlossen: thermische Doppelung als exakter Reset, Mischung der Stromreihenfolgen, Neustart im selben positiven Vakuumraum, der vorhandene Partner als positives Gedächtnis.",
      "Was bleibt: Ein global reversibler exakter Reset verlangt mehrere orthogonale Vakuumdarstellungen. Die entscheidende Probe liegt im gewichteten Originalschnitt.",
      "Gerechnet am 19.11.: Ein exakter Neustart, der alle Information behält, existiert genau dann, wenn die verborgene Seite mindestens r·d Dimensionen hat. Die Schranke ist scharf. Ausgeschlossen sind also nur die Fälle ohne Platz. P1 müsste dann ein lokaler Schnitt (Horizont) sein, nicht das ganze Vakuum."],
    req: { Z: "n", T: "o", P: "n", L: "o", W: "o", A: "n", O: "t" }, n: [], g: ["Reset", "Singularität", "Information"],
    l: [["p1", "dockt"], ["dyn-cyclic", "gleich"], ["gap-g3", "karte"]], src: ["semantic_explorer/research/program.json · problem:p1-singularity-seam-reset", "experiments/theory-contracts/landkarte-origin-scenarios-20261119/README.md"] });

  D({ id: "dyn-cyclic", grp: "Kosmische Lesarten", t: "Zyklische Deutung", st: "deutung",
    s: "Universum → Kollaps → Horizont als Naht → Attraktor als Gesetz → neues Universum.",
    b: ["Im Origin Papier ausdrücklich eine Deutung, keine Folge der Axiome. Ihre Kinematik ist maschinengeprüft, die Übergangsdynamik ist offen.",
      "Zwei Befunde begrenzen sie: Der bestehende Familienfixpunkt ist kein exakter Quantenzyklus, und Zustandsrückkehr ist kein hergeleiteter räumlicher Punkt."],
    req: { Z: "t", T: "o", P: "o", L: "o", W: "o", A: "t", O: "o" }, n: ["30", "2/3"], g: ["Zyklus", "Horizont", "Kosmos"],
    l: [["ext-origin", "karte"], ["attraktor", "dockt"], ["ext-ccc", "gleich"]], src: ["origin_theory.tex"] });

  D({ id: "dyn-oph", grp: "Vergleichsmodelle", t: "Kohärente OPH Randwege", st: "bedingt",
    s: "Eine winzige Grammatik aus Phasendrehung und Projektion auf den nächsten Flächenrand erzeugt die Lorentzgruppe.",
    b: ["In deklarierten endlichen Modellen: projektiver Abschluss PSL(2,C). Zwei kohärent zusammengeführte Rückwege geben diag(63,65)/128. Die Differenz liegt exakt zwischen −1/64 und +1/64, mit klassischer Wegaufzeichnung ist sie null.",
      "Offen: Eine gemeinsame positive Entwicklung muss Vakuummischung und geladene Rückkehr auf derselben ursprünglichen Zeit tragen."],
    req: { Z: "t", T: "t", P: "j", L: "t", W: "t", A: "n", O: "n" }, n: ["1/64", "60", "5/8", "25/64"], g: ["OPH", "Kohärenz", "Lorentz"],
    l: [["ext-oph", "karte"], ["quelle", "dockt"]], src: ["research/originlab/OPH-COHERENT-HISTORY-20261119.md"] });

  D({ id: "dyn-alphafield", grp: "Vergleichsmodelle", t: "α als Feld mit eigenem Potential", st: "neu",
    s: "Die Wurzel der Kubik ist ein Minimum des historischen Potentials U(α). Vielleicht ist α selbst dynamisch.",
    b: ["Die Kubik ist keine Laufgleichung (sie gäbe d(1/α)/dL = −15,48 statt +0,65). Die Lesart als Feld ist im Korpus benannt und nicht weiter geprüft."],
    req: { Z: "o", T: "o", P: "o", L: "o", W: "o", A: "o", O: "o" }, n: ["137"], g: ["Alpha", "Potential", "Test offen"],
    l: [["alpha", "dockt"]], src: ["experiments/theory-contracts/radical-global-hypotheses-20261118/README.md"] });

  window.TFPT.req = [
    ["Z", "Zustand", "Legt der Ansatz fest, in welchem Zustand die Welt ist?"],
    ["T", "Zeit", "Liefert er eine physische Zeit, nicht nur Schritte oder Perioden?"],
    ["P", "Positivität", "Sind Wahrscheinlichkeiten und Reflexion positiv?"],
    ["L", "Ort · 3+1D", "Entstehen lokale Orte und drei Raumrichtungen?"],
    ["W", "Wechselwirkung", "Trägt er echte gekoppelte Felder?"],
    ["A", "Ohne Zusatzwahl", "Kommt er ohne frei gewählte Zutat aus?"],
    ["O", "Originalquelle", "Ist er an die ursprüngliche P1 Quelle gebunden?"],
  ];

  /* ───────────── Vorhersagen (Werte aus predictions_frozen.json, Züge aus RADICAL20261118) ───────────── */
  window.TFPT.outputs = [
    { id: "ALPHA_INV", t: "Feinstrukturkonstante α⁻¹", f: "Wurzel von F(α) = 0", v: "137,0359992168", m: "137,035999177(21)", pull: 1.90, layer: "Kern", note: "Kein Fit: 5/4 ist Gitterzahl, der Fitwert wäre 1,2404 ± 0,0051. Nächster Test: genauere Einzelmessung." },
    { id: "LAMBDA_C", t: "Cabibbo Winkel λ_C", f: "√(φ₀(1−φ₀))", v: "0,224376", m: "0,2243(8)", pull: 0.10, layer: "Kern" },
    { id: "S23_CKM", t: "CKM s₂₃", f: "φ₀/(1+λ_C)", v: "0,043428", m: "0,0411(12)", pull: 1.94, layer: "Kern" },
    { id: "S13_CKM", t: "CKM s₁₃", f: "λ_C³/3", v: "0,003765", m: "0,00382(20)", pull: -0.28, layer: "Kern" },
    { id: "DELTA_CKM", t: "CKM Phase δ", f: "π/3 + 3λ_C²", v: "68,65°", m: "65,7(30)°", pull: 0.98, layer: "bedingt", note: "führende Ordnung" },
    { id: "SIN2_THETA12", t: "Sonnenwinkel sin²θ₁₂", f: "1/3 − φ₀/2", v: "0,306747", m: "0,308(12)", pull: -0.10, layer: "bedingt", note: "Vorhersage des Registers. Test: JUNO." },
    { id: "SIN2_THETA13", t: "Reaktorwinkel sin²θ₁₃", f: "φ₀·e^(−5/6)", v: "0,023108", m: "0,02215(57)", pull: 1.68, layer: "Kern" },
    { id: "OMEGA_B", t: "Baryonenanteil Ω_b", f: "(1 − 1/4π)·φ₀", v: "0,048941", m: "0,0493(9)", pull: -0.40, layer: "bedingt" },
    { id: "BETA", t: "Kosmische Doppelbrechung β", f: "φ₀/(4π)", v: "0,2424°", m: "0,20(8)°", pull: 0.53, layer: "bedingt", note: "Test: LiteBIRD, CMB-S4." },
    { id: "MMU_MTAU", t: "m_μ/m_τ", f: "(8/7)·φ₀", v: "0,060768", m: "0,05946", pct: "+2,20 %", layer: "Kern", note: "Quellverhältnis, Transfer zur Polmasse offen." },
    { id: "ME_MMU", t: "m_e/m_μ", f: "(12/7)·φ₀²", v: "0,0048467", m: "0,004836", pct: "+0,21 %", layer: "Kern", note: "Quellverhältnis, Transfer zur Polmasse offen." },
    { id: "MU_MD", t: "m_u/m_d", f: "55/117", v: "0,4701", m: "≈ 0,47", layer: "Kern" },
    { id: "MC_MS", t: "m_c/m_s", f: "(34/47)/φ₀", v: "13,605", m: "≈ 11 bis 14", layer: "Kern" },
    { id: "MT_MB", t: "m_t/m_b", f: "(3/26)/φ₀²", v: "40,81", m: "≈ 41", layer: "Kern" },
    { id: "MSCAL", t: "Skalaron M/M̄", f: "c₃^(7/2)", v: "1,2565·10⁻⁵", m: "Inflation, r", layer: "bedingt", note: "r > 0,01 widerlegt den R² Zweig." },
    { id: "RHOL", t: "Vakuumenergie ρ_Λ/M̄⁴", f: "3·e^(−2/α)/(4π²)", v: "7,125·10⁻¹²¹", m: "Planck", layer: "bedingt" },
    { id: "DCP_NU", t: "Neutrino CP Phase", f: "240° (Phasengitter)", v: "240°", m: "DUNE, Hyper-K", layer: "zugewiesen" },
    { id: "TH23", t: "sin²θ₂₃", f: "1/2 (Textur)", v: "0,5", m: "JUNO, DUNE", layer: "zugewiesen" },
    { id: "R", t: "Tensorverhältnis r", f: "12/N², N = 50 bis 60", v: "0,0033 bis 0,0048", m: "CMB-S4", layer: "Band" },
    { id: "NS", t: "Spektralindex n_s", f: "1 − 2/N", v: "0,960 bis 0,967", m: "Planck, DESI", layer: "Band" },
  ];
})();
