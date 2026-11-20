/* TFPT Landkarte · Teil 4: zweite Runde, zehn neue Hypothesen (19. November 2026).
   Gebaut aus den Urteilen der ersten Runde (HYPCHECK20261119) und aus einer parallelen Prüfung
   in einer anderen Sitzung. Deren Zahlen sind hier nur berichtet, nicht nachgerechnet.
   Der obere Teil ist der Vorschlag im Wortlaut. Der Block am Ende trägt das Urteil der Prüfung
   HYPCHECKB20261119 (experiments/theory-contracts/landkarte-hypothesen-runde2-20261119). Kein Ledgerwechsel. */
(function () {
  const T = window.TFPT;
  const SRC = ["tfpt_landkarte/checks/vorrechnung_runde2.py", "experiments/theory-contracts/landkarte-hypothesen-20261119/README.md"];
  const H = (o) => T.nodes.push(Object.assign({ k: "hyp", st: "neu", runde: 2, src: SRC }, o));

  H({ id: "hyp-nachbar", nr: 11, t: "Zeit ist die Richtung zwischen Viererschnitt und Fünferschnitt",
    s: "Die beiden Schnitte sind Nachbarknoten im E₈ Diagramm. Was sie gemeinsam haben, lässt genau eine Richtung frei.",
    b: ["Streicht man den Knoten mit Marke 4, bleibt D₅ ⊕ A₃, also Träger und Familie. Streicht man den Nachbarknoten mit Marke 5, bleibt A₄ ⊕ A₄, der Fünfertakt. Bisher waren das zwei getrennte Bilder.",
      "Die Schnittregel der ersten Runde sagt: Bei Einheitswindung verschwindet die Randphase nur am Fünferschnitt. Der Fünferschnitt ist der phasenfreie, der Viererschnitt der mit Ladungen.",
      "Streicht man beide Knoten, bleibt A₄ ⊕ A₃ vom Rang 7. Die achte Richtung ist ein einzelner Gittervektor der Norm 20 = 4 · 5.",
      "Hypothese: Diese eine Richtung ist die Zeit. Vierteluhr und Fünfertakt sind ihre zwei Ablesungen, 20 ist ihr gemeinsamer Takt."],
    pre: "Rest A₄ ⊕ A₃ hat Determinante 20, der senkrechte Gittervektor hat Norm 20.",
    why: ["Gemeinsame Schnittregel: Phase null bei Windung eins nur für p = 5", "Das Uhrengetriebe ist exakt, hat aber keine Zählrichtung", "Vierer und Fünferbild standen bisher unverbunden nebeneinander"],
    gives: ["T", "A"],
    test: "Die Richtung der Norm 20 mit Hyperladung und B − L im Träger vergleichen. Dann prüfen, ob die Drehung entlang ihr auf D₅ ⊕ A₃ als Vierteluhr und auf A₄ ⊕ A₄ als Fünfertakt wirkt.",
    kill: "Wenn die Richtung auf dem Träger nur eine schon vergebene Ladung ist oder keine der beiden Uhren erzeugt.",
    n: ["20", "4", "5", "7"], g: ["Zeit", "Clock", "Singularität", "McKay", "Gitter"],
    l: [["clocks", "testet"], ["e8", "testet"], ["pat-a4cut", "testet"], ["hyp-kac", "karte"], ["org-fivebeat", "testet"], ["dyn-clockgear", "testet"]] });

  H({ id: "hyp-beauville", nr: 12, t: "Die vier Marken sind die vier Ausnahmestellen einer Kurvenfamilie der Stufe 4",
    s: "Es gibt eine bekannte Familie von Tori über der Kugel, deren vier Ausnahmestellen genau harmonisch liegen. P1 wäre dann der Raum aller Torusformen mit einer Markierung der Stufe 4.",
    b: ["Die erste Runde fand eine Fläche mit den Fasern A₃ und D₅, aber nur drei Ausnahmestellen. Verlangt man zusätzlich alle Halbperioden, y² = x(x + t²)(x + 1), dann liegen die Ausnahmestellen bei 0, ∞, 1 und −1. Das ist die Lage der Marken.",
      "Die Fasern sind I₄, I₄, I₂, I₂, als Gitter A₃ ⊕ A₃ ⊕ A₁ ⊕ A₁ vom Rang 8. Das sitzt in D₅ ⊕ A₃, weil D₅ die Kette A₃ ⊕ A₁ ⊕ A₁ enthält. Die Breiten 4 + 4 + 2 + 2 geben 12.",
      "Schwierigkeit, die zugleich der Test ist: Die Marken wären von zwei Sorten. Die Vierteldrehung vertauscht die Sorten. Sie kann dann keine Symmetrie sein, sondern nur ein Übergang zur halbierten Kurve.",
      "Die Fläche ist seit Beauville (1982) bekannt, als eine von sechs mit genau vier Ausnahmefasern."],
    pre: "Diskriminante t⁴(t − 1)²(t + 1)², Doppelverhältnis −1, Punkt der Ordnung 4 bei x = t.",
    why: ["Ersatzfläche der ersten Runde hatte drei statt vier Stellen", "Markenlage braucht ein Element der Ordnung 4", "Kugel und quadratischer Torus sind exakt verbunden"],
    gives: ["A", "Z"],
    test: "Die Fläche durch einen Punkt der Ordnung 2 teilen und prüfen, ob sie dabei in sich übergeht und die zwei Markensorten tauscht. Dann die fünf Trägerfunktionen mit den Schnitten der Fläche vergleichen.",
    kill: "Wenn es diese Selbstabbildung nicht gibt. Dann widerspricht die Fläche der Vierteldrehung des Originals.",
    n: ["4", "12", "8", "2"], g: ["Marken", "Modul", "Singularität", "Auswahl", "Gitter"],
    l: [["viermarken", "testet"], ["taui", "testet"], ["d5", "testet"], ["a3", "testet"], ["hyp-milnor", "karte"], ["pat-modular", "karte"]] });

  H({ id: "hyp-genus", nr: 13, t: "Die Vier folgt daraus, dass nur der Torus eine Fouriertransformation hat",
    s: "Schließt die Lücke aus Hypothese 2: Stabilisatorzustände wählen das regelmäßige Vieleck, aber nicht die Zahl seiner Ecken.",
    b: ["Die Doppelüberlagerung der Kugel über n Marken hat Geschlecht (n − 2)/2. Nur vier Marken geben einen Torus.",
      "Nur ein Torus hat Verschiebungen, also Schieber und Uhr, also eine Fouriertransformation und einen Zustand kleinster Unschärfe. Sechs Marken gäben die Bolzakurve, zwölf eine Fläche vom Geschlecht fünf.",
      "Zusammen: Spinor (Doppelüberlagerung) plus Schieber und Uhr geben vier Marken. Stabilisatorzustand gibt das Quadrat. Von P2 blieben zwei Annahmen übrig, und keine nennt die Zahl vier."],
    pre: "Geschlecht 0, 1, 2, 3, 5 für 2, 4, 6, 8, 12 Marken.",
    why: ["Erste Runde: Auswahl der Lage exakt, Auswahl der Zahl fehlt", "Der Klebegenerator ist der Kommutator von Schieber und Uhr", "Laut Diagnose legen P1, P2 und E₈ die Markenlage nicht fest"],
    gives: ["A", "Z"],
    test: "Im Original nachsehen, ob Spinstruktur und Schieber und Uhr am Schnitt unabhängig von P2 vorkommen (Randprojektor, Klebegruppe).",
    kill: "Wenn das Original Schieber und Uhr erst aus den vier Marken gewinnt. Dann dreht sich das Argument im Kreis.",
    n: ["4", "1", "2", "5"], g: ["Marken", "Auswahl", "Fourier", "Qubit", "Ursprung"],
    l: [["p2", "testet"], ["viermarken", "testet"], ["hyp-ghz", "karte"], ["org-fourmark", "testet"], ["org-diagnose", "testet"]] });

  H({ id: "hyp-theta", nr: 14, t: "Die Raten 2/3 und 1/3 sind die Anteile der Spinsektoren am quadratischen Torus",
    s: "Die Gittersumme von E₈ zerfällt am Quadrat in drei Sektoren mit den Anteilen 1/6, 2/3 und 1/6.",
    b: ["Die E₈ Summe ist (θ₂⁸ + θ₃⁸ + θ₄⁸)/2. Am quadratischen Torus sind zwei der drei Teile gleich und der dritte viermal so groß. Die Anteile sind exakt 1/6, 2/3, 1/6.",
      "Die Raten des Attraktors sind 2/3 und 1/3. Die parallele Prüfung fand für E₈ ⊕ E₈ als gemischte Ladungsantwort genau 1/3. Der Randpunkt q = 2/9 der ersten Runde ist das Produkt beider Raten.",
      "Hypothese: Das ist dieselbe Zerlegung. Die Raten wären keine Parameter einer Bewegung, sondern die Gewichte der drei Randbedingungen, und nur am Quadrat rational.",
      "Vorsicht: Drei kleine Brüche treffen leicht zufällig zusammen."],
    pre: "Anteile 0,1667 · 0,6667 · 0,1667 auf 15 Stellen.",
    why: ["Raten 2/3 und 1/3 sind bisher gesetzt, nicht hergeleitet", "Parallele Prüfung: 1/3 unterscheidet E₈ ⊕ E₈ von D₁₆⁺", "q = 2/9 = (2/3)·(1/3) aus der ersten Runde"],
    gives: ["Z", "A", "P"],
    test: "Die Anteile am Sechseck und am Rechteck 2:1 ausrechnen und prüfen, ob die Vielfachheiten 6 und 6 im Spektrum {1, (2/3)⁶, (1/3)⁶} aus den Sektoren folgen.",
    kill: "Wenn die Vielfachheiten nicht folgen oder die Raten des Originals nicht von der Torusform abhängen dürfen.",
    n: ["2/3", "1/3", "1/6", "2/9", "8"], g: ["Attraktor", "Modul", "Zustand", "Gitter", "Fourier"],
    l: [["attraktor", "testet"], ["taui", "testet"], ["e8", "testet"], ["hyp-scharf", "karte"], ["dyn-markov", "testet"]] });

  H({ id: "hyp-phase", nr: 15, t: "Die Viertelphase der Naht ist der Abstand der beiden Seiten: (10 − 6)/16",
    s: "Die parallele Prüfung zeigt: Sechzehn gleiche Phasen ergeben eins, für eine Viertelphase braucht es eine Wirkung auf den geladenen Kanälen. Die liegt schon vor.",
    b: ["Der Träger D₅ steht auf Stufe 10 von sechzehn, die Familie A₃ auf Stufe 6. Die Summe 16 gibt die Phase eins, das ist E₈. Die Differenz 4 gibt exp(2πi · 4/16) = i.",
      "In der ersten Runde kam heraus: Die Klebegruppe wirkt als die Zahl i, und ihr Erzeuger ist der Kommutator von Schieber und Uhr. Hypothese: Das ist dasselbe i. Die Orientierung der Naht ist die Frage, welche Seite zuerst kommt.",
      "Dazu passt das Ordnungssignal der parallelen Prüfung: Es tritt bei allen vier Tripeln auf und verschwindet im Mittel ohne Marken. Ein Kommutator wechselt mit der Reihenfolge das Vorzeichen und mittelt sich weg."],
    pre: "exp(2πi·(10 − 6)/16) = i, exp(2πi·(10 + 6)/16) = 1.",
    why: ["Parallele Prüfung: einfache fermionische Restphase ausgeschlossen", "Parallele Prüfung: Ordnungssignal im Originalmodell v724 vorhanden", "Erste Runde: Klebegenerator = Kommutator, wirkt als i"],
    gives: ["T", "P"],
    test: "Die vier Ordnungssignale aus v724 mit der Umlaufrichtung des jeweiligen Tripels auf dem Markenkreis vergleichen. Vorhersage: gleicher Betrag, Vorzeichen folgt der Umlaufrichtung.",
    kill: "Wenn die vier Signale verschiedene Beträge haben oder das Vorzeichen nicht mit der Umlaufrichtung wechselt.",
    n: ["10", "6", "16", "4"], g: ["Phase", "Kitaev", "Zeit", "Verklebung", "Orientierung"],
    l: [["naht", "testet"], ["dyn-anomaly", "testet"], ["hyp-bell", "karte"], ["hyp-affin", "karte"], ["d5", "testet"], ["a3", "testet"]] });

  H({ id: "hyp-paarform", nr: 16, t: "Die sechs Paarübergänge geben, quadratisch gelesen, eine Zeit und drei Raumrichtungen",
    s: "Linear geht es nicht, das zeigt die parallele Prüfung (Rang höchstens drei). Die Summe über die sechs Paare ist aber von selbst eine Form mit Signatur (1,3).",
    b: ["Vier Alternativen mit Gewichten x₁ bis x₄. Die Form x₁x₂ + x₁x₃ + … + x₃x₄ hat einen positiven Eigenwert (alle vier gemeinsam) und drei negative (Unterschiede zwischen ihnen).",
      "Zeit wäre „alle vier zusammen“, Raum wären die Unterschiede. Der Lichtkegel ist die Bedingung, dass die Paarsumme null ist.",
      "Als antisymmetrische Paare geben die Sechs nur ein Feld vom Typ B_μν. Genau das fanden die Universalraum Runden. Eine Metrik braucht den symmetrischen Teil mit zehn Komponenten.",
      "Die parallele Prüfung ließ zwei Türen offen: eine antilineare Wirkung der Vierteluhr und quadratische Antworten mit D₅ Ladungen. Der Viererzyklus wirkt auf dieser Form als Drehspiegelung, also mit Spiegelung."],
    pre: "Eigenwerte 3 (einfach) und −1 (dreifach), Determinante −3.",
    why: ["Parallele Prüfung: lineare Abbildung der sechs Kanäle auf die Metrik hat höchstens Rang drei", "Universalraum Runden gaben nur B_μν, kein Spin 2", "Kein Ansatz liefert bisher 3 + 1 Dimensionen"],
    gives: ["L", "W"],
    test: "Die quadratischen Antworten des Originals auf den sechs Paaren auswerten: liegen sie mit gleichem Gewicht auf allen sechs, und erreichen sie mit D₅ Ladungen zehn unabhängige Komponenten?",
    kill: "Wenn die Gewichte der sechs Paare verschieden sind oder der Rang unter zehn bleibt.",
    n: ["6", "4", "3+1", "10", "3"], g: ["Raumzeit", "Lorentz", "Gravitation", "Familie", "Spin 2"],
    l: [["grav", "testet"], ["a3", "testet"], ["famop", "testet"], ["hyp-weyl", "karte"], ["dyn-lorentz", "testet"], ["ur-first", "testet"]] });

  H({ id: "hyp-cp", nr: 17, t: "Die messbare CP Verletzung ist der ungerade Teil unter Umkehr der Orientierung",
    s: "Die parallele Prüfung fand für die zwei Orientierungen verschiedene Beträge. Zerlegt man in ungerade und gerade, liegt der ungerade Teil beim Messwert.",
    b: ["Berichtet: J₊ = 3,3270 · 10⁻⁵ und J₋ = −2,7896 · 10⁻⁵. Eine reine CP Spiegelung müsste den Betrag erhalten und ist damit ausgeschlossen.",
      "Zerlegt: ungerader Teil 3,058 · 10⁻⁵, gerader Teil 0,269 · 10⁻⁵. Der Messwert liegt bei etwa 3,1 · 10⁻⁵ (vor einem Test nachschlagen).",
      "Erste Runde: Der feste Raum der modularen Wirkung ist zweidimensional, diagonal und gegendiagonal. Hypothese: Die zwei Orientierungen sind diese zwei Vektoren. Physikalisch zählt ihre Differenz, der gerade Teil ist eine Beimischung ohne CP Verletzung.",
      "Warnung: Die Zahlen sind berichtet, nicht hier nachgerechnet, und die Nähe zum Messwert fiel erst nachträglich auf."],
    pre: "(J₊ − J₋)/2 = 3,0583 · 10⁻⁵, (J₊ + J₋)/2 = 0,2687 · 10⁻⁵, Verhältnis 0,088.",
    why: ["Parallele Prüfung: Umkehr ändert den Betrag um 16 %", "Erste Runde: zwei feste Vektoren statt einem", "Die CP Phase ist eine der wenigen Größen mit freier Orientierung"],
    gives: ["W", "A"],
    test: "Vorher festlegen und dann rechnen: dieselbe Zerlegung für die CKM Phase und die Leptonphase. Der gerade Teil muss sich als reine Verschiebung der Mischungswinkel schreiben lassen.",
    kill: "Wenn der gerade Teil nicht wegtransformierbar ist oder die ungeraden Teile der anderen Größen neben den Messwerten liegen.",
    n: ["3,06", "16 %", "2"], g: ["Flavor", "CP", "Orientierung", "Verklebung"],
    l: [["flavor", "testet"], ["e8", "testet"], ["hyp-bell", "karte"], ["naht", "testet"]] });

  H({ id: "hyp-zelle", nr: 18, t: "TFPT ist ein Bruchteil einer OPH Zelle, keine Vereinigung von Zellen",
    s: "Die parallele Prüfung zeigt, dass Vergröbern von OPH nicht zu TFPT führt. Die Zahlen legen die umgekehrte Richtung nahe.",
    b: ["Berichtet: Der Kurzzeitquotient ist bei TFPT 1,125, bei OPH 7,6 · 7,9 · 7,975 für 240, 960 und 3840 Ports. Jede Vereinigung ganzer Zellen liegt bei mindestens 4.",
      "Die drei OPH Werte sind exakt 8 − 96/N. TFPT liegt mit 9/8 unter der Schranke 4. Es ist also kleiner als eine Zelle.",
      "240 ist die Zahl der Wurzeln von E₈. Hypothese: TFPT ist eine Zelle, geteilt durch ihre Symmetrie. Man fasst die Ports zusammen, die unter den Spiegelungen von D₅ und A₃ gleich sind. Das gibt vier Klassen mit 40, 12, 128 und 60 Ports."],
    pre: "8 − 96/N trifft alle drei berichteten Werte. 40 + 12 + 128 + 60 = 240.",
    why: ["Parallele Prüfung: Vergröberung nach ganzen Zellen scheitert", "OPH und TFPT teilen D₅ ⊕ A₃, aber nicht die Dynamik", "Erste Runde: Knotentafel mit 128 und 60 Wurzeln außerhalb des Teilgitters"],
    gives: ["L", "T"],
    test: "Die Zelle mit 240 Ports auf die vier Klassen zusammenfassen und denselben Quotienten rechnen. Vorhersage: 9/8.",
    kill: "Jeder andere Wert.",
    n: ["240", "9/8", "8", "96"], g: ["OPH", "Lokalität", "Gitter", "Dynamik"],
    l: [["ext-oph", "testet"], ["dyn-oph", "testet"], ["gap-oph", "testet"], ["e8", "testet"], ["hyp-kac", "karte"]] });

  H({ id: "hyp-pfeil", nr: 19, t: "Der Zeitpfeil ist ein Umlauf, der nie zurückkehrt",
    s: "In der ersten Runde fiel „Zeit als Umlauf“, weil die Umläufe unendliche Ordnung haben. Vielleicht ist genau das die richtige Eigenschaft.",
    b: ["Alle Uhren im Korpus (4, 5, 30) sind periodisch und haben deshalb keine Zählrichtung. Daran sind sie als Zeit gescheitert.",
      "Der Umlauf um eine Marke schert den Torus bei jedem Durchgang um vier weiter. Das ist eine Zahl, die nur wächst. Die Uhren wären die Reste dieser Zählung, die Zählung selbst wäre die Zeit.",
      "Die parallele Prüfung verlangt für einen Wechsel des Vakuumsektors bei endlicher Energie Anregungen, die auf großen Skalen beliebig weich werden. Ein stark gescherter Torus ist lang und dünn, seine tiefsten Anregungen werden weich."],
    pre: "Aus der ersten Runde: Umläufe T⁴, T und −T, alle von unendlicher Ordnung.",
    why: ["Erste Runde: kein Umlauf endlicher Ordnung", "Uhrengetriebe ohne Zählrichtung", "Parallele Prüfung: Sektorwechsel braucht weiche Moden, für P1 nicht hergeleitet"],
    gives: ["T"],
    test: "Am Torus mit Scherung 4n bei fester Fläche die tiefste Anregung ausrechnen und das Energiekriterium der parallelen Prüfung anwenden.",
    kill: "Wenn die Energielücke an der Marke nach unten beschränkt bleibt. Das flache Kissen des Originals hat Kegelspitzen und keinen dünnen Hals; das spricht dagegen.",
    n: ["4", "30", "∞"], g: ["Zeit", "Kosmologie", "Singularität", "Clock"],
    l: [["clocks", "testet"], ["kosmo", "testet"], ["hyp-milnor", "karte"], ["dyn-clockgear", "testet"], ["gap-g3", "karte"]] });

  H({ id: "hyp-norm", nr: 20, t: "Die absolute Normierung kommt aus dem Verkleben mit Drehung, nicht aus dem Zerschneiden",
    s: "Die parallele Prüfung zeigt: Zerlegen und Vergleichen lässt einen Faktor exp(−c · Fläche) frei. Verkleben mit Vierteldrehung könnte ihn festlegen.",
    b: ["Im Atlas steht: Am quadratischen Torus heben sich Schwingungsenergie −c₃ und Nullmode +c₃ genau auf, und nur dort.",
      "Die Fläche selbst ändert sich unter der Vierteldrehung nicht. Ein Flächenterm wäre also frei. Anders, wenn der Umfang feststeht: Dann ist Fläche gleich Umfang mal Länge, und der Term wirkt wie eine Verschiebung der Grundenergie. Die ist durch die Drehung festgelegt.",
      "Hypothese: Die Einheitswindung von P1 legt den Umfang fest, die Vierteldrehung dann c. Der Wert wäre c₃."],
    pre: "E₂ − 3/(π·Im τ) verschwindet bei τ = i.",
    why: ["Parallele Prüfung: reine Unterteilung legt die absolute Antwort nicht fest", "Balance −c₃ + c₃ = 0 nur am Quadrat", "h = 1 folgt aus einem Flussquant, falls die Gleichsetzung mit der Windung gilt"],
    gives: ["A", "Z"],
    test: "Die Familie exp(−c · Fläche) der parallelen Prüfung auf dem Torus mit festem Umfang auswerten und Verträglichkeit mit der Vierteldrehung verlangen. Vorhersage: genau ein c.",
    kill: "Wenn das Original den Umfang nicht festlegt. Dann bleibt c frei, wie die parallele Prüfung es fand.",
    n: ["1/(8π)", "2π", "1"], g: ["Alpha", "Normierung", "Modul", "Windung"],
    l: [["p1", "testet"], ["org-c3-modular", "testet"], ["alpha", "testet"], ["hyp-fluss", "karte"], ["gap-torus", "testet"]] });

  /* Parallele Prüfung in einer anderen Sitzung (andere zehn Hypothesen). Berichtet, hier nicht nachgerechnet. */
  T.parallel = [
    ["Raum aus vier Familienalternativen", "Zweikomponentenansatz scheitert an der Vierteluhr. Zweiter Sektor oder antilineare Wirkung bleibt möglich.", [16]],
    ["Absolute Normierung durch zwei Erzeugungswege", "Unterteilung allein reicht nicht, exp(−c · Fläche) bleibt frei.", [20]],
    ["Nahtorientierung aus fermionischer Restphase", "Einfache Fassung ausgeschlossen. Viertelphase braucht Wirkung auf geladenen Kanälen.", [15]],
    ["Nachbarschaft durch Wiederherstellbarkeit", "Der Fünfercode zeichnet keinen Nachbarn aus.", []],
    ["E₈ Orientierung als einzige CP Quelle", "Umkehr ändert den Betrag um 16 %. Keine reine CP Spiegelung.", [17]],
    ["Zeit aus der Reihenfolge von Zustandswechseln", "Ordnungssignal im Originalmodell v724 vorhanden, verschwindet im Mittel ohne Marken.", [15]],
    ["Gemeinsame Dynamik von OPH und TFPT", "Vergröberung nach ganzen Zellen scheitert: 1,125 gegen 7,6 bis 7,975.", [18]],
    ["Gitterauswahl durch markierte Amplituden", "E₈ ⊕ E₈ gibt 0 oder 1/3, D₁₆⁺ gibt 1/9 oder 7/9. Farbwurzel und schwache Wurzel geben 1/3.", [14]],
    ["Gravitation aus sechs Paarübergängen", "Lineare Abbildung hat höchstens Rang drei. Quadratische Antworten offen.", [16]],
    ["Kosmischer Übergang durch Vakuumwechsel", "Bei endlicher Energie nur mit weichen Moden, für P1 nicht hergeleitet.", [19]],
  ];

  /* ── Urteile der Prüfung HYPCHECKB20261119 ── */
  const C = "experiments/theory-contracts/landkarte-hypothesen-runde2-20261119/";
  const V = {
    "hyp-nachbar": { st: "aus", res: "Geometrie exakt, als Zeit gefallen", chk: "12/12", file: "check_h11_nachbar.py",
      found: ["Es gibt wirklich einen einzigen Kreis in E₈, der durch die Klebeerzeuger beider Schnitte läuft. Bei einem Viertel der Runde bleibt D₅ ⊕ A₃ stehen, bei einem Fünftel A₄ ⊕ A₄. Gemeinsamer Takt 20."],
      falls: ["Die Richtung ist eine schon bekannte Eichladung (die Ladung X von SO(10), Werte −5, −1, 3 auf dem Spinor).", "Auf den fünf Plätzen wirkt sie nur als Zahl. Die Vierteluhr des Originals trennt 3 von 2, der Fünfertakt verschiebt Plätze. Keine von beiden entsteht.", "Periodisch mit ganzzahligem Spektrum, also wieder eine Uhr ohne Zählrichtung."],
      opn: "" },
    "hyp-beauville": { st: "bedingt", res: "Trägt als Mathematik, Bindung an P1 offen", chk: "16/16", file: "check_h12_beauville.py",
      found: ["Die befürchtete Abbruchregel greift nicht: Teilt man die Kurve durch einen Punkt der Ordnung 2, erhält man dieselbe Familie an der gedrehten Stelle. Die Vierteldrehung der Marken ist eine Halbierung der Kurve.", "Die Halbdrehung ist eine echte Symmetrie. Das passt zur ersten Runde: Die Eichalgebra ist, was mit der Halbdrehung vertauscht.", "Die zwei Fixpunkte der Vierteldrehung, also die Pole neben den Marken, tragen genau den quadratischen Torus.", "Das Gitter A₃ ⊕ A₃ ⊕ A₁ ⊕ A₁ liegt in D₅ ⊕ A₃ und hat in E₈ den Quotienten Z₄ × Z₂, dieselbe Gruppe wie die Torsion der Fläche."],
      falls: ["Gleich sind die beiden Kurven erst nach Hinzunahme von i (Verdrehung mit −1)."],
      opn: "Dass P1 die Basis dieser Familie ist. Im Original sind die vier Marken gleichartig, hier gibt es zwei Sorten. Eine Abbildung von den fünf Trägerfunktionen auf die Fläche fand sich nicht." },
    "hyp-genus": { st: "belegt", res: "Richtig, steht aber schon im Original", chk: "11/11", file: "check_h13_genus.py",
      found: ["Das Original leitet die Vier schon so her: flache Kugel, lauter Halbdrehungen, Winkelsumme gibt vier Marken (Prüfung v216 läuft durch). Flach, Torus und Fouriertransformation sind dieselbe Bedingung.", "Schärfung: Es gibt genau vier flache Kugeln mit Kegelpunkten. Nur die mit lauter Halbdrehungen hat eine freie Markenlage. Ein Punkt der Ordnung 4 ließe keine Wahl."],
      falls: ["Keine neue Herleitung. Flachheit und Halbdrehung bleiben gesetzt.", "Schieber und Uhr kommen im Original aus der Klebegruppe, also von P2. Für diesen Teil dreht sich das Argument im Kreis."],
      opn: "" },
    "hyp-theta": { st: "aus", res: "Gefallen, eigene Abbruchregel greift", chk: "10/10", file: "check_h14_theta.py",
      found: ["Die Anteile 1/6, 2/3, 1/6 am Quadrat sind exakt.", "Die Gitterprobe der parallelen Prüfung ist nachgerechnet und stimmt: E₈ gibt immer 1/3, D₁₆ nur 1/9 oder 7/9. Der Grund ist einfach: 1/3 ist die Zahl jeder gleichmäßigen Verteilung."],
      falls: ["Die Vielfachheiten 6 und 6 folgen nicht aus drei Sektoren.", "Die Anteile hängen an der Torusform (am Rechteck 2:1 sind es 0,0004 · 0,515 · 0,485). Die Raten des Originals haben keinen Formparameter.", "Drei verschiedene Quellen für dieselben Brüche, keine Abbildung zwischen ihnen."],
      opn: "" },
    "hyp-phase": { st: "bedingt", res: "Kern exakt, Bindung an das Ordnungssignal gefallen", chk: "11/11", file: "check_h15_phase.py",
      found: ["Als Satz: Sind zwei Seiten verklebbar, dann ist das Verhältnis ihrer Phasen gleich dem Kommutator von Schieber und Uhr. Die zwei i der Hypothese sind zwingend dieselbe Zahl. Seitentausch gibt −i.", "Ein Ordnungssignal in v724 gibt es bei allen vier Tripeln."],
      falls: ["v724 ist reell. Der Teil, der von der Umlaufrichtung abhängt, ist bei allen vier Tripeln exakt null.", "Die vier Signale haben verschiedene Beträge (0,082 bis 0,157). Die Vorhersage trifft nicht."],
      opn: "Die Größe der parallelen Sitzung liegt nicht im Repository. Geprüft ist eine hier definierte Dreierspur." },
    "hyp-paarform": { st: "bedingt", res: "Form exakt, für die sechs fermionischen Kanäle gefallen", chk: "13/13", file: "check_h16_paarform.py",
      found: ["Die Paarsumme auf n Größen hat immer eine Zeit und n − 1 Raumrichtungen. Genau vier Größen geben 3 + 1.", "Die vier Alternativen selbst sind lichtartig. Ganzzahlig ist die Form eine hyperbolische Ebene plus ein Sechseckgitter.", "Der parallele Befund (linear höchstens Rang drei) ist nachgerechnet und stimmt."],
      falls: ["Mit voller Familiensymmetrie können Quadrate der sechs fermionischen Kanäle die Metrik gar nicht erreichen (Rang null).", "Diese sechs Kanäle tragen von sich aus eine Form mit Signatur (3,3)."],
      opn: "Ob das Original die sechs Paare gleich gewichtet. Aus einer Form allein folgt noch kein Ort." },
    "hyp-cp": { st: "aus", res: "Gefallen am zweiten Test", chk: "10/10", file: "check_h17_cp.py",
      found: ["Die berichteten Werte sind aus der Originalformel nachgerechnet und stimmen.", "Der ungerade Teil trifft J auf −0,16 σ (voll: +1,76 σ)."],
      falls: ["Die Phase selbst verfehlt mit der ungeraden Lesart den Messwert um 2,3 σ (voll: +1,0 σ). Zusammen passt es schlechter.", "Der gerade Teil ist messbar: Die zwei Orientierungen haben verschiedene |V_td|.", "Das Original erklärt den Überschuss bei J schon über die Skala von s₂₃."],
      opn: "" },
    "hyp-zelle": { st: "offen", res: "Nicht entscheidbar, Daten fehlen", chk: "10/10", file: "check_h18_zelle.py",
      found: ["Mit der Lesart Q = Zahl der Moden mal Σr²/(Σr)² geben die Raten 1/3, 1/3, 2/3 genau 9/8."],
      falls: ["Korrektur: Es sind fünf Klassen (40, 12, 64, 64, 60), nicht vier.", "Ersatzrechnung am Wurzelgraphen von E₈: Q = 1,014 statt 7,6, zusammengefasst 1,008 statt 9/8. Die OPH Zelle ist nicht der Wurzelgraph."],
      opn: "Das OPH Netz und die Definition des Quotienten liegen nicht im Repository." },
    "hyp-pfeil": { st: "aus", res: "Gefallen, eigene Abbruchregel greift", chk: "9/9", file: "check_h19_pfeil.py",
      found: ["Was bleibt: Uhr und Zählung unterscheiden sich genau wie endlich und unendlich viele Zustände. Ein Zeitpfeil braucht einen wachsenden Zustandsraum."],
      falls: ["Die Scherung ändert den Torus nicht. Mein Bild vom langen dünnen Torus war falsch.", "Auf jedem endlichen Zustandsraum hat der Umlauf endliche Ordnung.", "Das flache Kissen des Originals hat eine feste Lücke."],
      opn: "" },
    "hyp-norm": { st: "offen", res: "Halb: c ist festgelegt, aber auf null", chk: "9/9", file: "check_h20_norm.py",
      found: ["Mit festem Umfang lässt die Vierteldrehung genau ein c zu.", "Ohne absolute Länge hat ein Flächenkoeffizient gar keinen Platz."],
      falls: ["Der Wert ist null, nicht c₃. c₃ ist die Balance am Quadrat und steht so schon im Atlas."],
      opn: "Ob die absolute Antwort, die die parallele Prüfung meint, maßstabsfrei gemeint ist. Wenn nicht, bleibt c frei." },
  };
  T.nodes.forEach((n) => { const v = V[n.id]; if (v) { Object.assign(n, v); n.src = (n.src || []).concat([C + v.file, C + "README.md"]); } });
  T.parallelNote = "Nachgerechnet und bestätigt: Gitterprobe (0 oder 1/3 gegen 1/9 oder 7/9), Rang höchstens drei bei den sechs Paaren, die beiden J Werte. Nicht nachrechenbar: OPH Quotienten und das Energiekriterium.";
})();
