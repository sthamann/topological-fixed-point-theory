# Urteil: Weltgenerator und Realform (Exploration 27. November 2026)

Skripte: `w1_weltgenerator.py` (15/15 Prüfungen, etwa 15 s), `w2_bruecken.py` (24/24, etwa 3 s). Lane:
`research/originlab/whole_solution_search.py`, Abschnitt D, mit `test_whole_solution_search.py` (Klasse `SeedWorlds`,
9 Tests; zusammen mit den Nachbarsuiten 64 Tests grün). Ausgaben: `ausgabe_w1.txt` (sha256 19fefdcb…), `ausgabe_w2.txt`
(sha256 8a680257…), `w1_ergebnis.json`, `w2_ergebnis.json`. Nichts ist im Ledger, im Plan oder im Atlas befördert.

## Was gerechnet ist

**Weltgenerator (w1).** Ein Keim ist (Blätter k, Marken n, Träger N, Familien f, Dimension d), 16 464 Keime in den
Grenzen k ≤ 6, n ≤ 8, N ≤ 9, f ≤ 7, d ≤ 8. Jede Klausel wird aus dem Keim gerechnet, die Zahlen mit den
aufgezeichneten Gesetzformen und eingesetztem Keim.

| Frage | Ergebnis |
|---|---|
| Welche Keime schließen (flach, Träger und Familien aus dem Markendivisor, ℤ₄-Anomalie, Verklebung, eindeutige positive Teilung)? | nur (2, 4, 5, 3); die Dimension bleibt dabei frei |
| Was macht die ℤ₃-Welt (3, 3, 4, 2)? | scheitert an Anomalie, an beiden Verklebungen und an der Teilung |
| Regelsuche, 26 Prinzipien (7 Gegenprinzipien), 83 682 Regeln bis Länge 5 | 6692 Regeln wählen (2, 4, 5, 3, 4) eindeutig, die nächste Welt (d = 2) 3249; 399 minimale Regeln, 60 der Länge 3 |
| Zahlen über 336 Keime (k, N, f) | α trifft nur (2, 5, 3); genau ein Keim trifft alle vier Zahlen |
| Nachbarn von (2, 5, 3) | α⁻¹ = 34,18 (k = 1), 308,87 (k = 3), 156,11 (f = 2), 124,82 (f = 4), 119,63 (N = 6) |
| Verbünde aus 10 Lesarten, 386 Verbünde | keine vollständige Lösung; bester Abstand 7 |
| Rollen, die keine Lesart aus ihrem Grundobjekt liefert | Reflexionspositivität, Zeit, Kontinuum, Einsteindynamik, Flächennormierung, Restflavor |
| Einzelne Lesart, die den ganzen Keim erzwingt | nur „Spin 2 auf der Himmelskugel“ |
| Hebel | „ladder“ 3 Rollen, „equilibrium“ 2, „irreducible“ 1; alle drei zusammen senken den Abstand auf 1 (Restflavor) |
| Kann die Suche eine vollständige Lösung enthalten? | ja: eine eingepflanzte Lesart wird als vollständig erkannt, drei Gegenproben schlagen an |

**Brücken (w2).**

- B1: Index k², Markenzahl 2k/(k−1) und Liftordnung 2k sind nur für k = 2 gleich (alle 4). c₃(k) = 1/(4πk) ist nur
  für k = 2 die Einsteinkonstante 1/(8π).
- B2: Spannung (n−1)/(4n) durch Defizit 2π(n−1)/n ist für jedes n gleich 1/(8π). Für n = 2: 1/8 je Marke, vier Marken
  tragen 1/2 und schließen die Sphäre.
- B3: Zwistgewichte der fermionischen Replika; n = 2 gibt ±1/4, je 1/32, zusammen 1/16, Tausch der Ordnung 4.
  Einen Modul hat die Replikafläche nur für zwei Blätter und zwei Gebiete.
- B4: Affines E₈ als McKay-Graph: Der Schnitt am Knoten Sym³ (Marke 4) gibt A₃ mit den Dimensionen 1, 2, 3 und D₅
  mit 5, 6, 4, 2, 3. Kein anderer Knoten gibt D₅ + A₃.
- B5 (neu): Unter den 256 Involutionen exp(iπx) sind auf D₅ genau drei Typen trivial: kompakte E₈ mit su(4), E₈(8)
  mit kompakter su(4), E₈(−24) mit su(2,2). Kompakte Spin(10) mit nichtkompaktem A₃-Faktor gibt es also nur in
  E₈(−24), und der Faktor ist so(4,2). so(5,1), so(3,3) und su(3,1) kommen nicht vor. Nichtkompakt sind
  112 = 8 + 40 + 32 + 32; maximal kompakt ist E₇ × SU(2), die einzelne SU(2) liegt im A₃-Faktor. Die drei Involutionen
  vom Typ E₈(−24) sind die drei Paarungen der vier Gewichte der 4; der Stabilisator einer Paarung in S₄ hat Ordnung 8.
- B6: Schleife Dimension → Blätter → Marken → Familien → Dimension mit dem einzigen Fixpunkt d = 4.
- B7: Genau eine positive Teilung hat nur N = 5 (bis N = 40). D_N ⊕ A_f verklebt für N ≤ 16, f ≤ 12 nur bei (5, 3)
  und (13, 3).
- B8: Zahlen der Nachbarwelten wie in w1.

## Urteil

| Gegenstand | Urteil | Stufe |
|---|---|---|
| Weltgenerator als Suchraum für ganze Lösungen | brauchbar; eine vollständige Lösung ist darstellbar und wird erkannt | Werkzeug |
| „Dieselbe Vier“ als Eigenschaft von zwei Blättern (B1) | exakt | 1 (ein Stück, kein Fit) |
| Der Keim (2, 4, 5, 3) ist der einzige, der schließt | exakt in den Grenzen | 1 |
| Realform-Lesart: E₈(−24) ⊃ Spin(10) × SU(2,2) erzwingt 3+1 | Kandidat | 1, mit der Brücke zu T′ ein Stück von Stufe 4 |
| Geschlossene 2+1-Welt als Ursprung von c₃ | offen; die Zahl stimmt für jedes n, also keine Auswahl | keine |
| Drei Familien als drei Paarungen (drei Realformen) | Beobachtung, nicht geprüft | keine |
| Vollständige Lösung | keine; Abstand 7, mit drei benannten Aussagen 1 | – |

**Freie Parameter.** Keine Zahl ist angepasst. Gesetzt sind die Grenzen des Universums, der Klauselpool, die
verallgemeinerten Gesetzformen für fremde Keime (im Kopf von Abschnitt D genannt) und die Klassen der aus dem Bestand
gelesenen Rollen.

**Basislinie.** Sieben Gegenprinzipien in der Regelsuche; der Schluss ist über alle Keime gerechnet und hängt nicht
am Pool; die Zahlen sind an allen 336 Keimen gemessen. Die Zahl der Regeln, die eine Welt wählen, ist kein Beleg: 19
der 26 Prinzipien sind solche, die die aufgezeichnete Welt erfüllt. Belastbar sind der Schluss, die Ablationen der
minimalen Regeln und die Zahlen der Nachbarn.

**Stärkster Einwand.**

1. Gegen die Realform: In jeder Realform von E₈ ist die Materie bezüglich SL(2,ℂ) × Eichgruppe nicht chiral
   (Distler–Garibaldi 2009); (16,4) enthält 16 ⊗ (2_L ⊕ 2_R). Die Lesart trägt nur, wenn die Helizität die
   Leitergraduierung auf dem Twistorraum ist und die Orientierung der Naht eine Uhrklasse auswählt (Leitlinie
   „Chiralität aus Indizes“). Gerechnet ist bisher nur, dass die Realform die beiden Hälften trennt: 32 Wurzeln von
   (16,4) sind kompakt, 32 nichtkompakt.
2. Gegen die Brücke: Im Bestand (z79) verträgt sich die Markensymmetrie auf den Nahtlagen nur mit der euklidischen
   Realform. Mit kompakter Spin(10) gibt es die euklidische konforme Algebra so(5,1) in E₈ nicht (B5). Entweder ist
   die Realform der Nahtlagen nicht die Realform von A₃, dann fällt die Brücke, oder die euklidische Seite hat eine
   nichtkompakte Eichalgebra. Das ist nicht gerechnet.
3. Gegen den Rollenvertrag: 14 der 24 Rollen sind aus dem Bestand gelesene Status mit Quelle; die Suche beweist sie
   nicht. Die Leitergatter rechnen Positivität, Zeit und Kontinuum nur für das Leiterobjekt.
4. Gegen den Schluss: Träger N = n + 1 und Familien f = n − 1 sind als Sätze über den Markendivisor eingesetzt. Die
   physische Lesart (Eichwirkung, Familien als chiraler Index) bleibt eine eigene Frage.

**Erste frei gewählte Struktur.** Dass die physische Realform des geschlossenen Gitters nicht die kompakte ist. Die
Rechnung zeigt nur: Wenn sie nichtkompakt ist und die Eichgruppe kompakt bleibt, dann ist sie E₈(−24) mit SU(2,2).

## Abgleich mit Atlas und Bestand

- Reinheit zweier Bögen, Index-KMS, G₃₁, Gauß-Gitter, Twistorgerade, Himmelskugel: vorhanden; hier nur als Lesarten
  mit Quelle geführt.
- SU(2,2) entsteht im Bestand aus dem Nahtdublett mit Tomita-Partner (VOGELPERSPEKTIVE, Punkt 3); dass es die
  Realform des A₃-Faktors von E₈ ist, steht dort nicht. E₈(−24) kommt im Repositorium nicht vor.
- Der Kern ν = 1 und die SU(2,2)-Leiter sind im Bestand exakt unter T′. Der Hebel „ladder“ ist also T′, gelesen auf
  der Realform.
- Die sechs Rollen, die keine Lesart liefert, sind Abschnitt 5 von TFPT-STAND. Der Rollenvertrag aus den Grundlagen
  gibt denselben Rand wie die H6-Buchführung, aber mit drei Aussagen statt acht offenen Posten: T′ auf der Realform,
  Verschränkungsgleichgewicht (D), Irreduzibilität oder Verdrillung (G). C″ und E erscheinen nicht, weil Familien
  und Teilung auf der Strukturebene aus dem Keim folgen.

## Nächster entscheidender Test

Wirkt die aufgezeichnete Tomita-Konjugation des Dublettzustands (τ = (Ad J)⁺, VOGELPERSPEKTIVE Punkt 3) auf der 4 des
A₃-Faktors als eine der drei Involutionen aus B5, also als Paarung zwei gegen zwei, die mit Spin(10) vertauscht? Wenn
ja, ist T′ die Aussage „die Realform ist E₈(−24)“, und die Leitergatter gelten für das Original. Wenn nein, fällt die
Brücke zwischen Realform und Naht, und der Realformsatz bleibt ein Satz über E₈ ohne Bezug zum Original.

---

# Zweiter Durchgang (27. November): Bestand gelesen, Twistorbrücke, ℤ₄-Regel, Doppelbrechung

Skripte: `w3_twistorbruecke_chiralitaet.py` (34/34, 1 s), `w4_formsuche.py` (9/9, 3 s). Lane: Abschnitt E von
`research/originlab/whole_solution_search.py`, `research/originlab/ledger_map.py`, Tests `ShapeSearch` (4) und
`test_ledger_map.py` (4). Ausgaben: `ausgabe_w3.txt` (sha256 bcb25514…), `ausgabe_w4.txt` (sha256 03527353…).
w1 und w2 laufen unverändert (15/15, 24/24, gleiche Prüfsummen).

## Korrektur am ersten Durchgang

Der „nächste entscheidende Test“ oben war im Ledger bereits entschieden. WOIT.BETA3.DUALITY.01 (v565) zeigt: Uhr und
Woit-Struktur erzwingen auf dem Twistorraum der Naht eine hermitesche Form der Signatur (2,2). Der erste Durchgang
hatte nur `TFPT-STAND.md` gelesen, und der nennt 12 der 784 aktiven Ledger-Behauptungen. Richtig bleibt: Die Realform
E₈(−24) und die Gleichsetzung des Familienfaktors mit der Twistoralgebra stehen nirgends im Bestand. Die Ursache ist
behoben (`docs/research/LEDGER-KARTE.md`, Hinweis in `AGENTS.md`, `EXPLORATION.md` und `TFPT-STAND.md`).

## Was gerechnet ist

| Teil | Ergebnis |
|---|---|
| A, Paare von Realformen (diagonale Involutionen von O(16), Spinor über Cl(16)) | so(4,2) mit so(10), so(8,2) oder so(6,4); so(5,1) nur mit so(7,3); so(3,3) nur mit so(9,1) oder so(5,5); so(14,2) und so(10,6) kommen nicht vor |
| B, Blocksatz unabhängig nachgerechnet | Uhr: 16 → 8 Parameter; Minkowski-Zweig immer (2,2), Nullort auf der Nahtgerade ist der Nahtkreis; die diagonalen Mitglieder tragen die Paarungen {13\|24} und {14\|23}; die dritte Paarung {12\|34} liegt im euklidischen Zweig |
| C, Uhr und konforme Zeit | g = Vierteldrehung der konformen Zeit hat Ordnung 4 auf su(4), 8 auf E₈; g² = Cartan-Involution mal Klebeuhr; g⁴ = (−1)^F; E₈ enthält keine U(1), die mit D₅ ⊕ A₃ vertauscht |
| D, ℤ₄-Regel | Zentrum von Spin(10) = Zentrum von SU(2,2) in E₈(−24); das Zentrum wirkt auf der Helizitätsleiter h = C/2 als i^(C+2); erlaubt sind Klasse 0 mit h = ±1, Klasse 2 mit h = 0, die 16 mit genau einer Helizität, die 16bar mit der anderen; Skalare nur in 10, 120, 126 |
| E, Doppelbrechung | β = 0,24243503°; β/2 = c₃φ₀; Ω_b + β = φ₀; ACT und Planck PR4 zusammen 0,2415 ± 0,0614° (3,9σ von null, 0,02σ von TFPT); Frequenzexponent −0,35 ± 0,48, Faraday 3,4σ entfernt |
| Formsuche (w4) | 26 Rollen, 11 Lesarten, keine vollständige Lösung, Abstand 9; fünf Aussagen senken ihn auf 1 |

## Urteil

| Gegenstand | Urteil | Stufe |
|---|---|---|
| Brücke Realform ↔ Twistorbrücke: su(2,2) ist erzwungen, zwei der drei Paarungen sind die Minkowski-Mitglieder | exakt; die Gleichsetzung A₃ = Twistoralgebra bleibt Prämisse | 1 |
| z79 gegen kompakte Spin(10) | geklärt: beides zugleich geht in E₈ nicht; die euklidische Seite hätte so(7,3) | exakt |
| ℤ₄-Regel: eine chirale Generation aus der globalen Gestalt (Spin(10) × SU(2,2))/ℤ₄ | Kandidat; exakt als Satz über die Gruppe, bedingt als Aussage über Felder | 2 (Chiralität ohne neue Wahl außer einem Orientierungsbit) |
| Uhr μ₄ = Zentrum von SU(2,2) = Zentrum von Spin(10); halbe konforme Periode = Cartan-Involution mal Uhr | exakt | 1 |
| Doppelbrechung als Formbedingung | drei Bedingungen, siehe unten | Bedingung, kein Kandidat |
| Vollständige Lösung | keine; Abstand 9, mit fünf benannten Aussagen 1 | – |

## Was die Doppelbrechung über die Gestalt der Lösung sagt

1. Eine Drehung der Polarisation um β ist eine Phase 2β zwischen den Helizitäten ±1. Das ist ein Element der
   Helizitäts-U(1). E₈ enthält diese U(1) nicht, nur ihr ℤ₄, und das ℤ₄ wirkt auf Photonen trivial. β ist deshalb
   keine Größe der Kinematik, sondern eine Antwort des Zustands. Jede Lösung, deren Weltzustand das
   SU(2,2)-invariante Vakuum ist, gibt β = 0 und ist mit 3,9σ in Spannung.
2. β ist ungerade unter Parität. Das Vorzeichen ist dasselbe Bit wie die Orientierung in der ℤ₄-Regel, die
   festlegt, welche Helizität die 16 trägt. Eine vollständige Lösung muss das Vorzeichen von β an die
   Linkshändigkeit der schwachen Wechselwirkung binden. Das ist eine prüfbare Vorhersage der Form, die im Bestand
   fehlt: bisher ist nur der Betrag eingefroren.
3. β ist achromatisch und eine reine Zahl. Das passt zu einer festen Phase des Randes und nicht zu einem rollenden
   Feld, dessen Winkel von der Geschichte abhinge. In der Twistornormierung ist die Phase c₃φ₀.

Daraus die Spannung, die die Suche jetzt führt: Die kontinuierliche Zeit steht im Bestand nur im Vakuumsektor (Z72),
die Doppelbrechung verlangt einen Zustand, der nicht das invariante Vakuum ist. Eine vollständige Lösung muss beides
in einem Objekt zeigen: Die Darstellung ist die Leiter, der Weltzustand ist nicht ihr Vakuum.

## Stärkste Einwände nach dem zweiten Durchgang

1. **Familienfaktor = Raumzeit ist eine Umdeutung.** Im Ledger ist E₈ ganz intern (selbstduales Yang–Mills mit E₈
   auf PT/Γ), die Familien sind H₁(P¹∖μ₄) = ℤ³. Wenn A₃ die Twistoralgebra ist, sind die drei Familienrichtungen
   die Cartan-Richtungen der konformen Algebra. Ob die Flavoroperatoren (det 3, 4, 8, 20 auf ℤ³) das tragen, ist
   nicht gerechnet. Fällt das, fällt die Realform-Lesart.
2. **Die ℤ₄-Regel braucht eine Prämisse.** Felder müssen Darstellungen des Quotienten sein, masselose Felder
   Helizitätsleitern. Der Einwand von Distler und Garibaldi (keine chirale Generation in der Zerlegung der 248 unter
   SL(2,ℂ)) trifft diese Lesart nicht, weil die Chiralität aus dem Zentrum kommt; er zeigt aber, dass (16,4) keine
   Felder sein können.
3. **Zwei Minkowski-Mitglieder.** Für a·c > 0 ist die Uhr die Vierteldrehung der konformen Zeit, für a·c < 0 eine
   räumliche Vierteldrehung. Welches Mitglied vorliegt, ist ein freies Bit; WOIT.BETA1.GSO.01 führt genau diese
   Frage als unentschieden.
4. **Die Klassen der gelesenen Rollen bleiben Status aus dem Bestand.** Sechzehn der 26 Rollen rechnet die Suche
   nicht.

**Freie Parameter.** Keine Zahl ist angepasst. Gesetzt sind die Klassen der zwei neuen Rollen je Lesart und die
Texte der Aussagen `family_is_twistor` und `oriented_state`.

**Erste frei gewählte Struktur.** Dass der Familienfaktor A₃ die sl(4) des Twistorraums der Naht ist.

## Nächster entscheidender Test

Trägt die Lesart „Familienrichtungen = Cartan-Richtungen von su(2,2)“ die Flavoroperatoren? Zu prüfen: Wirken Q, K,
R, L aus dem Ledger (FLAV, Determinanten 3, 4, 8, 20, Produkt 1920 = |W(D₅)|) auf dem Wurzelgitter von A₃ so, dass
die Uhr die Vierteldrehung der konformen Zeit ist und die Paarung {13|24} erhalten bleibt? Wenn ja, ist
`family_is_twistor` ein Satz über vorhandene Objekte und die Chiralität folgt ohne neue Wahl. Wenn nein, ist die
Realform-Lesart verworfen, und die ℤ₄-Regel bleibt ein Satz über E₈(−24) ohne Bezug zur Theorie.

## Abgleich mit dem Bestand (zweiter Durchgang)

- Im Bestand, hier nur bestätigt: Signatur (2,2) aus Uhr und Woit-Struktur (WOIT.BETA3.DUALITY.01), die beiden
  Familien antilinearer Strukturen (WOIT.THETA.FREE.01), β = φ₀/(4π) und Ω_b + β = φ₀ (eingefroren).
- Im Bestand mit anderem Weg: Der Befund H6h48 vom 25. November gibt über die reelle Clifford-Struktur Cl(3,11)
  S⁺ = (S_L ⊗ 16) ⊕ (S_R ⊗ 16bar) und die Uhr als Zentrum von Spin(10). Das ist derselbe Inhalt wie die ℤ₄-Regel
  für Fermionen. Neu ist hier der Weg über das gemeinsame Zentrum in E₈(−24) und die Aussagen für Bosonen und
  Skalare. so(11,3) liegt in so(12,4) und damit in E₈(−24); beide Befunde vertragen sich (Signaturarithmetik, nicht
  weiter gerechnet).
- Nicht im Bestand: E₈(−24), die Paartabelle der Realformen, die Gleichsetzung A₃ = Twistoralgebra, die Lesart von β
  als Phase der Helizitäts-U(1) außerhalb von E₈ und die Bindung des Vorzeichens von β an die Chiralität.

---

# Dritter Durchgang (27. November): Flavortest der Realform-Lesart

Skript: `w5_flavortest_familie_twistor.py` (50/50, 1 s), Ausgabe `ausgabe_w5.txt` (sha256 e7c9d9cc…),
`w5_ergebnis.json`. Lane: `family_clock_gates()` in `research/originlab/whole_solution_search.py`, Test
`ShapeSearch.test_family_clock_gates`. w1 bis w4 laufen unverändert (gleiche Prüfsummen).

## Was geprüft wurde

Der zweite Durchgang hatte den Test so gestellt: Tragen die Flavoroperatoren des Ledgers die Lesart
„Familienrichtungen = Cartan-Richtungen von su(2,2)“, mit der Uhr als Vierteldrehung der konformen Zeit und
erhaltener Paarung {13|24}? Wenn nein, sei die Realform-Lesart verworfen.

Eingaben, alle aus dem Bestand: die Drehung der vier Marken auf H₁(P¹∖μ₄) (GATE.QGEOM.02, QGEO.STATE.01), die
Transportmatrizen Q, R, K = R + QΣ, L = R + Q(1 + Σ) mit Σ = diag(1, −1, −1) (FLAV.KQ.01), das ganzzahlige Deck
G = T_A Σ (v98, v141, v146), die Formen ω_k = z^(k−1)dz/(z⁴−1) mit Charakter i^k (FLAV.QGEO.04, FLAV.QGEO.07), die
Twistoruhr ρ₄ = diag(i, 1, i, 1) (WOIT.BETA3.DUALITY.01), das Wörterbuch 4 = H⁰(O(1)) ⊗ ℂ², Familien = H⁰(K(D))
(Befund h6h33). Dass H₁ = ℤ³ das Gitter des A₃-Faktors von E₈ ist, steht in SEL.THM.01.

## Was gerechnet ist

| Teil | Ergebnis |
|---|---|
| A, Familienuhr | Die Drehung der Marken ist ein Coxeter-Element von W(A₃) = S₄, Eigenwerte i^m mit m = 1, 2, 3 |
| B, Familienraum als Cartan von sl(4), Uhr als Weyl-Element | Ad-Spektren (Vielfachheit von 1, i, −1, −i): ρ₄ (7, 4, 0, 4), Coxeter-Klasse (3, 4, 4, 4). Nicht konjugiert, auch nicht mit dem äußeren Automorphismus. Die Quadrate sind dasselbe Element θ = diag(−1, 1, −1, 1), Typ (2,2). Coxeter-Vertreter Z = ρ₄ · J mit J = diag(1, 1, −1, −1); θ, J, θJ sind die drei Paarungen. Auf dem Vektor 6 von SO(4,2): ρ₄ dreht die zeitartige Ebene um π/2, Z zusätzlich eine raumartige Ebene um π |
| C, Realform aus Uhr² | In W(A₃) haben nur Vierer-Zykel ein Quadrat, das eine nichttriviale Involution ist; Signatur (1, −1, −1). Unter den fünf Realformen von A₃ hat nur su(2,2) diese Cartan-Involution auf der maximal zerfallenden Cartan. Mit trivialer Wirkung auf D₅: Cartan-Involution von E₈(−24), 45 + 7 + 20 + 64 = 136. Zentralisator von Uhr² ist D₄ = W(C₂); eingeschränkte Wurzeln C₂ mit Vielfachheiten 1 und 2, reeller Rang 2 |
| D, ganzzahlige D₄ des Ledgers | T = (−1)^(Q₊) = [[−1,0,0],[0,1,0],[0,4,−1]] ist ganzzahlig. {1, Σ, T, ΣT} ist eine Vierergruppe, über ℤ diagonal in der Basis e₁, (0,1,2), e₃ (Q₊ = 3, 2, 1). Das Deck des Ledgers hat G² = ΣT; seine selbstkonjugierte Gerade ist ℤe₃ mit Q₊ = 1. Eine zweite ganzzahlige Uhr r mit r² = T existiert (eindeutig bis auf r⁻¹) |
| D, Q₊ gegen den Charaktergrad | Im Charakterraum kippt jede Spiegelung von D₄ den Grad: s N s⁻¹ = 4 − N. Σ ist eine Spiegelung und vertauscht mit Q₊. In beiden ganzzahligen D₄ vertauscht genau die Vierergruppe mit Q₊, die Uhr nicht |
| D, Gittertyp | B₁ ⊕ E hat im Transportgitter Index 1 (für G und für r), im Punktierungsgitter A₃* (v70) Index 2, im Wurzelgitter A₃ Index 2. G ist über ℚ zur Uhr von v70 konjugiert, über ℤ nicht. Das Gitter zwischen A₃ und A₃* (Gewichte von Λ²ℂ⁴) hat Index 1 |
| D, Q, K, R, L | Einfache Eigenwerte; rationale Kommutanten ℚ × ℚ(√13), ℚ × ℚ(√3), zwei kubische Körper; keiner enthält i. Keiner der vier vertauscht mit G, r, Σ oder T. R hat ein komplexes Eigenwertpaar (Diskriminante −7996) |
| E, Familien im Wörterbuch des Bestands | ω_k = v_k eingesetzt in q, mit v_k = z^(k−1)∂_z (Möbius-Algebra der Nahtgeraden) und q = dz²/(z⁴−1); q ist unter z → iz und z → 1/z ungerade. In ℂ² ⊗ ℂ² ist ρ₄ = 1 ⊗ diag(i, 1); auf 1 ⊗ sl(2) wirkt Ad(ρ₄) mit (1, i, −i); die Familienuhr ist −Ad(ρ₄). Mit der Form der Signatur (2,2) ist diese sl(2) die su(1,1) = so(2,1): eine kompakte Richtung (L₀ = ω₂, erzeugt die Uhr), zwei nichtkompakte |

## Urteil

Nach dem Wortlaut des Tests aus dem zweiten Durchgang ist die Antwort **nein**:

1. Als Cartan von sl(4) gelesen ist die Familienuhr nicht die Vierteldrehung der konformen Zeit. Sie liegt in der
   Coxeter-Klasse, ρ₄ ist zentral in der maximal kompakten Untergruppe.
2. Q, K, R, L erhalten weder eine Uhr noch die Paarung.

Die Aussage `family_is_twistor` in der Fassung von w4, „die drei Familienrichtungen sind die Cartan-Richtungen von
su(2,2) mit der Uhr ρ₄“, ist damit verworfen.

| Gegenstand | Urteil | Stufe |
|---|---|---|
| Familienrichtungen = Cartan von su(2,2) mit Uhr ρ₄ | widerlegt, exakt | – |
| Im Wörterbuch des Bestands sind die Familien die so(2,1) der Nahtgeraden in so(4,2); Familienuhr = −Ad(ρ₄); die mittlere Familie ist die Richtung der konformen Zeit | exakt als Rechnung; setzt h6h33 und FLAV.QGEO.04 zusammen | 1 (Umformulierung) |
| Quadrat der Familienuhr hat die Signatur (1, −1, −1): su(2,2) als einzige Realform von A₃, E₈(−24) mit kompakter Spin(10) | exakt; folgt allein aus der Ordnung 4 der Uhr | 1 |
| Q₊ ist in der ganzzahligen D₄ nicht der Charaktergrad | exakt; schärft den offenen Rest von GATE.QGEO (Modulidentifikation) | Beobachtung zur Prüfung |
| Deck G und Homologie-Uhr sind über ℤ nicht konjugiert; Typ des Zwischengitters | exakt | Beobachtung, ohne Folgerung |
| Q, K, R, L aus der Lesart | nicht erklärt | offen (`flavour_rest`) |
| ℤ₄-Regel für Chiralität | unverändert: exakt als Satz über die Gruppe, bedingt als Aussage über Felder | 2 |

**Warum die Realform-Lesart nicht ganz fällt.** Der zweite Punkt des Tests war falsch gestellt: Auch das Deck des
Ledgers vertauscht nicht mit Q. Er entscheidet deshalb nicht zwischen Lesarten. Der erste Punkt trifft die
Gleichsetzung „Familienraum = Cartan“. Der Bestand selbst liest die Familien anders (H⁰(K(D)), Möbius-Algebra), und
in dieser Lesart stimmt die Uhr: Sie ist die Vierteldrehung der konformen Zeit, bis auf das Vorzeichen von q.

**Was das wert ist.** Wenig mehr als vorher. Die Signatur (1, −1, −1) ist keine unabhängige Bestätigung: jedes
Element der Ordnung 4 in W(A₃) liefert sie (Basislinie: die einzige der fünf Klassen von S₄ mit Ordnung 4). Teil E
setzt zwei Aussagen des Bestands zusammen und rechnet das Vorzeichen aus; neu ist daran die Einbettung in su(2,2)
mit einer kompakten und zwei nichtkompakten Richtungen. Eine Korrektur nach einem gescheiterten Test ist schwächer
als ein bestandener Test.

**Die Spannung, die bleibt.** Der Bestand liest den dreidimensionalen Familienraum auf zwei Arten: als
Cartan-Gitter von A₃ mit der Uhr als Coxeter-Element (SEL.THM.01, QGEO.STATE.01) und als Möbius-Algebra der
Nahtgeraden mit der Uhr −Ad(ρ₄) (h6h33, FLAV.QGEO.04). Als Teilräume von sl(4) sind das verschiedene Dinge: der
erste ist abelsch, der zweite nicht, und die Uhren liegen in verschiedenen Klassen. Für E₈(−24) muss der A₃-Faktor
von E₈ die konforme Algebra sein; das verlangt die erste Lesart für das Gitter und die zweite für die Uhr.

**Beobachtung zu GATE.QGEO.** FLAV.QGEO.07 nennt nach dem Abgleich der vollen D₄ „keine prüfbare Abweichung“. Auf
der Ebene der Operatoren gibt es eine: Der Charaktergrad N vertauscht mit der Uhr und wird von jeder Spiegelung
gekippt; Q₊ vertauscht mit der Spiegelung Σ und nicht mit dem Deck G. Im ganzzahligen Modell trägt die
selbstkonjugierte Gerade Q₊ = 1, in der Kohomologie den Grad 2. Die Identifikation von Q₊ mit dem Charaktergrad
gilt also für das Spektrum, nicht für den Operator. Das ist keine Korrektur am Ledger (er führt die kanonische
Abbildung als offenen Rest), aber eine Bedingung an jede solche Abbildung. Nicht befördert, zur Prüfung.

**Freie Parameter.** Keine Zahl ist angepasst. Gewählt ist der Vertreter Z = diag(i, 1, −i, −1) der Coxeter-Klasse,
der mit ρ₄ vertauscht; „Z = ρ₄ · J“ gilt nur für diesen. Die hermitesche Form in Teil E ist die diagonale mit den
Blöcken {e₁, e₃} und {e₂, e₄} des Ledgers.

**Stärkste Einwände.**

1. „Cartan-Involution = Uhr²“ ist die Regel der Lesart, kein Satz des Bestands.
2. Teil E benutzt, dass der Twistorraum von WOIT.BETA3.DUALITY.01 das ℂ² ⊗ ℂ² von h6h33 ist und dass ρ₄ auf dem
   Faktor H⁰(O(1)) wirkt. Die Eigenwerte passen (i, 1, i, 1), geprüft ist die Gleichsetzung der Räume nicht.
3. Die Lesart liefert keine der Zahlen 3, 4, 8, 20 und keinen der Operatoren Q, K, R, L.

## Nächster entscheidender Test

Welche der beiden Lesarten des Familienraums ist der A₃-Faktor von E₈? Zu prüfen an SEL.THM.01: Die Verklebung
D₅ + A₃ zu E₈ braucht das Gitter A₃ mit Diskriminantengruppe ℤ₄. Wenn die Möbius-Algebra der Nahtgeraden mit dem
Vorzeichencharakter von q ein Gitter dieses Typs trägt (Perioden der ω_k über die Zykel um die Marken, mit der Form
der Verklebung), dann sind beide Lesarten ein Objekt und die Realform ist an die Twistorbrücke gebunden. Der
Befund D spricht dagegen, dass es das Punktierungsgitter selbst ist: Das Deck des Ledgers hat den Typ des
Zwischengitters. Wenn kein solches Gitter existiert, bleibt su(2,2) eine Aussage über die Familienuhr ohne Brücke
zu E₈, und die ℤ₄-Regel verliert ihre Begründung aus dem Bestand.

## Abgleich mit dem Bestand (dritter Durchgang)

- Im Bestand, hier nur benutzt: Uhr als Coxeter-Element von W(A₃) (QGEO.STATE.01), Matrix der Uhr auf der
  Homologie (GATE.QGEOM.02), Spec(Q₊) = Exponenten von A₃ (FLAV.KQ.01, FLAV.QGEO.02), Charaktere und Möbius-D₄ auf
  den ω_k, ganzzahliges Modell ⟨G, T_A, Σ⟩ und die Klasse von Σ (FLAV.QGEO.04, FLAV.QGEO.07, v98, v146),
  Exponentendualität (GATE.UWALL.13), H₁ → A₃ (SEL.THM.01), Klassengrößen 52, 64, 60, 64 (GATE.METRIC.09),
  Familien = H⁰(K(D)) = holomorphe Vektorfelder der Naht und 4 = H⁰(O(1)) ⊗ ℂ² (h6h33).
- Erst nach der ersten Fassung dieses Durchgangs gefunden: das ganzzahlige D₄-Modell des Ledgers. Die erste
  Fassung hatte eine eigene ganzzahlige Uhr als „die“ Uhr des Transportgitters ausgegeben; sie ist die zweite Wahl
  r² = T neben dem Deck G² = ΣT des Ledgers und steht jetzt so im Skript.
- Nicht im Bestand gefunden (Ledger-Karte nach SU(2,2), SO(4,2), Coxeter element, W(A3), Lambda^2, exponent duality,
  sheet parity): der Vergleich der Uhrklassen in sl(4), die Realform aus dem Quadrat der Familienuhr, D₄ als
  eingeschränkte Weyl-Gruppe, die Einbettung der Familien als so(2,1) in su(2,2), der Unterschied von Q₊ und
  Charaktergrad als Operatoren, der ganzzahlige Typ des Decks.
