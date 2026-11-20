# Urteil

Gerechnet am 20. November 2026 mit `.venv/bin/python` (sympy 1.14, numpy). Alle neun Skripte
enden ohne Fehler; die Ausgaben liegen als `ausgabe_*.txt` daneben. Erste Runde: A1 bis A6.
Fortsetzung am selben Tag: B1 bis B4, weiter unten. Auf ausdrücklichen Wunsch ist das Ergebnis
der Fortsetzung als offener Forschungsvertrag im Ledger (E8.STANDARDFORM.CHART.01, kein Skript,
keine Statusänderung an einem bestehenden Eintrag) und als Befund im Atlas eingetragen. PLAN
und progress sind unberührt.

## Übersicht

| Nr | Urteil | Stufe | Erste freie Wahl |
|---|---|---|---|
| A1 kleinstes Verbundsystem, mit B1 bis B4 | **Kandidat** | 2 für P2 unter benannter Lesart, 1 für Teile von P1 | Ganzzahligkeit (Gitter der Stabilisatorzustände) |
| A2 perfekter Code | **Kandidat** | 1 | Markierung einer der zehn 3 + 2 Teilungen |
| A3 Anomalie mit μ₄ | **Kandidat** | 1 | Dimension 4 oder die Zahl 16 |
| A4 Arithmetik über 𝔽₂ | Zahlen **verworfen**, Struktur **offen** | 0 | Kurve und Erweiterungsgrad |
| A5 universelle Optimalität | **offen** | 1 nur mit gesetzter Dimension 8 | Dimension 8 |
| A6 selbstdualer Code | **Kandidat**, als Abschluss in TFPT bekannt | 1 | Zerlegung in D₅ ⊕ A₃ |
| A7 Quantenzellautomat | **offen**, nicht gerechnet | – | Raumdimension des Gitters |
| A8 Kausalmengen | **offen**, nicht gerechnet | – | Dimension der Einbettung |

Kein Ansatz erreicht Stufe 2 ohne Vorbehalt, keiner Stufe 3 oder 4. Eine vollständige Lösung
liegt nicht vor.

## A1: kleinstes Verbundsystem

**Befund (exakt, `a1_zwei_qubits.py`).**
- Alle 6 vollständigen Sätze von fünf unverzerrten Basen bestehen aus 3 Produktbasen und
  2 verschränkten; alle 6 Sätze von fünf antikommutierenden Paulis aus 3 nichtlokalen und
  2 lokalen. Die lokale Symmetrie wirkt als S₃ auf dem Tripel und S₂ auf dem Paar.
- Die 16 geraden Monome der fünf Erzeuger sind eine Basis von M₄(ℂ). Mit der spurfreien
  Klassenfunktion der Teilung entstehen ν^c, u^c, Q, e^c, L, d^c mit Y = 0, −2/3, 1/6, 1,
  −1/2, 1/3; ΣY, ΣY³, SU(3)²Y und SU(2)²Y verschwinden.
- γ⁰ und γ⁵ sind die beiden lokalen Erzeuger, γ¹, γ², γ³ die drei nichtlokalen; Cl(1,3) gilt.
  Drehungen sind lokal, Boosts verschränkend.
- Cliffordgruppe der Ordnung 46080 mit Zentrum μ₄; invariantes Gitter mit 240 = 16 + 96 + 128
  kürzesten Vektoren; L/(1+i)L hat 15 Klassen zu 16 Wurzeln, bijektiv zu den 15 Basen; die fünf
  Erzeuger geben 5 × 48 = 144 + 96; Hopfbild zehn Ecken ±4e_a mit Faser 24-Zell.
- Standardform: J Δ^(1/2) a Ω = a* Ω für alle Matrixeinheiten; modulare Ladungen 0, 0, +1, −1.

**Basislinie.** Für zwei p-Niveau-Systeme ist die Teilung (p²+1; p+1, p²−p): 5 = 3 + 2 nur bei
p = 2, sonst 10 = 4 + 6, 26 = 6 + 20, 50 = 8 + 42. Innerhalb von zwei Qubits tragen 6 von 6
Spreads und 6 von 6 Ovoiden die Teilung. Die Zahl 5 allein ist in diesem System billig (auch
der spurfreie symmetrische Teil des Korrelationstensors hat Dimension 5); gewertet wird nur die
erzwungene Teilung zusammen mit der Graduierung 1 + 10 + 5. Die Hyperladungen sind kein
zusätzlicher Treffer: mit markierter 3 + 2 Teilung liegt Y bis auf die Skala fest.

**Stärkster Einwand.** Dieselben drei nichtlokalen Erzeuger werden als Farbplätze und als
Raumrichtungen gelesen. Das ist eine Struktur mit zwei physischen Lesarten und keine
Herleitung, dass innerer Träger und Raumzeit dieselben fünf sind. Die kontinuierliche
Eichgruppe folgt nicht; natürlich ist nur S₃ × S₂ beziehungsweise SO(3) × SO(2).

**Trägt / liefert nicht / bricht (Stand erste Runde).** Trägt 16 mit Hyperladungen, μ₄, E₈
Wurzeln, 5 × 48. Liefert nicht drei Familien, D₅ ⊕ A₃ mit ℤ₄ Glue, c₃, φ₀, 41, α, Massen,
Metrik, Gravitation. Bricht kein vorhandenes Ergebnis.

**Korrektur nach dem vertieften Abgleich.** Die 16 mit Hyperladungen ist die Identität aus
v774 und das Gitter die Struktur aus v783; beides ist nicht neu. Die Lesart der 16 als
Materietabelle auf Gitterebene ist durch v775 ausgeschlossen. Was A1 darüber hinaus beiträgt,
steht im Abschnitt „Fortsetzung“.

## A2: perfekter Code

**Befund (`a2_perfekter_code.py`).** Mit k = 1 ist (D, n) = (2, 5) die einzige Lösung für
D ≤ 64, n ≤ 2000. Der Code hat A = (1,0,0,0,15,0), B = (1,0,0,30,15,18); die 15 Einzelfehler
haben 15 verschiedene Syndrome; die fünf Syndromgeraden sind ein Spread mit 𝔽₄ Struktur
(genau zwei simultane Drehungen) und für genau 3 symplektische Formen lagrangesch. ZZZZZ ist
logisch, der gerade Sektor hat einen Zustand je Syndrom. Die Symmetrie wirkt als S₅ auf den
Plätzen und transitiv auf den zehn Paaren.

**Basislinie.** Die Familie n = D² + 1 mit vier Prüfsymbolen existiert für jedes D (10, 17,
26, …); k = 1 und damit n = 5 gibt es nur bei D = 2. Ein Treffer in 63 geprüften Alphabeten.

**Stärkster Einwand.** Die 16 zerfällt hier als 1 + 15 (Syndrome), bei TFPT als 1 + 10 + 5;
die Schwelle 3 + 2 markiert keine bestimmte Teilung. Warum die Natur perfekt und
nichtentartet codieren sollte, ist nicht begründet.

## A3: Anomalie mit μ₄

**Befund (`a3_anomalie_eta.py`).** η(RP^(2k−1)) mit μ₄ Twist: d = 4 gibt −1/16, Ladung 3 das
Gegenteil. Auf Λ^p ℂ⁵ ist X = 5 − 2p für gerades p immer 1 mod 4. 16 · η = −1, 15 · η = −15/16.

**Basislinie.** Ordnung 16 tritt mit μ₄ Twist nur bei d = 4 auf (d = 8 gibt 64, Literatur 32);
ohne Twist auch bei d = 6 (Literatur ℤ₁₆ für Pin⁻). Ungerade ℤ_n bis 9 geben 9, 5, 7, 27.

**Stärkster Einwand.** Das Ergebnis steht in der Literatur (García-Etxebarria und Montero
2018; Wang und Wen). Es verknüpft 16, μ₄ und 3 + 1, erzeugt aber weder 5 noch 3 + 2 noch die
Familienzahl, und eine der beiden Größen d = 4 und 16 bleibt Eingabe. Für d ≥ 6 gibt die
einfache Fixpunktformel die doppelte Ordnung der Bordismusgruppe; dort ist sie nur Schranke.

## A4: Arithmetik über 𝔽₂

**Befund (`a4_endliche_koerper.py`).** ℙ¹(𝔽₄) hat drei rationale Punkte und ein konjugiertes
Paar; über 𝔽₂ gibt es fünf Isogenieklassen elliptischer Kurven, drei supersingulär, zwei
gewöhnlich. Die Kurve mit Frobenius 1 + i hat 5 Punkte über 𝔽₄ und 41 über 𝔽₃₂.

**Basislinie.** 12 Treffer auf die zehn Zielzahlen bei erwarteten 7,6; für Ziele über 16 ein
Treffer (41) bei erwarteten 0,74. 68 = q(q²+1) ist einer von sieben Werten unter 800 und
5 + 63 eine Identität der Formel. Nichts liegt über der Basislinie.

**Urteil.** Zahlen verworfen. Offen bleibt nur die Strukturaussage, dass ℙ¹(𝔽₄) derselbe
Spread ist wie in A2; dort ist sie bereits gewertet.

## A5 bis A8

- **A5** liefert E₈ nur, wenn Dimension 8 gesetzt ist; die Sätze gelten ebenso für 24. Offen.
- **A6** (`a6_selbstdual.py`): keine doppelt geraden selbstdualen Codes für n = 2, 4, 6; genau
  30 für n = 8, alle äquivalent (|Aut| = 1344); Konstruktion A gibt 240. Ohne Wahl bis E₈. In
  TFPT als Abschluss vorhanden; als Ursprung fehlt der Schritt zu D₅ ⊕ A₃. 30 als Zahl nicht
  über der Basislinie.
- **A7, A8** nicht gerechnet. A7 wäre die natürliche Dynamik zu A1 (vier Komponenten), setzt
  aber die Raumdimension. A8 hat keinen Strukturtreffer in Sicht.

## Fortsetzung: die drei entscheidenden Tests (B1 bis B4)

Gerechnet am selben Tag mit `.venv/bin/python`: `b1_klebegleichung.py` (8 Prüfungen),
`b2_woerterbuch.py` (22), `b4_fuenf_basen.py` (11), dazu `b3_gegenproben.py` (8 Gegenproben,
alle feuern). Alle vier enden mit Rückgabewert 0; Ausgaben in `ausgabe_b*.txt`.

**Urteil: Kandidat. Stufe 2 für P2 unter einer benannten Lesart, Stufe 1 für Teile von P1.
Keine vollständige Lösung.**

### Was exakt gezeigt ist

1. **Wörterbuch (B2, Teil 1).** Die 240 Wurzeln der Standardform von M₂(ℂ), geschrieben in
   Paulikomponenten des Operators a in a Ω, sind Koordinate für Koordinate die 240 Wurzeln der
   TFPT Karte (v634, v689). Dabei ist: Familienpaare = X, Y, Z Komponenten; Ankerpaar =
   Einskomponente mal i, also die Zustandsgerade ℂΩ; J = Multiplikation mit i; σ = die
   Clifforddrehung X → Y → Z, verwirklicht als u ⊗ ū; die äußere antiunitäre Spiegelung aus
   v783 = die modulare Spiegelung a Ω ↦ a* Ω = diag(+,−,+,−,+,−,−,+).
2. **Fünf Basen (B4).** In dieser Karte ist σ lokal. σ lässt genau drei der 15 Messkontexte
   fest, alle verschränkt, darunter die Bellbasis, und sie ist die Ankerklasse. Drei der sechs
   Sätze von fünf unverzerrten Basen sind σ-invariant; in jedem dreht σ die 3 Produktbasen und
   lässt die 2 verschränkten fest. Klasse in L/(1+i)L = Kontext, und ħ = 0 genau dann, wenn zwei
   Kontexte einen Pauli teilen. Der Ledger führt die fünf Trägerplätze als ħ-Ovoid der
   σ-erhaltenen Form q* (NORMALFORM.CFIN.01). Damit sind die fünf Plätze fünf unverzerrte Basen
   und die Markierung 3 + 2 ist die Teilung Produkt / verschränkt. Die fünf 48er Pakete sind die
   σ-Bahnen der Kontexte: 3 Bahnen aus Produktbasen, 1 Bahn aus verschränkten, 1 festes Tripel.
3. **Lage von D₅ ⊕ A₃ (B2, Teile 2 bis 4).** Von den 7560 Lagen von A₃ in E₈ sind 18
   σ-invariant. Sie bilden genau zwei Bahnen unter dem Zentralisator von σ (Ordnung 72): 12 ohne
   Bezug zu J und 6, in denen D₅ eine komplexe Gerade enthält. In allen 18 dreht σ drei der fünf
   D₅ Achsen und lässt zwei fest; σ-Invarianz allein erzwingt also 3 + 2. In den 6 sind die zwei
   festen Achsen eine der drei σ-festen komplexen Geraden, und die drei gedrehten sind die
   J-Partner der A₃ Achsen. Die beiden σ-stabilen Koordinatenteilungen aus v775 gehören dazu;
   die Arbeitskarte aus v634 (D₅ auf den Koordinaten 0 bis 4) ist ein anderer Typ und nicht
   σ-invariant.
4. **Klebegleichung (B1).** n = 4, m = 5 ist die einzige gemeinsame Lösung von m = n + 1,
   n² = 2^(m−1) und m/4 + (n−1)/n = 2; aus (i) und (iii) folgt n² − 3n − 4 = 0. Mit n = k² ist
   k = 2. Das Gitter für (5, 4) hat 240 = 52 + 60 + 64 + 64 Wurzeln nach Klebeklasse.
5. **Familien und Anker (B2, Teile 5 und 6).** Unter σ ist 4 = 1 + 3 und 16 = 4 feste + 4
   Dreierbahnen. Die Normen von 1, i, 1 + i auf der Zustandsgeraden sind (1, 1, 2), Summe 4:
   der Anker, auf den AX.P1.01 im Ledger reduziert ist.

### Basislinie

- Wörterbuch: Die naheliegenden Alternativen scheitern klar (Ankerachse 1 statt i: 112 von
  240; Rechenbasis: 0 von 240; einseitige Drehung: 240 von 240 Abweichungen). Dass das gaußsche
  E₈ Gitter überhaupt eine Zwei-Qubit-Lesart hat, ist nicht neu (v783) und zählt nicht.
- Lage: 18 von 7560 Lagen sind σ-invariant, 6 von 7560 mit komplexer Gerade.
- Klebegleichung: Suche über n, m bis 400; die Gitterbedingung allein gibt eine unendliche
  Reihe, erst die zweite Forderung wählt aus. 5/4 + 3/4 = 2 selbst ist bekannt (v1).
- 4 = 1 + 3 und (1, 1, 2) sind kleine Zahlen; sie zählen nur, weil sie an der Stelle stehen,
  an der TFPT sie benutzt. Für sich liegen sie nicht über der Basislinie.

### Erste freie Wahl

Die Ganzzahligkeit: Benutzt wird das ℤ[i]-Gitter, das die Stabilisatorzustände aufspannen
(gleichwertig: die Cliffordgruppe als Normalisator der Pauligruppe). Die Standardform von
M₂(ℂ) allein gibt ℂ⁴ mit Ω, a ↦ a* und i, aber kein Gitter. Zweite Wahl: die Lesart „Platz =
Basis einer vollständigen unverzerrten Beschreibung“.

### Stärkster Einwand

Das Wörterbuch ist eine exakte Umbenennung eines Objekts, das TFPT schon hat. Es ersetzt die
getrennten Eingaben Deck, Orientierung, Anker und Länge 8 der Kette aus SEAM.CODE.TYPEII.01
durch ein Objekt, beweist aber nicht, dass die Naht dieses Objekt ist; die Nahtaxiome A1 bis A3
bleiben daneben stehen. Die Klebegleichung ist eine Eindeutigkeit innerhalb einer gewählten
Familie von Gittern. Der gefeuerte Negativbefund ARF.ROOTCLASS.01 (v775) gilt unverändert:
Klassen sind keine reinen Materiezustände. Die Hyperladungstabelle aus A1 ist die Identität
aus v774 mit genau dieser eingezäunten Lesart und wird hier nicht als Materieaussage gewertet.

### Trägt / liefert nicht / bricht

Trägt: fünf Plätze mit markierter Teilung 3 + 2, Deck σ, Orientierung J, Anker, D₅ ⊕ A₃ mit
Kleber, 4 = 1 + 3. Liefert nicht: Ganzzahligkeit, 2π in c₃, nichttriviale modulare Zeit (für
den Spurzustand trivial; jedes andere Gewicht wäre ein freier Parameter), kontinuierliche
Eichgruppe, die Coxeteruhr, φ₀, 41, α, Massen, Kontinuum, Gravitation. Bricht kein vorhandenes
Ergebnis.

## Abgleich mit dem Korpus (Atlas, Ledger und Skriptköpfe; keine vollständige Durchsicht)

Die erste Runde hat den Korpus unterschätzt. Korrigierter Stand:

- **Im Korpus vorhanden:** E8.CLIFFORD2Q.01 (v783: die 60 gaußschen Geraden sind die 60
  Stabilisatorstrahlen, G31 ist bis auf Index 2 die Zwei-Qubit-Cliffordgruppe, Spreads = Sätze
  unverzerrter Basen); E8.PROJHAMMING.01 (v752); ARF.SPINORCOMPILER.01 (v774: 16 = 1 + 5 + 10
  auf L/(1+i)L, Paritätslift auf C_even(5), Hyperladungstabelle, σ als 3 + 2 auf den fünf
  Plätzen); NORMALFORM.CFIN.01 (v845); SEAM.CLIFFORD.MODULAR_S.01 (v798);
  SEAM.CODE.TYPEII.01 (v799: P2 auf die Nahtaxiome A1 bis A3 eingeengt); die beiden σ-stabilen
  Koordinatenteilungen (v775); q(D₅) + q(A₃) = 2 (v1). Damit sind A1 Teil 2 (Hyperladungen),
  A1 Teil 4 (Gitter) und A6 keine neuen Befunde.
- **Negativbefund mit passenden Prämissen:** ARF.ROOTCLASS.01 (v775) schließt aus, dass eine
  Klasse von L/(1+i)L genau einen Gewichtszustand der 16 trägt, für jede Einbettung von
  D₅ ⊕ A₃. Das trifft jede Lesart von A1 Teil 2 als Materietabelle auf Gitterebene. Es trifft
  nicht das Wörterbuch, die Lageklassifikation und die Klebegleichung; diese machen keine
  Aussage über Materiezustände.
- **In den gelesenen Quellen nicht gefunden:** die Lesart der Karte als Standardform von M₂(ℂ)
  (Paulikomponenten, Ankerpaar = Zustandsgerade, σ = u ⊗ ū, Spiegelung = a ↦ a*); 3 + 2 als
  Produkt / verschränkt; die Klassifikation aller 7560 Lagen; die drei Forderungen der
  Klebegleichung; die Diophantische Eindeutigkeit des perfekten Codes (A2); die ℤ₁₆ Anomalie
  auf RP⁵ (A3). Gelesen wurden die Köpfe von v634, v689, v752, v774, v775, v783, v798, v845
  und die zugehörigen Ledgerzeilen. Fehlende Treffer sind kein Beleg, dass es keinen Eintrag gibt.
- **Vorläufer im Atlas:** `finding:mixed-window-pauli-five-point-isometry-20261116` hat bereits
  fünf feste spurfreie Cliffordrichtungen mit der Teilung 2 + 3 nach Fensterstütze (Geltung:
  Stütze, Algebra, unmarkierte Quadratik). Die Teilung „2 lokal + 3 nichtlokal“ aus A1 Teil 1
  ist damit für antikommutierende Fünfer nicht neu; neu bleibt die duale Aussage für die fünf
  Basen und ihr Anschluss an σ.
- **Einschränkung mit passenden Prämissen:**
  `finding:tfptfollowup20261110-local-tomography-deficit-and-coxeter-obstruction-to-gaussian-i`
  zeigt, dass das gaußsche J mit σ kommutiert, aber nicht mit dem nativen Coxeterelement. Die
  Coxeteruhr (Ordnung 30) ist deshalb keine Symmetrie der Standardformkarte; die Karte trägt
  μ₄, σ und J, aber nicht die Coxeterstruktur.
- **Parallele Explorationen vom selben Tag (Remote-Zweige, nicht auf main):**
  `tfpt/exploration-20261120-foreign-origins` verwirft den Spurzustand als Zeitanfang (trivialer
  modularer Fluss in M₂, exakter Gegentest), hält fest, dass endliche Algebra den Faktor 1/(8π)
  nicht erzeugen kann, wertet den perfekten Code mit demselben Einwand wie A2 und den Typ-II-Code
  als bekannten Kontrollweg wie A6, und zeigt mit einer Gegenkontrolle, dass ein E₈ aus einem
  beliebigen 𝔽₂⁴ Syndromraum generisch ist. Das trifft A6 und die bloße Aussage „E₈ entsteht“,
  nicht das koordinatengenaue Wörterbuch (Gegenproben G1 bis G3).
  `exploration/2026-11-20-vertauschung-und-auslese` erreicht fünf markierte Richtungen aus einem
  Zopfmodell über 𝔽₄ (Kandidat, Stufe 1); die Ansätze überschneiden sich nicht.
- **Weitere Ausschlüsse:** `finding:identical-overlapping-native-five-codes-have-no-common-state`
  betrifft starre identische Codes auf überlappenden Registern; A2 benutzt einen einzelnen
  Block. `finding:p1-z4-cut-eta-equals-d5-a3-glue-spin-weights-pattern-20261119` rechnet η auf
  S³/ℤ_p; A3 rechnet auf RP⁵ mit anderer Fragestellung. Beide treffen nicht.

## Stand der drei Tests aus der ersten Runde

1. **Wörterbuch: ja**, in der dualen Form. Die fünf Plätze sind fünf unverzerrte Basen (nicht
   fünf antikommutierende Paulis), und die Markierung 3 + 2 ist Produkt / verschränkt (B4).
2. **Glue: ja** als Eindeutigkeit von (4, 5) unter drei Forderungen (B1). Offen bleibt, warum
   die Zerlegung in D₅ ⊕ A₃ physisch ausgezeichnet ist.
3. **Familien: teilweise.** 4 = 1 + 3 unter σ ist exakt (B2); die Dreiheit ist die der drei
   Pauliachsen eines Qubits. Massen und Mischungen folgen daraus nicht.

## Kleinste entscheidende Tests für eine Folgesitzung

1. **Ganzzahligkeit.** Gibt es eine Forderung an die Standardform (zum Beispiel: die
   Beschreibung ist unter ihrer eigenen Symmetriegruppe abgeschlossen und diskret), die das
   Gitter der Stabilisatorzustände ohne Wahl erzwingt? Abbruch, wenn schon für ein einzelnes
   Qubit mehrere Gitter die Forderung erfüllen.
2. **Zeit.** Der Spurzustand hat trivialen modularen Fluss. Legt eine der vorhandenen TFPT
   Größen (φ₀) das Gewicht des Zustands fest, und bleibt das Gitter dabei erhalten?
3. **Naht.** Folgen die Nahtaxiome A1 bis A3 aus der Standardform, oder sind sie unabhängig?

## Wildcard

Sp(4,2) ≅ S₆ ist die einzige symmetrische Gruppe mit äußerem Automorphismus. Er vertauscht in
der Zwei-Qubit-Geometrie Punkte mit Geraden, also die sechs Sätze von Erzeugern mit den sechs
Sätzen von Messbasen. These: Die Welt ist der Fixpunkt dieser Vertauschung von Erzeugen und
Messen; 6 = 5 + 1 sind fünf Plätze und der Beobachter (die sechsbeinige perfekte
Tensorform des Codes aus A2). Abwegig, weil eine Gruppenkuriosität Physik begründen soll;
nicht auszuschließen, weil A1 und A2 genau diese beiden dualen Fünfen erzwingen. Nicht gerechnet.
