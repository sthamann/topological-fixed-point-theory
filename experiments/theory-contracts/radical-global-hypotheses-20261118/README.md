# Verrennen wir uns? Sechs Gesamthypothesen mit Gegenprobe

18. November 2026 · F1 und F2 · Kennung RADICAL20261118

Geltungsbereich: neue Gesamthypothesen zur Frage, was TFPT insgesamt ist, mit je
einer ausgeführten exakten oder numerischen Probe. Kein Ledgerwechsel, keine
Beförderung, keine physische Herleitung. Gelesener Stand: `progress.json`
Revision 405, Atlas `9b394151bbd090f3eead2480`. Die Messwerte in Probe B sind
aus dem Gedächtnis zitierte Literaturwerte (CODATA 2022, PDG 2024, NuFIT 6.0,
Planck 2018, ACT DR6) und wurden in diesem Lauf nicht gegen die Quellen
nachgeschlagen.

Reproduktion: `python3 -B check.py` (Exit 0, 20 Prüfungen bestanden, Ausgabe in
`run.log`, Daten in `results.json`).

## Ergebnis in vier Sätzen

1. Empirisch ist TFPT nicht widerlegt: kein Zug über 2 σ, ausgenommen die
   Leptonquellverhältnisse, die 0,2 % und 2,2 % neben den Polmassen liegen.
2. Die zehnstellige Übereinstimmung von α trägt allein wenig Beweislast. In
   einer Nachbarschaftsgrammatik derselben Gleichung sind 0,40 Zufallstreffer
   zu erwarten, gefunden wird genau der TFPT Treffer.
3. Neu: Die Normierung c₃ = 1/(8π) ist exakt E₂(i)/24, die Ladungsenergie pro
   Cartanrichtung jedes selbstdualen Gitters am selbstdualen Punkt τ = i. Das
   verklebte E₈ trifft c₃, das unverklebte D₅ + A₃ liefert 0,289 c₃.
4. Neu: Die absolute Konstante der Stromzweipunktfunktion von (E₈)₁ auf dem
   Torus ist −E₆(τ)/(12 E₄(τ)). Sie verschwindet genau bei τ = i. Dort ist die
   gewichtete reflektierte Zweifensterform ohne Referenzwahl gleich −℘(w; i).

Punkt 3 und 4 ergeben zusammen einen konkreten Kandidaten für das fehlende
Auswahlprinzip (G1). Die Gesamtlösung ist damit nicht gefunden. Der Anschluss
an die vier R113 Marken ist nicht gerechnet.

## Hypothesen und Verdikte

| Kennung | Gesamthypothese | Probe | Verdikt |
| --- | --- | --- | --- |
| G0 | Kompression statt Physik: Die Treffer sind nachträgliche Passung. | B, C | Weder bestätigt noch widerlegt. α allein entscheidet nicht. |
| G1 | Selbstdualer Punkt statt Quelle: Das fehlende Auswahlprinzip ist Stationarität am selbstdualen Modul. | D1 bis D11, Nachtrag | Als Aussage über Energien exakt. Als Produzent der Ortsfaktoren der R113 Marken nicht gestützt, siehe Nachtrag. |
| G2 | Familie gleich Kanal: Die drei Familien sind die drei Paarungen der vier Marken (s, t, u). | E1 bis E3 | Gruppenstruktur exakt. In der Flavormatrix R nicht sichtbar. |
| G3 | Horizont zuerst: Die Naht ist die Bifurkationsfläche, der Produzent ist der thermofeldverdoppelte Zustand. | D6 | Physische Lesart von G1. Keine eigene Probe. |
| G4 | Trägerskala gleich Skalaronskala. | F1 | Größenordnung passt, Faktor 2,97 bei einer Schleife. Als exakte Gleichheit ungeprüft. |
| G5 | Vier Familiengewichte sind physisch, drei sind eine bedingte Erwartung. | Bestand | Nicht neu gerechnet (ORIGINDIAG20261116, E8WORD20261116). |

## G0 Kompression statt Physik

**Probe B, Datenabgleich.** Züge (TFPT minus Messung, in σ):

| Größe | TFPT | Messung | Zug |
| --- | --- | --- | --- |
| α⁻¹ (CODATA 2022) | 137,0359992168 | 137,035999177(21) | +1,90 |
| α⁻¹ (Rb 2020) | | 137,035999206(11) | +0,99 |
| α⁻¹ (Cs 2018) | | 137,035999046(27) | +6,33 |
| λ_C gegen \|V_us\| | 0,224376 | 0,2243(8) | +0,10 |
| s₂₃ gegen \|V_cb\| | 0,043428 | 0,0411(12) | +1,94 |
| s₁₃ gegen \|V_ub\| | 0,003765 | 0,00382(20) | −0,28 |
| δ_CKM | 68,65° | 65,7(30)° | +0,98 |
| sin²θ₁₂ | 0,306747 | 0,308(12) | −0,10 |
| sin²θ₁₃ | 0,023108 | 0,02215(57) | +1,68 |
| Ω_b | 0,048941 | 0,0493(9) | −0,40 |
| β (ACT DR6) | 0,2424° | 0,20(8)° | +0,53 |

Die Cs Messung widerspricht auch der Rb Messung und zählt nicht als Befund
gegen TFPT. Die Leptonquellverhältnisse (12/7)φ₀² und (8/7)φ₀ liegen +0,21 %
und +2,20 % neben m_e/m_μ und m_μ/m_τ der Polmassen. Der Korpus führt sie als
Quellverhältnisse mit offenem Transfer. Dieser Transfer ist die schärfste
empirische Stelle, weil die Messgenauigkeit dort um Größenordnungen besser ist.

**Probe C, Nachbarschaftsgrammatik für α.** Gleichung
`α³ − t c₃³ α² − K c₃⁶ log(1/φ_s) = 0` mit `φ_s = 1/(kπ) + q(1−q)^(−p)`,
`q = m c₃⁴ exp(−eα)`. TFPT: t = 2, K = 164/5, k = 6, m = 48, e = 2, p = 5/4.
Variiert: t in 0..4, k in 1..12, m in 1..96, e in 0..4, p in Vierteln bis 2, K
rational mit Nenner bis 8. Das sind 259200 Kombinationen. Im Fenster CODATA
2022 ± 2 σ liegt genau ein Treffer, der TFPT Satz. Die Zufallserwartung ist
0,40.

Folgerung: Die Stellenzahl von α ist in dieser Grammatik kein Beleg. Die
Beweislast liegt in den getrennten Herleitungen von 164/5 = (4/5)·41, 48, 5/4
und den Faktoren 2. Der Korpus nennt für die gemeinsame deklarierte Grammatik
aller Ausgaben eine weit kleinere Nullwahrscheinlichkeit (AUDIT.NULLMC.01).
Beide Aussagen gelten für ihre jeweilige Grammatik und widersprechen sich
nicht.

## G1 Selbstdualer Punkt statt Quelle

ORIGINDIAG20261116 hielt fest, dass P1, P2 und E₈ die vier Marken nicht
festlegen und ein Auswahlprinzip fehlt. G1 schlägt als Prinzip vor: Die Naht
steht am selbstdualen Modul, an dem der Vakuumcharakter stationär ist und keine
absolute Referenz gebraucht wird.

Exakt geprüft (D1 bis D11):

- E₂(i) = 3/π und E₆(i) = 0. Damit c₃ = 1/(8π) = E₂(i)/24, und der Baumterm
  von φ₀ ist 1/(6π) = E₂(i)/18 = (4/3) c₃.
- Die Energie des Vakuumcharakters χ = E₄/η⁸ ist q∂_q log χ = −E₆/(3E₄). Sie
  verschwindet bei τ = i. Dort ist χ(i) = 12 und die Entropie log 12.
- Die Gitterladungsenergie ⟨p²/2⟩ pro Cartanrichtung bei τ = i ist für E₈
  exakt c₃. Die Oszillatoren tragen −c₃ pro Richtung. Für D₈ ergibt sich
  0,562 c₃, für das unverklebte D₅ + A₃ 0,289 c₃, für Z⁸ wieder c₃.
- Die Stromzweipunktfunktion auf dem Torus ist
  ⟨J^a(w) J^b(0)⟩ = δ^{ab} (−℘(w; τ) − E₆(τ)/(12 E₄(τ))). Die Konstante ist für
  drei Moduln gegen die direkte Gitterrechnung geprüft und verschwindet bei
  τ = i. Die Modensumme stimmt dort mit einer unabhängigen Gittersumme von −℘
  auf 7·10⁻⁹ überein.
- Der antipodale Zweipunktkern ist auf der Ebene −1/4 (so in E8WORD20261116
  verwendet). Auf dem Torus bei τ = i ist er −(ϖ/2π)² = −0,174150 mit der
  lemniskatischen Konstante ϖ. Der Unterschied beginnt mit 2c₃ = 1/(4π).

Was daraus folgt:

- P1 und der E₈ Abschluss hängen an einer Bedingung. Die Normierung c₃ pro
  Richtung gilt genau für selbstduale Gitter. Die Z₄ Verklebung von D₅ + A₃ ist
  dafür nötig. Die Geradheit, also die Abgrenzung gegen Z⁸, kommt aus der T
  Phase und nicht aus c₃.
- τ = i ist der einzige Modul mit verschwindender absoluter Konstante. E₆ hat
  nach der Valenzformel im Fundamentalbereich nur diese Nullstelle. Das ist ein
  Auswahlprinzip für den harmonischen Viermarkenpunkt, das die Marken nicht als
  Ziel einsetzt.
- Für die Cartanströme ist die gewichtete reflektierte Zweifensterform ohne
  Wahl bestimmt: Q₀ = −℘(w; i). Ihre Werte an den drei Halbperioden stehen im
  Verhältnis 1 : 0 : −1. Die vierte Marke ist der Pol.

Grenzen:

- Die Rechnung betrifft den holomorphen (E₈)₁ Torus und Cartanströme. Die
  Gleichsetzung des reflektierten Paars mit diesem Torus ist eine Annahme.
- Der Korpus leitet c₃ auch als 1/(|Z₂|·2π·χ(S²)) her (SEAM.KHALF.01). Dass
  beide Ausdrücke gleich sind, ist eine Konsistenz und beweist die
  Identifikation nicht.
- Kein Wörterbuch zum skalaren Higgsmarker, kein Feldmaß für wechselwirkende
  Felder, keine physische Zeitkarte. Die Bedingungen (c) und (d) aus
  ORIGINDIAG20261116 sind nicht erfüllt.

Entscheidender nächster Test: Die algebraischen Spurfaktoren aus
E8WORD20261116 (ZY/Z0 = 1/4, verbundene Hyperladungsantwort −1/32 und −1/72)
sind Nullmodenspuren und ändern sich auf dem Torus nicht. G1 ändert den
Ortsfaktor, von −1/4 zu −0,174150 am Gegenpunkt. Zu rechnen ist die Torus
Zwei und Dreipunktfunktion der geladenen Wurzelströme bei τ = i an den Marken
1, i, −1 und der Vergleich mit einer Ausgabe, die in G1 nicht eingeht. Der
erste Kandidat ist die Schranke (0, 1/(8π)] der antipodalen Antworten aus der
lokalen Deformation. Liefert der Torus sie, ist Bedingung (d) erfüllt. Liefert
er einen anderen Wert, ist G1 als Produzent der R113 Marken ausgeschlossen.

## G2 Familie gleich Kanal

W(A₃) = S₄ permutiert die vier Marken. Auf den drei Paarungen wirkt S₃ mit Kern
V₄ (E2). Am harmonischen Punkt hat jeder der drei Doppelverhältniswerte −1, 1/2,
2 den Stabilisator Z₂, am Dreierpunkt Z₃, generisch trivial (E1). τ = i würde
also eine Zerlegung 2 + 1 der Familien mit Austauschsymmetrie zweier Familien
liefern, passend zur zugewiesenen Textur θ₂₃ = 45°.

Gegenprobe E3: Die Flavormatrix R mit den Zeilen (1,3,0), (1,5,2), (2,5,3) hat
unter Zeilen und Spaltenpermutationen nur die Identität als Symmetrie. Die
Austauschsymmetrie ist im Quarksektor nicht vorhanden. Verdikt: als Herkunft
der Quarkhierarchie nicht gestützt, als Herkunft der Leptontextur ungeprüft.

## G3 Horizont zuerst

Physische Lesart von G1. Das reflektierte Paar ist ein thermofeldverdoppelter
Zustand bei selbstdualer Temperatur. Zustand, absolute Referenz und Zeitkarte
(modularer Fluss mit β = 2π) kämen dann aus einer Wahl. Der Preis: 3+1
Dimensionen und ein Horizont werden vorausgesetzt, Aufgabe E1 kehrt sich um.
Keine eigene Probe in diesem Lauf.

## G4 Trägerskala gleich Skalaronskala

Mit den Einschleifenkoeffizienten des Standardmodells erreicht sin²θ_W den Wert
3/8 bei 1,03·10¹³ GeV. Die Skalaronmasse c₃^(7/2) M̄ ist 3,06·10¹³ GeV. Das
Verhältnis ist 2,97. Die starke Kopplung trifft sich dort nicht mit den beiden
anderen. Verdikt: gleiche Größenordnung, keine belegte Gleichheit. Entscheidend
wäre die Zweischleifenrechnung mit Schwellen.

## Was das für den Plan heißt

F1 und F2 bleiben in Bearbeitung. Neu für F2: G1 ist ein benanntes
Auswahlprinzip mit einer einzigen zusätzlichen Annahme (reflektiertes Paar
gleich Torus bei stationärem Modul). Der nächste Schritt ist der oben genannte
Test an den geladenen Wurzelströmen.

## Nachtrag RADICALCHARGED20261118: geladene Wurzelströme bei τ = i

Reproduktion: `python3 -B check_charged.py` (Exit 0, 8 Prüfungen bestanden,
Ausgabe in `run_charged.log`, Daten in `results_charged.json`).

Gerechnet auf zwei unabhängigen Wegen, Gitterspur im D₈ Chart mit Primformen
und elliptische Funktionen aus den Ward Identitäten, für drei Moduln:

- Die Zweipunktfunktion der geladenen Wurzelströme ist dieselbe wie die der
  Cartanströme, −℘(w; τ) − E₆/(12 E₄). Das ist die E₈ Invarianz der Spur.
- Die Dreipunktfunktion ist bis auf eine feste Phase
  −(G³ + ½ Σ ℘′)/3 − κ G mit G = ζ(w₁₂) + ζ(w₂₃) + ζ(w₃₁) und
  κ = E₆/(12 E₄). Bei τ = i ist κ = 0.
- An den Marken 1, i, −1 bei τ = i: Zweipunktwert am Gegenpunkt
  −(ϖ/2π)² = −0,174150, Dreipunktbetrag (1 + √2)(ϖ/2π)³ = 0,175454. Die Ebene
  liefert −1/4 und 1/4. Die Spurfaktoren aus E8WORD20261116 bleiben unverändert.

Korrektur des angekündigten Tests: Die Schranke (0, 1/(8π)] aus
ORIGINLABH0INTERPRETATION20261116 und ORIGINLABSPHERICALPAIRING20261116 betrifft
eine Gaußkovarianz vom Gewicht 1/2 auf der runden Sphäre. Ein Stromkorrelator
ist damit nicht vergleichbar. Der Vergleich war so, wie er oben angekündigt
ist, falsch gestellt. Der gleichartige Vergleich ist der Kern vom Gewicht 1/2
auf dem Torus. Er liefert am Gegenpunkt für die drei geraden Spinstrukturen 0,
θ₄(i)²/2 und θ₃(i)²/2, in Einheiten des Sphärenwerts 0, 0,8346 c₃ und
1,1803 c₃. Der Sphärenwert c₃ wird von keiner Spinstruktur getroffen.

Verdikt:

- Bedingung (d) ist im einzigen gleichartigen Vergleich nicht erfüllt. G1 ist
  als Produzent der Ortsfaktoren der R113 Marken nicht gestützt.
- Alle Punktwerte am selbstdualen Torus enthalten die lemniskatische Konstante
  ϖ. Die TFPT Ausgaben enthalten nur π. Falls der selbstduale Torus die Quelle
  ist, können die Auslesungen nur Energien und Nullmodengrößen sein, keine
  Punktwerte. Auf dieser Ebene bleiben c₃ = E₂(i)/24 und das Verschwinden der
  absoluten Konstante exakt.
- Kein Ausschluss von τ = i als Auswahlprinzip. Kein Ausschluss anderer
  Träger, etwa der Sphäre mit Gewicht 1/2, auf der der Bestand c₃ am Gegenpunkt
  bereits findet.

Offene Frage, die daraus entsteht: Auf der runden Sphäre ist c₃ ein Punktwert
(Gegenpunkt, Gewicht 1/2), auf dem selbstdualen Torus eine Energie pro
Cartanrichtung. Ob beide Ausdrücke dieselbe Größe sind, entscheidet eine
Rechnung, die Sphäre und Torus über dieselbe Viermarkenüberlagerung verbindet.

## Nachtrag RADICALBRIDGE20261118: Sphäre und Torus

Reproduktion: `python3 -B check_bridge.py` (Exit 0, 9 Prüfungen bestanden).

Frage: Ist c₃ als Punktwert auf der Sphäre und c₃ als Energie auf dem Torus
dieselbe Größe?

- Lemma: Gilt f(1/t) = t^k f(t), dann ist die Energie am Fixpunkt t = 1 gleich
  k/(4π). Eine Gitterrichtung hat Gewicht 1/2. Daraus folgt c₃ = (1/2)/(4π)
  ohne jede Normierungswahl, für jedes selbstduale Gitter, auch im Rang 24.
- Exakte Zerlegung für E₈ pro Cartanrichtung: Energie = c₃/Im τ plus ein
  modularer Rest, der bei τ = i verschwindet. Dazu die Summenregel
  E(t) + E(1/t)/t² = 1/(4πt).
- Auf der Sphäre ist der Kern 1/(4π(1 + r)) unter der Kelvin Inversion
  kovariant mit Gewicht 1. Das Gewicht legt nur die logarithmische Ableitung
  fest, also den Fluss c₃/2 = 1/(16π), den der Bestand als k führt
  (SEAM.KHALF.01). Der Wert c₃ selbst braucht die Normierung 4π und den Radius
  eins. Ein zweiter Kern gleichen Gewichts hat am Fixpunkt √2 c₃.
- Die paarweise Zuordnung Sehne = 2 Im τ scheitert: Sehnen (2, √2, √2) gegen
  (2, 2, 1).

Verdikt: Beide c₃ sitzen am Fixpunkt einer Inversion vom Gewicht 1/2. Der
Sphärenwert folgt nicht aus dem Toruswert und umgekehrt. Die Torusform ist die
stärkere, weil sie keine Normierung wählt.

## Nachtrag RADICALHIST20261118: frühe Fassungen und historischer Austausch

Reproduktion: `python3 -B check_historical.py` (Exit 0, 8 Prüfungen bestanden).

Gelesen: `_archive/paper-latex/old` (Theory-Root V1.0 vom 27. September 2025,
Unified Field Equation, ImprovedCubic, Eliminating K) und die anonymisierten
Abschriften in `historical-dynamics-source-recovery-20261103`. Die rohen
Nachrichten des Austauschs liegen außerhalb des Repositorys und wurden nicht
gelesen.

- Die Schließung von φ wurde dreimal geändert. Gegen CODATA 2022 liefert
  Baum plus δ den Zug +23918 σ, δ(1 − 2α) den Zug −418 σ, δ exp(−2α) den Zug
  −242 σ und die heutige Form q(1 − q)^(−5/4) den Zug +1,9 σ. Jede Stufe rückte
  näher an den bekannten Messwert. Zusammen mit Probe C stärkt das G0. Ob der
  Exponent 5/4 vor dem Vergleich mit dem Messwert hergeleitet war, ist hier
  nicht geprüft.
- Herkunft der 41: 8 b₁ = 164/5 mit b₁ = 41/10, dem Einschleifenkoeffizienten
  der Hyperladung für drei Familien und ein Higgs.
- Die Kubik ist keine Laufgleichung. Sie liefert d(1/α)/dL = −15,48, ein
  Einschleifenlauf gäbe +0,65.
- In der Variable y = 1/(64π²α) lautet sie 8 b₁ L y³ + y/(4π) = 1. Die Wurzel
  ist ein Minimum des historischen Potentials U(α). Lesart: α als Feld mit
  eigenem Potential. Nicht weiter geprüft.
- Die historische k Kubik (V25-C) gibt 2πk = 1,0019 φ₀. Das geschlossene Paar
  φ = 2πk verfehlt α um mehr als 1000 σ und ist ausgeschlossen.
- ξ = c₃/φ_Baum = 3/4 und 12π φ_Baum = 2.

## Nachtrag RADICALCHRONO20261118: Chronologie des Koeffizienten 5/4

Reproduktion: `python3 -B check_chrono.py` (Exit 0, 9 Prüfungen bestanden).

Quelle: die datierten Fassungen unter `documents/archive--paper-latex--*`. Das
Git Protokoll beginnt erst am 27. April 2026 und enthält 5/4 von Anfang an,
die Reihenfolge stammt deshalb aus den archivierten Fassungen.

| Fassung | Datum | Zweidefektterm | α⁻¹ | Zug gegen CODATA 2022 |
|---|---|---|---|---|
| 2.2 | 21. Januar 2026 | keiner | 137,03599410 | −242 σ |
| 2.4 | 22. Januar 2026 | δ_top², als Diagnose bezeichnet | 137,03599819 | −46,9 σ |
| 2.5 | 27. Januar 2026 | (g/4) δ_top² mit g = 5 | 137,035999216 | +1,9 σ |
| 3.1 | später | resummiert, q(1 − q)^(−5/4) | 137,0359992168 | +1,9 σ |

- Der Messwert 137,035999177 steht schon in den Fassungen 2.2 und 2.4. Der
  Koeffizient 5/4 erscheint fünf Tage, nachdem der Koeffizient 1 den Messwert
  um 47 σ verfehlt hatte.
- Fassung 2.5 rechnet den Wert von δ₂ aus, der den Messwert genau trifft
  (Zeile „targeted match"). In Einheiten von δ_top² ist das 1,2404 ± 0,0051.
  Der Koeffizient 5/4 liegt 1,9 σ darüber. Fassung 2.5 nennt den alten
  Koeffizienten 1 nachträglich „g = 4" und den neuen „Derived, no fits".
- Suchraum: Unter den 18 Brüchen mit Nenner bis 4 im Bereich (0, 3] trifft nur
  5/4 innerhalb von 2 σ, die Zufallserwartung ist 0,12. Bei Nennern bis 12 sind
  es 138 Brüche, Erwartung 0,94, Treffer wieder nur 5/4. Die Nachbarn 6/5, 4/3,
  9/8 und 3/2 verfehlen um mehr als 7 σ.
- Die Resummation zu (1 − q)^(−5/4) kam später und verschiebt den Zug nur von
  1,865 auf 1,897 σ. Sie ist an α nicht prüfbar.

Urteil: Der Koeffizient 5/4 war nicht vor dem Vergleich mit dem Messwert
hergeleitet. Er wurde eingeführt, als der Messwert und der nötige Wert von δ₂
bekannt waren. Das beweist keine Anpassung: 5/4 ist ein einfacher Bruch mit
einer Strukturbegründung (g = 5), und ein Zufallstreffer in der engen Auswahl
hat etwa 12 % Wahrscheinlichkeit. Es heißt aber, dass die Übereinstimmung auf
1,9 σ kein unabhängiger Test der Theorie ist. Ohne nachträglich gewählte
Bausteine bleibt die Fassung 1.0 mit einer Abweichung von 3,7 ppm.

Nicht gelesen: das in Fassung 2.5 zitierte Prüfskript
`defect_partition_g5_audit`. Im gelesenen Satz wird der Nenner 4 gesetzt und
nicht hergeleitet. Die rohen Nachrichten des Austauschs liegen außerhalb des
Repositorys.

Nächster Schritt: Das Zweidefektgewicht 5/4 ohne Blick auf Alpha herleiten: aus der tatsächlichen Zweidefekt Zustandssumme auf dem orientierbaren Doppelüberzug den Koeffizienten von q^2 rechnen und mit 5/4 vergleichen. Alpha selbst kann das nicht mehr prüfen (die Resummation verschiebt den Zug nur um 0,03 sigma).

## Nachtrag RADICALTWODEFECT20261118: Zweidefektgewicht ohne Blick auf α

Reproduktion: `python3 -B check_twodefect.py` (Exit 0, 9 Prüfungen bestanden).

Gelesen: `documents/tfpt_1_architecture_e8/sections/007-source.tex` und
`verification/v341_alpha_quillen.py`. Eine Zweidefekt Zustandssumme ist in den
Quellen nicht definiert. Die heutige Begründung ist die Gleichsetzung
5/4 = q(D₅), die Norm des Klebevektors von D₅. Das Skript v341 nennt die
Herleitung des Funktionals selbst offen. Das in Fassung 2.5 zitierte Prüfskript
`defect_partition_g5_audit` ist im Repository nicht vorhanden.

Vorab festgelegt wurden zwei Listen, erst danach kam der Vergleich.

- Neun übliche Zählungen eines Defektgases (verdünntes Gas, Bose, Fermi, g
  Moden, geordnete und ungeordnete Paare, Doppelüberzug) geben für den
  Koeffizienten von q² die Werte −1, 0, 1/2, 3/5, 4/5, 1, 2, 3. Keine gibt 5/4.
  Ausschlussstatistik bräuchte den Parameter −1/8, außerhalb des zulässigen
  Bereichs.
- Zweipunktfunktionen der Primärfelder: ⟨V(1) V†(q)⟩ = (1 − q)^(−2h). Der
  Spinor von so(10)₁ hat h = 5/8, das Fundamentalfeld von su(4)₁ hat h = 3/8,
  zusammen 1. Damit ist q(1 − q)^(−5/4) genau q mal die Zweipunktfunktion des
  D₅ Spinorfelds, in allen Ordnungen (1, 5/4, 45/32, 195/128).
- Vergleich mit CODATA: Exponent 3/4 gibt −95,5 σ, 1 gibt −46,8 σ, 5/4 gibt
  +1,9 σ, 2 (ganzer E₈ Klebevektor) gibt +148 σ. Von vier nativen Exponenten
  trifft einer.
- α legt nur eine Zahl fest. Die abgebrochene Reihe, (1 − q)^(−5/4) und das
  Produkt über alle n unterscheiden sich um weniger als 0,1 σ.
- Der Wert q = 1,186·10⁻⁴ ist nicht das Doppelverhältnis der vier Marken (1/2).

Urteil: Das Gewicht 5/4 folgt aus keiner Zählung von Defektpaaren. Die
Begründung „g = 5 Holonomieentartung" der Fassung 2.5 ist damit nicht
gestützt. Es gibt ein natives Objekt, das die ganze Funktion liefert: die
Zweipunktfunktion des D₅ Spinorfelds. Das ist eine Benennung und keine
Herleitung, solange zwei Dinge fehlen: warum nur der D₅ Faktor des
Klebevektors eingeht und warum der Einsetzungsabstand gleich q ist.

Arbeitshypothese G6: Die Rückkopplung ist die Zweipunktfunktion des Feldes,
das die Z₄ Verklebung trägt, eingeschränkt auf den Träger D₅.

## Nachtrag RADICALSEPARATION20261118: Ist q ein Einsetzungsabstand?

Reproduktion: `python3 -B check_separation.py` (Exit 0, 9 Prüfungen bestanden).

Gelesen: `verification/v1033_charged_disorder.py` (dort h = 5/8 + 3/8 = 1 für
den Klebevektor) und der Abschnitt zur elektromagnetischen Schließung in
`documents/tfpt_1_architecture_e8/sections/007-source.tex`.

- Auf dem Zylinder (Naht als Kreis mal Zeit) lautet die Zweipunktfunktion
  eines Felds vom Gewicht h: x^h (1 − x)^(−2h) mit x = exp(−t).
- Das D₅ Spinorfeld allein gibt q^(5/8) (1 − q)^(−5/4) und verfehlt den
  Messwert um mehr als 10⁷ σ. Der ganze Klebevektor gibt q (1 − q)^(−2) und
  verfehlt um +148 σ.
- Die Formel ist ein gemischtes Objekt: führende Potenz 1 = 5/8 + 3/8 vom
  ganzen Klebevektor, Exponent 5/4 nur von D₅. Gleichwertig: Zweipunktfunktion
  des ganzen Klebevektors, geteilt durch den Faktor (1 − q)^(−3/4) von A₃.
  Eine Regel, die den A₃ Faktor entfernt, steht in keiner gelesenen Quelle.
- Zweite, gleich kurze Lesart: q (1 − q)^(−5/4) = 4 q d/dq (1 − q)^(−1/4),
  also die Ableitung nach α der Zweipunktfunktion eines Felds vom Gewicht 1/8.
  Die Lesart ist nicht eindeutig.
- Als Abstand müsste q die Zeit t = 2α + t₀ mit t₀ = log(1/δ_top) = 9,02548
  bedeuten. α legt t₀ auf 6·10⁻⁷ fest. Keine von 19 vorab benannten nativen
  Zeiten (2π, Vielfache von ϖ, Logarithmen von 240, 248, 1920, 12 und andere)
  trifft; die nächste liegt 0,0255 daneben, das sind mehr als 40000 σ.
- δ_top = 48 c₃⁴ = 3/(256 π⁴) = 1/(2⁸ vol S⁷). Das ist ein Gewicht (vier
  Kernwerte mal 48 oder eine Dichte auf der Einheitssphäre), kein Abstand.
- Bei gleicher Zeit geben die Markenabstände feste Zahlen ohne α (0,420 am
  Gegenpunkt, 0,648 am Viertelkreis).

Urteil: G6 ist in der strengen Form nicht gestützt. q ist ein Gewicht und kein
Einsetzungsabstand, und keine native Zweipunktfunktion auf dem Zylinder gibt
die Formel. Die Benennung 5/4 = 2h des D₅ Spinors bleibt eine exakte
Gleichheit ohne Mechanismus. Siehe die Korrektur nach dem Atlasabgleich weiter unten.

## Nachtrag RADICALPROSPECT20261118: Prospektive Testgröße

Reproduktion: `python3 -B check_prospective.py` (Exit 0, 8 Prüfungen bestanden).

Gesichtet: Atlas (Suche nach prospektiven Vorhersagen), `verification/v283_live_killtest_scorecard.py`
und die aktuellen Paper. α wird dort nur mit den CODATA Mittelwerten 2018 und
2022 verglichen. Die drei einzelnen Bestimmungen und das magnetische Moment
des Elektrons stehen in keiner Quelle des Repositorys. Die Messwerte wurden am
18. November 2026 per Websuche bestätigt.

| Bestimmung | α⁻¹ | TFPT minus Messwert |
|---|---|---|
| Rb Rückstoß, Morel et al. 2020 | 137,035999206(11) | +0,99 σ |
| Elektron g−2 mit QED, Fan et al. 2023 | 137,035999166(15) | +3,39 σ |
| Cs Rückstoß, Parker et al. 2018 | 137,035999046(27) | +6,33 σ |
| CODATA 2022 (Mittel, Fehler aufgeweitet) | 137,035999177(21) | +1,90 σ |

- Die drei Bestimmungen widersprechen einander (Rb gegen Cs 5,5 σ, g−2 gegen
  Cs 3,9 σ, Rb gegen g−2 2,2 σ). Ihr gewichtetes Mittel ist 137,0359991778(84)
  mit χ² = 31 bei 2 Freiheitsgraden. Der CODATA Fehler ist deshalb aufgeweitet.
- Die eingefrorene Wurzel 137,0359992168407 liegt beim Rb Wert, 3,4 σ über dem
  g−2 Wert und 6,3 σ über dem Cs Wert.
- Mit Standard QED folgt aus dem TFPT α ein magnetisches Moment
  a_e = 0,00115965218016. Gemessen ist 0,00115965218059(13), also 4,3·10⁻¹³
  mehr (3,3 σ der Messung allein).
- Der Koeffizient 5/4 liefert nicht das angepasste Mittel: Der CODATA optimale
  Koeffizient 1,2404 gäbe 137,035999177, die Wurzel liegt 4,0·10⁻⁸ darüber.
- Andere Ablesungen benutzen den ungedämpften Seed φ₀. Er weicht relativ um
  3,2·10⁻⁵ vom Nahtwert der α Gleichung ab. Keine andere Ablesung wird mit
  dieser Genauigkeit verglichen, dort gibt es keinen Test.

Testgröße: α⁻¹ selbst gegen die nächste Einzelbestimmung. Bei einer künftigen
Unsicherheit von 3·10⁻⁹ trennt die Vorhersage sich von den heutigen
Zentralwerten um 3,6 σ (Rb), 16,9 σ (g−2) und 56,9 σ (Cs). Regel: Eine
Bestimmung mit Unsicherheit bis 5·10⁻⁹, die mehr als 5 σ von 137,0359992168
abweicht, widerlegt die Gleichung samt Exponent 5/4. Übereinstimmung innerhalb
von 2 σ bei dieser Genauigkeit wäre die erste unabhängige Stütze.

Grenze: Die heutige Nähe zum Rb Wert ist keine unabhängige Stütze, weil der
Koeffizient mit Blick auf CODATA gewählt wurde und Rb das CODATA Mittel
mitbestimmt. Der Test verlangt neue Messungen. Die Rechnung für a_e benutzt
nur die führende Empfindlichkeit.

## Korrektur nach Atlasabgleich (18. November 2026)

Beim Eintragen in den Atlas (`problem:alpha-fixed-input-root`) zeigte sich,
dass zwei Aussagen der Nachträge RADICALTWODEFECT und RADICALSEPARATION zu
weit gingen. Der Atlas führt bereits:

- `finding:marked-e8-sector-grade-response-realizes-five-quarters`: In der
  markierten E₈ Quelle gilt exakt
  G_u(q₅, q₃) = [(1 − u)/(1 − u q₅)]^(5/4) [(1 − u)/(1 − u q₃)]^(3/4). Die
  getrennte D₅ Auslese gibt (1 − u q)^(−5/4), die gemeinsame Zeitdiagonale den
  Exponenten 2. Das ist dieselbe Struktur wie das hier gefundene gemischte
  Objekt und war vorher bekannt. Offen sind dort u = 48 c₃⁴, q = exp(−2α) und
  die physische Präparation.
- `finding:alpha-archive-term-provenance-and-master-functional-scope`: Die
  späteren Fassungen begründen 5/4 als B γ = (3/2)(5/6) und als Antwort eines
  Rang eins Determinanten.

Richtig ist deshalb: Der Exponent 5/4 hat eine exakte Herleitung in einem
festgelegten Modell (getrennte D₅ Auslese). Nicht hergeleitet sind die Wahl
dieser Auslese statt der gemeinsamen Zeit, das Gewicht 48 c₃⁴ und die Kopplung
exp(−2α). Die Chronologie bleibt: Der Wert 5/4 kam nach dem Vergleich mit dem
Messwert, die Modellherleitung später. Der Satz, es gebe keine nachträglich
gelungene Herleitung, ist zurückgezogen.

## Nachtrag RADICALNORMAL20261118: Die Gleichung in der Variable g

Reproduktion: `python3 -B check_normal.py` (Exit 0, 8 Prüfungen bestanden).

- Mit g = α/(2 c₃²) = 32π² α = e²/c₃ lautet die Gleichung exakt
  g² (g − c₃) = b₁ log(1/φ). Die Wurzel ist g = 2,30470. Die Faktoren 8 und 2
  der üblichen Schreibweise verschwinden. Ihre Deutungen als Rang von E₈ und
  als Ordnung von Z₂ hängen an der Wahl der Variable, ebenso die Potenz c₃⁶.
- Die ganze Gleichung braucht fünf Zahlen: c₃, b₁ = 41/10, 4/3 (Basis
  φ = (4/3) c₃), 48 und 5/4.
- Leiter gegen CODATA 2022: g = (b₁ log 6π)^(1/3) gibt +5567 ppm, der Zusatz
  c₃/3 gibt −219 ppm, die exakte Kubik mit φ = 1/(6π) gibt −252 ppm, mit
  48 c₃⁴ sind es +3,67 ppm, mit exp(−2α) −0,037 ppm, mit (1 − q)^(−5/4)
  +0,0003 ppm.
- Empfindlichkeit: d ln α / d ln b₁ = 0,331 und d ln α / d ln c₃ = 1,89.
  Andere Feldinhalte geben 135,9 (zwei Higgs), 138,2 (kein Higgs), 124,8 (vier
  Familien) und 156,1 (zwei Familien).
- Alles unterhalb von 3,7 ppm steckt in einer Zahl: Ohne Dämpfung bräuchte
  CODATA das Gewicht 47,31150(3) c₃⁴ statt 48 c₃⁴. Dämpfung und Resummation
  liefern zusammen 47,31155.

Einordnung: Die Kubik mit c₃, b₁ und der Basis 1/(6π) trifft auf 252 ppm, mit
dem Gewicht 48 auf 3,7 ppm. Diese Stufen standen schon in Fassung 1.0. Die
letzten vier Größenordnungen hängen an der einen Zahl 47,3115 und damit an
Dämpfung und Exponent.

Atlas: Die Ergebnisse RADICALHIST, RADICALCHRONO, RADICALTWODEFECT,
RADICALSEPARATION und RADICALPROSPECT stehen seit dem 18. November 2026 als
vier Befunde und eine Strategie in `problem:alpha-fixed-input-root`
(Atlasrevision c564aadd0af506c8f73fdc34). Die Torusergebnisse RADICAL,
RADICALCHARGED und RADICALBRIDGE sowie RADICALNORMAL sind nicht eingetragen.

## Nachtrag RADICALTHREE20261118: Die drei offenen Stücke

Reproduktion: `python3 -B check_three.py` (Exit 0, 15 Prüfungen bestanden).

### Stück 1: Form der Kubik

- Die Urform der Fassung 1.0 lautet α³ − A α² − A c₃² κ = 0 mit
  A = 1/(256π³) und κ = (b₁/2π) log(1/φ₀) = 1,9147. κ ist der Betrag, um den
  die inverse Hyperladungskopplung in einer Schleife über das Skalenverhältnis
  1/φ₀ läuft. Das ist dieselbe Gleichung wie heute.
- Fassung 2.7 nennt die Koeffizienten des Potentials selbst „postulated
  input" in einem erklärten Schema. Eine Herleitung steht in keiner Quelle.
- Suchraum auf dem Stand der Fassung 1.0, ohne spätere Zutaten: 8064
  Nachbargleichungen (Potenz 2 bis 4, Vorfaktor des zweiten Terms 0, 1, 2, 4,
  Potenz von c₃ von 2 bis 8, acht Zweierpotenzen, drei Potenzen des
  Logarithmus, zwei Normierungen von b₁, zwei Werte von φ). 26 liegen
  innerhalb von 1 % des Messwerts. Die Zufallserwartung für einen Treffer
  innerhalb von 3,7 ppm ist 0,0096. Genau eine Gleichung trifft, die Kubik
  der Fassung 1.0.

Urteil: Die Form ist postuliert. Auf dem Stand vor allen späteren Zutaten ist
ihr Treffer bei 3,7 ppm in diesem Suchraum aber kein Zufallsprodukt
(Erwartung etwa 1 %). Grenze: Den Suchraum habe ich in Kenntnis der Form
gewählt, und ob vor Fassung 1.0 Varianten probiert wurden, ist nicht geprüft.

### Stück 2: Regel für die D₅ Auslese

- Jede Cartanrichtung einer Spinorwurzel von E₈ trägt genau 1/4 zur Norm
  bei. Liest die Rückkopplung n Richtungen, ist der Exponent n/4. Das „g/4"
  der Fassung 2.5 ist n = 5, ihr „g = 4" ist n = 4.
- Züge gegen CODATA 2022 für n = 0 bis 8: −242, −193, −144, −96, −47, +1,9,
  +51, +99, +148 σ. Nur n = 5 trifft.
- Benannte Möglichkeiten: nur elektrische Ladung (n = 1), Isospin und
  Hyperladung (n = 2), Eichrang des Standardmodells (n = 4), Rang des Trägers
  D₅ mit B−L (n = 5), ganz E₈ (n = 8).
- Die Trägerzahl 5 stand schon vor dem Koeffizienten 5/4 fest: Fassung 2.4
  setzt γ(0) = 5/6 aus SU(5).

Urteil: Die Regel lässt sich auf eine ganze Zahl bringen. Die Rückkopplung
liest alle fünf Richtungen des Trägers einschließlich B−L und keine der drei
Familienrichtungen. Der Eichrang 4 des Standardmodells ist mit 47 σ
ausgeschlossen. Warum B−L mitzählt, ist nicht hergeleitet.

### Stück 3: Gewicht 48 c₃⁴ und Dämpfung

- 48 = 16 · 3: ein Spinor von D₅ mit 16 Zuständen mal drei der vier
  Gewichte von A₃. So steht es auch im vorhandenen Atlasbefund.
- exp(−2α) = [exp(−e² c₃)]⁴ = exp(−α |β|²) mit der Wurzelnorm |β|² = 2,
  und 2α = g/(16π²). Damit ist q = 48 [c₃ exp(−e² c₃)]⁴.
- Die Dämpfung wirkt nur auf den Vierfachterm: Dämpft man auch die Basis
  (4/3) c₃, verfehlt die Wurzel um mehr als 10⁶ σ.
- Die Dämpfung ist für alle 48 Zustände gleich: Eine von der elektrischen
  Ladung oder der Hyperladung abhängige Dämpfung verfehlt um mehr als
  10⁴ σ.
- Dämpfung exp(−s α) mit s = 3/4, 5/4, 2, 4: nur s = 2 trifft. Die Dämpfung
  benutzt die ganze Wurzelnorm, der Oszillatorfaktor nur den Anteil in D₅.
- Das Gewicht 48 c₃⁴ steht schon in den Fassungen 1.0 und 2.2.

Urteil: Gewicht, Dämpfung und Exponent lassen sich in eine Aussage fassen:
φ − (4/3) c₃ = x^h (1 − x)^(−2 h₅) mit x = 48 c₃⁴ exp(−2α), h = 1 und
h₅ = 5/8. Die Dämpfung ist keine elektrische Selbstenergie der einzelnen
Zustände. Nicht hergeleitet sind die Zeit 2α und die Potenz c₃⁴.

### Stand nach den drei Stücken

Hergeleitet ist keines der drei Stücke. Eingegrenzt sind sie auf: eine
postulierte Kubik mit Suchraumerwartung 0,0096, eine ganze Zahl n = 5 und
eine Zeit 2α. Offen bleibt ein Operator, der die Zeit am ganzen Klebevektor
und die Oszillatorgrade nur an den fünf Trägerrichtungen liest.
