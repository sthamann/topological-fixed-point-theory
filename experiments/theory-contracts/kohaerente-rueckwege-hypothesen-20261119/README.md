# RUECKWEGE20261119: Prüfung von dreizehn Hypothesen zu kohärenten Rückwegen

Datum 19. November 2026. Atlasrevision `1823e3ca9e529dc86c8e8c0a`. Lokaler Stand gleich `origin/main`.

Geprüft wurden die sieben Hypothesen des ersten eingereichten Textes (H1 bis H7) und die sechs des
zweiten (A bis F). Jede Hypothese bekam die kleinste Rechnung, die sie stützen oder treffen kann.
Sieben Skripte, 89 Prüfzeilen, alle bestanden (darunter einige reine Urteilszeilen ohne eigene Rechnung). Nichts befördert, kein Ledgerwechsel, kein Commit.

## Urteile

| Hypothese | Urteil | Kern der Rechnung |
|---|---|---|
| H1 neutrale Vierwurzelkohärenz trägt die Rückwirkung | **Teil ausgeschlossen, Teil bestätigt** | Die Kohärenz ist der einzige neutrale Zweierblock der Faser. Die 600 Farb und Schwachmarken sehen sie nicht. |
| H1/C ein Zwischenoperator, verschiedene Kopplungen | **im erklärten Menü ausgeschlossen** | 0 von 32 795 Kopplungspaaren treffen Farbe und Schwach gemeinsam. |
| H2 geschützte Stromauslesung für exp(−2α) | **bedingt, kein neuer Herkunftsschritt** | Gaußform folgt im Vakuum von selbst. Die Varianz liefert die Stromalgebra nicht. |
| H3 ein gemeinsames Quellenfunktional | **im endlichen Chart ohne Auswahlkraft** | Beim flachen Zustand ist die Prüfung identisch erfüllt. |
| H4 Prozess mit Gedächtnis | **als Methode gestützt, Originaltest nicht ausführbar** | Geordnete Sonden heben den Rang von 12 auf 32. Zwei Prozesse bleiben trotzdem übrig. |
| H5 P1 als Fixpunkt der Verfeinerung | **nicht entschieden** | Im Vergleichsmodell hängt die Einbettung von der Energie ab. |
| H6 dritte Raumrichtung aus Rückwegen | **Modellrechnung bestätigt, Prüfkriterium zu schwach** | 0,8σ und Rang drei stimmen. Eine Kontrolle ohne σ_z hat denselben Rang und keinen Kegel. |
| H7 kausales Wirkungsprinzip | **als Ableitung aus dem Bestand gescheitert, als Erweiterung offen** | Mit positiven Formen gibt es keine raumartige Trennung. |
| A 60er Uhr auf zwei E₈ Kopien | **bestätigt und erweitert** | χ_U = Φ₆₀ exakt. Alle 16 Galoiskonjugationen sind Gitterisometrien. |
| B flacher sichtbarer Zustand mit Gedächtnis | **bedingt** | Modell stimmt. Die aufgezeichneten Rückkehrwerte sehen das Gedächtnis, wenn sie auf dem ganzen Raum gebildet werden. |
| C vierte Familie als virtueller Rückkehrraum | **Zählung bestätigt, Mechanismus im Menü ausgeschlossen** | 480 von 1440 Übergängen treten aus (ein Drittel). Gemeinsamer Zwischenoperator siehe H1/C. |
| D drei Raumrichtungen aus drei Paarungen | **in der gestellten Form nicht gestützt** | Die drei Paarungen geben drei orthogonale Drehachsen, keine vertauschenden Verschiebungen. |
| E P1 als Grenzfläche zweier Phasen | **algebraisch eingegrenzt** | Auf E₈ ⊕ E₈(−1) dreht U das Vorzeichen der Form um. Der Spiegel ist der Tausch, nicht U. |
| F Schnittmenge der Bedingungen wählt die Quelle | **Eindeutigkeit im endlichen Chart widerlegt** | Zwei verschiedene Prozesse erfüllen alle geprüften Bedingungen. |

## Was neu ist

**1. Die neutrale Kohärenz hat eine Adresse** (`check_coherence.py`, 13/13).
Unter Standardmodell mal Familien SU(3) zerfallen die 128 Zustände in 24 Besetzungsklassen. Genau ein Typ kommt
doppelt vor: das vollständige Singulett, auf den Zuständen 0 (alles leer) und 255 (alles voll). Die Kommutante hat
Dimension 26 = 22 + 2². Der Block verbindet das rechtshändige Neutrino der ausgelassenen Familie mit seinem
Konjugierten. Die Kohärenz ändert die Trägerzahl um 5 und die Familienzahl um 3. Nimmt man die Trägerzahl als
Ladung hinzu, verschwindet sie (Kommutante 24).

- Kürzestes Wort: vier Wurzeln mit allen geraden Wurzeln, fünf mit den 60 gemischten allein.
- Alle 600 Marken ZY und ZUY (60 Ströme, fünf äußere Auswahlen) bleiben bei eingeschalteter Kohärenz gleich.
  Sie kann den Unterschied zwischen Farb und Schwachantwort daher nicht tragen.
- Die Komplementpaarung macht sie reell. Frei bleibt ein Vorzeichen. Das passt zu den aufgezeichneten Zählern
  +1/2 und −1/2. Positivität begrenzt nur den Betrag.

**2. Ein gemeinsamer Zwischenoperator scheitert im Menü** (`check_common_middle.py`, 6/6).
Ansatz: W_Typ(t) = δ + Σ_β G(t + β), β aus einer Kopplungsmenge K_Typ, G ≥ 0 Funktion der Besetzungsklasse,
δ ≥ 0 gemeinsam. Kopplungsmengen: Vereinigungen von höchstens vier der zwölf Eichklassen gerader Wurzeln,
Einheitskopplung. Kontrollen: zwei getrennte Blöcke treffen alles; ein gemeinsamer Block verfehlt auf der vollen
Spur um genau 7/528 (wie aufgezeichnet), bei allen fünf Auswahlen zugleich um 0,0204. Ergebnis: kein Paar trifft,
kleinster Rest 6,2 · 10⁻⁴. Gleitkommarechnung mit HiGHS, kein exakter Beweis.

**3. Die 60er Uhr trägt ihre ganze Galoisgruppe** (`check_clock60.py`, 17/17).
U = [[0, C], [1, 0]] hat χ_U = Φ₃₀(x²) = Φ₆₀ und erhält E₈ ⊕ E₈. Für jede der 16 Einheiten k modulo 60 gibt es eine
Gitterisometrie G_k = diag(g, C^((k−1)/2) g) mit G_k U G_k⁻¹ = U^k. G₇ wirkt mit Ordnung vier auf den Fünfertakt.
Die zweite Gegenprobe der Hypothese ist damit bestanden. Grenzen:

- Kein D₅ + A₃ Rahmen bleibt unter C^k erhalten außer bei C¹⁵ = −1. Fünfertakt und Dreiertakt bewegen jede feste
  Teilung in Träger und Familie.
- Die Verdopplung ist allgemein: aus den Uhren 16, 20, 24, 30 werden primitive Uhren 32, 40, 48, 60. Nichts wählt 60.
- Auf E₈ ⊕ E₈(−1) gilt Uᵀ G U = −G. Keine ungerade Potenz ist eine Involution, U¹⁵ quadriert zu −1.
  Die Tauschform auf E₈ ⊕ E₈ hat Signatur (8, 8).

**4. Die Alpha Dämpfung bleibt eine Normierungsfrage** (`check_current_readout.py`, 9/9).
Im Vakuum liegt der Vektor im Fockmodul der Ladung null; jeder verschmierte Strom ist exakt gaußisch und der
Nullmodenanteil ist 1. Die Varianz ist β² Σ n |f_n|², also eine Zahl der Testfunktion. Verlangt ist 1/(2π).
Kanonische Normierung gäbe exp(−4πα), 122 000 σ daneben. Bei τ = i trägt der Oszillator coth(π): −87,9 σ.
Die reine Nullmode bei τ = i gibt +141,6 σ (aufgezeichnet +141,7). H2 fällt damit auf den schon geführten
bedingten Befund zum Fourier selbstdualen Grundzustand zurück.

**5. Das Gedächtnismodell und der Satz zur flachen Dichte** (`check_memory_functional.py`, 13/13).
ρ = (1 + ε Z⊗Z)/4 ist sichtbar flach, die Zweizeitantwort ist cos 2κt + iε sin 2κt. Der Rückkehrwert mit Referenz
Tr ρ² auf dem erweiterten Raum ist (1 − ε²)/(1 + ε²). Allgemein gilt 1 − Rückkehr = ‖[ρ, X]‖²/(2 Tr ρ²). Im Chart
128 × 2 ändern alle drei getesteten Kopplungen die Anker 1/4. Hypothese B hält nur, wenn die Rückkehrwerte allein
auf der sichtbaren Algebra definiert sind.

**6. Der Kegel** (`check_cone.py`, 17/17).
Die Zahlen der Vergleichsrechnung stimmen exakt, die Händigkeit ist umgekehrt (V_x V_y V_z = −0,512 i). Mit
B = g + v sin k_z ohne σ_z hat (1/2) tr V_i V_j ebenfalls Rang drei, die Dispersion ist aber
ω = −0,8 k_z ± 0,8 |k_⊥|. Der Test muss {V_i, V_j} = 2 G_ij verlangen. Der einzige neutrale Zweierblock der lokalen
Faser ist der Block aus Punkt 1; seine σ_x und σ_y verletzen die Trägerzahl.

**7. Marken und Geometrie** (`check_geometry.py`, 14/14).
Die drei Paarungen sind drei spurlose Möbiusmatrizen, die für jedes Doppelverhältnis paarweise antikommutieren.
Sie spannen sl(2, ℂ) = so(3, 1) auf und haben einen gemeinsamen Fixpunkt. Das ist ein Lorentzrahmen für beliebige
vier Punkte, keine Auswahl und keine Verschiebung. Für kausale Fermionsysteme: bei positiver Spinform ist das
Spektrum jeder geschlossenen Kette reell und nichtnegativ (2000 Zufallsfälle), bei skalaren Auswertungen ebenfalls
reell. Erst das Paar aus Punkt und Spiegelpunkt mit der Tauschform gibt beide Typen (2630 raumartig, 1370 zeitartig
von 4000). Kern und Maß sind dabei gewählt.

## Nachtrag CARTANDEPHASE20261119: ein gemeinsamer Prozess statt eines gemeinsamen Operators

`check_dephasing.py`, 16/16. Ausgangspunkt ist das Hindernis aus Punkt 2: Kein gemeinsamer positiver Operator
zwischen den Strömen trifft Farbe und Schwach. Ein gemeinsamer vollständig positiver Prozess trifft beide exakt.

- **Der Prozess.** Dephasierung durch ein Cartanelement h: L_h(A) = h A h − (h² A + A h²)/2, eine Zeit für alle
  Ströme. Jeder gemischte Strom trägt unter der Hyperladung eine feste Ladung q = 2 a_i.
- **Die Marken.** Tr(h X* h X) − (früh + spät)/2 = −(q²/2) · ZY in allen 300 Fällen (60 Ströme, fünf Auswahlen).
  Das sind genau die aufgezeichneten Mischwerte −1/72 (Farbe) und −1/32 (Schwach) und die Ziele der beiden
  angepassten Blöcke. Bei endlicher Zeit fällt die Rückkehr eines Farbstroms wie exp(−τ/18), eines schwachen wie
  exp(−τ/8). Der Unterschied ist das Ladungsquadrat 1/9 gegen 1/4, kein zweiter Prozess.
- **Verträglich mit der flachen Dichte.** Der Prozess ist unital, der flache Zustand bleibt stehen, Besetzungen
  ändern sich nie. Er gibt trotzdem eine gerichtete Zeit.
- **Die neutrale Kohärenz.** Unter Dephasierung durch alle Eich Cartanelemente (Farbe, Schwach, Hyperladung,
  Familien SU(3)) überlebt von allen Kohärenzen der 128 Zustände genau eine: die zwischen 0 und 255.
  Trägerzahl löscht sie mit Rate 25/2, der volle Cartantorus mit Rate 4.
- **Anschluss an Alpha.** Dephasierung über den vollen Cartantorus dämpft alle 112 Wurzelströme gleich mit
  exp(−τ): die Rate ist |β|²/2, das konforme Gewicht 1 (Nullmodenanteil von L₀). Mit τ = 2α ist das der
  aufgezeichnete zustandsunabhängige Faktor exp(−2α). Der Hyperladungsanteil am Gewicht ist 1/15 (Farbe) und
  3/20 (Schwach).

Grenzen. Die Gleichung für die Marken ist eine Identität: die Zielwerte waren als zweiseitige Cartanwerte
definiert. Neu ist, dass ein einziger vollständig positiver Prozess sie trägt und die beiden Blöcke mit
gegenläufiger Ordnung überflüssig macht. Nicht hergeleitet: warum die Quelle dephasiert, welche Cartanrichtungen
und der Wert τ = 2α. Die Zeile zur Varianz ist Buchführung der schon festgelegten Konvention. Kein originaler
P1 Kern.

## Nachtrag 2 · KERNDEPHASE20261119: Liefert der originale Randkern die Dephasierung?

`check_kernel_dephasing.py`, 21/21. Gelesen: Randkernquelle Abschnitt 001 (`B_Σ`, Calderón Projektor,
`χ_seed = 1/λ₁⁺`, `H_int` als „analytic realization choice“).

Lesart: Ein Strom mit Gewichtsdifferenz Δw trägt die Phase exp(i e Δw·φ). φ hat in jeder inneren Richtung den
Kern als Kovarianz, e² = 4πα. Das Mittel über φ ist ein Zufallsunitärkanal mit Faktor exp(−τ |PΔw|²/2).

- **Die Zeit: ja, bedingt.** Die positive Hälfte von |B|⁻¹ hat auf der Einheits S² je Kanal das Gewicht 1/(2π).
  Ein Kanal gibt τ₁ = e²/(2π) = 2α exakt. Der Saatkanal (λ₁⁺ = 1, Entartung 2, Fläche 4π) ist genau dieser eine Kanal.
- **Halbgruppe nur auf S².** Weil das Spektrum flach ist, dephasieren N Kanäle mit exp(−N τ₁ |Δw|²/2). Die Kanalzahl
  ist die Zeit. Das gilt nur auf S² (S¹ bis S⁷ geprüft). Die Lindbladform aus Nachtrag 1 ist also das flache Spektrum.
- **Die Richtungen: alle acht, gleich stark.** Der Kern wirkt als B ⊗ 1 auf den inneren Marken. Das gibt den vollen
  Cartantorus und damit denselben Faktor exp(−τ) für alle 112 Wurzelströme, mit τ = 2α also α⁻¹ = 137,0359992168.
- **Jede engere Auswahl fällt.** Nur Hyperladung: Raten 1/15 und 3/20, z = +22547 und +20522. Nur elektrische Ladung
  auf einem Feld: höchstens exp(−α), z = +12042. Sechs Eich Cartanrichtungen: Raten 1/3 bis 1, also ungleich.
  Zwei Kanäle oder beide Vorzeichen: z = −23816. Quantenschwingung bei τ = i: z = −87,9.
- **Die Hyperladungsrichtung liefert der Kern nicht.** Der Generator des Kernprozesses gibt für alle 60 gemischten
  Ströme denselben Wert −ZY, ohne Unterschied zwischen Farbe und Schwach. Die Werte −1/72 und −1/32 sind sein Anteil
  entlang der Hyperladung: Anteil (1/15, 3/20) mal |h|² = 5/6 mal ZY. Diese Richtung kommt aus der Marke U = h² der
  Auslese, nicht aus dem Kern.
- **Folgen.** Die neutrale Kohärenz 0 gegen 255 ist unter dem Kernprozess nicht geschützt: Faktor exp(−8α) = 0,9433.
  Die Marken sind Generatorwerte bei τ = 0. Nach der endlichen Kernzeit gelesen wären sie um 1,45 Prozent kleiner.
  Marken und Alpha Dämpfung sind also zwei verschiedene Auslesen desselben Prozesses.

Urteil: **halb.** Zeit und Gleichverteilung über alle Richtungen kommen aus dem Kern, die Hyperladungsrichtung nicht.
Annahmen der Lesart, die nicht im Original stehen: die Kopplung e² = 4πα der Cartanphase an das Kernfeld, die
Einheitskugel (bei Radius R ist das Gewicht 1/(2πR)) und die runde S². Die Zahl 2α ist dieselbe wie in P1KERN20261119
und zählt nicht als zweiter Beleg. Nichts befördert.

## Nachtrag 3 · COULOMBKOPPLUNG20261119: Folgt die Kopplung e² = 4πα aus dem Original?

`check_coulomb_coupling.py`, 19/19. Gelesen: Papier v28, Lemma „Defect Sector Reweighting“; aktuelle EM Quelle,
Abschnitt 001; Randkernquelle, Verdopplungssatz.

- **Was im Original steht.** Der Exponent ist dort α · J mit J = d · SF, d = 2 Blätter des Doppelüberzugs und
  SF = 1 als Einheit des Spektralflusses. Die Einheit je Blatt ist dort gesetzt („We identify“), nicht gerechnet.
  Eine Kopplung e oder eine Feldkovarianz kommt nicht vor.
- **Die Kopplung ist das Coulombgesetz des Kerns.** Der skalare Vergleichskern bei κ = 1 ist exakt die Greenfunktion
  des Laplace Operators im R^d, eingeschränkt auf die Einheitskugel (d = 3 bis 7 geprüft). In d = 3 ist das
  1/(4π|x−y|) mit Legendre Reihe Σ P_l: jeder Multipol hat Koeffizient 1. Damit ist κ = 1 das Gaußsche Gesetz, und
  e² mal Kanalgewicht ist α, in jedem Einheitensystem. Das 4π der Kopplung hebt das 4π des Kerns. Je Kanal und Blatt
  ist der Exponent α · |β|²/2 = α · 1. Das ist genau die Einheit J₁ = SF = 1 des Originals.
- **Die Zwei.** Zwei Zählungen geben k = 2: zwei Blätter mit je einem skalaren Kanal, oder ein Blatt mit den zwei
  Spinorkomponenten der positiven Hälfte. Beides zusammen gäbe k = 4 (z = −23816). Das Original stützt, dass nicht
  multipliziert wird: Die Calderón Hälfte ist die spiegelgerade Hälfte der verdoppelten Daten, also eine Kopie.
  Welche der beiden Zählungen gemeint ist, sagt das Original nicht, und α kann es nicht trennen.
- **Radius.** In Saateinheiten ist das Gewicht skalenfrei: χ_seed² · Entartung / Fläche = 1/(2π) für jeden Radius.
  Dass die Kugel in Einheiten ihrer eigenen Saat gemessen wird, bleibt Lesart.
- **Sektoren.** Das Gesetz k = |Δw|² gibt auf den 128 Gewichten die Sektoren 2, 4, 6, 8. Das passt zur Bemerkung des
  Originals, dass SF ≠ 1 ein anderer Sektor mit anderem k ist.
- **Mehrere Defekte: neuer trennender Test.** Die festgelegte Funktion q(1−q)^(−5/4) dämpft n Defekte mit
  exp(−2nα): z = +1,9. Kohärente Addition derselben Wurzel im selben Kanal (exp(−2n²α)) gibt z = −5,1. Ein Paar,
  das zu einer Wurzel verschmilzt (exp(−2α) bei n = 2), gibt z = +5,5. Nur unabhängige Addition hält. Der Kern
  liefert sie für Defekte in verschiedenen Kanälen oder für zueinander senkrechte Wurzeln.

Urteil: **bedingt ja.** e² = 4πα ist keine eigene Eingabe mehr, sobald das Kernfeld ein Cartan Eichpotential mit
Gaußscher Normierung in drei Raumrichtungen ist. Diese Gleichsetzung steht nicht im Original; sie ist die eine
verbleibende Annahme und ersetzt drei (Kopplung, κ, Einheit je Blatt). Die 5σ Aussagen hängen an CODATA 2022 und
an der festgelegten Funktion. Nichts befördert.

## Nachtrag 4 · KERNFELD20261119: Legt das Original das Kernfeld fest? Und was legt die Coulomb Lesart sonst fest?

`check_kernel_field.py`, 20/20. Gelesen: Randkernquelle Abschnitte 001 und 002, EM Quelle 001, Trägerquelle 006.

**A. Das Original.** Nein. Der primitive Kern enthält dort ausdrücklich „no carrier split, no hypercharge operator,
no family geometry“. Eichfelder kommen erst später als „carrier-compatible bosonic fluctuations“. Das Kernfeld als
Eichpotential ist also keine Angabe des Kerns. Eine Brücke gibt es: Die Norm der Hyperladungsrichtung im Chart,
|h|² = 5/6, ist die Trägernorm γ = tr_E Y² des Originals.

**B. Eine Varianz für alle Eingaben.** Mit dem Gewicht w = 1/(4π) je Komponente (das gibt exp(−2α)):
- c₃ = w/2 und δ_top = 48 c₃⁴ = 3 w⁴.
- In einem Menü von 32 Gaußschen Lesarten trifft δ_top genau zweimal. Nur eine benutzt dieselbe Varianz wie die
  Dämpfung: **drei Paarungen von vier Marken**, jedes Paar einmal auf jeder der zwei Komponenten verbunden (w · w).
- Die Zählung ist scharf: 2 Paarungen (komplexes Feld) und 9 (unabhängig auf beiden Komponenten) liegen weit
  außerhalb, 3 liegt bei z = +1,9. Die Paarung muss auf beiden Komponenten dieselbe sein.
- Nicht erklärt: φ_base = (2/3) w. Kein Gaußsches Objekt des Menüs gibt 2/3.

**C. Die Lage der Marken.** Lässt man die vier Marken über denselben Kern 1/(4π d) wechselwirken:
- Auf dem Nahtkreis ist das Quadrat (Marken bei μ₄) das einzige Gleichgewicht (60 Zufallsstarts).
- Frei auf der Kugel ist es das regelmäßige Tetraeder, mit tieferer Energie. Das Quadrat ist dort ein Sattel.
- Quadrat: j = 1728, quadratischer Torus. Tetraeder: j = 0, Sechseck. Das sind genau die zwei Formen aus
  P1KERN20261119 (Sechseck als echtes Extremum, Quadrat nur mit Spiegelung).
- Damit wäre die Markenlage keine freie Eingabe mehr: Spiegelung plus Abstoßung über den Kern gibt μ₄.
  SWEEP20261119 hatte mit Spiegelung allein nur ein Rechteck mit freiem Seitenverhältnis.
- Muster, keine Herleitung: Die Paarungssumme des Coulombkerns am Quadrat ist exakt 5/4 (Rechteck:
  1/4 + 1/sin²(2θ)). Für den flachen Saatkanal ist dieselbe Summe 3. 5/4 ist auch der festgelegte Exponent.
  Kein Schritt macht aus einer Paarungssumme den Exponenten von (1 − q).

Urteil: **Im Original nein, als Lesart tragfähig.** Neue Annahmen der Lesart: Die vier Marken tragen gleiche
Ladung und stoßen sich über den Kern ab; die Paarung ist auf beiden Komponenten dieselbe. Offen bleiben φ_base,
der Exponent 5/4, die Kubik mit 41 und die Hyperladungsrichtung. Eine vollständige Herleitung ist das nicht.

## Grenzen

- Alle Rechnungen stehen im endlichen D₈ Chart oder in erklärten Vergleichsmodellen. Der originale P1 Kern, die
  physischen Higgs und Yukawamarken und eine regulierte Randquelle wurden nicht ausgewertet.
- Die Suche in Punkt 2 und die Kommutantenränge in Punkt 1 sind Gleitkommarechnungen. Die Zufallsproben in Punkt 7
  sind kein Beweis; die Aussage über positive Formen folgt aber aus der Ähnlichkeit zu einer positiven Matrix.
- Die im eingereichten Text genannten Arbeiten (Wassner und andere, Finster und Fischer, Dittrich, Mason und Tuite)
  wurden nicht nachgelesen.
- H4, H5 und F sind als Ganzes nicht durch eine Rechnung entscheidbar. Für H4 fehlt ein Originalmodell mit zwei
  beschränkten geladenen Eingriffen, für H5 eine regulierte Randquelle mit zwei verschiedenen Verfeinerungen.

## Wiederholen

```bash
cd experiments/theory-contracts/kohaerente-rueckwege-hypothesen-20261119 && for s in clock60 coherence common_middle current_readout memory_functional cone geometry dephasing kernel_dephasing coulomb_coupling kernel_field; do python3 -B check_$s.py; done
```

Laufzeit etwa 80 Sekunden, davon 50 für `common_middle`.
