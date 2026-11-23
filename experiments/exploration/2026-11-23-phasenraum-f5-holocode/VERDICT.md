# Urteil

Gerechnet am 23. November 2026 mit `.venv/bin/python` (numpy, scipy, sympy). Fünf Skripte,
alle mit Rückgabewert 0; Ausgaben in `ausgabe_*.txt`. PLAN, progress, Ledger und Atlas sind
unberührt. Ausgangspunkt war der konsolidierte Stand in `docs/research/TFPT-STAND.md`
(Abschnitt „Befunde 20.–23. November“): Die offenen Eingaben laufen auf „woher die Vier“ und
„warum minimal“ hinaus.

## Übersicht

| Nr | Urteil | Stufe | Erste freie Wahl |
|---|---|---|---|
| W Weil-Darstellung von SL(2, 𝔽₅) | **Kandidat** | 2 für P2 unter benannter Lesart | Auswahlregel „fermionischer Teil = treues SU(2) Dublett“; Lesart „fünf Plätze = Ortsbasis eines Fünfersystems“ |
| W2 Clifford-Gruppe eines Fünfersystems | Uhrenrelation **bestätigt**, sechs Fasern **verworfen** | – | – |
| W3 die 16 unter 2I | **bestätigt**, exakt | – | – |
| H minimaler holografischer Code | **offen** | 1, an der Basislinie | Abstand 3; Holografie |
| P Poincaré-Sphäre | bereits in TFPT (v219), kein neuer Weg | – | – |
| T c₃ als Bekenstein-Viertel mal Unruh-Periode | bereits in TFPT, kein neuer Weg | – | – |
| J Jordanalgebra J₃(𝕆), R Renormierung | **offen**, nicht gerechnet | – | siehe IDEAS.md |
| F Familientest zu W (Nachtrag) | Kopienweg **verworfen**; Wörterbuch über 𝔽₅ **offen** | 1 für die Zahl 3 | siehe Abschnitt F |

Keine vollständige Lösung. Was sich ändert: Die Vier braucht in W keine Eingabe mehr, und die
Fünf folgt aus einer einzigen Auswahlregel statt aus d = 4, ℤ₄-Anomalie und Minimalität.

## W: Weil-Darstellung von SL(2, 𝔽_p) (`w1_weil_darstellung.py`)

**Befund (Gruppen exakt erzeugt, Matrizen numerisch auf 10⁻⁹).**
- Für p = 3, 5, 7, 11, 13 erzeugen Fouriertransformation S und Chirp T genau SL(2, p) der
  Ordnung p(p²−1); die Darstellung ist linear.
- Das zentrale −I wirkt als ± Parität und zerlegt ℂ^p in einen bosonischen (−I ↦ +1) und einen
  fermionischen Teil (−I ↦ −1).
- **Nur bei p = 5** ist der fermionische Teil ein treues SU(2) Dublett (Dimension 2, Kern
  trivial, det = 1, Frobenius–Schur −1) **und** der bosonische Teil eine treue reelle
  3D-Drehung von SL(2,5)/±1 = A₅ (Dimension 3, Kern {±I}, det = 1, Frobenius–Schur +1).
  Bei p = 3 hat das Dublett det ≠ 1 und der bosonische Teil Dimension 1; ab p = 7 ist der
  fermionische Teil größer als 2. Beide Weil-Darstellungen bei p = 5 (additiver Charakter mit
  Rest oder Nichtrest) geben dasselbe Bild.
- SL(2,5) ≅ 2I (binäre Ikosaedergruppe), also dieselbe Gruppe, deren 120 Einheitsquaternionen
  E₈ erzeugen und deren McKay-Graph Ê₈ mit Coxeterzahl 30 ist (in TFPT bekannt).
- **Hyperladung erzwungen.** Der Kommutant der ganzen Darstellung auf ℂ⁵ hat Dimension 2. Die
  einzige spurfreie Größe darin ist Y ∝ I − 5P (P = S² = Parität), mit Skala 1/12 also −1/3 auf
  dem bosonischen Dreier und +1/2 auf dem fermionischen Zweier. Das ist TFPTs Y_par.
- **Vierteldrehung ohne Eingabe.** S ist die Fouriertransformation des diskreten Phasenraums,
  S⁴ = I für jedes p, S² = Parität. Ihr Spektrum ist {1, −1, −1, i, −i}, also bis auf das globale
  Vorzeichen das Spektrum {1, 1, −1, i, −i} von TFPTs kohärentem Vierteloperator R₅
  (Farbe 1, 1, −1; schwach ±i).
- **Kein Spin-1 Vektor des Dubletts.** Der bosonische Dreier ist nicht Sym² des Dubletts, sondern
  Sym² des galoiskonjugierten Dubletts (√5 → −√5), also das andere Ikosaedertriplett 3′.

**Basislinie.** Die Teilung ((p+1)/2, (p−1)/2) ist bei jedem p die balancierteste; das Verhältnis
−2/3 der Hyperladungen folgt aus 3 + 2 und Spurfreiheit. Informativ ist allein, dass eine
einfache physikalische Regel p = 5 auswählt und dass dann der fermionische Teil das schwache
Dublett ist. Zwei natürliche Regeln („Spinorteil ist ein Dublett“, „Bosonteil ist eine reelle
3D-Drehung“) wählen beide p = 5; einfache Gegenregeln („Bosonteil eindimensional“,
„Spinorteil dreidimensional“) wählen p = 3 bzw. keine Primzahl. Ein strukturelles Ergebnis,
kein Zahlentreffer.

**Stärkster Einwand.** Es entstehen endliche Gruppen, keine Eichgruppen: 2I ⊂ SO(3) × SU(2) und
der Kommutant U(1) × U(1), nicht SU(3) × SU(2). Der Dreier ist ein reelles Ikosaedertriplett,
kein komplexes Farbtriplett. Und die Auswahlregel ist eine Wahl: Warum soll der fermionische
Teil eines endlichen Oszillators das schwache Dublett sein?

**Abgleich mit Atlas und Prüfkatalog.** Die Permutation x ↦ 2x auf den fünf Plätzen, Y_par =
(I − 5G²)/12 und PGL(2,5) als Uhrengruppe stehen seit dem 11. November bzw. Oktober in TFPT
(`TFPT-CLOCK-LADUNG-TRANSFER-20261111.txt`, Uhrenbündel 20261028); 2I, E₈ und McKay seit v219.
Neu sind die Weil-Darstellung als gemeinsamer Ursprung, die Auswahlregel für p = 5, der
Kommutant als Begründung von Y und die Zerlegung der 16 (W3). **Benannte Identifikations-
änderung:** W übernimmt Y_par, also die Teilung in gerade und ungerade Funktionen auf 𝔽₅, nicht
die ursprüngliche Ladung, die auf markierten Platzkoordinaten diagonal ist. Laut Befund vom
11. November erhält der dort vorgeschlagene Operator H* zwar Y_par, aber nicht die markierte
Ladung (‖[H*, Y_original]‖² = 19/12 bei der Standardaufteilung, 17/6 bei den mit G²
verträglichen Aufteilungen); der Wechsel ändert also die gemeinsame Wirkung. W behauptet,
dass die physische Farb/Schwach-Basis die Paritätsbasis des Fünfersystems ist; das ist zu
prüfen, nicht zu übergehen.

**Trägt / liefert nicht / bricht.** Trägt: fünf Plätze, 3 + 2 mit kanonischer Markierung, Y,
μ₄ (als Fourier-Vierteldrehung), die 16 mit Hyperladungen, 2I und damit E₈ und h = 30.
Liefert nicht: drei Familien, D₅ ⊕ A₃ mit ℤ₄ Glue, c₃, φ₀, 41, α, Eichgruppen, Raumzeit,
Gravitation. Bricht: nichts, ersetzt aber die markierte Koordinatenladung durch Y_par.

## W2: Clifford-Gruppe eines Fünfersystems (`w2_qudit5_clifford.py`)

- Die Dilatationen x ↦ 2x und x ↦ 3x liegen bis auf das Vorzeichen −1 im Weil-Bild; die
  Parität mit Vorzeichen +1.
- Verschiebung X (TFPT C₅) und Dilatation G erfüllen G X G⁻¹ = X², G⁴ = I und erzeugen
  AGL(1,5) = F₂₀. **TFPTs Uhrenrelation ist die affine Gruppe des Körpers 𝔽₅**, also der
  affine Teil der Clifford-Gruppe. Y kommutiert mit SL(2,5) und G, nicht mit X.
- **Verworfen:** Die sechs gegenseitig unverzerrten Basen eines Fünfersystems werden von der
  Clifford-Gruppe nur als PSL(2,5) ≅ A₅ (Ordnung 60, Stabilisator 10) permutiert, auch mit
  antiunitärer Konjugation. TFPTs sechs Fasern tragen PGL(2,5) mit Stabilisator F₂₀. Die
  Lesart „sechs Fasern = sechs Basen“ fehlt ein Element mit Nichtrest-Multiplikator.

## W3: die 16 unter 2I (`w3_sechzehn_unter_2I.py`)

Auf Λ^gerade ℂ⁵ = 16 ist **jeder Hyperladungssektor eine irreduzible Darstellung von 2I**:
ν^c 1, u^c 3′, Q 6, e^c 1, L 2, d^c 3′. Q ist die sechsdimensionale irreduzible Darstellung
(der Verzweigungsknoten von Ê₈). Das liegt am galoiskonjugierten Dreier: Wäre der Dreier der
Spin-1 Vektor des Dubletts, zerfiele Q als Spin 1 ⊗ Spin ½ = 2 ⊕ 4. Damit spricht W3 gegen die
Lesart „drei Farbplätze = drei Raumrichtungen“ (der stärkste Einwand gegen A1 vom 20.11.) und
für einen inneren Dreier.

## H: minimaler holografischer Code (`h1_holocode.py`)

**Befund.** [[5,1,3]] ist ein perfekter Sechsbein-Tensor (alle 20 Dreierreduktionen = I/8);
vier Qubits erreichen bei der Reinheitssumme nur 1 statt 0,75 (kein perfekter Vierbein-Tensor).
Quanten-Singleton gibt n = 5 für eine logische Einheit mit Abstand 3; das kleinste hyperbolische
q für Fünfecke ist 4. Die Dreiecksgruppe (2,4,5) hat als kleinstes glattes endliches Bild S₅
(Bring-Fläche, Geschlecht 4, 30 Ecken, 60 Kanten, 24 Fünfecke). In allen 120 erzeugenden
Paaren ist a eine Transposition (markiert 3 + 2 der Plätze), b² hat auf ℂ⁵ die Eigenräume 3 + 2.

**Basislinie.** Andere Signaturen liefern ebenso TFPT-nahe Gruppen (kleinste gefundene Bilder
in transitiven Darstellungen bis Grad 8): (2,3,7) → PSL(2,7) ≅ GL(3,2) (Hamming und Fano),
(2,5,5) → A₅, (2,4,6) und (2,3,8) → Ordnung 48. Für (2,4,5) ist S₅ beweisbar minimal: Gruppen der
Ordnung 40 oder 80 haben keine glatte (2,4,5)-Erzeugung (Sylow-Argument). Mit Abstand 2 statt 3 entsteht
(2,3,7) statt (2,4,5). S₅ ist also kein seltener Treffer.

**Stärkster Einwand.** Die Uhrenrelation G C₅ G⁻¹ = C₅² wird in keinem der 120 Paare erfüllt
(b normalisiert ⟨ab⟩ nie). Der Innenraum ist zweidimensional, liefert also 2 + 1, nicht 3 + 1.
Die Gravitationsseite (Ryu–Takayanagi) ist nicht gerechnet.

## P, T, J, R

P und T stehen bereits in TFPT (Poincaré-Sphäre und McKay in v219; Einstein 8π = 1/c₃ mit
Bekenstein 1/4 = 1/|μ₄|). J und R sind nicht gerechnet; ihre kleinsten Tests stehen in IDEAS.md.

## Was jetzt noch Eingabe ist (aus Sicht von W)

1. Die Auswahlregel für p = 5 (fermionischer Teil = treues SU(2) Dublett).
2. Die Lesart „fünf Plätze = Ortsbasis“ samt Übernahme von Y_par statt der markierten Ladung.
3. Familienzahl 3, D₅ ⊕ A₃, c₃, φ₀ und alles Dynamische wie bisher.

## F: Familientest (Nachtrag am selben Tag, `f1_familien.py`)

**A. Mehrere Fünfersysteme (der angekündigte Test).** Weil-Darstellung von Sp(2n, 5) auf
ℂ^(5ⁿ), zerlegt nach dem Howe-Paar aus diagonalem 2I und O(n, 5); O(n,5) vertauscht exakt mit
der Weil-Wirkung.
- n = 1: Dublett einmal, Partner das Vorzeichen von O(1) (das ist W).
- n = 2 (bosonisch 13, fermionisch 12): **kein** Dublett, weder 2 noch 2′.
- n = 3 (bosonisch 63, fermionisch 62): Dublett 2 **fünfmal**, sein Multiplizitätsraum ist eine
  irreduzible Fünferdarstellung von O(3,5) ≅ PGL(2,5) × ℤ₂; 2′ kommt nicht vor.

**Urteil A: verworfen.** Die Abbruchregel greift: Weder n = 2 noch n = 3 liefert drei Kopien
oder ein Triplett, und D₅ ⊕ A₃ erscheint nicht. Nebenbefund: Bei n = 3 trägt das Dublett
wieder eine Fünf. Dass der bosonische Teil bei n = 3 die Dimension 63 hat (TFPT: 68 = 5 + 63
Gedächtnisrichtungen), ist ein ungeprüfter Zufall; 63 = 2⁶ − 1 ist ebenso naheliegend.

**B. Wörterbuch über 𝔽₅.** Exakt geprüft:
- Die Nullstellen von z⁴ − 1 in 𝔽₅ sind **ganz 𝔽₅^×**: die vier Marken sind die Einheiten.
- Polynome vom Grad ≤ 4 sind über 𝔽₅ genau alle Funktionen auf 𝔽₅ (3125 von 3125): die fünf
  Schnitte von O(D) entsprechen den fünf Punkten der affinen Geraden.
- Die Dilatation x ↦ 2x (TFPTs G) hat auf ℂ^𝔽₅ den invarianten Teil span(δ₀, 1_{𝔽₅^×}) und die
  drei nichttrivialen Dirichlet-Charaktere mod 5 mit Eigenwerten −i, −1, +i. Das sind genau die
  μ₄-Charaktere {−1, i, −i}, die der Atlas dem Familienraum H¹ zuschreibt (Befund vom 23.11.).
- Die sechs Fasern sind P¹(𝔽₅) (in TFPT seit Oktober bekannt), ihr Punktstabilisator
  AGL(1,5) = F₂₀ ist die Uhrengruppe ⟨C₅, G⟩.
- Gauss-Summen: |g(χ)|² = 5 für alle drei, die Legendre-Summe ist √5 (goldener Schnitt, Ikosiane).
- Auswahl: „Die Einheitengruppe ist genau μ₄“ gilt nur für p = 5. Das ist eine zweite,
  unabhängige Regel neben „Spinorteil = Dublett“; beide wählen p = 5.

Damit ist das Zahlengerüst 3 Familien, 4 Marken, 5 Plätze, 6 Fasern dasselbe wie p − 2, p − 1,
p, p + 1 auf der projektiven Geraden über 𝔽₅, und die Familiencharaktere sind die nicht-
trivialen Einheitencharaktere.

**Basislinie zu B.** In TFPT sind 3, 4, 5 schon n − 1, n, n + 1 für n = 4 Marken
(Riemann–Roch). Das Wörterbuch ersetzt die Frage „warum vier Marken“ durch „warum 𝔽₅“; neu ist
die Lesart der Marken als Einheiten und der Familiencharaktere als Dirichlet-Charaktere, nicht
die Zahlen selbst.

**C. Familien als getwistete Sektoren.** Induziert man auf den sechs Punkten von P¹(𝔽₅) mit
jedem Charakter von 𝔽₅^×, entstehen: trivial 1 ⊕ 5 (Steinberg), Legendre 3 ⊕ 3′ (bosonisch),
die beiden quartischen Charaktere jeweils dieselbe irreduzible fermionische 6.
Die drei Familiensektoren sind also **nicht gleich**: Unter der Galoiswirkung i ↦ −i zerfallen
die Familiencharaktere in 1 + 2 (Legendre reell, quartisches Paar konjugiert).

**Urteil F.** Kopienweg verworfen. Das Wörterbuch ist **offen**, Stufe 1 für die Zahl 3 und die
Familiencharaktere, aber ohne Mechanismus für drei gleiche Kopien der 16. Stärkster Einwand:
In TFPT sind Familien eine eigene Kohomologie H¹, hier liegen die Familiencharaktere im selben
ℂ⁵ wie die Plätze; und die drei Sektoren sind verschieden statt gleich. Die Zerlegung 1 + 2 der
Familien unter Galois ist eine prüfbare Vorhersage gegen die TFPT-Flavorstruktur, kein Befund.

## Kleinste entscheidende Tests für eine Folgesitzung

1. **Familien als Kohomologie über 𝔽₅.** Die étale Kohomologie H¹ von P¹ ohne μ₄ über dem
   algebraischen Abschluss von 𝔽₅ mit der μ₄-Wirkung berechnen und mit dem TFPT Familienraum
   (Charaktere {−1, i, −i}, Rang 3) vergleichen. Abbruch, wenn Rang oder Charaktere abweichen.
2. **Galois 1 + 2 gegen Flavor.** Prüfen, ob eine TFPT Flavorgröße (Leptonkoeffizienten
   16/7, 4/3, 7/6, Transportlängen) eine Familie gegen zwei auszeichnet, und ob es die zum
   Legendre-Charakter passende ist.
3. **Sechs Fasern.** Welches natürliche Element mit Nichtrest-Multiplikator (etwa J aus der
   Standardform kombiniert mit einer Dilatation) erweitert A₅ auf PGL(2,5) mit Stabilisator F₂₀?
4. **Y_par gegen die markierte Ladung.** Den Befund ‖[H*, Y]‖² = 19/12 gegen 17/6 mit der Weil-
   Lesart neu bewerten: Welche Dynamik ist mit der Paritätsbasis verträglich?
5. **Fockspur-Anschluss.** Wählt die positive Fockspur über Λ ℂ^(𝔽₅) mit der Weil-induzierten
   quadratischen Quelle (Exploration F1 vom 23.11.) dieselbe Orientierung Σ ∝ −6Y?
