# Urteil

## Aktueller Stand: H6h6 — Ursprung der Naht (`z7_ursprung_naht.py`)

```
        K1 … K15 K16 K17 K18 K19    P  L  F  O  E
   H6h5  P … P   L   L   L   E    14  4  0  0  1
   H6h6  P … P   L   L   L   L    14  5  0  0  0
```
K19 (Existenz der Naht) ist neu im Katalog; bei allen früheren Hypothesen war die Naht Eingabe. In
H6h6 ist nichts mehr als Eingabe gesetzt; alles steht unter benannten Prämissen.

**Hypothese G.** Die Naht ist der verdrillte Defekt eines komplexen Paarungsmassenfelds
m e^{iθ} ψχ zwischen Weyl-Fermionen ψ und Partnern χ. Um den Defekt werden ψ und χ mit der
Vierteldrehung J: ψ → χ → −ψ vertauscht, und die Massenphase windet um π. Eindeutig ist das, weil J die
Masse ψχ umpolt.

**Gerechnet (alle sieben Prüfungen bestehen).** 3D Gitter, Defektlinie entlang z, Bogoliubov-de-Gennes
mit Weyl-Kinetik und unter J invarianter Wilson-Majoranamasse gegen die Verdoppler. Kern- und
Randmoden sind über den Kernprojektor getrennt, die Chiralität kommt aus ⟨∂H/∂k_z⟩.

| Defekt | chirale Majorana-Weyl Moden am Kern |
|---|---|
| verdrillt (J, Phase windet um π), ein Paar | **1**, Geschwindigkeit +1,00 |
| verdrillt, zwei Paare | **2** |
| gewöhnlicher Wirbel (ohne J, Windung 1), ein Paar | 2 |
| nur J, ohne Massenwindung | 0 |

Der äußere Rand trägt jeweils dieselbe Zahl mit umgekehrter Richtung. Das Ergebnis ist robust bei
L = 14 und 18 und für Massen von 0,4 bis 0,9.

**Was daraus für die Naht folgt.**
- Kodimension 2: Eine Defektlinie im Raum ist eine Fläche in der Raumzeit.
- Geschlossen und flach heißt Torus (Gauß-Bonnet, χ = 0).
- J² = −1 auf Fermionen, also wirkt J auf allem Bosonischen mit Ordnung 2. Das sind die **zwei
  Blätter**; auf Fermionen ist es **μ₄**, wie die gelifte Vierteldrehung (g̃⁴ = −1, Z1).
- Der Verzweigungsschnitt ist der Ort, an dem die halbe Phase springt. Seine Lage ist eine Eichwahl;
  das ist der unbeobachtbare Schnitt aus F (Z5).
- Die chiralen Moden laufen mit Lichtgeschwindigkeit. Ihre Zahl ist die Zahl der Paare, deren Masse am
  Defekt windet; **16 Paare geben 16 Majorana-Weyl Moden, also c = 8**.
- Gegenprobe: Die 16 Verdoppler eines naiven Feldes liefern die Moden nicht. Ohne Wilson-Glied haben
  sie abwechselnde Chiralität, und der verdrillte Defekt trägt dann keine Mode. Gebraucht werden 16
  echte Paare.

**Korrektur zur Gabel (Z5).** Z5 hat nur die Operatorlesart (H|ψ⟩ = 0) ausgeschlossen. c = 24 mit der
Energiebedingung im Mittel ist damit nicht ausgeschlossen. Z3 und Z4 gelten für E₈³ genauso: Dort ist
nur i stationär, und die Energie ist null. Die Wahl c = 8 stützt sich auf TFPTs Familien in einem E₈
und auf C′ (Familien als Spinstruktursektoren). Unter G ist c = 8 × (Zahl der 16er-Sätze, deren Masse
windet).
- Windet ein Satz, ist c = 8 und F nötig.
- Winden drei, ist c = 24. Dann verschwindet die Rahmungsphase ohne F, und drei ist die kleinste Zahl
  mit 8N ∈ 24ℤ. Die Familiencharaktere {−1, i, −i} kämen dann aber nicht mehr aus der Naht, denn J
  wirkt innerhalb der Paare und vertauscht keine Sätze.

Die Gabel ist damit eine Frage nach der Zahl der windenden 16er-Sätze.

**Offen.**
- Die Identität der 16 Paare. Sie sind im Volumen schwer (Masse m) und daher nicht die leichten
  Fermionen des Standardmodells.
- Warum genau ein (oder drei) Satz windet.
- Warum der Grundzustand einen geschlossenen Defekt enthält.

## Davor: H6h5 — Gabel entschieden, Kontinuum (`z5_gabel.py`, `z6_kontinuum.py`)

```
        K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 K11 K12 K13 K14 K15 K16 K17 K18    P  L  F  O  E
   H6h4  P  P  P  P  P  P  P  P  P  P   P   P   L   P   P   O   L   L    14  3  0  1  0
   H6h5  P  P  P  P  P  P  P  P  P  P   P   P   L   P   P   L   L   L    14  4  0  0  0
```
Erstmals ist kein Kriterium mehr offen. Vier stehen unter benannter Lesart: Gravitation (D),
Kontinuum (physisches Gitter, perturbativ), Anomalie (F) und Neutrinomischung (TFPT-Kanäle).

**Gabel c = 8 gegen c = 24: entschieden für c = 8 mit F.** Alle fünf Prüfungen bestehen.
1. Die Operatorlesart H|ψ⟩ = 0 widerspricht H6 selbst. Der KMS-Zustand bei β = L hat Var(L₀) = 0,2426,
   ist also kein Eigenzustand. Unter H|ψ⟩ = 0 wäre e^{−βH} auf allen physischen Zuständen eins, und β,
   die Kragenbreite L/2 und die Zeitneutralität aus Z3 und Z4 wären unbestimmt (K6 fiele).
2. Ein E₈ trägt genau drei Spinstruktursektoren, die Terme θ₂⁸, θ₃⁸, θ₄⁸ der GSO-Summe. S hält (A,A)
   fest und tauscht AP mit PA; T tauscht (A,A) mit PA, mit der gemeinsamen Phase e^{−2πi/3}. Bis auf
   diese Phase wirkt SL(2,ℤ) als S₃ = SL(2, ℤ/2) auf genau drei Objekten, passend zu C′. E₈³ hätte neun
   Sektoren. TFPT bettet die Familien als SU(3)-Triplett in ein E₈ ein: 248 = 78 + 8 + 2·27·3.
3. **Neu: Die 3D Wand ist genau bei c ∈ 8ℤ unsichtbar.** Verschiebt man sie über ein 4D Gebiet, ändert
   sich die Phase um exp(2πi c σ/8) (Anomaliepolynom (c/24)p₁, ∫p₁ = 3σ). Weil die GSO-Summe die Naht
   spinstrukturfrei macht, zählen alle orientierten 4-Mannigfaltigkeiten, und c = 8 ist der kleinste
   Wert mit trivialer Phase. Das ist dieselbe 16 wie in Z1. Lesart: Die zwei Blätter von P1 sind
   entlang eines 3D Schnitts verklebt, der an der Naht endet. Trägt der Schnitt die invertierbare E₈
   Phase, ist seine Lage unbeobachtbar wie ein Dirac-String. F ist dann keine zusätzliche Struktur,
   sondern der Verzweigungsschnitt der zwei Blätter.

**Kontinuum (K16).** Das Gitter ist physisch (Z2: a ≥ 3 ℓ_P). Verlangt ist daher kein Limes a → 0,
sondern Universalität im Infraroten. Alle acht Prüfungen bestehen.
- **Lorentz ohne Feinabstimmung.** Auf dem isotropen Gitter ist die Lichtgeschwindigkeit für jede Art
  und jede Richtung exakt eins. Die Abweichung ist ~ p² (Dimension 6) und richtungsabhängig, zum
  Beispiel −1/3 entlang der Achse und −2/9 entlang der Raumdiagonale für das naive Fermion. Die
  Ein-Schleifen-Vakuumpolarisation ist entlang Zeit- und Raumachse gleich (−0,252432). Entlang der
  Diagonale weicht sie um 5·10⁻⁴ ab und skaliert mit q² (Verhältnis 0,260 gegen 0,25).
- **Gegenprobe mit ausgezeichneter Zeit.** Mit einem eigenen Zeitsprung ν je Art bekommt jede Art ihre
  eigene Geschwindigkeit, und das induzierte Photon liegt dazwischen: Zwei Arten mit ν = 1 und 0,9
  geben c_γ² = 1,122 gegen 1 und 1,235. Das ist Lorentzverletzung bei Dimension 4, nur durch
  Feinabstimmung zu beheben. Schon eine einzelne Art bekommt einen nicht-universellen Gitteranteil von
  etwa 2 %. **Die Gleichwertigkeit der vier Achsen, die die Kritik als Hindernis für Lorentz las, ist
  genau die Bedingung für Lorentz ohne Feinabstimmung.** Die Vierteldrehung zwischen Zeit- und
  Raumachse verbietet jedes ν ≠ 1 in jeder Ordnung.
- **Reflexionspositivität.** Naive Fermionen haben die Transfereigenwerte ±e^{−E}. T ist nicht positiv,
  T² schon (G(t+2)/G(t) = e^{−2E} exakt). Die OS-Rekonstruktion läuft über Doppelschritte; der
  Zeitverdoppler ist eine eigene Art.
- **Größenordnungen.** Quadratische Lorentzverletzung setzt bei etwa 2–7·10¹⁸ GeV ein; bei 100 TeV ist
  δv ~ 10⁻²⁷. Die Beobachtungsschranken liegen viele Größenordnungen darunter. Der Landau-Pol von
  U(1)_Y liegt bei etwa 10⁴¹ GeV, weit über der Gitterskala. Weil das Gitter physisch ist, ist die
  Trivialität von U(1) kein Hindernis; die nichtabelschen Faktoren sind asymptotisch frei.

**Nicht gezeigt.**
- Nichtperturbative Universalität (Confinement, Higgs-Sektor).
- Die chirale Gitterformulierung. Naive Fermionen haben 16 Geschmäcker in Paaren entgegengesetzter
  Chiralität (Nielsen–Ninomiya). Ein Weg aus der Literatur, hier nicht gerechnet: die Spiegelpartner
  symmetrisch massiv machen. Das geht nur, wenn die ℤ₁₆-Anomalie einer Generation verschwindet, also
  mit 16 Weyl-Fermionen je Generation einschließlich ν_R.

**Prämissen von H6h5:** A\*, C′, D, E und F. F ist jetzt an K1 gebunden: Es ist die E₈ Phase auf dem
Verzweigungsschnitt.

## Davor: Stabilität der Lage τ = i (`z4_stabilitaet.py`)

**Ergebnis: Die Lage ist stabil. Der Sattel aus Z3 ist keine Instabilität**, sondern die gewöhnliche
Struktur zweier konjugierter Paare. Alle neun Prüfungen bestehen.

In Rekonstruktion A (Raum L, Zeit β) ist τ = x + iy. Die beiden Richtungen des Sattels bedeuten
Verschiedenes.
- **Verdrillung x** koppelt an den Impuls: Z = Tr e^{−βH} e^{iθP}. Bei positivem Spektrum gilt
  |Z(x + iy)| ≤ Z(iy), mit Gleichheit nur bei ganzem x. Verdrillung null ist für jede unitäre Naht ein
  echtes Maximum von |Z|. Für E₈ gilt |Z(x+i)|/Z(i) = 0,953 bei x = 0,1 und 0,373 bei x = ½, also
  F_xx = +9,58: **stabil**. Die „eingefrorene Verdrillung“, die frühere Befunde als Zusatzannahme
  brauchten, ist damit automatisch.
- **Streckung y = β/L** ist die inverse Temperatur. log Z ist in β konvex, denn die zweite Ableitung
  ist die Energievarianz. Für E₈ gilt bei y = 1 exakt F_yy = −(2π)²·Var(L₀) = −9,58, mit
  Var(L₀) = 0,2426. Längs der Temperatur ist −log Z daher **für jede unitäre Naht** konkav. Das ist
  positive Wärmekapazität, keine Instabilität. Ein Minimum über die Temperatur gibt es grundsätzlich
  nicht; die Temperatur ist eine Legendre-Variable.
- **Mikrokanonisch:** Eine hintergrundfreie, geschlossene Naht hat kein äußeres Wärmebad. Ihre
  Temperatur wird durch die Energiebedingung festgelegt, nicht gewählt. S(E) = min_β[log Z + βE] ist
  konkav (stabil). Bei E = 0 liegt genau β = L, also τ = i, mit S(0) = log χ(i) = log 12. Dort gilt
  ⟨L₀⟩ = c/24 = 1/3 exakt: Casimir- und Anregungsenergie heben sich auf.

Die Lage τ = i ist damit der stabile Gleichgewichtszustand einer geschlossenen Naht bei Energie null
und Impuls null: stationär (Z3) und stabil (Z4).

**Neu sichtbar: die Gabel F gegen F′ als Frage nach der Energiebedingung.** Bei c = 8 liegen die
Niveaus bei n − 1/3; einen Zustand mit Energie genau null gibt es nicht. E = 0 gilt also nur im
Mittel des thermischen Zustands. Soll die Bedingung als Operatorgleichung gelten (H|ψ⟩ = 0,
Wheeler-DeWitt-artig), muss c durch 24 teilbar sein. Bei E₈³ (Z = j) hat das Nullniveau 744 = 3 × 248
Zustände. Das ist dieselbe Bedingung wie die T-Invarianz aus Z1. Die Lesart im Mittel passt zu TFPT,
das ein E₈ mit 248 = 78 + 8 + 2·27·3 nutzt (SEAM.EQUIV.02). Entscheiden muss das eine eigene Prüfung.

**Prüfstand.** Die Bilanz bleibt 14 P, 3 L, 0 F, 1 O. K2 trägt jetzt „stationär und stabil“.

**Offen.**
- K16: wechselwirkender Kontinuumslimes, Reflexionspositivität und SO(4).
- a/ℓ_P mit vollem Feldinhalt.
- Die Energiebedingung im Mittel (F, c = 8) oder als Operator (F′, c = 24).
- Warum es überhaupt eine geschlossene 2D Naht mit chiralen Fermionen gibt.

## Davor: H6h4 — Prämisse A hergeleitet (`z3_hintergrundfreiheit.py`)

```
        K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 K11 K12 K13 K14 K15 K16 K17 K18    P  L  F  O  E
   H6h3  P  P  P  P  P  P  P  P  P  P   P   P   L   P   P   O   L   L    14  3  0  1  0
   H6h4  P  P  P  P  P  P  P  P  P  P   P   P   L   P   P   O   L   L    14  3  0  1  0
```
Die Bilanz ist gleich; geändert hat sich die Grundlage. Die setzende Prämisse A („der orientierte
Tausch der Zyklen ist eine Symmetrie“) ist ersetzt durch:

**A\*: Die Naht hat keine Hintergrundgeometrie.** Ihre Basiswechsel (große Diffeomorphismen) sind
Symmetrien, und ihre Form ist nicht vorgegeben, sondern erfüllt ihre Bewegungsgleichung: Die Antwort
des Zustands auf jede Deformation der komplexen Struktur verschwindet, ⟨T⟩ = 0. Das ist die
Bedingung einer dynamischen Nahtmetrik, wie die Virasoro-Bedingung.

**Was daraus folgt (alle sieben Prüfungen bestanden).**
- Für die E₈-Naht gilt exakt d log χ/dτ = −(2πi/3)·E₆/E₄ (geprüft auf 10⁻²¹). ⟨T⟩ = 0 heißt also
  E₆(τ) = 0. Nach der Valenzformel hat E₆ im ganzen Fundamentalbereich genau eine Nullstelle, die
  halbe am elliptischen Punkt i. **Die Form der Naht ist damit eindeutig τ = i**, ohne
  Orientierungsvoraussetzung und ohne eingefrorene Verdrillung. Am Sechseck ist χ = 0; das ist keine
  zulässige Form.
- Dort ist die Gesamtenergie in beiden Rekonstruktionen exakt null: Die thermische Energie hebt die
  Casimirenergie auf. Bei β/L = 0,8 und 1,25 vertauschen sich die Energien (+2,40 und −1,54). Beide
  Rekonstruktionen sind am stationären Punkt dasselbe System (β = L). **Zeitneutralität ist damit
  eine Folge, keine Prämisse.** Der Stabilisator von i ist μ₄; daraus folgen wie bisher die
  zyklischen Marken, die Viertelzählung und c₃.
- **Verbot (neu):** Für eine chirale Naht ist log|Z| harmonisch (Hesse +9,58 / −9,58, Spur 0). Es gibt
  daher nirgends ein lokales Minimum. Energieminimierung kann die Form grundsätzlich nicht wählen; nur
  Stationarität kann es. Das erklärt, warum frühere Minimierungsansätze (KERNSIX, KERNSYM) bei τ = i
  einen Sattel fanden.
- **Gegenproben.** Bei der Dirac-Naht mit fester Spinstruktur (A,A) ist i im Modulraum der
  Thetagruppe ebenfalls die einzige stationäre Form. Bei E₈³ (Z = j) ist es nur i. Beim freien Boson
  und beim Leech-Gitter ist zusätzlich das Sechseck ρ stationär; dort entscheidet erst die
  Spinstruktur (Z1: am Sechseck bleibt nur PP mit Nullmoden). Die Auswahl ist also kein Artefakt
  jeder beliebigen Naht.

**Was davon schon bekannt war.** Dass die Energie des E₈-Charakters genau bei τ = i verschwindet und
dass E₆ dort seine einzige Nullstelle hat, steht im Atlas (Befund vom 18.11., „c₃ ist die
Gitterladungsenergie …“). Ebenso ist dort festgehalten, dass τ = i ein Sattel ist. Neu ist:
- diese Nullstelle als Bewegungsgleichung einer hintergrundfreien Naht zu lesen, die den gesetzten
  Fixpunkt ersetzt;
- das Verbot für Minimierung über die Harmonizität;
- die Gegenproben.

**Prämissen von H6h4.**
- **A\*** (hintergrundfreie Naht) ersetzt A.
- C′, D, E, F bleiben wie in H6h3.
- A\* trägt zweimal: Die großen Diffeomorphismen geben die GSO-Summe und damit E₈ mit T-Phase (daher
  F). Die Bewegungsgleichung gibt τ = i.

**Was A\* nicht leistet.**
- (Nachtrag: Z4 zeigt, dass der Sattel keine Instabilität ist; siehe oben.) Die Lage ist stationär, aber nicht stabil: Sattel, instabil gegen Streckung. Wenn die Form wirklich
  dynamisch ist, braucht die Stabilität einen weiteren Grund. Kandidaten sind die 3D E₈ Phase (F),
  die Ankopplung an die 4D Geometrie oder das Modulintegral, in dem ein Sattel der maßgebliche
  Punkt ist.
- Dass es überhaupt eine geschlossene 2D Naht mit chiralen Fermionen gibt, bleibt Ausgangspunkt.

## Davor: H6h3 — externe Kritik geprüft (`z1_kritik_pruefung.py`, `z2_verschraenkungsdichte.py`)

```
        K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 K11 K12 K13 K14 K15 K16 K17 K18    P  L  F  O  E
   H6h2  P  P  P  P  P  P  P  P  P  P   P   P   L   P   P   O   O   L    14  2  0  2  0
   H6h3  P  P  P  P  P  P  P  P  P  P   P   P   L   P   P   O   L   L    14  3  0  1  0
```
K17 (Rahmungsanomalie der Naht) und K18 (Neutrinomischung gegen Daten) sind neu im Katalog, damit
die Punkte der Kritik sichtbar bleiben. H6h2 wird rückwirkend ehrlich bewertet: K17 offen.

**Die Kritik Punkt für Punkt.** Alle 16 Rechenprüfungen in `z1` bestehen.
1. *Kovarianz wählt kein Quadrat.* Richtig für die Zustandssumme: Z_AA(−1/τ) = Z_AA(τ) gilt auch bei
   τ = 2i. Prinzip A verlangt aber mehr, nämlich dasselbe thermische System in beiden
   Rekonstruktionen. Bei τ = 2i sind es verschiedene Systeme (β/L = 2 gegen 1/2, Besetzung der
   tiefsten Mode 0,0019 gegen 0,172). A ist also die Fixpunktbedingung. Die von der Kritik
   vorgeschlagene Verschärfung (echte Drehung der Ordnung größer zwei, Spinstruktur ohne Nullmoden
   erhalten) ist dieselbe Bedingung. Drehungen der Ordnung > 2 gibt es nur bei i und ρ; bei ρ
   bleibt nur PP erhalten, das Nullmoden hat; also τ = i mit (A,A). Das deckt sich mit dem
   Atlasbefund vom 16.11. (nur j = 1728 erlaubt Ordnung vier auf den Marken). Offen bleibt, warum
   die Natur den Fixpunkt realisiert: Das ist Prämisse A selbst.
2. *Spinorwirkung der Vierteldrehung.* Gerechnet mit der gelifteten Drehung g̃ um den Viertelpunkt:
   g̃ vertauscht mit dem Dirac, g̃⁴ = −1, es gibt keinen invarianten Sektor, und der naive Projektor
   (1/4)Σg̃ᵏ ist keiner. Die Sektoren sind die vier Eigenwerte exp(iπ(2m+1)/4). Weil kein Impuls in
   (ℤ+½)² fix ist, trägt jeder Sektor **exakt ein Viertel jedes Eigenraums**, ohne
   Fixpunktkorrektur. Damit ist c₃ = 1/(8π) je Spinorsektor (Wärmespur 0,039788735773). Die
   Viertelzählung überlebt die fermionische Drehphase.
3. *Woher das Familienvorzeichen?* Neu und exakt: Die drei geraden Spinstrukturen der Naht sind die
   nichttrivialen Charaktere der Markengruppe E[2]. Die zyklische Markendrehung g(z) = iz + ½ ist
   affin. Ihr Verschiebungsteil wirkt auf diesen Charakteren über die Weilpaarung als Vorzeichen,
   und das Spektrum ist **genau {−1, i, −i}**: (A,A) mit −1, AP und PA mit ±i. Kein Vorzeichen ist
   eingesetzt. Eine Drehung um eine Marke statt um den Viertelpunkt gäbe {1, 1, −1}. Dasselbe g auf
   H¹ der vier Marken gibt dasselbe Spektrum. Für alle 24 Elemente von S₄ = ⟨R, C⟩ = Aut(Q₈) gilt
   Charakter auf Im H = Signum × Standard = Λ²(Standard). Die chiralen 2-Formen der Raumzeit sind
   also das Λ² des Familienraums, und das Signum ist dessen Determinante. Das ist genau der
   orientierte Viermarkenzustand, den die Kritik vermutet. Prämisse C wird damit zu
   **C′: Familien = gerade Spinstrukturen der Naht, chirale 2-Formen = ihr Λ²**.
4. *E₈ und die T-Phase.* Richtig: Die Vierteldrehung braucht nur S, und S ist für jedes N exakt.
   T-Kovarianz der GSO-Summe verlangt N ≡ 0 mod 16, strenge T-Invarianz N ≡ 0 mod 48. E₈ (N = 16)
   ist daher nur bis auf exp(−2πi/3) basisunabhängig. Es gibt zwei konsistente Ergänzungen:
   **(F)** eine invertierbare 3D E₈ Phase mit Rand c = 8. Ihr T²-Hilbertraum ist eindimensional,
   also bricht die Wahl der Füllung die Vierteldrehung nicht. Oder (F′) drei identische E₈ Blöcke,
   c = 24, mit Zustandssumme χ_E8³ = j (χ_E8(i) = 12, j(i) = 1728). H6h3 nimmt F als neue Prämisse;
   F′ würde TFPTs c = 8 in der Gravitation ändern. Nebenbefund: arg exp(−2πi/3) = 240° ist TFPTs
   δ_PMNS. Die Zahl ist gleich; ein Zusammenhang ist nicht hergeleitet.
5. *Lorentzsignatur.* Richtig: Die Gleichwertigkeit der vier Achsen ist eine euklidische Aussage.
   H6 nutzt sie nur so. Die lorentzsche Zeit entsteht durch OS-Rekonstruktion entlang der Achse, die
   der Zustand wählt. Dafür braucht es Reflexionspositivität der Gitterwirkung und die
   Wiederherstellung von SO(4) im Limes. Beides gehört zu K16 und ist nicht gezeigt.
6. *Eichdynamik.* Verbindung, Krümmung und induzierte Wirkung sind in `y1` gerechnet. Die gemeinsame
   Normierung kommt aus denselben Fermionschleifen bei vollständigen Multipletts. Schwellen durch
   ungleiche Massen sind nicht behandelt.
7. *Gravitation braucht η mikroskopisch (Gegentest D, `z2`).* Gerechnet: Ein naives Diracfeld auf
   ℤ³ trägt je Gitterfläche κ = 0,78657 nat (konvergiert in Ring und Querimpulsen; nicht π/4),
   also 0,0983 je Geschmack. η = 1/4 entspricht 2,54 Dirac-Geschmäckern. Drei Generationen mit ν_R
   geben schon ohne Eichfelder κ ≥ 2,36. Unter D ist daher a/ℓ_P = 2√κ ≥ 3,07 (naiv mit Verdopplern
   8,7). **Korrektur:** Das Lipschitzgitter sitzt nicht in Planckeinheiten. Unter D ist η = 1/4 die
   Definition der Planckeinheit, und a/ℓ_P wird eine Vorhersage, sobald der volle Feldinhalt samt
   Kontaktterm der Eichfelder feststeht (Sakharovs Artenproblem in Gitterform).

**θ₁₃ neu eingeordnet.** TFPT hat θ₁₃ bereits als eigenen Trägerkanal: sin²θ₁₃ = φ₀·exp(−tr Y²) mit
tr Y² = 5/6 über der Hyperladung der fünf Plätze (FLAV.TH13.01). In H6 ist Y der Generator, den der
lokale 3+2-Split (Prämisse E) liefert. H6 gibt die führende Ordnung (tribimaximal aus S₄), TFPT die
Korrekturen: sin²θ₁₂ = 0,3067, sin²θ₂₃ = 1/2, sin²θ₁₃ = 0,0231. Laut TFPTs eigener Prüfliste liegt das
etwa 2σ von den Daten; JUNO entscheidet. θ₁₃ ist damit keine Lücke von H6 mehr, sondern eine
geerbte Vorhersage (K18 unter Lesart).

**Der gemeinsame Nahtoperator (Vorschlag 5 der Kritik).** Ein einziger Operator, die geliftete
Vierteldrehung g̃ um den Viertelpunkt, leistet zwei der drei geforderten Eigenschaften exakt.
- Er permutiert die vier Marken zyklisch und gibt c₃ = 1/(8π) je Spinorsektor, mit Spinphasen.
- Auf den geraden Spinstrukturen gibt er die Familien mit {−1, i, −i}, ohne eingesetztes Vorzeichen.

Die dritte Eigenschaft, E₈ mit vollständig behandelter Anomalie, gilt für g̃ selbst (S ist exakt).
Für vollständige Basisunabhängigkeit braucht sie die Ergänzung F.

**Prämissen von H6h3.**
- **A:** Der orientierte Tausch der Zyklen ist eine Symmetrie (Fixpunkt statt Kovarianz).
- **C′:** Familien = gerade Spinstrukturen der Naht; die chiralen 2-Formen sind ihr Λ².
- **D:** Gravitation als Zustandsgleichung (Jacobson).
- **E:** Der 3+2-Split ist lokal.
- **F:** Die Rahmungsanomalie wird durch eine invertierbare 3D E₈ Phase kompensiert.

**Offen.**
- K16: wechselwirkender Kontinuumslimes, Reflexionspositivität und SO(4).
- a/ℓ_P mit vollem Feldinhalt.
- F gegen F′.
- Warum der Fixpunkt.
- Der Weg von der T-Phase zur CP-Phase, falls es einen gibt.

## Davor: H6h2 — Eichdynamik und Prämissen (`y1_eichdynamik_praemissen.py`)

```
        K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 K11 K12 K13 K14 K15 K16    P  L  F  O  E
   H6h   P  P  P  P  P  P  P  P  P  P   P   P   L   P   O   O    13  1  0  2  0
   H6h2  P  P  P  P  P  P  P  P  P  P   P   P   L   P   P   O    14  1  0  1  0
```

**Eichdynamik (K15) hergeleitet unter einer kleinen Prämisse E** (der 3+2-Split ist lokal):
- Der örtlich variierende Split definiert auf jeder Gitterkante eine Paralleltranslation
  U = polar(V_x* V_y); sie ist unitär und eichkovariant (U ↦ h_x* U h_y).
- Fermionschleifen auf dem Lipschitzgitter induzieren den Yang–Mills-Term. Mit Fermionmasse als
  Infrarotregler und im Schwachfeldbereich (F/M² ≤ 0,12) ist die Wirkung **quadratisch** in F
  (Fluss 1 und 2 gleich auf 2 %), **positiv**, **eichinvariant**, und die Kopplung **läuft
  logarithmisch** mit d(1/g²)/d log(1/M) = 0,287 gegen 0,270 im Kontinuum (16 Diracarten).
  Ohne Masse setzt das Feld selbst die Infrarotskala (F² log F), wie erwartet.
- Eine Generation hat für SU(3), SU(2) und U(1) (GUT-normiert) dieselbe Dynkin-Summe 2. Die
  induzierten Kopplungen sind daher gleich, **sin²θ_W = 3/8 an der Gitterskala**; TFPT hat 3/8 an
  der Trägerskala.

**Prämissen geprüft.**
- **A** (keine ausgezeichneten Zyklen) ist tragend. Ohne A wären τ, c₃ (1/(2π) oder 1/(4π)), die
  Kragenbreite, die Fermionenzahl der Naht, E₈ und die Spinstruktur frei.
- **B entfällt.** Drei chirale Mittelformen gibt es nur in d = 4. Mit der 2D Naht ist die
  Kodimension 2, und A auf allen vier Zyklen macht die Normalebene gleichwertig. B folgt also aus
  A, C und der Naht.
- **C gegen Daten.** Die Familiengruppe S₄ (Ordnung 24) mit Restsymmetrien Z₃ für die geladenen
  Leptonen und Z₂ × Z₂ für Majorana-Neutrinos gibt genau die **tribimaximale Mischung**:
  sin²θ₁₂ = 1/3 (TFPTs Leitterm; mit −φ₀/2 sind es 0,3067 gegen gemessen 0,307), sin²θ₂₃ = 1/2
  (gemessen 0,47 bis 0,57), aber **sin²θ₁₃ = 0 gegen gemessen 0,022**. C besteht zwei von drei
  Winkeln in führender Ordnung; θ₁₃ verlangt eine Brechung, die H6 noch nicht liefert. Die
  Majorana-Forderung passt zu TFPTs Kandidaten vom 18. November.
- **D** bleibt der Rahmen (Jacobson wie TFPT). Entscheidender Gegentest für später: die
  Flächendichte der Verschränkung des Gittervakuums ausrechnen; sie legt fest, ob die Gitterskala
  die Planckskala ist (2π/η = 8π = 1/c₃).
- **E** (neu, mild): der Split ist ein lokales Feld.

**Was jetzt noch offen ist.**
1. **K16:** Der freie Sektor hat die renormierbare Struktur (logarithmischer Lauf). Ein
   wechselwirkender Kontinuumslimes ist nicht gezeigt.
2. **θ₁₃:** Die S₄-Restsymmetrie gibt null; der gemessene Wert verlangt eine Korrektur. Das ist
   jetzt die schärfste Datenprüfung von H6.
3. **D:** Gravitation nur im Jacobson-Rahmen; Gegentest oben.

## Davor: H6 gehärtet (H6h)

Gleicher Prüfstand, Katalog um K15 (Eichdynamik in 4D) und K16 (Kontinuumslimes) erweitert,
damit nichts Fehlendes verdeckt wird.

```
        K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 K11 K12 K13 K14 K15 K16    P  L  F  O  E
   H6    P  P  L  P  P  P  L  P  L  L   P   L   O   P   O   O     8  5  0  3  0
   H6h   P  P  P  P  P  P  P  P  P  P   P   P   L   P   O   O    13  1  0  2  0
```

**Was die Härtung leistet.** Die fünf Lesarten von H6 sind durch Rechnungen ersetzt, die nur noch
auf vier offen genannten Prämissen stehen:
- **A** Keine ausgezeichneten Zyklen (zeitneutral und basisunabhängig).
- **B** Die Naht ist eine Kodimension-2-Fläche mit gleichwertiger Normalebene.
- **C** Familien = chirale Formen mittlerer Dimension der Raumzeit.
- **D** Gravitation als Zustandsgleichung (Jacobson, wie TFPT).

Die Ersetzungen im Einzelnen:
1. **GSO, 16 Fermionen und E₈ aus A.** Eine chirale Fermion-Naht, die von der Wahl der
   Basiszyklen nicht abhängt, ist nur für N ≡ 0 mod 16 modular konsistent (S exakt, T bis auf eine
   Phase; geprüft N = 8, 16, 24, 32, 48). Das kleinste N = 16 gibt c = 8 und das Gitter D₈⁺ = E₈.
   Der (P,P) Sektor trägt die 8 flachen Nullmoden des E₈ Worts. Die frühere Lesart „B ist über
   Spinstrukturen summiert“ ist damit hergeleitet.
2. **Familien, ihre Charaktere und das Drehzentrum aus C.** σ und die Vierteldrehung erzeugen auf
   den drei chiralen 2-Formen eine Gruppe der Ordnung 24 mit dem Charakter Standard × Signum von
   S₄. Das ist genau die Darstellung, die TFPT für die Familienmonodromie aufgezeichnet hat (v117,
   „signgetwistete Standarddarstellung von S₄“). H¹ der Marken ist die Standarddarstellung (Satz).
   Der früher „gewählte“ q-Charakter ist also das Signum von S₄. Die beiden Lesarten derselben
   Familien stimmen nur bei zyklischer Wirkung auf den Marken überein; das erzwingt das
   Drehzentrum auf dem Viertelpunkt ({−1, i, −i} gegen {−1, 1, 1}).
3. **3 + 1 aus B und C.** Drei chirale Mittelformen gibt es nur in d = 4 (d = 8: 35, d = 12: 462).
   Naht plus gleichwertige Normalebene ist ℤ[i] ⊕ ℤ[i]·j, also das Lipschitz-Gitter.
4. **Gravitation unter D mit Eingaben aus H6.** TFPTs Herleitung braucht c = 8, |μ₄| = 4 =
   |ℤ₂|·χ und 2π/η = 8π = 1/c₃. H6 liefert alle drei; c = 8 ist jetzt abgeleitet (N = 16), nicht
   gesetzt, und stimmt mit TFPTs g_car + N_fam = 8 überein.

**Neue Folgerung, zu prüfen.** Eine chirale Naht mit c = 8 hat die Rahmungsanomalie
exp(−2πi c/24) = exp(−2πi/3). Sie verlangt einen Zufluss aus einem 3D Volumen, also eine
E₈-artige Wand, die an der Naht endet.

**Was noch fehlt (ehrlich).**
- **K15 Eichdynamik:** Die Erhaltungsgruppe (SU(3) × SU(2) × U(1))/ℤ₆ folgt aus 3 + 2 (TFPT,
  exakt). Yang–Mills-Felder in 4D folgen aus H6 nicht.
- **K16 Kontinuumslimes:** Das Lipschitzgitter sitzt in Planckeinheiten; ein wechselwirkender
  Skalenlimes ist nicht gezeigt.
- **Prämissen B, C, D** sind physikalische Postulate, keine Herleitungen. Die quantitativen
  Ergebnisse (α, φ₀, Massen) kommen wie bisher aus dem TFPT Compiler mit dessen eigenen
  bedingten Schritten.

**Nächste Schritte am Ganzen (kein Nebenschauplatz).** K15 auf dem Lipschitzgitter angehen:
TFPTs Befund „Grassmannkrümmung spannt s(u(3) ⊕ u(2))“ als Gitter-Eichfeld auf den Verbindungen
von ℤ⁴ formulieren und prüfen, ob die Wirkung Yang–Mills ergibt. Danach K16 und für jede Prämisse
einen entscheidenden Gegentest. Nach jedem Schritt läuft der ganze Prüfstand neu.

## Vorher: vollständige Hypothesen auf einem festen Prüfstand (`v1_vollhypothesen.py`)

Arbeitsweise seit dem Hinweis des Nutzers: Jede Hypothese deckt die ganze Lösung ab und wird
gegen denselben Katalog K1…K14 geprüft. Ein Fehlschlag beendet die Prüfung nicht.

Katalog: K1 vier Marken mit zwei Blättern · K2 quadratische Lage τ = i · K3 μ₄ permutiert die
Marken zyklisch · K4 Spinstruktur (A,A) · K5 c₃ = 1/(8π) · K6 Zustand und Zeit mit festgelegter
Breite · K7 flacher Sektor für das E₈ Wort · K8 P2 über P1 · K9 drei Familien als gleiche Kopien ·
K10 Familiencharaktere {−1, i, −i} · K11 E₈ · K12 3+1 · K13 Gravitation · K14 keine freien
Parameter.

```
        K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 K11 K12 K13 K14    P  L  F  O  E
   H1    E  P  L  P  P  P  F  P  O  P   E   F   F   P     7  1  3  1  2
   H2    P  P  L  P  P  F  F  P  P  F   P   L   O   P     8  2  3  1  0
   H3    E  E  E  P  O  O  F  E  E  E   E   E   P   F     2  0  2  2  8
   H4    E  P  L  P  P  P  L  P  O  P   P   F   F   P     8  2  2  1  1
   H5    P  P  L  P  P  F  L  P  P  F   P   L   O   P     8  3  2  1  0
   H6    P  P  L  P  P  P  L  P  L  L   P   L   O   P     8  5  0  1  0
```
P hergeleitet, L unter benannter Lesart, F gescheitert, O offen, E Eingabe.

**H6 ist die erste vollständige Hypothese ohne Fehlschlag.** In Worten:
- **Raumzeit** ist das ganzzahlige Gitter der größten assoziativen Divisionsalgebra, die
  Lipschitz-Quaternionen ℤ⁴ = ℤ⟨1, i, j, k⟩. Ihre Einheiten ±1, ±i, ±j, ±k machen alle vier
  Achsen gleichwertig (keine vorab ausgezeichnete Zeit). Die Fermionen sind auf allen vier
  Zyklen antiperiodisch (invariant, ohne Nullmoden). 3 + 1 ist Re ⊕ Im.
- **Die Naht** ist eine Gaußsche Ebene ℤ[i] darin: quadratisch, τ = i, Spinstruktur (A,A). Die
  π-Drehung −1 gibt zwei Blätter und vier Marken. Die Vierteldrehung um einen Viertelpunkt
  permutiert die Marken zyklisch; c₃ = 2/(4π·4) = 1/(8π); die Kragenbreite ist ℓ = L/2.
- **P2** folgt aus der Naht wie in TFPT (P2.PARTITION.01).
- **E₈** ist dasselbe Gitter zweimal: aus ℤ⁸ = (ℤ⁴)² mit Halbzahlkleber (112 + 128 = 240) und als
  GSO-Summe von 16 Majorana-Weyl-Fermionen auf der Naht (120 + 128 = 248). Der (P,P) Sektor hat
  genau 8 komplexe Nullmoden, flach für jede Breite; dort lebt das E₈ Wort. Das löst die Gabel
  aus Abschnitt N.
- **Drei Familien** sind die drei selbstdualen 2-Formen (die drei Gaußschen Geraden ℤ[i], ℤ[j],
  ℤ[k]), von σ zyklisch permutiert, also gleiche Kopien. Die Vierteldrehung wirkt auf ihnen mit
  {1, i, −i}; zusammen mit dem Vorzeichencharakter von q, den TFPT auch für die Trägermoden
  nutzt, ergibt das genau {−1, i, −i}.
- **Gravitation**: 4D lokal; die Einsteingleichung mit Koeffizient 1/c₃ liefert nur TFPTs
  bedingte Herleitung (GRAV.NONLINEAR.01).

**Die fünf Lesarten von H6 (dort kann sie noch falsch sein).**
1. K3: Keine Marke ist ausgezeichnet, also liegt das Drehzentrum auf einem Viertelpunkt.
2. K7: Der Randoperator B ist der Naht-Dirac summiert über die Spinstrukturen (GSO).
3. K9: Familien = selbstduale 2-Formen.
4. K10: Die Familien tragen den Vorzeichencharakter von q. Diese Lesart ist nach Ansicht des Ziels
   gewählt; derselbe Charakter steht aber schon in TFPTs Trägermoden.
5. K12: Raumzeit = größte assoziative Divisionsalgebra (Frobenius).

**Ausschlüsse dieses Prüfstands (wiederverwendbar).**
- *Dreierdrehungen verbieten zeitneutrale invertierbare Fermionen.* Das hexagonale 2D-Gitter und
  das Hurwitz-4D-Gitter haben keine einheiteninvariante Spinstruktur ohne Nullmoden; das
  quadratische und das hyperkubische (Lipschitz) haben sie. Deshalb μ₄ und nicht μ₆, und deshalb
  Lipschitz statt Hurwitz.
- *Tensorprodukt-Spektraltripel* (H3) verlieren den flachen Sektor: mit invertierbarem
  geometrischem Dirac gibt es keine Nullmoden.
- *Naht allein* (H1, H4) liefert kein 3+1 und keine Gravitation.

**Nächster Schritt (für H6 als Ganzes).** Jede der fünf Lesarten durch eine Herleitung oder einen
entscheidenden Test ersetzen und danach den ganzen Prüfstand neu laufen lassen. Dazu die
Gravitation (K13) auf der Lipschitz-Raumzeit selbst prüfen, statt sie aus TFPT zu übernehmen.

---

Gerechnet am 23. November 2026 mit `.venv/bin/python` (numpy, scipy, sympy, mpmath). Zehn Skripte,
alle mit Rückgabewert 0; Ausgaben in `ausgabe_*.txt`. PLAN, progress, Ledger und Atlas sind
unberührt. Ausgangspunkt war der konsolidierte Stand in `docs/research/TFPT-STAND.md`
(Abschnitt „Befunde 20.–23. November“): Die offenen Eingaben laufen auf „woher die Vier“ und
„warum minimal“ hinaus.

## Übersicht

| Nr | Urteil | Stufe | Erste freie Wahl |
|---|---|---|---|
| **N Zeitneutrale Fermionnaht (Ursprung von P1)** | **Kandidat** | 1 für P1: Lage der Marken, Viertelzählung in c₃, Zustand; nicht 2, weil Naht und zwei Blätter Eingabe bleiben | Prinzip „keine vorab ausgezeichnete Zeitrichtung“ |
| **G Gaußscher Torus ℂ/ℤ[i]** | **Kandidat, korrigiert**: P2-Teil ist Doppelherleitung (P2 folgt aus P1) | 1 für P1 (Lage der Marken, c₃) | Startklasse CM-Torus; Auswahl über eine P2-Regel |
| W Weil-Darstellung von SL(2, 𝔽₅) | **Kandidat** | 2 für P2 unter benannter Lesart | Auswahlregel „fermionischer Teil = treues SU(2) Dublett“; Lesart „fünf Plätze = Ortsbasis eines Fünfersystems“ |
| W2 Clifford-Gruppe eines Fünfersystems | Uhrenrelation **bestätigt**, sechs Fasern **verworfen** | – | – |
| W3 die 16 unter 2I | **bestätigt**, exakt | – | – |
| H minimaler holografischer Code | **offen** | 1, an der Basislinie | Abstand 3; Holografie |
| P Poincaré-Sphäre | bereits in TFPT (v219), kein neuer Weg | – | – |
| T c₃ als Bekenstein-Viertel mal Unruh-Periode | bereits in TFPT, kein neuer Weg | – | – |
| J Jordanalgebra J₃(𝕆), R Renormierung | **offen**, nicht gerechnet | – | siehe IDEAS.md |
| F Familientest zu W (Nachtrag) | Kopienweg **verworfen**; Wörterbuch über 𝔽₅ **offen** | 1 für die Zahl 3 | siehe Abschnitt F |
| K Familien als étale Kohomologie (Nachtrag) | Test **bestanden**, aber nicht unabhängig; Frobeniusgleichheit **neu** | 1 | siehe Abschnitt K |

Keine vollständige Lösung (Stufe 4 fehlt: physische Zeit und Zustand, 3 + 1, Kontinuum,
Gravitation).

**Korrektur (Hinweis des Nutzers, mit Quellen geprüft).** P2 folgt in TFPT bereits aus P1: Aus der
markierten Viermarkennaht folgen Fünferträger und Teilung 3 + 2 exakt (P2.PARTITION.01,
AX.P2.01). W, F und der P2-Teil von G sind daher eine zweite Herleitung von schon Hergeleitetem,
kein neuer Ursprung. Was von G bleibt, ist die P1-Seite: die quadratische Lage der Marken und
c₃ aus demselben Torus. Dort wählt G den Torus aber über eine P2-Regel (Spinordublett im
Weil-Oszillator); eine Begründung von der P1-Seite selbst fehlt. Daran setzt Abschnitt N an.

## N: Ursprung von P1 — die zeitneutrale Fermionnaht (`n1_zeitneutrale_naht.py`)

**Prinzip.** Die Naht ist eine geschlossene, orientierte, flache 2D-Fläche mit freien Fermionen,
und sie zeichnet **keine Zeitrichtung** aus: Rekonstruiert man nach Osterwalder–Schrader entlang
jedes ihrer beiden Basiszyklen ein Quantensystem, entsteht dasselbe orientierte thermische
System. Jede Rekonstruktion soll ein echter KMS-Zustand sein (Fermionen in der Zeit
antiperiodisch) mit invertierbarem Hamiltonoperator (keine Nullmoden). Zeit kommt dann aus dem
Zustand, nicht aus der Geometrie.

**Befund (exakt bis auf die Gittersuche).**
- Rekonstruktion A hat (β/L, Drall/L) = (Im τ, Re τ), die orientierte Vertauschung B hat
  (Im τ/|τ|², −Re τ/|τ|²). Gleichheit erzwingt |τ| = 1 und Re τ = 0, also **τ = i**. Zusammen mit
  KMS und Invertierbarkeit bleibt im ganzen Fundamentalbereich nur **τ = i mit Spinstruktur
  (A, A)**. Lässt man Parität zu, bleibt die ganze rhombische Familie |τ| = 1; die Orientierung
  von P1 ist also tragend. Am hexagonalen Punkt wechselt der Drall das Vorzeichen.
- Die Zeitneutralitätsgruppe ist die Vierteldrehung S, Ordnung 4 = μ₄. Kein Impuls in (ℤ+½)²
  ist unter ihr fix, alle Bahnen haben Länge 4. Die bisher gewählte **Viertelzählung ist damit
  exakt die Zählung zeitneutraler Bahnen**, und c₃ = (Diracrang 2)/(4π · |μ₄|) = **1/(8π)** ohne
  Wahl einer Fläche (s·Bahnzahl → 0,039788735773). Mit Parität wäre es 1/(4π), ohne Neutralität
  1/(2π).
- **Zustand und Zeit.** KMS-Zustand der Diracfermionen bei selbstdualer Temperatur β = L,
  Einteilchenenergien 2π|n + ½|/L, Lücke π/L. Zeit ist jeder der beiden Zyklen; die zwei Wahlen
  sind durch μ₄ äquivalent. Die Zustandssumme am Fixpunkt ist exakt Z(i) = |θ₃(i)/η(i)|² = 2.

**Was neu ist.** Bekannt waren: τ = i als einziger symmetrieerzwungener kritischer Punkt (KERNSYM,
dort ist die Symmetrie τ ↦ −1/τ vorausgesetzt und τ = i nur Sattel), die Wahl (−,−) durch
Viertelinvarianz und fehlende Nullmoden, und 1/(8π) aus einer gewählten Viertelzählung (alles
23.11. oder früher). N ersetzt die zusätzlich genannte quadratische Lage und die gewählte
Viertelzählung durch **eine** physikalische Forderung und liefert dazu einen Zustand mit Zeit.
Energieminimierung leistet das nicht (KERNSYM: Sattel).

**Stufe.** Kandidat, Stufe 1 für P1: Lage der Marken (μ₄, τ = i), Spinstruktur, Viertelzählung
und Normierung c₃, KMS-Zustand. Nicht Stufe 2, weil die Naht selbst (flach, 2D, Fermionen) und
die zwei Blätter mit den vier Marken Eingabe bleiben. Die Blattinvolution der Marken ist eine
π-Drehung um eine Marke; ihr Zentrum liegt um eine halbe Periode neben dem Zentrum der
Vierteldrehung, sie ist also nicht einfach μ₄².

**Basislinie.** Die Rechnung ist elementare Gittergeometrie; die Stärke liegt im Prinzip. Die
Zahlen 4 und 1/(8π) folgen ohne Anpassung; die Gegenproben (Parität zugelassen, keine
Neutralität) geben andere Werte.

**Stärkster Einwand.** „Keine vorab ausgezeichnete Zeit“ ist eine Prämisse (Nähe zur
Thermal-Time-Hypothese von Connes und Rovelli), keine Herleitung. Und ob der KMS-Zustand bei
β = L die operative Saat des Originals ist (Fockzustand der Calderón-Polarisation, Erzeuger
|B_Σ|), ist nicht geprüft.

### N gegen die operative Saat des Originals (`c1_vergleich_zustand.py`)

**Wie das Original den Zustand baut** (KERNSIX20261120 Schritte 1, 3, 4; P1BIND20261120). Die Naht
Σ ist der Raum, die Zeit läuft normal dazu. Kragen D = γ_n(∂_n + B), Polarisation = Calderón-
Projektor P_C, auf den zwei Normalhälften flach. Der Randzustand eines Kragens der Breite ℓ ist
rein; auf einer Fläche ist er thermisch für B bei β = 2ℓ. Breite null gibt den flachen
Kernzustand C = ½, aus dem die aufgezeichneten Zahlen des E₈ Worts stammen. **ℓ ist im Original
frei.** Nullmoden von B bleiben für jede Breite flach; dort lebt das E₈ Wort.

Der Vergleich verlangt, N auf die volle euklidische Geometrie zu heben: thermischer Kreis 2ℓ mal
Naht T²(L, L), Fermionen überall antiperiodisch. N in seiner ersten Form (Zeit innerhalb der Naht)
war eine andere Schnittrichtung derselben Geometrie.

**Befund.**
1. Die Rekonstruktion entlang des thermischen Kreises und die entlang eines Nahtzyklus sind genau
   bei **2ℓ = L** dasselbe System (Spektren β·λ identisch; Abweichung bei ℓ = 0,25…1,0 sonst
   8,7 bis 36). Zeitneutralität **legt die offene Breite des Originals fest**.
2. Damit wird die im Original freie Größe T/λ₁ = 1/(2ℓλ₁) zu **1/(√2 π) = 0,2251** (λ₁ = √2 π/L).
   Die Verschränkung der beiden Kragenflächen ist endlich, **S = 0,5154**, ohne Abschneiden;
   der flache Kernzustand hätte log 2 je Mode.
3. Bei (A, A) hat das geometrische B **keine Nullmoden** (nur (P, P) hat zwei).
4. Lage der Marken: Die mit λ₁ normierte Spektralwirkung hat am Quadrat in der rechteckigen
   Familie ein Minimum, unter Verdrillung einen Sattel (−0,33 bei Re τ = ±0,1), wie KERNSIX
   Schritt 4. N schließt die Verdrillung durch die orientierte Zeitneutralität aus und braucht
   keine Minimierung.

**Urteil: gemischt, mit scharfer Gabel.**
- *Übereinstimmung:* Der zeitneutrale Zustand ist genau der Kragenzustand des Originals (rein,
  auf jeder Fläche thermisch für B), und N liefert die im Original fehlende Breite ℓ = L/2.
- *Abweichung:* Die aufgezeichneten E₈ Zahlen stammen aus dem flachen Zustand der Breite null.
  Den schließt N für jeden Mode aus, der an das geometrische B koppelt, und bei (A, A) gibt es
  keine geometrischen Nullmoden. Die Abbruchregel greift damit für den flachen Sektor, nicht für
  den Zustand als Ganzes.
- *Gabel:* Beide passen nur zusammen, wenn das E₈ Wort auf Nullmoden eines **inneren, nicht an
  die Nahtgeometrie gekoppelten Teils** von B sitzt (B = B_geom ⊕ B_int statt Tensorprodukt).
  Das Original lässt das offen: P1BIND liest B als den endlichen D_F, KERNSIX modelliert B als
  Dirac auf der Naht. Diese Frage entscheidet, ob N trägt.

**Kleinster nächster Test (Kern).** Im Original nachlesen, ob B_Σ die geometrische Naht und den
inneren Teil als direkte Summe oder als Tensorprodukt trägt, und ob die Nullmoden von D_F auch bei
eingeschalteter Nahtgeometrie Nullmoden bleiben. Tensorprodukt mit invertierbarem geometrischem
Teil: N ist mit dem E₈ Wort unverträglich. Direkte Summe: N trägt und fixiert die Breite.
Parallel verfolgt eine andere Sitzung „P1 aus positiver, verlustfreier Nahtverklebung“; beide
Wege zielen auf dieselbe Lücke und sollten vor einer Atlasübernahme abgeglichen werden.

## G: Ein Objekt für P1 und P2 — der Gaußsche Torus ℂ/ℤ[i] (`g1_gauss_torus.py`)

**Prinzip.** Startklasse sind flache 2-Tori, deren Gitter ein Ganzheitsring O_K ist (Tori mit
komplexer Multiplikation, Klassenzahl 1, neun Fälle). Auf einem Torsionsstück des Torus lebt ein
endlicher Oszillator, die Weil-Darstellung auf den Funktionen dieses Stücks. Eine einzige
physikalische Forderung: Auf dem Torsionsstück, auf dem die Symmetrien des Torus als volle
Einheitengruppe eines Primkörpers wirken, hat der Oszillator einen Spin-½-Dublett als
fermionischen Teil.

**Kette.** Spin-½-Dublett ⇒ p = 5 (eindeutig, W) ⇒ Einheitengruppe der Ordnung p − 1 = 4 ⇒ ℤ[i],
der quadratische Torus. Aus diesem einen Torus folgt:
- **P2.** Plätze = Funktionen auf der Torsionsgeraden E[2−i] ≅ 𝔽₅. 3 + 2 = bosonischer Dreier ⊕
  fermionisches Dublett. Y ist die einzige spurfreie Größe im Kommutanten (−1/3, +1/2). Die 16
  trägt die Hyperladungen, jeder Sektor ist unter 2I irreduzibel (W3). SL(2,5) = 2I führt zu E₈
  (Ikosiane, McKay, h = 30).
- **μ₄ und Uhren.** Die komplexe Multiplikation mit i ist auf der 5-Torsion genau das Weyl-Element
  [[0,−1],[1,0]]; ihr Weil-Bild ist die Fourier-Vierteldrehung, also TFPTs R₅ bis auf das
  Vorzeichen. Auf E[2−i] wirkt sie als t ↦ 2t = G, die Torsionstranslation ist C₅, und
  G C₅ G⁻¹ = C₅² gilt exakt.
- **P1.** Randnormierung c = Kovolumen/(2π · |Einheiten|) = 1/(8π) = c₃ (Weylterm der
  Dirac-Wärmespur pro Einheitenbahn, numerisch auf 12 Stellen). Die einzige invariante
  Spinstruktur ohne Nullmoden ist (−,−). Der Modul τ = i ist der Fixpunkt, den KERNSYM als
  symmetrieerzwungen kennt. Die vier Marken sind die Bilder der 2-Torsion unter der
  Doppelüberlagerung y² = z⁴ − 1, mit Windung +1 des quadratischen Differentials (Atlas 23.11.).
  Die Orientierung kommt aus der komplexen Struktur, die Reflexion z ↦ z̄ erhält ℤ[i].
- **Familien.** H¹ von P¹ ohne die Marken: Rang 3, Charaktere {−1, i, −i}, über 𝔽₅
  frobeniusgleich (K), Residuenbild der Trägermoden.

**Auswahltabelle (alle neun CM-Tori der Klassenzahl 1).**
- (A) Volle Einheitenwirkung auf einem Primkörper: p = 5 für ℤ[i], 7 für ℤ[ω], 3 für d = 2 und 11,
  sonst kein p.
- (B) Dublett: nur bei p = 5.
- (C) Invariante invertierbare Spinstruktur: alle außer ℤ[ω].

Nur d = 1 erfüllt alle drei Regeln. Die c-Werte der neun Tori liegen zwischen 0,023 und 0,51;
nur ℤ[i] gibt 1/(8π).

**Neu gegenüber TFPT.** Die Einzelteile waren bekannt: quadratischer Torus mit μ₄ (τ = i,
j = 1728), c₃ aus der Viertelzählung (23.11., bedingt), 5 = (2+i)(2−i) als Zahlenbeobachtung
(ARITH.TRICHOTOMY.01), x ↦ 2x und Y_par, PGL(2,5), 2I und E₈. Neu ist die Brücke: P2s fünf
Plätze sind ein Torsionsstück desselben Torus, dessen Wärmespur P1s c₃ gibt. Und eine einzige
Forderung (Spin-½-Dublett) wählt diesen Torus unter allen CM-Tori aus. Damit beantwortet G beide
Leitfragen der Konsolidierung: **Die Vier ist die Einheitengruppe ℤ[i]^× des Torus**, erzwungen
durch das Dublett. **Die Minimalität ist das erste Torsionsstück**, auf dem die Einheiten voll
wirken.

**Stufe.** Stufe 3 unter benannten Lesarten: P1 (Normierung, Spinstruktur, Modul, Marken,
Orientierung) und P2 (fünf Plätze, 3 + 2, Y) kommen aus demselben Objekt, ohne angepasste
Konstante; die vorhandenen Ergebnisse trägt der TFPT Compiler aus P1 und P2. Nicht Stufe 4:
physische Zeit und Zustand, 3 + 1, Kontinuum und Gravitation liefert G nicht.

**Benannte Lesarten (erste freie Wahlen).**
1. Startklasse: Der Rand ist ein flacher Torus mit Ring-Gitter.
2. Plätze = Schrödinger-Modell der Weil-Darstellung auf dem ersten Torsionsstück mit voller
   Einheitenwirkung.
3. c₃ = Weylterm der Dirac-Wärmespur, gezählt pro Einheitenbahn (in TFPT bisher bedingt).
4. Die Paritätsbasis (Y_par) ist die physische Farb/Schwach-Basis. Das ist eine Identifikations-
   änderung gegenüber der markierten Ladung, siehe W.

**Basislinie.** Die Auswahl ist strukturell, kein Zahlentreffer. Effektiv entscheidet
„Einheitengruppe der Ordnung 4“, und das Dublett erzwingt sie. Der Wert 1/(8π) ist keine
unabhängige Vorhersage, weil die Viertelzählung am quadratischen Torus formuliert wurde.
Unabhängig ist, dass die Regel von der P2-Seite denselben Torus wählt. Gegenproben: ℤ[ω] (sechs
Einheiten) scheitert an zwei Regeln, die Tori mit nur ±1 an der Dublettregel.

**Stärkster Einwand.** Die Lesarten 2 und 3 sind gewählt. SU(3) × SU(2) entsteht nicht, nur 2I
und U(1) × U(1). Die Dynamik fehlt ganz. G erklärt die Eingaben, nicht die Physik dahinter.

**Kleinster nächster Test (Stufe 4, kein Nebenschauplatz).** Zeit und Zustand aus demselben
Torus: Die Vierteldrehung vertauscht die beiden Zyklen (τ ↦ −1/τ hat τ = i als Fixpunkt).
Prüfen, ob der Weil-Oszillator auf E[5] zusammen mit der Torus-Wärmespur eine reflexionspositive
Halbgruppe mit ausgezeichneter Zeitrichtung definiert, deren KMS-Zustand die Normierung c₃
trägt. Abbruch, wenn Zeit oder Zustand eine zusätzliche freie Wahl brauchen.

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

## K: Familien als étale Kohomologie (Nachtrag, `k1_kohomologie.py`)

Test 1 aus Abschnitt F: H¹ von X = P¹ ohne die vier Marken über 𝔽̄_p, mit μ₄-Wirkung, gegen den
TFPT Familienraum (Rang 3, Charaktere {−1, i, −i}). Abbruch bei Abweichung.

**Befund (Punktzählung in explizit konstruierten Körpern 𝔽_{p⁴}, Lefschetz).**
- |X(𝔽_{p^k})| ergibt tr(Frobᵏ | H¹_c). Für p = 5 und 13 (p ≡ 1 mod 4) hat der Frobenius auf
  H¹_c die Eigenwerte **{1, 1, 1}**, auf H¹ also dreimal p (Gewicht 2). Für p = 3 und 7
  (p ≡ 3 mod 4) sind es {1, 1, −1}.
- Getwistete Zählungen bei p = 5: #Fix(Frob ∘ g_a) = 2, 6, 6, 6 für a = 1, 2, 4, 3. Die
  Fourierzerlegung über μ₄ gibt: trivialer Charakter **0**, jeder der drei nichttrivialen
  Charaktere **genau einmal** mit Frobenius 1. Dasselbe bei p = 13.
- **Residuenbrücke.** Die Residuen der TFPT Trägermoden e_m = z^m/(z⁴ − 1) an den vier Marken
  bilden die fünf Moden modulo Konstanten **surjektiv auf den Familienraum** ab (Rang 3). Der
  Kern ist span(e₃, e₄ − e₀), die Markenverschiebung ζ ↦ iζ hat auf dem Bild die Eigenwerte
  −1, i, −i. Über 𝔽₅ ist das Gegenstück die Einschränkung Fun(𝔽₅) → Fun(𝔽₅^×)/Konstanten mit
  Kern span(δ₀, 1).

**Urteil K: Test bestanden, aber nicht unabhängig.** Rang 3 und die Charaktere {−1, i, −i} folgen
schon aus der komplexen Geometrie von P¹ ohne vier Punkte; über 𝔽̄₅ ändert sich daran nichts.
**Neu** sind zwei Dinge:
1. **Arithmetisch gleiche Familien.** Genau wenn i in 𝔽_p liegt (p ≡ 1 mod 4, kleinstes p = 5),
   wirkt der Frobenius auf alle drei Familienklassen gleich. Die Familien unterscheiden sich
   dann nur durch ihre μ₄-Ladung (i, −1, −i), sonst durch nichts. Bei p ≡ 3 mod 4 trennt der
   Frobenius 2 + 1.
2. **Familien = Trägermoden modulo Invarianten.** Die Residuen machen aus den fünf Trägermoden
   genau den Familienraum; zwei Moden tragen keinen Familienanteil. Das beantwortet den Einwand
   aus F teilweise: Die Familien liegen nicht „im selben ℂ⁵“, sie sind dessen Residuenbild in
   H¹.

Weiter ohne Mechanismus: drei gleiche Kopien der 16. Der Befund legt eine Lesart nahe (Materie =
Familienraum ⊗ 16, die drei Kopien arithmetisch gleich und nur durch eine ℤ₄ Familienladung 1, 2, 3
unterschieden), aber diese Lesart ist gewählt, nicht hergeleitet.

## Kleinste entscheidende Tests für eine Folgesitzung

1. **ℤ₄ Familienladung gegen Flavor.** Wenn die Familien die μ₄-Ladungen i, −1, −i tragen, ist
   die Massenmatrix in Potenzen nach Ladungssumme geordnet. Prüfen, ob TFPTs Transportlängen L
   und Potenzen K (oder die Leptonkoeffizienten 16/7, 4/3, 7/6) mit Ladungen 1, 2, 3 mod 4
   verträglich sind. Abbruch, wenn keine Zuordnung der Ladungen die Ordnung der K trägt.
2. **Galois 1 + 2 gegen Flavor.** Prüfen, ob eine TFPT Flavorgröße (Leptonkoeffizienten
   16/7, 4/3, 7/6, Transportlängen) eine Familie gegen zwei auszeichnet, und ob es die zum
   Legendre-Charakter passende ist.
3. **Sechs Fasern.** Welches natürliche Element mit Nichtrest-Multiplikator (etwa J aus der
   Standardform kombiniert mit einer Dilatation) erweitert A₅ auf PGL(2,5) mit Stabilisator F₂₀?
4. **Y_par gegen die markierte Ladung.** Den Befund ‖[H*, Y]‖² = 19/12 gegen 17/6 mit der Weil-
   Lesart neu bewerten: Welche Dynamik ist mit der Paritätsbasis verträglich?
5. **Fockspur-Anschluss.** Wählt die positive Fockspur über Λ ℂ^(𝔽₅) mit der Weil-induzierten
   quadratischen Quelle (Exploration F1 vom 23.11.) dieselbe Orientierung Σ ∝ −6Y?
