# KERNSYM20261120: Spiegelung des Originals und erzwungener kritischer Punkt

Datum: 20. November 2026. Theorievertrag mit begrenzten Urteilen. Nichts wird befördert, kein
Ledger und kein Paper wird geändert. Anschluss an KERNORIG20261120
(`../original-profile-time-modulus-20261120`), dessen drei offene Punkte hier bearbeitet werden:
Vorzeichen des Determinantenterms, Spiegelung Θ gegen die modulare Konjugation, Residuenterm und
andere Spektralfunktionale am elliptischen Überlagerer.

Reproduktion (14 Prüfungen, alle bestanden):

```bash
cd experiments/theory-contracts/original-reflection-and-forced-critical-point-20261120
for s in a_reflection b_forced; do ../../../.venv/bin/python -B check_$s.py; done
```

## Was im Original steht und was nicht

Pfade relativ zu `documents/archive--tfpt-45--source_extracts--`.

- `04_qft_source/sections/008-…tex`, Zeilen 140–253: Form (F,G) = ⟨ΘF,G⟩ auf Observablen
  positiver Zeit, fermionischer Kern „as a reflected quadratic form".
- In den gelesenen Abschnitten steht weder eine Formel für Θ als Operator noch die fermionische
  Wirkung in Z_adm. Das Vorzeichen des Determinantenterms lässt sich am Text deshalb nicht
  entscheiden.
- `03_em_flavor_source/sections/001-…tex`, Zeilen 124–146 und 276–290: Der quartische Teil der
  Alpha Gleichung ist als Gegenterm durch seine Ableitung definiert; die Gleichung wird als
  „logarithmic scale closure" bezeichnet. Eine Herleitung der Potenzen steht dort nicht.

## A: Spiegelung (`check_a_reflection.py`, 8/8)

Fockraum mit drei Moden je Normalhälfte, Kernzustand Ω.

- Die Lagrange Relation G P_C G = 1 − P_C macht die antilineare Feldabbildung a(u) → a*(Gu) zu
  einer Symmetrie der Zweipunktfunktion. Mit jeder anderen relativen Phase zwischen den Hälften
  ist sie es nicht.
- Sie wird durch einen antiunitären Operator Θ realisiert, der Ω festhält.
- Θ bildet die Felder einer Normalhälfte auf dieselbe Hälfte ab. Der Spiegel J aus P1BIND bildet
  a(e) auf a*(Ue) der anderen Hälfte ab. Θ und J sind verschiedene Operatoren.
- Als ordnungserhaltende Abbildung gibt Θ keine positive Form (Wert −1/4 für F = a₁a₂). Mit der
  Ordnungsumkehr der fermionischen Rekonstruktion ist Θ(F) = F* auf einer Hälfte, und
  ⟨Ω, F*GΩ⟩ = Tr(F*G)/2^d.

Urteil: Die Zeitspiegelung bei Zeit null ist nicht die modulare Konjugation. Sie ist die zweite
antiunitäre Symmetrie des Kernzustands, und es gibt sie genau wegen der Lagrange Relation. Auf
jeder Normalhälfte stimmen die reflektierte Form des Originals und die Spiegelform aus P1BIND
überein, weil der Zustand dort die Spur ist. Die Annahme „zwei Seiten = zwei Normalhälften" wird
für die Form nicht gebraucht. Offen: Θ ist hier aus der Symmetrie konstruiert, nicht aus einer
Formel des Originals abgelesen.

## B: Erzwungener kritischer Punkt (`check_b_forced.py`, 6/6)

Wärmesummen H(τ, s) = Σ exp(−s λ²) über das Dirac Spektrum des flachen Torus der Fläche eins.
Jedes Spektralfunktional (Zeta Familie, Residuum, endlicher Teil, Determinante mit Masse) ist ein
Integral von H über s.

- Antiperiodische Spinstruktur: Der Gradient von H verschwindet bei τ = i für jedes s. τ = i ist
  Fixpunkt von τ → −1/τ, und diese Abbildung erhält die Spinstruktur.
- Am hexagonalen Punkt verschwindet er nicht. Bei periodischer Spinstruktur verschwindet er an
  beiden Punkten.
- Für jedes s ist H am Quadrat entlang der Spiegelfamilie am größten und quer dazu am kleinsten,
  und entlang der Spiegelfamilie monoton.
- Für kleine s ist s·H = 1/(2π) = 4c₃ je Flächeneinheit, unabhängig von τ: Der Residuenterm
  trägt keinen Modul.

Urteil: Das Quadrat ist kritischer Punkt jedes Spektralfunktionals der antiperiodischen
Spinstruktur, erzwungen durch Symmetrie, und es ist der einzige so erzwungene Punkt. Für jedes
Funktional mit positivem Gewicht in s ist es ein Sattel fester Orientierung ohne weiteren
kritischen Punkt auf der Spiegelfamilie. Ein Vorzeichenwechsel vertauscht nur die Richtungen.
Damit sind drei offene Punkte aus KERNORIG für die Lage der Marken unerheblich: das Vorzeichen
des Determinantenterms, die Masse μ und der Residuenterm bei s = 1. Es bleiben die zwei
Voraussetzungen: elliptischer Überlagerer und antiperiodische Spinstruktur. Die Referenz D_ref
ist nicht nachgebildet; sie hängt nach Definition nicht vom Modul ab.

## Stand der Kette

1. Randzustand = Fockzustand der Polarisation: bedingt auf die Reflexionspositivität des
   fermionischen Kerns, die das Original voraussetzt.
2. Spiegelform Tr(F*G)/2^d: aus dem Kernzustand, jetzt auf zwei Wegen (modulare Konjugation und
   Zeitspiegelung).
3. Zeit: rekonstruierte Zeit mit Erzeuger |B_Σ| und Lücke λ₁.
4. Marken bei μ₄: erzwungener kritischer Punkt; bedingt auf Überlagerer und Spinstruktur.
5. Wörterbuch zum E₈ Wort: Zustand, Form, Marker folgen; die Aufteilung 5 + 3 der acht Nullmoden
   ist Gegenstand von NULLMODEN20261120.
6. Kubik der Alpha Gleichung: im Original gesetzt, hier nicht hergeleitet.
7. Gravitation: nicht hergeleitet.

F1 bleibt in Bearbeitung.
