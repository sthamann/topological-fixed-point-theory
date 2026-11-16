# F1: Zwei elliptische Doppeldecker als GENUS-FIVE-Test für eine Vorstruktur größer als E8

**16. November 2026. Exakte topologische Konstruktion und unabhängige Rootkontrollen; Literaturtheorem der höheren Thetafunktionen. Physische F1 NICHT gelöst.**

## 1. Ausgangspunkt: zwei unabhängige Tori vor dem postulierten P1

Setze ohne E8- oder SM-Vorgabe zwei verschiedene glatte elliptische Doppelüberlagerungen derselben abstrakten Riemannschen Sphäre an:

```text
E_a: y^2 = f_4(x),   E_b: z^2 = g_4(x)
C:   y^2 = f_4(x), z^2 = g_4(x).
```

`f_4` und `g_4` sind Polynome vom Grad 4 mit jeweils vier einfachen Nullstellen. Die hier prüfbare Originalgeometrie benutzt z.B. die **disjunkten** Branch Sets
`B_a={-4,-3,-2,-1}`, `B_b={1,2,3,4}`. Beide `E_a` und `E_b` haben Genus eins. Das normierte Faserprodukt `C` ist zusammenhängend und besitzt die `(Z2)^2` Deckgruppe.

Ist `k=|B_a intersect B_b|<4`, so besitzt die gesamte 4:1-Abbildung `C->P1` genau `8-k` verschiedene einfache Verzweigungswerte. Über jedem von ihnen tragen zwei der vier Blätter Verzweigungsindex zwei. Riemann-Hurwitz ergibt

```text
2g(C)-2 = 4*(-2) + 2*(8-k) = 8-2k
g(C) = 5-k.
```

Unabhängig bestätigt die Projektion `C->E_a`, eine 2:1-Abbildung, genau `2(4-k)` Verzweigungspunkte. Weil `g(E_a)=1` ist, gilt `2g(C)-2 = 2(4-k)`, dieselbe Formel.

| k | g(C) | Parent-Auslese |
|---|---:|---|
| 0 | **5** | Erster interessanter Genus-Sektor |
| 1 | 4 | Schottky-Genus-4-Auslöschung |
| 2 | 3 | gleiche Gitter-Thetas |
| 3 | 2 | gleiche Gitter-Thetas |
| 4 | nicht verbunden | zwei gleiche elliptische Komponenten |

Wichtig: Die Konstruktion zweier *unabhängiger* elliptischer Doppelüberlagerungen ist nicht dasselbe wie eine unramifizierte 2:1-Überlagerung E->E oder der Orientierungsdoppeldecker einer Kleinschen Flasche.

## 2. Zwei mögliche Rank-16-Parentgitter sind auf einer einzelnen Torusfläche UNSICHTBAR verschieden

Die zwei geraden positiv definiten unimodularen Rang-16-Gitter sind
`L_A=E8+E8` und `L_B=D16+`. Beide besitzen 480 Norm-zwei-Wurzeln und 16 Cartanrichtungen, also 496 affine Stromdimensionen. Dennoch sind sie nicht isomorph:
- `E8+E8`: Der nichtorthogonale Wurzelgraph zerfällt in **240+240**.
- `D16+`: Der Wurzelgraph ist **verbunden mit 480**; die Wurzeln erzeugen das D16-Untergitter von Index zwei, nicht das ganze D16+.

Diese Unterschiede wurden durch **wirkliche exakte ganzzahlige Rootvektoren** reproduziert. E8 wird als 112 Koordinatenwurzeln plus 128 Halbspinorwurzeln konstruiert; beide Eltern besitzen 480 Norm-zwei-Wurzeln, nicht aber denselben Wurzelgraphen.

Die **Genus-eins** Thetafunktion ist dennoch bei beiden durch Modulformen von Gewicht acht identisch: `Theta_L=E4^2`. Die Literatur zeigt stärker, dass die Gitter-Thetas **auf Jacobimatrizen echter Riemannscher Flächen bis Genus 4** gleich bleiben: `Theta_(E8+E8)^(g) = Theta_(D16+)^(g)` für `g<=3` auf dem ganzen Siegelraum; für `g=4` verschwindet ihre Differenz `J_8` auf dem Jacobian locus (Schottky-Igusa).

**Bei Genus 5** dagegen ist die Differenz

```text
F_5(Omega) = Theta_(E8+E8)^(5)(Omega) - Theta_(D16+)^(5)(Omega)
```

auf echten Jacobis nicht identisch null. Der veröffentlichte Satz von Grushevsky und Salvati Manni präzisiert, dass ihr Nullstellenort auf dem glatten Genus-5-Jacobian locus genau den TRIGONALEN Kurven entspricht.

## 3. Der konstruktive Test schließt den trigonalen Ausfall aus

Für DISJUNKTE `B_a,B_b` hat `C` Genus 5 und einen 2:1-Morphismus `C->E_a`. Sie ist somit **bielliptisch**. Angenommen, sie hätte gleichzeitig eine 3:1-Abbildung auf `P1` (trigonal). Weil die Grade zwei und drei teilerfremd sind, wäre der kombinierte Morphismus birational auf das Bild in `E_a x P1`. Castelnuovo-Severi liefert dann

```text
g(C) <= 2*g(E_a) + 3*g(P1) + (2-1)*(3-1)
     = 2*1 + 3*0 + 2 = 4,
```

im Widerspruch zum exakt konstruierten `g(C)=5`.

Für den vollständigen Schluss aus dem Brill Noether Divisor ist auch
Hyperelliptizität auszuschließen: Eine zusätzliche 2:1-Abbildung auf P1
könnte nicht durch denselben Grad-zwei-Zwischenkörper wie `C->E_a`
faktorisieren, weil die beiden Zielkurven Genus null beziehungsweise eins
haben. Beide Abbildungen wären daher unabhängig. Castelnuovo-Severi ergäbe
`g(C)<=2*1+2*0+(2-1)*(2-1)=3`, erneut ein Widerspruch. Somit hat C auch
kein Grad-drei-Büschel mit Basispunkt aus einer hyperelliptischen Struktur.
Der trigonale Abschluss im Nullstellensatz ist damit vollständig ausgeschlossen.

**Konditionales mathematisches Ergebnis:** Wenn die beiden Parentoptionen tatsächlich durch ihre gewöhnlichen Rang-16-Gitter-Thetapartitionen auf der Jacobischen Periodenmatrix des hier konstruierten glatten `C` ausgelesen werden, dann sind sie bei `k=0` unterscheidbar:

```text
F_5(Omega_C) != 0.
```

Bei `k=1` entsteht dagegen Genus 4 und der mathematische Genus-4-Jacobian-Thetatest gibt `F_4(Omega_C)=0`. **Das gemeinsame Ausgangsobjekt aus zwei Torusquellen erhält hier einen scharfen topologischen An/Aus-Schalter für die Elternunterscheidung.**

Wichtige Präzisierung: Der nichtverschwindende Genus-5-Thetawert wird hier aus einem PUBLIZIERTEN SATZ und der bewiesenen Nichttrigonalität geschlossen; er wurde NICHT numerisch für die konkreten Polynome durch Berechnung der kompletten Periodenmatrix und der Rang-16-Siegel-Thetas verifiziert. Die exakten eigenen Rechentests decken Verzweigungen, Genus, Castelnuovo-Severi und Wurzelgraph ab. Ein numerischer Wertebeleg ist eine weitere Testaufgabe.

## 4. Warum das F1 und P1 noch NICHT löst

Das wäre eine ungewöhnlich konkrete mathematische Möglichkeit, einen vor P1 liegenden größeren Ursprung zu unterscheiden. Es beweist ABER weder, dass die ursprüngliche TFPT Quelle überhaupt ein `c=16` Parentgitter trägt, noch, dass die beiden Torusüberlagerungen physisch entstehen, noch, dass das realisierte Gitter durch diese Observable ausgewählt wird.

Vier weitere unabhängige Verpflichtungen bleiben:
1. **P1 Lesekopf:** Eine von der Originalquelle festgelegte Auslese `C->P1` mit echter positiver `Q0/P_adm` Wirkung, nicht eine beliebig angenommene Theta-Auswertung.
2. **E8_1:** Die naive diagonale Summe zweier `E8_1` Ströme hat `k=2`; eine einzelne physisch selektierte E8_1 Kopie ist eine extra Operation. `D16+` enthält dagegen E8 nicht automatisch als gewählte Kac-Moody-Algebra vom Niveau eins.
3. **Dynamik/Zeitrichtung:** Weder die genus-5 Geometrie noch die Schottky-Form liefert den physikalischen OS-positiven Zeitgenerator oder 3+1D Lokalität. Eine zugrunde liegende positive gemeinsame Quellform muss die tatsächlichen Zustände und ihre Zeitentwicklung erzeugen.
4. **Originale R113 Antwort:** `whole_source_bundle=null`. Unmarkiertes absolutes `Z0`, sämtliche markierten `ZU,ZY,ZUY`, original `P_adm`, Higgs/Gaugesektor, Spinororientation und wirkliche zusammengesetzte 2-Stufen-Projektion fehlen.

Die mathematische hier gezeigte Alternative nimmt sogar ZWEI ursprüngliche elliptische Deckquellen als Startannahme; sie ist nicht eine ab initio Auswahl dieser Geometrie. Auch die beobachtete Parent-Unterscheidbarkeit ist keine Quellenwahrscheinlichkeit für A oder B.

## 5. Reproduktion und Quellen

Ausführbar mit Python-Standardbibliothek:
`cd research/originlab && python -S -m unittest -v test_f1_double_torus_genus5_parent`.

**6 von 6 endlichen Tests lokal bestanden**. Code:
[f1_double_torus_genus5_parent.py](f1_double_torus_genus5_parent.py) und
[test_f1_double_torus_genus5_parent.py](test_f1_double_torus_genus5_parent.py).

Wissenschaftliche Grundlagen:
- Henriksson, Kakkar, McPeak, *Classical codes and chiral CFTs at higher genus*, JHEP (2022) 159: https://link.springer.com/article/10.1007/JHEP05%282022%29159
- Grushevsky, Salvati Manni, *The superstring cosmological constant and the Schottky form in genus 5*, American Journal of Mathematics 133 (2011): https://arxiv.org/abs/0809.1391
- Cris Poor: https://faculty.fordham.edu/poor/publications.html
- Bielliptic Genus ≥5 ist nicht trigonal, Castelnuovo-Severi: https://link.springer.com/article/10.1007/s10231-025-01643-3
- Ursprüngliche TFPT Q0 Quelle und einzelne Energieform: [F1-DIRECT-SINGLE-ORIGINAL-SOURCE-20261116.md](F1-DIRECT-SINGLE-ORIGINAL-SOURCE-20261116.md)
- Erster spektraler Deck Lesekopf: [F1-P1-TORUS-DECK-SPECTRAL-READER-20261116.md](F1-P1-TORUS-DECK-SPECTRAL-READER-20261116.md)

**Status: exakte nichttriviale Vergleichsgeometrie, keine Auswahl des physischen Elternzustands, F1 OFFEN.**
