# P1BIND20261120 · Welchen Marker, welchen Schnitt und welchen Spiegel liest der originale P1 Kern?

Datum: 20. November 2026. Atlasrevision beim Lesen `b5b809c43458554296b079da`. Stand des Repos gleich `origin/main`.
Auftrag: die eine benannte Annahme aus FOURMARKS20261118 (die ursprüngliche Randform ist ein Vektorzustand auf zwei
Seiten der Naht mit dem aufgezeichneten Spiegel) am originalen Randkern prüfen und Marker, Schnitt und Spiegel dort
festmachen. Drei Skripte, 54 von 54 Prüfungen. Kein Ledgerwechsel, keine Beförderung. `physical_solution = false`.
F1 bleibt in Bearbeitung.

## Was im Original steht

Quelle: `documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/001-…tex`, Zeilen 42 bis 177.

- Das einseitige Randdatum ist `(A₊, H₊, D₊, J, Γ, B_Σ)`. Im Kragen gilt `D₊ = γ_n(∂_n + B_Σ)`.
- Die Polarisation ist der Calderón Projektor `P_C` des inneren Randproblems, `ι_C = 2P_C − 1`.
- Die Verdopplung `M₊ ∪ M₊^op` liefert die Deckinvolution `τ`. Satz (iii) schreibt `ι_C = R_τ`.
- Die Saatantwort ist `1/λ₁⁺(|B_Σ|)`.
- Das Korollar danach verbietet, geometrische Spiegelung, Polarisation und Trägersplit wörtlich gleichzusetzen.

Schon im Atlas: Für die vollständige Dirac Cauchyprojektion gilt `N C N* = 1 − C` (Lagrange Lemma, mehrere Befunde vom
9. und 10. November). Das wird hier benutzt, nicht neu behauptet. Neu ist, was daraus für Zustand, Spiegel und Marken folgt.

## Ergebnis in vier Sätzen

1. **Zustand.** Die originale Polarisation ist auf den zwei Normalhälften (Eigenräume von `G = −iγ_n`) genau flach:
   `P_C = ½ [[1, U*], [U, 1]]` mit unitärem `U`. Jede Hälfte für sich trägt den Spurzustand, beide zusammen sind rein und
   voll korreliert. Das gilt für jedes Innere und jede Schnittfläche, nicht nur für die runde Kugel.
2. **Spiegel.** Der Fockzustand dieser Polarisation hat auf einer Normalhälfte die reflektierte Form `Tr(F†G)/2^d` mit
   Referenz 1. Sein Spiegel ist die modulare Konjugation `J` mit `J² = 1`. Auf den Feldern wirkt sie als
   `J a(e) J = (−1)^N · a*(U e)`: Seitentausch durch den Kreuzblock von `ι_C` plus Adjunktion. Das ist der einfache Spiegel
   der Sitzung vom 18. November. Die Realstruktur des Datums wird dafür nicht gebraucht. Der Seitentausch ohne Adjunktion
   ist nicht positiv.
3. **Marker und Schnitt.** Der Kern enthält genau zwei antikommutierende Involutionen, `G` und `ι_C`. Das Paar senkrecht
   zur Polarisation ist `(U, Y) = (G, iGι_C)`, bis auf eine Drehung um `ι_C`. Auf jeder Linie der Seite `η` gilt
   `Z_0 = 1`, `Z_U = Z_Y = 0`, `Z_UY = iη`, `ABAB = −1`. Die gemischte Marke ist die Polarisation selbst: `G·Y = i·ι_C`.
   Im Produktkragen ist `Y` die Phase des tangentialen Teils `γ_n B_Σ` von `D₊`, bis auf das Vorzeichen.
4. **Der aufgezeichnete native Operator ist ein Fall davon.** Mit `G` = aufgezeichnete Graduierung, `B` = `D_F`,
   `ι_C` = `sign(D_F)` ist die aufgezeichnete Richtung des Potentials exakt `−iG·sign(D_F)`, jede Graduierungshälfte
   trägt die flache Dichte 1/48, und `λ₁⁺` ist nicht entartet: die Saatlinie des Kerns ist eine der 12 Kappenlinien.

## Rechnungen

### 1. Einteilchenebene (`check_polarization.py`, 22/22)

- Symmetrie von `D₊` verlangt `γ_n B + B γ_n = 0`. Die Randform `⟨γ_n u, v⟩` ist dann entlang der Normalen erhalten; mit
  einem geraden Anteil in `B` nicht (Fluss 17,2 im Beispiel).
- Geschichtetes Inneres (40 Schichten, jede ungerade, kein Produkt): die Cauchydaten der abklingenden Lösungen sind ein
  Projektor vom halben Rang, verschieden vom Spektralprojektor des Randoperators (Abstand 0,50), und exakt lagrangesch
  (Defekt 7·10⁻¹⁵).
- Daraus: beide Diagonalblöcke sind ½, der Kreuzblock ist unitär durch 2. Kontrolle: 200 zufällige Polarisationen vom
  halben Rang sind nie flach (kleinster Defekt 0,22). Umkehrung: jeder Projektor, der auf beiden Hälften flach ist, ist
  lagrangesch. Die positive Hälfte eines Operators mit geradem Anteil ist nicht flach (Defekt 0,32).
- Die Cauchydaten der gespiegelten Seite sind `γ_n Ran(P_C) = ker(P_C)`.
- Marken auf 40 zufälligen Linien je Seite und unter Drehung von `(U, Y)`: Abweichung höchstens 8·10⁻¹⁵.
- Für jeden Marker `c`, der mit `G` vertauscht, auf dem flachen Zustand einer Seite: `Z_c = ⟨c⟩`, `Z_Y = 0`,
  `Z_cY = iη⟨Gc⟩`. Das ist das Gesetz der Runde 45 (gemischte Marke gleich chiraler Mittelwert des Markers), jetzt als
  Folge der Kernstruktur.
- Modell mit dem Spektrum der runden S² (Eigenwerte ±(k+1), Vielfachheit 2(k+1)): `λ₁⁺ = 1`, positive Hälfte flach.

### 2. Vielteilchenebene (`check_fock.py`, 23/23)

Je `d = 1, 2, 3` Moden pro Normalhälfte, zufälliges unitäres `U`, Jordan Wigner.

- Der Fockzustand hat die Zweipunktfunktion `P_C`. Seine Einschränkung auf die erste Hälfte ist die Dichte `1/2^d`, auf
  die zweite Hälfte die Spur auf allen `4^d` Monomen.
- `Ω` ist zyklisch für die Algebra einer Hälfte. `J(XΩ) = X*Ω` ist antiunitär, `J² = 1`, `JΩ = Ω`. Weil `J` schon
  antiunitär ist, ist der Modularoperator 1: der Zustand hat keinen eigenen Fluss.
- `⟨Ω, (JFJ) G Ω⟩ = Tr(F*G)/2^d`, Referenz 1. `J A₊ J` vertauscht mit `A₊`.
- `J a(e_k) J = (−1)^N a*(L e_k)` mit `L = U` bis auf eine gemeinsame Phase. Die drei anderen Formen (ohne Parität, ohne
  Adjunktion) haben Restfehler von Ordnung 1.
- Der Seitentausch ohne Adjunktion (`a → a`, antilinear fortgesetzt) nimmt auf hermiteschen `F` beide Vorzeichen an
  (−7,7 bis +6,3 bei `d = 2`).

### 3. Der aufgezeichnete native Operator (`check_native.py`, 9/9)

Unveränderter Assembler `native-diagnostic-v1`, 96 Fasern.

- `D_F` ist selbstadjungiert, invertierbar (kleinster Betrag 0,0852) und ungerade unter der Graduierung (48 + 48).
- `sign(D_F)` antikommutiert mit der Graduierung; die positive Hälfte ist auf jeder Graduierungshälfte exakt ½.
- `sign(−i g D_F)`, die aufgezeichnete Richtung des Potentials, ist exakt `−i g sign(D_F)`.
- 12 eindeutige Eigenlinien (sechs Werte, beide Vorzeichen: 0,0852; 1,1676; 1,3931; 2,2901; 4,1741; 5,5684).
  `λ₁⁺ = 0,0852` ist nicht entartet.
- Die aufgezeichnete Realstruktur `K` lässt die Polarisation fest und tauscht die Graduierungshälften. Unter ihr ist der
  Schnitt gerade und der Marker ungerade. Das passt zur Rollenteilung der Runde 43: `K` hält den Zustand fest, der
  einfache Spiegel reflektiert.

## Was sich an F1 ändert

Vorher (Runde 47): F1 ist geschlossen unter der Annahme, die ursprüngliche Randform sei ein Vektorzustand auf zwei
Seiten mit dem aufgezeichneten Spiegel.

Jetzt folgt diese Form aus dem Originaldatum unter drei schwächeren, benannten Voraussetzungen:

- **(a)** Der Randzustand ist der Fockzustand der originalen Polarisation `P_C` (die übliche Bedeutung von
  „Polarisation“; im Original steht das Wort, nicht dieser Satz).
- **(b)** Die „zwei Seiten“ der Form sind die zwei Normalhälften `G = ±1`, nicht die Cauchyräume der zwei geometrischen
  Seiten. Auf den Cauchyräumen ist derselbe Zustand ein Produkt (gefüllt und leer).
- **(c)** Die Cauchydaten sind lagrangesch (Standardsatz für Diracoperatoren mit eindeutiger Fortsetzung; im
  geschichteten Modell exakt geprüft, im Atlas als Lemma geführt).

Unter (a) bis (c) sind Form, Referenz, Spiegel, Marker, Schnitt und Vorzeichen bestimmt. Die 60 aufgezeichneten
Rückkehrwerte, die in Runde 28 und 40 die flache Dichte erzwangen, werden dafür nicht mehr gebraucht; die Flachheit
kommt aus der Polarisation.

## Hinweis zum Originaltext

Satz (iii) schreibt `Ran(P_C) = ker(R_τ − 1)`. Liest man `R_τ` als punktweisen Lift der Spiegelung auf Randdaten, dann
ist `R_τ` bis auf eine Phase `γ_n`, also `G`. `G` tauscht die Cauchyräume der zwei Seiten, und ein Eigenprojektor von `G`
ist nie lagrangesch. In dieser Lesart kann (iii) nicht gelten. Der Satz stimmt, wenn `R_τ` als die Involution gelesen
wird, die auf inneren Daten +1 und auf Daten der gespiegelten Seite −1 ist; dann ist er die Definition von `ι_C`. Das
Korollar des Originals (drei Involutionen nicht gleichsetzen) passt zu dieser zweiten Lesart. Am Papier ist nichts
geändert.

## Was nicht gezeigt ist

- Voraussetzung (a). Sie ist eine Lesart des Wortes Polarisation.
- Das Wörterbuch zum E₈ Wort. Der primitive Kern enthält keinen Trägersplit und keine Hyperladung; Higgsquadrat und
  Hyperladung können deshalb auf dieser Stufe nicht Marker und Schnitt sein. Wie die Marken des E₈ Worts nach dem
  Trägersatz aus `G` und `ι_C` entstehen, ist offen.
- Die Zahlenwerte der gewichteten R113 Form (1/4 auf 128, 1/6 auf 96 Zuständen). Sie hängen am E₈ Wort.
- Eine Zeit mit Länge. Der Zustand ist auf jeder Hälfte die Spur, sein Modularoperator ist 1. Aus dem Kernzustand allein
  folgt keine Zeit; das bestätigt „Dauern sind Reihenfolgen“ (Runde 40) von der Quelle her.
- Die Lage der vier Marken und der Wert c₃.
- Die Gleichsetzung der 96 nativen Fasern mit Randmoden. Geprüft ist, dass der native Operator dieselbe Struktur hat.

## Nächster Test

Voraussetzung (a) am Original entscheiden: die primitive Wirkung enthält `tr f(D_rel/χ) − tr f(D_ref/χ)` mit
`D_ref`, `D_rel` als geradem und ungeradem Teil unter der Deckinvolution. Zu rechnen ist, ob der Zustand, den diese
relative Wirkung auf den Randdaten erzeugt, der Fockzustand von `P_C` ist, und ob `D_rel` auf den Randdaten der Schnitt
`Y` ist. Fällt das, bleibt (a) eine Annahme.

## Wiederholen

```
cd experiments/theory-contracts/p1-kernel-marker-cut-mirror-20261120
for s in polarization fock; do ../../../.venv/bin/python -B check_$s.py; done
../../../.venv/bin/python -B check_native.py ../../..
```
