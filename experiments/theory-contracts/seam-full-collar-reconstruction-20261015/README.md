# Vollständige Nahtdaten: Operator, Zustand und Zeit gemeinsam rekonstruieren

**Ergebnis:** Eine falsche Identifikation im archivierten Ursprungsbeweis ist konkret
korrigiert. Für einen gegebenen exakten Produktkragen liefern die vollständigen Daten
beider Randflächen den vorzeichenbehafteten Operator eindeutig zurück. Unter der
zusätzlichen geladenen quasifreien CAR-Realisierung folgen daraus auch Kovarianz
und Zeitwirkung. Die physische Auswahl dieser Quelle aus P1/P2 ist damit nicht bewiesen.

Der Unterschied ist anschaulich: Eine einzelne Vorzeichenprojektion sagt, welche
Moden in welche Richtung gehören. Die vollständige Verbindung beider Randflächen
enthält zusätzlich, wie stark jede Mode übertragen wird und auf welche Felder sie
wirkt. Diese Information darf vor der Rekonstruktion nicht weggeworfen werden.

## Was korrigiert wird

In `_archive/tfpt-45/source_extracts/01_boundary_kernel_source.tex:130–136,154–166`
wird die geometrische Deckspiegelung mit `2 P_C − I` gleichgesetzt. Ein expliziter
kompakter Dirac-Kragen widerlegt diese allgemeine Identifikation. Die Spiegelung,
der Normal-Cliffordoperator und die Cauchy-Polarisierung werden nun getrennt geführt.
Eine algebraische Vierteldrehung ist keine Herleitung der geladenen Viertelholonomie.
Das Archiv bleibt als Originalbeleg unverändert.

Außerdem behält die ursprüngliche Zulässigkeitsprojektion nur den positiven
Polarisationsraum. Wird danach dieselbe Polarisation auf den behaltenen Träger
eingeschränkt, ist sie überall `+1`: Der negative Dreierblock der späteren
`3+2`-Markierung ist bereits entfernt. Der Beweis in `PROOF.txt`, Abschnitt 8,
zeigt das ohne numerische Näherung. Die interne Trägermarkierung muss deshalb vor
der Zustandskompression mit ihrer tatsächlichen Einbettung erhalten bleiben.
Auch `2G−I` auf `Ran(G)` wäre wieder identisch `+1`; der neue Graph darf nicht
als Vorwand für dieselbe Verwechslung dienen.

## Was konstruktiv folgt

Für die tatsächliche Produktgleichung `(∂r+B)ψ=0` und bekannte Länge `ell` lautet
die vollständige Randrelation `v=T u`, `T=exp(−ell B)`. Ihr orthogonaler Projektor ist

```text
G = [ A      Cross ]       A     = (1 + exp(−2 ell B))⁻¹
    [ Cross  1−A   ]       Cross = (2 cosh(ell B))⁻¹
```

Daraus folgt `B = log(A/(1−A))/(2 ell)`. Der Checker erhält kein `B` als Eingabe.
Bei gegebenem geladenem quasifreiem Feldzustand mit Kovarianz `G` ist die
Kovarianz der zweiten Fläche `n=1−A`; sie liefert denselben Generator über
`log((1−n)/n)=2 ell B`. Die ganze Operatorwirkung wird im gegebenen Rahmen
rekonstruiert, nicht bloß eine Eigenwertliste.

Die Logit-Inversion war im Atlas bereits bekannt. Neu in diesem TFPT-Arbeitsweg
sind die konkrete Verbindung zur vollständigen Kragen-Randrelation, die geprüfte
Korrektur des geometrischen Operators und der ausführbare inverse Datenpfad.

## Ausführen

Vom TFPT-Repository aus:

```sh
python3 experiments/theory-contracts/seam-full-collar-reconstruction-20261015/checker.py check
python3 experiments/theory-contracts/seam-full-collar-reconstruction-20261015/checker.py verify
python3 experiments/theory-contracts/seam-full-collar-reconstruction-20261015/checker.py reconstruct --kernel experiments/theory-contracts/seam-full-collar-reconstruction-20261015/input-example.json --output /tmp/tfpt-collar-reconstruction.json
python3 experiments/theory-contracts/seam-full-collar-reconstruction-20261015/checker.py audit-carrier --kernel experiments/theory-contracts/seam-full-collar-reconstruction-20261015/input-example.json --output /tmp/tfpt-carrier-audit.json
```

Die Eingabe enthält eine rationale Projektormatrix, zwei benannte Randflächen mit
Modenlisten, die ausdrückliche Deklaration
`fiber_identification: "ordered_shared_orthonormal_frame"` und die unabhängig
gegebene positive Normallänge. Ohne den gemeinsamen orthonormalen Faserrahmen
ist die behauptete Selbstadjungiertheit zwischen zwei Flächen nicht definiert.
Der Prüfer lehnt
reine Vorzeichenprojektoren, nichtpositive Transfers und ungültige Projektionen ab.
Das Beispiel ist ausdrücklich eine Diagnosematrix, kein neu gefundener P1-Kern.
Alle Rechnungen sind klein und exakt; GPU-Suche bringt für diese Inversion keinen Vorteil.

Der Träger-Audit prüft eine gelieferte isometrische Einbettung `W` und interne
Involution `E` gegen die tatsächliche Gleichung `(2G−I)W=WE`. Wenn `GW=W`,
erzwingt dies `E=I`; ein negatives Etikett wird mit diesem Grund abgewiesen.
Ein korrektes gegebenes Wörterbuch wird als geprüft, ausdrücklich nicht als
aus P1/P2 hergeleitet ausgegeben.

## Was als Nächstes tatsächlich gebraucht wird

Die ursprüngliche Quelle muss beide Randspuren samt Vergleichsabbildung und Skala
liefern. Aus einem schon gegebenen `D_+` darf man den endlichen Kragen ausschneiden;
seine Rückrechnung liefert aber keine unabhängige Auswahl von `D_+`.
Die vorhandenen P1-Definitionen geben den Operator beziehungsweise den operativen
Seed bereits vor. Der neue Weg verhindert, dass zusätzlich zur Quellenwahl noch
Zeitgenerator und Zustand unabhängig geraten werden.

Eine ursprüngliche nichtkommutierende Markierung benötigt die lokale Algebra oder
Windungswirkung; sie folgt nicht aus einer bloßen Funktion des Spektrums von `B`.
P2-Träger, E8-Erweiterung, physisches Ereignisgesetz und T1–T8 werden hier nicht als
geschlossen verbucht. Der nächste Nachweis muss diese tatsächlichen Operatoren
aus derselben Quelle liefern, nicht weitere spektral passende Hilfsmodelle.

## Belege und Grenzen

- `PROOF.txt`: analytischer Beweis, kompakter Gegenfall, Definitionsbereiche und Inputs.
- `checker.py` / `results.json`: exakte Rechnungen, Eingabeprüfung und Replay.
- `lean/FullCollar.lean` / `lean-receipt.json`: kernel-geprüfte rationale Graphalgebra
  und antikommutierende Involutionen; keine Formalisierung der unbeschränkten Analysis.
- `contract_index.json`: semantische Verknüpfung mit der bisherigen Transferrekonstruktion
  und dem offenen gemeinsamen Quellenproblem.

Eine unabhängige Prüfung des vollständigen Beweises bestätigte die Operatorwirkung,
den kompakten Gegenfall, die Abschlussdomäne und die bedingte CAR/KMS-Aussage.
