# Alpha Herkunft nachgeprüft, drei Ursprungsideen und Twistor/Hopf getestet

Datum: 19.11.2026. Ergebniskennung: ORIGINIDEAS20261119. Stand des Repos: origin/main 9c695949.
Kein Ledgerwechsel, keine Papieränderung, kein Commit. physical_solution=false.

## Frage

1. Stimmt die Aussage vom 18.11. (RADICALCHRONO), der Koeffizient 5/4 der Alpha Gleichung sei erst nach dem Vergleich mit dem Messwert gesetzt worden?
2. Tragen drei neue Ursprungsideen (Reed Muller Code, Kopplungsdualität, Gittermittel) eine Rechnung?
3. Liefern Hopf Hebung und Twistorfaserung die Konstanten der Alpha Gleichung?

## Ergebnis in einem Satz je Teil

| Teil | Skript | Checks | Urteil |
|---|---|---|---|
| Alpha Herkunft | check_alpha_provenance.py | 13/13 | Aussage vom 18.11. zurückgenommen: kein Fit |
| Reed Muller | check_reed_muller.py | 18/18 | exakte Struktur, teilweise tragfähig |
| Kopplungsdualität | check_duality.py | 7/7 | als Flusssumme ausgeschlossen |
| Gittermittel | check_debye_waller.py | 8/8 | strenge Lesart ausgeschlossen (+142 sigma), Gaußanteil trifft exakt |
| Hopf und Twistor | check_hopf_twistor.py | 12/12 | drei Volumentreffer, drei Fehlstellen, nicht eindeutig |

## 1. Alpha Herkunft (ALPHAPROV20261119)

Quellen: documents/archive--paper-latex--* (Fassungen 2.1 bis 2.8), Atlas finding:alpha-archive-term-provenance-and-master-functional-scope und finding:marked-e8-sector-grade-response-realizes-five-quarters.

Abfolge der archivierten Fassungen, alle Werte nachgerechnet:

| Fassung | Datum | Rückkopplung | alpha^-1 | Abstand zu CODATA 2022 |
|---|---|---|---|---|
| 2.1 | 21.01.2026 | delta (1 - 2 alpha) | 137.035990390 | -418 sigma |
| 2.2 | 21.01.2026 | delta e^(-2 alpha) | 137.035994102 | -242 sigma |
| 2.4 | 22.01.2026 | + delta^2, ausdrücklich „diagnostic refinement“ | 137.03599819 | -46,8 sigma |
| 2.5 | 27.01.2026 | + (g/4) delta^2, g = 5 | 137.035999216 | +1,9 sigma |
| 2.7 | 03.02.2026 | wie 2.5, mit Begründung der Zählung | 137.035999216 | +1,9 sigma |

Was die Rechnung zeigt:

* Es wurde kein stetiger Parameter angepasst. Der Wert, den ein Fit liefern würde, ist 1,2404 ± 0,0051. Der Wert 5/4 liegt 1,9 sigma daneben. Ein Fit hätte die Abweichung auf null gebracht.
* 5/4 ist eine Strukturzahl: 5/4 = q(D5) = (3/2)(5/6) = g/4 mit g = 5. Die Zahl 5/6 steht als diskreter Abschluss schon in Fassung 2.1, also vor jedem zweiten Ordnungsterm.
* Fassung 2.7 hält schriftlich fest, g sei „prior to any alpha comparison“ durch die SU(5) Holonomieentartung festgelegt.
* Seit Oktober 2026 gibt es eine exakte Modellherleitung des Exponenten 5/4 als getrennte D5 Gradantwort.

Was richtig bleibt und was ich falsch gemacht habe: Der Term (g/4) delta^2 erscheint im Archiv erstmals in Fassung 2.5, fünf Tage nach der Fassung mit Koeffizient 1. Aus dieser zeitlichen Reihenfolge habe ich am 18.11. auf eine Anpassung an den Messwert geschlossen. Das war ein Fehlschluss. Reihenfolge ist kein Beleg für einen Fit, und die Zahlen sprechen dagegen.

Ehrliche Grenze: Unter den Exponenten n/4 liegt nur n = 5 innerhalb von 2 sigma (n = 4: -46,8; n = 6: +50,6). Die Messung kann also benachbarte ganze Zählungen unterscheiden. Ob die Zählung g = 5 vor oder nach dem Blick auf die Zahl aufgeschrieben wurde, lässt sich aus dem Archiv allein nicht beweisen; die Quelle sagt „vorher“, und nichts im Archiv widerspricht dem. Das zitierte Skript defect_partition_g5_audit liegt nicht im Repo. Die Angabe, xi_tree = 3/4 stehe schon in der Fassung vom September 2025, habe ich nicht verifizieren können.

## 2. Reed Muller (RMFAMILY20261119)

* Die Pascalgleichung 2^g = g^2 + g + 2 ist gleichbedeutend mit dim RM(2,g) = 2^(g-1), also mit der Selbstdualität des Codes zweiter Ordnung. Sie gilt nur für g = 5.
* RM(2,5) ist ein doppelt gerader selbstdualer Code [32,16,8]. RM(1,3) ist der Hammingcode [8,4,4], aus dem E8 entsteht.
* Die Monome vom Grad höchstens 2 in fünf Variablen tragen genau die Ladungen 1 + 5 + 10. Das passt zu einer Familie erst nach der Komplementabbildung auf Grad eins (5 gegen 5 quer).
* Die Symmetrie AGL(3,2) hat Ordnung 1344. Der Stabilisator einer Dreiermenge hat Ordnung 24, zerlegt acht Punkte in 3 + 1 + 4 und wirkt als volle S4 auf dem Viererblock. Alle 56 Dreiermengen sind gleichwertig.
* Teilung der Wurzeln mit lauter halben Koordinaten in 5 + 3 gibt die Normen 5/4 und 3/4.
* Fehlstelle: Die Verdopplung 16 = 7 + 4 + 4 + 1 gibt nicht 1 + 10 + 5.

Urteil: eine exakte neue Lesart von g = 5, nicht aus dem Bestand abgeleitet. Sie erklärt die Zahl fünf, nicht die Zahl drei. Die Wahl der Dreiermenge bleibt eine Wahl.

## 3. Kopplungsdualität (COUPLINGMODULUS20261119)

Vorab festgelegte Konvention: e^2 = 4 pi alpha. Das Maxwellgewicht eines Flusses m ist exp(-(pi/(2 alpha)) m·m). Ein Faktor exp(-2/alpha) verlangt m·m = 4/pi, der duale Faktor exp(-2 alpha) verlangt n·n = 1/pi. Kein ganzzahliges Gitter hat diese Normen, und keine ganzzahlige Ladungsskalierung repariert das.

Urteil: als Paar einer Thetafunktion ausgeschlossen. Übrig bleibt nur die formale Vertauschung alpha gegen 1/alpha. Die 2 in exp(-2 alpha) ist in den Quellen eine Zählung (Überlagerungsgrad mal Spektralfluss), keine Flussnorm.

## 4. Gittermittel (DEBYEWALLER20261119)

Vorab festgelegt: Phase exp(i e beta·p), e^2 = 4 pi alpha, beta eine E8 Wurzel, zwei Blätter. Exaktes Objekt: phi(z) = E_{4,1}(i,z)/E_4(i), z = e/(2 pi).

* Zweites Moment: <p_i p_j> = delta_ij/(4 pi) = 2 c3, exakt. Daraus folgt pro Blatt exp(-alpha), mit zwei Blättern exp(-2 alpha). Der aufgezeichnete Dämpfungsfaktor ist genau dieser Gaußanteil.
* Selbstdualität erzwingt log phi = -pi z^2 + a4 z^4 + a8 z^8, ohne z^6. a4 = 7,87803000538, a8 = -15,9591.
* Das volle Gittermittel verschiebt die Wurzel auf alpha^-1 = 137.036002152, also +141,7 sigma gegen CODATA 2022 (aufgezeichnet: +1,9).
* Das zweite Moment 1/(4 pi) hat schon das einfache Gitter Z bei tau = i. Es wählt E8 nicht aus.

Urteil: Die strenge Lesart (Mittel über die Impulsgitterpunkte) ist mit 142 sigma ausgeschlossen. Tragfähig bleibt nur: exp(-2 alpha) ist eine reine Gaußantwort mit Varianz 2 c3 pro Richtung und zwei Blättern. Welcher Zustand exakt gaußisch ist, bleibt offen. Der Test war vorab festgelegt und parameterfrei, das Ergebnis ist negativ.

## 5. Hopf Hebung und Twistor (HOPFTWISTOR20261119)

Formel: w_k = 1/(2^(2k) vol S^(2k-1)) = (k-1)! 2^(k-1) c3^k. Das ist eine reine Umschreibung von c3 = 1/(8 pi).

| k | Sphäre | Hopf Faserung | w_k | aufgezeichnete Größe | Befund |
|---|---|---|---|---|---|
| 1 | S1 | reell | c3 | c3 | Treffer |
| 2 | S3 | komplex | 2 c3^2 | Skala in g = alpha/(2 c3^2) | Treffer |
| 3 | S5 | keine | 8 c3^3 | 2 c3^3 | Faktor 4 daneben |
| 4 | S7 | quaternionisch | 48 c3^4 | delta_top | Treffer |
| 6 | S11 | keine | 3840 c3^6 | (4/5)·41 c3^6 | daneben |
| 8 | S15 | oktonionisch | 645120 c3^8 | (5/4) delta_top^2 = 2880 c3^8 | Faktor 224 daneben |

* Die Twistorfaserung CP1 → CP3 → S4 ist dieselbe Kette S1 ⊂ S3 ⊂ S7. Damit delta_top = c3 · (2 c3) · (24 c3^2): Kreis, Twistorlinie, Basis S4.
* Nicht eindeutig: 1/(2^11 vol B8) und 1/(2^11 vol CP4) geben dieselbe Zahl delta_top.
* Das Instanton auf dem quaternionischen Hopfbündel hat Gewicht exp(-2 pi/alpha), rund 10^-374. Es spielt numerisch keine Rolle.
* Die Kette wählt weder g = 5 noch drei Familien noch 5/4 aus.

Urteil: drei exakte Volumentreffer bei den Hopfsphären der reellen, komplexen und quaternionischen Zahlen, drei Fehlstellen, keine Herleitung. Bezeichnungsebene.

Stand der vorhandenen Twistor und Hopf Verträge im Ledger, unverändert:

* CELEST.SEAM.01: offen. WP1 bis WP5 ausgeführt. v492 heute erneut 36/36.
* WOIT.OS.TWISTOR.01: offen. Freie Reflexionspositivität hält (v519 heute erneut 34/34); die wechselwirkende Algebra ist nicht gebaut.
* SEAM.EQUIV.TWISTOR.01: offen.
* SEAM.BULK4D.RECON.01: offen. Das Ledger hält fest, dass Hopf Hebung, 4D Kegel und Einsteinmetrik (v957, heute erneut 17/17) diesen Vertrag nicht schließen.

Die Hopf Hebung ist also als Kinematik zertifiziert. Was fehlt, ist unverändert die Dynamik: ein wechselwirkender reflexionspositiver Zustand und die Rekonstruktion des 4D Volumens.

## Nicht geprüft

* Messwerte von Alpha sind aus dem Gedächtnis zitiert (CODATA 2022, Rb 2020).
* Atlas: eingetragen am 19.11.2026 (Revision ed949b8a…), siehe contract_index.json.
* Die Archivfassungen wurden nach Stichworten durchsucht, nicht Zeile für Zeile gelesen.

## Wiederholen

    cd experiments/theory-contracts/origin-ideas-alpha-recheck-20261119
    for s in alpha_provenance reed_muller duality debye_waller hopf_twistor; do python3 -B check_$s.py; done
