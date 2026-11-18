# Warum die fundamentale Lösung fehlt, und fünf neue Ursprungsproben

16. November 2026 · F1 und F2 · Kennung ORIGINDIAG20261116

Geltungsbereich: Befund über den Suchprozess aus `progress.json` Revision 399
und exakte Gitterrechnungen im D8 Chart von E8. Kein Ledgerwechsel, keine
physische Herleitung, keine Beförderung. Atlasstand beim Lesen
`92ba6b8310357b488a3076f8`.

## 1. Befund: Die Suche wartet auf einen Produzenten, den die Eingaben nicht enthalten

1. Von 246 Einträgen `latest_*` in `progress.json` nennen 55 dieselbe offene
   Stelle: die ursprüngliche gewichtete Form, die absolute Leerreferenz oder die
   physische Zeitkarte. Der direkte P1 Aufruf vom 12. November liefert
   `original_full_source_operation_not_supplied`, null Originalatome und 27
   kinematische Diagnosekonfigurationen. Der autonome Lauf vom 16. November
   beendet 355 Jobs und findet null passende Originalpakete.
2. Alle 52 abgegrenzten F1 Abschlüsse liegen auf der Verbraucherseite. Sie
   sagen, wie eine gelieferte Quelle zu lesen wäre. Keiner erzeugt sie.
3. Die Unterbestimmtheit ist im Bestand bereits mehrfach exakt belegt:
   isospektrale Verklebungen mit gleichen Einzelmarken, zwei unimodulare
   Abschlüsse über derselben D8 Stufe, das v252/v88 Paar mit gleichem Spektrum
   und verschiedenem CP Vorzeichen, die lokale Deformation mit allen antipodalen
   Antworten in (0, 1/(8π)]. Daraus folgt: P1, P2 und E8 legen die verlangten
   vier Marken nicht fest. Eine weitere Prüfung der Originalquelle kann eine
   Information nicht finden, die dort nicht steht.
4. Die beiden Hälften der Theorie leben auf verschiedenen Algebren. Die
   diskrete Schließung benutzt die E8 Stromalgebra, die Randantwort den 96
   dimensionalen Spektraltripelraum. Probe T1 unten zeigt, dass eine E8
   kovariante Bindung auf drei Familien nicht existiert.
5. Das F1 Erfolgskriterium ist in endlich vielen Schritten nicht entscheidbar.
   F2, die Frage nach den wirklich unabhängigen Annahmen, hat Priorität P2 und
   keinen Owner, obwohl genau sie beantwortet, was fehlt.

Folgerung: Es fehlt ein Auswahlprinzip, also eine zusätzliche Annahme. Sie
lässt sich nicht aus dem Original herleiten. Die sinnvolle Frage lautet, welche
kleinste zusätzliche Regel die vier Marken festlegt und ob sie die vorhandenen
Ausgaben reproduziert.

## 2. Hypothesen und exakte Proben

Reproduktion: `python3 -B check.py` (Exit 0, alle Prüfungen bestanden, Ausgabe
in `run.log`, Daten in `results.json`).

| Kennung | Hypothese | Probe | Verdikt |
| --- | --- | --- | --- |
| H1 | Die R113 Bindung ist ein Typfehler: E8 geladene Wörter haben auf dem Raum dreier Familien keinen Platz. | T1 | Bestätigt im Bereich E8 kovarianter Bindungen |
| H2 | Die Naht ist die Spiegelwand eines 3+1D reflexionssymmetrischen Bulks. Beide Seiten tragen je c = 4 mit drei Fermionen (D4). | T2, T3 | Gitterseitig möglich, mit der 3+2 Hyperladung unverträglich |
| H3 | Die absolute Leerreferenz der einzelnen Naht ist keine Zahl. | T4 | Bestätigt, bekannte Tatsache mit neuer Folge für R113 |
| H4 | Eine goldene (ikosische) Struktur liefert auf E8 eine kanonische 1+3 Zerlegung, die zum 5+3 Schnitt passt. | T5 | Ausgeschlossen für D5 und A3 einzeln |
| H5 | Fixpunkt statt Quelle: Die vier Marken folgen aus Positivität, Komposition und Clocksymmetrie allein. | Bestand | Nicht neu gerechnet. Die Gegenbeispiele aus Abschnitt 1 Punkt 3 widerlegen Eindeutigkeit in den geprüften Klassen. |

### T1 Typgrenze

Im D8 Chart zerfallen die 240 Wurzeln in 40 (D5), 12 (A3), 60 (10,6) und 128
Materiewurzeln (16,4) + (16bar,4bar). Full96 entspricht drei der vier
Familiengewichte.

- Die Summe zweier Materiewurzeln ist nie eine Materiewurzel. Materieströme
  führen aus dem Materieraum heraus.
- Für jede der vier Wahlen der ausgelassenen Familie erhält keiner der 60
  (10,6) Ströme den Raum dreier Familien. 960 Matrixelemente bleiben innen,
  480 treten in die vierte Familie aus. Die Kompression auf Full96 ist daher
  kein Homomorphismus der Stromalgebra.
- Den Raum erhalten genau die 40 D5 Wurzeln und die 6 Wurzeln der A2 in A3.
- Die drei Familienkoordinaten sind die drei perfekten Paarungen der vier
  Familiengewichte.
- e8 besitzt keine nichttriviale Darstellung unterhalb der Dimension 248.

Folge: Der native R113 Aufruf kann ein E8 geladenes Wort mit (10,6)
Mittelstrom auf Full96 grundsätzlich nicht auswerten. Der kleinste Raum, der
das geordnete Wort samt kompensierenden Enden trägt, enthält alle vier
Familiengewichte. Der Ausschluss gilt für E8 kovariante Bindungen. Zusätzliche,
frei gewählte Operatoren auf Full96 erfasst er nicht.

### T2 und T3 Spiegelwand

Der Stabilisator von D5 + A3 in W(E8) hat Ordnung 92160. Die Probe klassifiziert
alle Involutionstypen nach Minusrang, Fix und Antifixgitter (35 Klassen).

- Minusrang eins haben nur Wurzelspiegelungen (E7 | A1).
- Vier Klassen zerlegen E8 in D4 | D4. In der reinen Vorzeichenform sind vier
  Trägerkoordinaten negiert, fest bleiben eine Trägerkoordinate und die drei
  Familienkoordinaten.
- D4 hat c = 4 mod 8 und drei nichttriviale Klassen mit Spin 1/2. D5 und A3
  haben c = 5 und c = 3 mod 8 mit Spins 5/8 und 3/8.
- Unter den 18 Involutionen mit Minusrang vier, die die Hyperladung
  (−1/3, −1/3, −1/3, 1/2, 1/2) bis aufs Vorzeichen erhalten, hat jede den Typ
  A1⁴ | A1⁴. Keine zerlegt in D4 | D4.

Verdikt: Die symmetrische Drei Fermionen Spiegelwand ist mit dem P2 Schnitt
3+2 unverträglich. Der Ausschluss betrifft Gitterinvolutionen in W(E8), die
D5 + A3 und ±Y erhalten. Eine antilineare Wirkung außerhalb der Weylgruppe und
ein Bulk ohne Reflexionssymmetrie bleiben ungeprüft.

### T4 Leerreferenz

Der Vakuumcharakter E4/η⁸ der einzelnen Naht hat χ(i) = 12 und χ(ρ) = 0. Unter
τ → τ + 1 erhält er die Phase exp(−2πi/3), unter τ → −1/τ ist er invariant.
Das Betragsquadrat ist invariant. Eine absolute Referenz ist nur für das
reflektierte Paar definiert. Am Dreierpunkt verschwindet die unreflektierte
Referenz, eine Division durch sie ist dort unzulässig. Die Tatsache ist im
Bestand bekannt (CELEST.WP5E.ALPHA.01). Neu ist die Folge für R113: absolute
Werte auf reflektierte Verhältnisse beschränken.

### T5 Goldene Struktur

Ein selbstadjungierter ganzzahliger Operator F mit F² = F + 1 existiert nur auf
Gittern geraden Rangs. D5 und A3 haben Rang 5 und 3. Für den Fünfertakt g der
Trägerkoordinaten erfüllt F = 1 + g + g⁻¹ die Gleichung genau auf dem A4
Komplement. Das Fixgitter von g hat Rang 4, Determinante 5 und Wurzelsystem A4.
Es besteht aus der Familie A3 und der Trägerspur.

## 3. Was für eine Lösung zählen würde

Eine Regel zählt, wenn sie (a) auf dem Raum mit allen vier Familiengewichten
formuliert ist, (b) nur reflektierte Verhältnisse benutzt, (c) die vier Marken
eindeutig festlegt und (d) eine vorhandene Ausgabe vorhersagt, die nicht in
ihre Formulierung eingeht. Keine der fünf Hypothesen erfüllt (c) und (d).
Die fundamentale Lösung ist mit diesem Lauf nicht gefunden.

Nächster entscheidender Schritt: das geordnete E8 Wort auf dem 128
dimensionalen Materieraum mit allen vier Familiengewichten und reflektierter
Referenz auswerten und prüfen, ob die Projektion auf drei Familien als
bedingte Erwartung (statt als Kompression) die vier Marken festlegt.
