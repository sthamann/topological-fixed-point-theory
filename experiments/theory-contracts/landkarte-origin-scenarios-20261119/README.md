# Zustand, Uhren, Neustart: drei Ursprungsproben und Ansätze aus der Landkarte

Datum: 19.11.2026. Ergebniskennung: LANDKARTE20261119. Vier Skripte, 39/39 Checks.
Kein Ledgerwechsel, keine Papieränderung, kein Commit. physical_solution=false.

## Frage

Die Landkarte (tfpt_landkarte) zeigt in ihrer Dynamiktabelle: Kein Ansatz liefert zugleich den Zustand (Z), kommt ohne Zusatzwahl aus (A) und ist an die Originalquelle gebunden (O). Geprüft werden drei Ansätze, die genau dort ansetzen, und eine nie probierte empirische Stelle.

| Teil | Skript | Checks | Urteil |
|---|---|---|---|
| Zustand | check_vacuum_state.py | 11/11 | exakte Lesart gefunden, bedingt |
| Uhren | check_clock_rank.py | 13/13 | exakter Struktursatz |
| Neustart an der Singularität | check_singularity_reset.py | 9/9 | Ausschlüsse eingegrenzt, physische Probe offen |
| Leptonen | check_lepton_transfer.py | 6/6 | QED Lauf ausgeschlossen, Koide halbiert die Lücke |

## 1. Zustand (VACUUMSTATE20261119)

Ausgangslage aus DEBYEWALLER20261119: exp(-2 alpha) ist exakt gaußisch, ein einzelnes E8 Gitter liefert das nicht.

* **Viele unabhängige Keime.** Mit der gemeinsamen Koordinate P = Summe p_k / Wurzel N fällt die Abweichung wie 1/N: N = 1 gibt +141,6 sigma, N = 240 gibt +2,48, N = 1920 gibt +1,97. Innerhalb 2 sigma liegt der Wert ab N = 1360, innerhalb 3 sigma ab N = 127. Das ist eine Familie mit einem Parameter. Die Daten geben nur eine Untergrenze.
* **Grundzustand der Vierteldrehung.** Mit dem Fourierkern exp(-2 pi i x p) (also h = 1) ist psi0 = 2^(1/4) exp(-pi x^2) seine eigene Fouriertransformierte. Es gilt <x^2> = <p^2> = 1/(4 pi) = 2 c3. Die Nullpunktsenergie je Quadratur ist c3. Der Erwartungswert des Verschiebeoperators exp(i e beta·(x+p)) mit Wurzelnorm 2 ist exakt exp(-2 alpha). Alle höheren Kumulanten sind null, die Wurzel bleibt 137,0359992168.
* **Gegenprobe.** Der nächste unter der Fouriertransformation feste Zustand (n = 4) gibt exp(-2 alpha) L_4(4 alpha) und liegt um mehr als 10^5 sigma daneben. Nur der Grundzustand passt.
* **Passt zur Landkarte.** Gaußische Momente enthalten nur pi. Die Punktwerte am Torus enthalten die lemniskatische Konstante, die Ausgaben nicht.

Urteil: Der aufgezeichnete Dämpfungsfaktor ist exakt die Antwort des tiefsten Zustands, den die Vierteldrehung festhält. Damit wären τ = i, die Viertel-Uhr, c3 und die Dämpfung eine einzige Aussage.

Grenzen: Die Konvention h = 1 und das Mitzählen beider Quadraturen habe ich in Kenntnis des Zielwerts gewählt. Dieselbe Zahl kam vorher als „Gittervarianz mal zwei Blätter“ heraus. Die Fouriertransformation hält viele Zustände fest; dass der tiefste gilt, ist eine Zusatzannahme (kleinste Unschärfe). Warum die Originalquelle in diesem Zustand sitzt, ist nicht hergeleitet. Das zweite Moment wählt E8 nicht aus.

## 2. Uhren (CLOCKRANK20261119)

* Eine Uhr der Ordnung n braucht ein Gitter mit Mindestrang. In der Ebene gehen nur 1, 2, 3, 4, 6. Der Fünfertakt braucht Rang 4.
* Die Vier erzwingt das Quadrat, die Drei erzwänge das Sechseck. Beide auf einer Fläche bräuchten Ordnung 12 und Rang 4. Die Picardkurve y^3 = x^4 - 1 trägt genau diese Ordnung 12.
* Primitive Uhren vom Rang genau 8 sind die Ordnungen 15, 16, 20, 24, 30. In W(E8) kommen 15, 20, 24, 30 jeweils als eine primitive Uhr vor: Je zwei der Takte 3, 4 (bzw. 8), 5 verschmelzen in E8.
* Ein Element der Ordnung 60 gibt es in W(E8) nicht (Stichprobe von 60000 Wörtern findet genau die bekannte Ordnungsliste; Beweis: Aut(A4) = S5 x Z2 hat kein Element der Ordnung 12, und der andere Fall scheitert am Grad). Als eine Uhr bräuchten 3, 4, 5 zusammen Rang 16.
* Deshalb wirkt die Vier auf die Fünf als Getriebe: Sie ist die Galoisgruppe (Z/5)^x = Z4, und G C5 G^-1 = C5^2 ist die Frobeniusgruppe F20 in S5 = W(A4).
* Zahlenfolgen: Am Sechseck wäre die Gitterenergie (2/Wurzel 3) c3 und alpha^-1 = 104,42. Kohärente Windung 2 gibt 37,30, Index 8 gibt 512,23.

Urteil: exakter Struktursatz. Er erklärt, warum das Getriebe nicht kommutiert und warum die Fläche quadratisch sein muss, sobald die Vier auf ihr wirkt. Er leitet nicht her, warum die Vier auf der Fläche wirkt und die Drei nicht. Die Liste der Ordnungen in W(E8) ist per Stichprobe bestätigt, nicht durch Klassenaufzählung.

## 3. Neustart an der Singularität (SINGRESET20261119)

* Die aufgelöste Singularität C^2/Z4 trägt die Schnittform -A3 (Rang 3), Eulerzahl 4, Rand S^3/Z4 mit Verschlingungsform 3/4 und Signatur-Eta 1/2. Das sind die Zahlen der Familienseite und der vier Marken.
* Ein exakter Neustart V mit V*AV = omega(A) 1 für alle sichtbaren A, der als Isometrie alle alte Information behält, existiert genau dann, wenn die verborgene Seite mindestens r·d Dimensionen hat (r Rang der Referenz, d Größe des neuen Codes). Mit r·d gelingt er, mit r·d - 1 scheitert er. Die im Atlas genannte Schranke ist also scharf.
* Die vier im Atlas ausgeschlossenen Fälle (thermische Doppelung bei endlichem Regulator, selbes positives Vakuum, endlicher Kommutant) sind genau die Fälle ohne Platz. Eine echte kleinere sichtbare Algebra ist dort ausdrücklich nicht ausgeschlossen.

Urteil: P1 als Neustart ist nur möglich, wenn P1 ein lokaler Schnitt mit unbegrenzter Rückseite ist (Horizont), nicht das ganze Vakuum. Das ist die Karte G3 „Horizont zuerst“, die noch keine eigene Probe hat. Dass lokale Algebren eines konformen Netzes solche Isometrien für jedes d besitzen, ist ein Standardsatz und hier zitiert, nicht nachgerechnet. Die entscheidende physische Probe im gewichteten Originalschnitt bleibt unberechnet, weil das Originalpaket fehlt.

## 4. Leptonen (LEPTONTRANSFER20261119)

* Die Lücken der Landkarte sind reproduziert: +2,20 % (mu/tau) und +0,21 % (e/mu).
* Einschleifen QED Lauf auf eine gemeinsame Skala gibt -0,98 % und -1,86 %. Falsches Vorzeichen: ausgeschlossen als Transfer.
* Die Quellverhältnisse erfüllen Koide nicht (Q = 0,66446), die Polmassen schon.
* Nimmt man Koide als Schnittstelle, gibt das Quellverhältnis e/mu allein m_mu/m_tau auf -0,02 %. Die zwei Lücken werden zu einer. Die Vorhersage m_tau = 1777,30 MeV liegt etwa 4 sigma zu hoch, getrieben von der Lücke 0,21 %.
* Nur als Zahl vermerkt: Die Lücke mu/tau ist gleich 3 alpha innerhalb von 0,7 bis 1,7 sigma von m_tau. Die Lücke e/mu (0,29 alpha) hat keine solche Form. Das ist ein Zufallskandidat, kein Ergebnis.

## Was das für die Lösung heißt

Aus der Landkarte und diesen Proben ergibt sich eine zusammenhängende Arbeitshypothese, nicht bewiesen:

1. Die Fläche ist quadratisch, weil die Vier auf ihr wirkt (Uhrensatz).
2. Die Vierteldrehung ist die Fouriertransformation. Ihr Grundzustand ist der Zustand der Quelle. Er gibt c3 als Nullpunktsenergie und exp(-2 alpha) exakt.
3. P1 ist ein lokaler Schnitt mit unbegrenzter Rückseite. Nur dann ist ein Neustart möglich.

Der nächste entscheidende Test: den Grundzustand an die Originalquelle binden, also prüfen, ob der originale P1 Kern (Randkern mit Windung eins) als Zweipunktfunktion genau die Kovarianz 1/(4 pi) je Quadratur hat. Die Landkarte nennt dafür die Karte „Native Operatorgleichheit C_nat“ mit Priorität A.

## Nicht geprüft

* Messwerte (CODATA 2022, Leptonmassen) aus dem Gedächtnis.
* Ordnungsliste von W(E8) per Stichprobe.
* Typ III Aussage zitiert.
* Bindung an die Originalquelle fehlt bei allen vier Teilen.

## Wiederholen

    cd experiments/theory-contracts/landkarte-origin-scenarios-20261119
    for s in vacuum_state clock_rank singularity_reset lepton_transfer; do python3 -B check_$s.py; done
