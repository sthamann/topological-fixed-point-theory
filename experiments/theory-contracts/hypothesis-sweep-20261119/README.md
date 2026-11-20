# SWEEP20261119 · Reihenprüfer: viele Ursprungshypothesen gleichzeitig

Datum: 19. November 2026. Auftrag: nicht eine Hypothese nach der anderen, sondern viele gleichzeitig
und automatisch bewerten. Kein Ledgerwechsel, keine Beförderung, kein Commit.

## Was das Werkzeug tut

`sweep.py` erzeugt Hypothesen familienweise, sagt für jede die festen Zielgrößen der Theorie voraus
und bewertet automatisch: exakter Treffer je Zielgröße und, wo die Alpha Gleichung berührt ist,
der Abstand in σ. Ausgabe: `results_sweep.json` und die lesbare Rangliste `RANKING.md`.

Zielgrößen: c₃ = 1/(8π), Dämpfung exp(−2α), Exponent 5/4, φ_tree/c₃ = 4/3, δ_top/c₃⁴ = 48,
Flächenuhr der Ordnung 4, zweite Uhr der Ordnung 5, Doppelverhältnis der vier Marken 1/2,
Koeffizienten 2 und 164/5 der Kubik.

Eigene Hypothesen ohne Code: Datei `hypotheses_extra.json` anlegen, zum Beispiel
`[{"id": "X:mein_Test", "text": "Dämpfung 3", "wahl": 1, "pred": {"k": 3}}]`, dann `python3 -B sweep.py`.
Neue Familien: eine Funktion `fam_*()` ergänzen. Laufzeit etwa 2 Sekunden.

## Ergebnis des ersten Laufs: 141 Hypothesen, 7 Familien (`check_sweep.py`, 11/11)

| Familie | geprüft | trägt | fällt | wählt nichts |
|---|---|---|---|---|
| A Schnittgeometrie | 26 | 4 | 19 | 3 |
| B Zustand (Dämpfung) | 16 | 5 | 11 | 0 |
| C Exponent (Gitterklasse) | 25 | 1 | 24 | 0 |
| D Singularität (ADE Schnitte) | 46 | 1 | 45 | 0 |
| E Auswahl der vier Marken | 11 | 4 | 5 | 2 |
| F Zahlen φ_tree, δ_top | 95000 Kombinationen | 1 | – | – |
| G Koeffizienten der Kubik | 16 | 2 | 14 | 0 |

Die wichtigsten Aussagen:

- **Zahlen sind im Suchraum einzig.** Unter 95000 Kombinationen φ_tree = (a/b)·c₃, δ_top = N·c₃^m
  (a, b aus 1 bis 8, N bis 1000, m von 2 bis 6) liegt nur (4/3, 48, 4) innerhalb 3σ.
- **Exponent.** Von 25 Gitterklassen (D₂ bis D₈, A₁ bis A₇, E₆, E₇) trägt nur die Spinorklasse von D₅.
- **Kubik.** Nur Koeffizient 2 und 8·b₁ mit b₁ = 41/10 (Standardmodell) tragen. Ohne Higgs, zwei Higgs
  und MSSM fallen.
- **Singularität.** Von 46 ADE Schnitten gibt nur Z₅ = A₄ die Fünferuhr.
- **Marken (die gestellte Frage).** Die Spiegelung des Kragens allein wählt nichts: sie erlaubt
  rechteckige und rhombische Tori. Liegen die Marken als vier Ecken auf der Naht (Doppel eines flachen
  Rechtecks), ist der Torus rechteckig, aber das Seitenverhältnis bleibt frei. Eine transitive Wirkung
  durch die Kleinsche Vierergruppe gibt es bei jeder Markenlage, sie wählt ebenfalls nichts.
  Erst ein Element der Ordnung 4 erzwingt das Quadrat. Die Klebegruppe von D₅ ist zyklisch (Z₄) und
  leistet das; bei D₄ oder D₆ wäre sie Z₂ × Z₂ und leistete es nicht.
- **Entartung.** c₃ allein trennt „Quadrat mit chiraler Nullmode“ nicht von „Rechteck 2:1 mit voller
  Nullmode“. Erst die Uhr oder die Markenlage trennt.

## Grenzen

- „trägt“ heißt: kein Widerspruch zu den festen Zielgrößen. Es ist kein Beweis und keine Herleitung.
- Die Zielgrößen sind die Zahlen der Theorie selbst. Der Lauf misst Einzigkeit im erklärten Suchraum,
  nicht Wahrheit. Familie F hat genau den erklärten Suchraum und keinen größeren.
- Nicht gezeigt: dass die Klebegruppe im Original tatsächlich auf den vier Marken wirkt.
- Viele Einträge wiederholen Rechnungen aus P1KERN20261119 und LANDKARTE20261119; neu sind die
  Familien E, F, G und die gemeinsame automatische Bewertung.
- Messwert CODATA 2022 aus dem Gedächtnis.

## Wiederholen

```
cd experiments/theory-contracts/hypothesis-sweep-20261119
python3 -B check_sweep.py
```
