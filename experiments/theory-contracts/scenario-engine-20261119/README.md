# Szenariosuche S42 (SZENARIO20261119)

Stand 19. November 2026. Lesart und Suchwerkzeug, kein Beweis, keine Beförderung, kein Ledgerwechsel.

**Abgelöst am selben Tag:** Die Suche ist jetzt Teil von OriginLab (`research/originlab/structure_scenarios.py`,
Methode S42, Ergebnis SZENARIO20261119B) mit Fingerabdruck aus dem Atlas. Dieser Ordner bleibt als Beleg des
ersten Laufs mit handgewähltem Fingerabdruck stehen und wird nicht weitergeführt.

## Was anders ist

Der Reihenprüfer S41 bewertet Hypothesen, die jemand vorher hingeschrieben hat. Diese Suche schreibt keine
Hypothesen hin. Sie erzeugt Objekte der Mathematik mit all ihren berechneten Merkmalen und fragt, welche Objekte
und welche Verbindungen von Objekten den Fingerabdruck der Theorie von selbst tragen. Sie setzt keinen
bestimmten Engpass voraus und benutzt weder R113 noch eine andere einzelne Rechnung als Tor.

## Ablauf

1. `fingerprint.json`: 36 Merkmale der Theorie (Zahlen, Brüche, Typen) aus den Karten der Landkarte. Frei erweiterbar.
2. `generators.py`, `gen_mehr.py`: Erzeuger. Jeder liefert Objekte mit Merkmalen `(art, wert, wofür)` und kennt den
   Fingerabdruck nicht. Zurzeit 554 Objekte aus 17 Familien: Teilgitter von E8, E7, E6; Spiegelungsgruppen;
   Drehgruppen; Linsenräume; Modulkurven; elliptische Flächen; Codes und Quantencodes; Qubitsysteme; Anyonen,
   Klebegruppen und Kitaev Stufen; Kegelflächen; Polytope; Darstellungen; Modulformen; Randtheorien.
3. Gewicht je Merkmal: log2(N / (1 + Zahl der Objekte, die es tragen)). Eine 4 zählt wenig, eine 46080 viel.
4. Zufallsvergleich: 1000 falsche Fingerabdrücke aus Werten, die im Vorrat wirklich vorkommen und ähnlich groß
   sind. `p gesamt` ist der Anteil der falschen Fingerabdrücke, bei denen das beste Objekt (Paar, Tripel) des
   ganzen Vorrats mindestens so gut abschneidet. Das Herumsuchen ist damit eingerechnet.
5. Szenarien: Paare und Tripel aus verschiedenen Familien, die ein seltenes Merkmal teilen (Brücke).
6. Jedes Szenario bekommt von selbst Annahme, Prüfung, Fallbedingung, die Liste der fehlenden Merkmale und
   einen Vermerk, ob das Objekt in Landkarte, RESULTS oder SEARCH-METHODS schon vorkommt.
7. Zwei Läufe: über alles, und über den Rest, den der Anker D5 + A3 in E8 nicht trägt.

Aufruf: `python3 -B suche.py --fakes 1000 --seed 1` (drei Sekunden je 200 Vergleiche). Ausgabe `ERGEBNIS.md`,
`results.json`. Prüfung: `python3 -B check_engine.py` (12 von 12).

Erweitern: eine Funktion `gen_*()` in einer Datei `gen_*.py` ablegen, sie wird von selbst geladen. Oder Merkmale in
`fingerprint.json` ergänzen. Jede Erweiterung ändert Gewichte und Zufallsvergleich mit.

## Ergebnis

| Lauf | Bestes | Anteil am Gewicht | Zufall im Mittel | p gesamt |
|---|---|---|---|---|
| alles, einzeln | D5 + A3 in E8 | 32 % (57,5 von 178,1) | 23,7 | 0,001 |
| alles, Paar | dazu Darstellungen von SO(10) | 39 % | 32,7 | 0,001 |
| alles, Tripel | dazu A4 + A4 in E8 und G20 | 44 % | 42,3 | 0,001 |
| Rest, einzeln | E8 selbst | 16 % (19,0 von 120,6) | 16,2 | 0,20 |
| Rest, Paar | Gosset Polytop und W(E8) | 20 % | 22,2 | 0,32 |
| Rest, Tripel | E8, Gosset Polytop, W(E8) | 26 % | 28,9 | 0,27 |

Lesart in drei Sätzen:

1. **Der Anker ist kein Zufall.** Unter 554 Objekten steht D5 + A3 in E8 allein vorn, weit über dem Zufall. Das ist
   teils eingebaut (der Fingerabdruck stammt zum Teil von dort) und deshalb kein neuer Beleg, aber eine Kontrolle,
   dass die Maschine misst, was sie soll. Zweiter ist A4 + A3 in E8 (p 0,002).
2. **Der Rest hat in diesem Vorrat keinen gemeinsamen Träger.** Was der Anker nicht trägt (41, 48, 96, 1920, 46080,
   4/3, 2/3, 1/3, 4/5, 5/3, 8/7, 12/7, 7/2 und weitere), wird von keinem Objekt, Paar oder Tripel besser getragen
   als ein gleich großer Zufallssatz. Kein Szenario der Liste ist ein Fund.
3. **Neu aufgetauchte Spuren, alle schwach.** Parafermionen Z5 haben zentrale Ladung 8/7 und hängen über Z5 an
   A4 + A4 (Klebenorm 4/5); A6 hat Klebenorm 12/7 und Dimension 48. Beides kam ohne Vorgabe heraus und kommt in
   Landkarte und Plan bisher nicht vor. Dagegen steht: der Atlas führt 8/7 als Verhältnis zweier Zahlen 4 und 7/2,
   nicht als eine Ladung. Solange keine Abbildung das verbindet, ist es dieselbe Zahl in zwei Bedeutungen
   (p 0,33).

## Grenzen

- Gleiche Zahl ist nicht gleiche Sache. Die Suche zeigt Kandidaten, die Prüfung bleibt Handarbeit oder Vertrag.
- Der Vorrat entscheidet, was gefunden werden kann. Reelle Größen (c₃ = 1/(8π), α) und Verhältnisse zweier
  Merkmale eines Objekts fehlen noch; der Atlasbefund zu 8/7 zeigt, dass Verhältnisse der nächste Erzeuger sein sollten.
- Der Fingerabdruck ist eine Auswahl aus der Landkarte, keine Liste des Ledgers.
- Die Brücken sind Zahlengleichheiten zwischen Familien; viele sind zufällig (40 Wurzeln von D5 gegen 40 Spiegelungen von G20).

## Nächste Ausbaustufen

1. Erzeuger für Verhältnisse und Differenzen innerhalb eines Objekts.
2. Erzeuger für reelle Invarianten (Volumina, Regulatoren, Werte von L Funktionen) mit Toleranz.
3. Fingerabdruck aus dem Atlas statt aus der Landkarte ziehen.
4. Größerer Vorrat: Niemeier Gitter, K3 Fasern, Hecke Gruppen, weitere Codes.
