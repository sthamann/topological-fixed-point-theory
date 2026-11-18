# Bedingte Quellenkontrolle: geladene Viertelstruktur

Die vorhandene einchiral geladene Teilchen und Loch Quelle hat bei der
ausdrücklich gesetzten Länge ℓ = 4 log 4 und α = 1/4 auf ihrem relativen
Fockträger die Antwortwerte 1/4 und 1/64. Die eingegebene Beziehung lautet
ω(γ₀Tγ₀) − 16 ω(γ₂Tγ₂) = 0. Sie stammt aus dem gegebenen Kontrollfall,
nicht aus einer neu behaupteten Auswahl durch P1/P2.

[source-control.json](source-control.json) enthält nur diese Beziehung und
zwei ursprüngliche Vakuummomente. Der Lauf mit m=2 und Grad 2 erreicht SCS
`optimal`. Sein unbestimmter numerischer Moment lässt sich mit der vorhandenen
Nennergrenze 64 nicht rekonstruieren. Deshalb sind seine beiden
Vervollständigungen unentschieden.

Eine getrennte verstärkte Kontrolle speist zehn ursprüngliche Quellmomente bis
Grad 4 und die Normierung ein. Der endliche Vervollständiger baut einen exakt
geprüften vierdimensionalen Kern. Diese Kontrolle verwendet bereits die beiden
einzelnen Quellwerte und erbringt keine unabhängige Auswahl der Quelle.
Nur die zehn eingespeisten Quellmomente sind gemeinsam; eine Übereinstimmung
aller Wörter bis Grad 4 wird nicht behauptet.

Für das Wort γ₀γ₂Tγ₂γ₀ ergibt die ursprüngliche Quellmatrix exakt 1/256,
der gewonnene exakte Vervollständigungskern dagegen 0. Beide Werte wurden
separat aus den jeweiligen Original bzw. Zertifikatsmatrizen neu berechnet und
mit ihren gespeicherten Angaben verglichen. Somit bestimmen diese Eingaben den
gesamten ursprünglichen Kern nicht.

Die Suchclock ist eine interne CAR Rotation. Die geladene ursprüngliche
Viertelwirkung hat entgegengesetzte Teilchen und Loch Phasen. Die erfolgreiche
Kontrolle identifiziert diese beiden Clockoperatoren nicht und liefert keine
Herleitung von Raumzeit, dem vollständigen Higgs oder E8 Träger, Z[0] oder
einem physikalisch ausgewählten Quellzustand.

Aus dem Projektverzeichnis mit installiertem `requirements.txt`:

```sh
python charged-clock-relation-control-20261102/reproduce.py --out /tmp/tfpt-charged-control-new
```

Der Ausgabeordner muss neu sein. Ergebnis und benutzte Daten des hier
ausgeführten Laufs stehen unter [run-20261102-verified](run-20261102-verified/).

## Vollständige Erzeugerrelationen der festen freien Kontrollquelle

Das mitgelieferte Programm `generator_relations.py` setzt den vollständigen
gegebenen Einteilchentransfer `diag(1/4, 1/64)` voraus und übersetzt beide
Erzeugerrelationen `T c_i† = a_i c_i† T` in 32 Majorana Wortrelationen für
alle vier Bra und Ket Basiszustände. In der bestehenden Grad 3 Suche werden
nur diese Relationen und die beiden Vakuummomente eingespeist; die Paarantwort
`1/256` ist dort keine zusätzlich eingegebene Zielzahl. Der bisherige exakte
`x=0` Kandidat wird aus seinem gespeicherten Zertifikat geladen.

Der ursprüngliche freie Transfer erfüllt alle Relationen. Der Gegenkandidat
verletzt zwei. Die exakte symbolische Lösung der vollen Matrizenrelationen
auf der gegebenen vierdimensionalen Fockdarstellung ist eindeutig
`diag(1, 1/4, 1/64, 1/256)`. Die numerische Suche gibt jeweils optimal einen
Minimal und Maximalwert nahe `1/256`; numerische Optimierung ist kein exaktes
Zertifikat. Da die Wortgleichungen für ganzzahlige Koeffizienten mit dem
doppelten Erzeugeroperator arbeiten, sind die gespeicherten Residuen doppelt
so groß wie die Residuen der unskalierten Operatorgleichungen.

```sh
uv run --quiet --with cvxpy --with sympy python -B -OO charged-clock-relation-control-20261102/generator_relations.py
```

Die freie Operatorrelation ist hier eine ursprüngliche Kontrolle **vor** einer
gemeinsamen Integration dynamischer Felder. Die Aussage begründet keine
Gaußförmigkeit der integrierten Quelle und keine physische Auswahl von α,
Länge, Zustand, absolutem Z[0] oder E8 Einbettung. Die Suchebene Grad 3 ist eine
Momentenrelaxation; die Eindeutigkeit wird in der genannten Fockdarstellung
exakt bewiesen, nicht für alle möglichen abstrakten Quellen behauptet.
