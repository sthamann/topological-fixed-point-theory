# Rekursive Maschine aus TFPT und Kohonen

Die Maschine läuft als Python Prozess mit einer lokalen Bedienoberfläche:

**[Maschine öffnen](http://127.0.0.1:8769)**

60 Zellen erzeugen gemeinsam eine Antwort. Ein lernender Speicher verarbeitet diese Antwort. Das Ergebnis fließt in die Zellen zurück und wird zur nächsten Eingabe. Es kommen während des Laufs keine externen Trainingsdaten hinzu.

## Bedienen und erneut starten

In der Oberfläche lassen sich einzelne Schritte und längere Läufe ausführen, Störungen einspeisen und Zustände sichern. Startwert, Lernrate, Rückkopplung und die Vorzeichen der Kopplung lassen sich für einen neuen Lauf ändern.

Zum erneuten Starten im Terminal:

```sh
cd /Users/stefanhamann/Documents/Codex/2026-10-26/checl-mal-was-teuvo-kohonen-herausgefundwn
./run.sh
```

Ein gespeichertes System lässt sich einschließlich Speicher und Zufallszustand fortsetzen:

```sh
./run.sh --resume runs/autosave.npz
```

Ohne Oberfläche:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 machine.py --steps 10000 --save runs/my_run.npz
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 machine.py --resume runs/my_run.npz --steps 10000 --save runs/continued.npz
```

Python 3 und NumPy genügen. Die getestete Umgebung verwendet NumPy 2.4.2. `native_operators.npz` enthält die Eingabematrizen, sodass der laufende Rechner das ursprüngliche Repository nicht benötigt. `native_data.py` kann diese Datei aus den angegebenen Originalquellen neu erzeugen.

## Was tatsächlich eingebaut ist

Die ursprünglichen 240 E8 Wurzeln ergeben unter der vorhandenen komplexen Struktur 60 Strahlen. Die normierten Tensoren ihrer vierten Potenz werden in einer orthonormalen Basis von `Sym^4(C^4)` dargestellt. Alle 35 komplexen Komponenten bleiben im Signal und im Speicher erhalten.

Der kollektive Kopplungskern ist die ursprüngliche Gram Matrix

\[
K_{ab}=\langle z_a,z_b\rangle^4.
\]

Seine Eigenwerte sind `3` mit Vielfachheit 5, `1.5` mit Vielfachheit 30 und `0` mit Vielfachheit 25. Die ursprüngliche Basis des Fünfercodes dient der Auswertung. Sie wird weder als Zielantwort noch als Projektor in der Lernregel verwendet. Der fünfdimensionale Spektralvorsprung steckt bereits im gewählten ursprünglichen Kern.

Die Herkunft der Matrizen und Prüfsummen der verwendeten Quelldateien stehen in [native_manifest.json](native_manifest.json). Der maximale numerische Fehler der kontrollierten ursprünglichen Momentidentität beträgt `1.33e-15`.

Der Kohonen Anteil verwendet Wettbewerb um den ähnlichsten Speicher und eine Anpassung seiner Nachbarn. Als Primärquelle dient [Kohonen und Honkela, Kohonen Network, 2007](https://www.scholarpedia.org/article/Kohonen_networks). Der Artikel beschreibt die Lernidee; die konkrete rekursive Verbindung mit den TFPT Matrizen ist die hier untersuchte Versuchsregel.

## Die gewählte Rückkopplung

`X` enthält 60 Zeilen mit jeweils 35 Komponenten. `W` enthält die 60 Speicher. `N` normiert jede von null verschiedene Signalzeile auf Länge 1. Eine Nullantwort bleibt null.

1. Die Zellen erzeugen die kollektive Antwort `Y = N(((1-c) I + c K/3) X)`.
2. Für jede Antwort `Y_a` gewinnt der Speicher mit dem kleinsten euklidischen Abstand. Komplexe Komponenten entsprechen dabei 70 reellen Komponenten.
3. Der Gewinner und seine Nachbarn erhalten gewichtete Beiträge dieser Antworten. Aus diesen Beiträgen wird je Speicher ein Mittelwert `T_i` gebildet.
4. Der Speicher lernt mit `W'_i = (1-eta) W_i + eta T_i`, sofern er Beiträge erhalten hat.
5. Die neue Eingabe lautet `X'_a = N((1-beta) Y_a + beta W'_{winner(a)})`.

Die erste zusätzlich gewählte Dynamik ist die Mischung von Eigenzustand und kollektiver Antwort mit anschließender Normierung. Hinzu kommen die Lernrate, das Nachbarschaftsgewicht und die Rückkopplung des Speichers. Diese Regeln sind hier explizite Modellannahmen. Aus den TFPT Originalprinzipien werden sie in diesem Experiment nicht hergeleitet.

Standardwerte: `c=0.5`, `eta=0.15`, `beta=0.25`. Der Gewinner erhält Gewicht 1, seine Nachbarn mit ursprünglichem Strahlüberlappungsquadrat `1/2` erhalten Gewicht `0.25`. Die Speicher starten als zufällig umgeordnete Kopien der Signale. Zum Start kommt normiertes komplexes Rauschen mit Stärke `0.05` hinzu. Der Vergleich ohne dieses Rauschen wurde ebenfalls ausgeführt.

Die ursprüngliche Matrix `K` bleibt fest. Speicher, Gewinnerzuordnung und Signale verändern sich im Lauf. Dadurch verändert sich die wirksame Antwort der Maschine auf ihren eigenen Zustand.

## Was beim Ausführen passiert ist

Die folgenden Zahlen stammen aus tatsächlich ausgeführten Läufen mit Startwert 7 nach jeweils 10.000 Schritten:

| Versuchsregel | Effektiver Rang | Anteil im Fünfercode | Mittlere Ähnlichkeit |
|---|---:|---:|---:|
| Kopplung, Lernen und Rückkopplung | 3.794 | 98.07 % | 0.263 |
| Rückkopplung aus | 5.000 | 99.93 % | 0.186 |
| Speicher eingefroren | 18.336 | 63.89 % | 0.071 |
| Wechselwirkung aus | 1.000 | 7.35 % | 1.000 |
| Vorzeichen durch Beträge ersetzt | 1.000 | 10.46 % | 1.000 |
| Ohne anfängliches Rauschen | 3.794 | 98.14 % | 0.263 |

Der effektive Rang ist `exp(-sum p log p)`, wobei `p` das normierte Spektrum von `X* X` bezeichnet. Er misst die Verteilung über unabhängige Signalrichtungen und ist keine ganzzahlige Raumdimension. Die mittlere Ähnlichkeit ist der Mittelwert der Betragsquadrate der Signalüberlappungen zwischen verschiedenen Zellen.

Mit 20 Startwerten von 0 bis 19 und jeweils 5.000 Schritten lag der effektive Rang zwischen **3.46 und 5.00**, der Anteil im Fünfercode zwischen **97.95 und 98.45 %**. Bei ausgeschalteter Wechselwirkung sowie im Vergleich mit entfernten Vorzeichen gleichen sich im untersuchten Lauf alle Signalrichtungen an. Mit der ursprünglichen Kopplung bleiben mehrere Richtungen erhalten.

Das Lernen verändert das Ergebnis gegenüber dem Vergleich ohne Rückkopplung. Es liefert bei diesen Regeln keine startunabhängige eindeutige Antwort. Die Unterschiede zwischen den Startwerten gehören zum beobachteten Ergebnis.

## Die Struktur bewegt sich weiter

Der Lauf mit Startwert 7 wurde bis **100.000 Schritte** fortgesetzt. Der effektive Rang beträgt dann `3.794077`, der Anteil im Fünfercode `98.1094 %`. Die Signaländerung je Schritt bleibt bei ungefähr `0.000377119`. Stabile Gesamtwerte bedeuten hier also keinen unbewegten Zustand.

Für den beobachteten Zustand wurde zusätzlich eine lineare Antwortmatrix aus den tatsächlichen Gewinnern und Normierungsfaktoren rekonstruiert. Bei festgehaltenen Faktoren reproduziert sie die nächsten 1.000 Schritte der vollständigen Maschine mit einem relativen Fehler von `5.48e-14`. Sie besitzt zwei komplexe Eigenwertpaare auf dem numerischen Einheitskreis. Ihre Winkel entsprechen ungefähr 11.201 und 64.540 Schritten pro Umlauf. Das sind Schwingungsanteile der gewählten rekursiven Lernregel, gemessen in Iterationen.

Diese Matrix ist eine eingefrorene Antwortdarstellung am untersuchten Lauf, nicht der vollständige Jacobi Operator der nichtlinearen Maschine und kein Beweis für sämtliche Startwerte. Die Zahlen werden nicht als Takte einer physikalischen Zeit eingesetzt.

Nach einer Störung mit Stärke `0.1` fiel die mittlere Signalabweichung gegenüber dem ungestörten Vergleich von `0.09891` unmittelbar nach der Störung auf `0.01182` nach zehn Schritten. Eine kleine bleibende Verschiebung blieb erhalten. Der effektive Rang kehrte dabei in die Nähe von `3.794` zurück.

## Ergebnisse und Prüfung

* [machine.py](machine.py): tatsächlicher Ablauf, Zustände und Wiederaufnahme.
* [server.py](server.py) und [index.html](index.html): lokale interaktive Bedienung.
* [experiments.py](experiments.py): Vergleichsläufe, 20 Starts und Störungstest.
* [inspect_dynamics.py](inspect_dynamics.py): Erklärung und Prüfung der fortlaufenden Bewegung.
* [runs/experiment_results.json](runs/experiment_results.json): vollständige gemessene Vergleichsdaten.
* [runs/long_run.json](runs/long_run.json): Messdaten bis 100.000 Schritte und eingefrorene Antwortmoden.

Erneut ausführen:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 experiments.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 inspect_dynamics.py
```

Die Prüfung umfasst die ursprünglichen Eigenwerte, Signalnormen, einen ursächlichen Kontrolllauf mit ausgeschalteter Rückkopplung sowie Speichern, frisches Laden, Stören und Fortsetzen. Letzteres reproduziert in derselben Laufzeitumgebung die Signale, Speicher und Zufallszustände bitgenau. Auch die Bedienfolge aus 1.000 Schritten, Sichern, Stören und Zurückladen wurde im Browser ausgeführt und am sichtbaren Zustand geprüft.

**Ergebnisstatus: ausgeführte, endlich getestete Modellstudie mit offen ausgewiesener Rückkopplungsregel.**

## Weitere Untersuchung vom 26. Oktober

Die [fundamentale Fortsetzung](FUNDAMENTALE_FORTSETZUNG.md) enthält die inzwischen ausgeführte Prüfung der Codeausrichtung und der Lernraten. Sie präzisiert den hohen Codeanteil: Die Update Regel allein ist unter allen gemeinsamen unitären Drehungen des Signalraums kovariant. Die ursprüngliche Ausrichtung von Eingang, Kopplung und Auslesung liefert daher eine wesentliche Voraussetzung. Die Schwingungen hängen von den gewählten Raten ab. [Vollständige Messdaten](runs/fundamental_audit.json).

Der anschließend ausgeführte [Quellanschluss mit geometrischer Rückwirkung](QUELLANSCHLUSS_UND_RUECKWIRKUNG.md) verwendet den tatsächlichen vorhandenen Diracoperator. Für eine kontrollierte Änderung seiner Geometrie wurden die Kopplungen und die Gedächtnisrekursion aus diesem Operator berechnet. Die ursprüngliche Kohonen Maschine bleibt als Vergleich erhalten.

Die [Untersuchung der Auswahl dieser Rückwirkung](AUSWAHL_DER_RUECKWIRKUNG.md) leitet die gekoppelte Variation des ursprünglichen Referenzpaares her. Ein exakter Gegenfall zeigt, dass gleiche endliche Operatorantworten verschiedene Geometriekräfte zulassen. Der bereits festgelegte skalare Nahttransfer behält seine eindeutige Ableitung. Der nächste Anschluss ist seine ursprüngliche geometrische Einfügung samt vollständiger Quellantwort. Reproduktion: `python3 check_primitive_selection.py`.

Die [geometrische Antwort mit Gedächtnis](GEOMETRISCHE_ANTWORT_MIT_GEDAECHTNIS.md) transportiert den alten Nahtdefekt und den ursprünglichen See auf derselben Geschichte. Die Maschine berechnet jetzt dessen geometrische Antwort. Ein Impulsversuch zeigt verschiedene Feldantworten bei exakt derselben abschließenden Geometrie und demselben relativen Indexfaktor. Unmittelbarer und zeitverzögerter Anteil sind hergeleitet und gegen die direkte Entwicklung geprüft. Reproduktion: `python3 check_geometric_sector_response.py`.

Die [Formgeometrie und ihre Zulässigkeitsprüfung](RANDKOPPLUNG_AUS_FORMGEOMETRIE.md) enthalten die abgeleitete Randantwort, den neutralen Vierfeldanteil und den Ausschluss dieser Modulusvariation im vorhandenen festen Compiler.

Die [berechnete Verbindungsantwort](VERBINDUNGSANTWORT_UND_QUELLENEINFUEGUNG.md) enthält den tatsächlichen elliptischen Determinanten, beide Umlaufphasen und den Vergleich mit der ursprünglichen Besetzungseinfügung.
