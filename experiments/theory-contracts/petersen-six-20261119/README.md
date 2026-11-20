# SECHS20261119: Woher kommt die Sechs im Transfer?

Datum: 19. November 2026. Anlass: In der Restliste der Szenariosuche (SZENARIO20261119C) blieben 1/729 und 64/729
unerklärt. Das sind (1/3)^6 und (2/3)^6. Frage: Legt das gefundene Objekt, die Teilung (3+2)+3, die Sechs fest?

Aufruf: `../../../.venv/bin/python -B check_six.py` (20 von 20 Prüfungen, exakt, `results_six.json`).

## Was schon im Original steht

`verification/v808_petersen_sixthroot.py` enthält den Kern bereits: Die zehn Paare aus fünf Plätzen (drei Farbe, zwei
schwach) bilden den Petersengraphen. Um das Paar der zwei schwachen Plätze liegen Schalen mit 1, 3 und 6 Paaren. Die
äußere Schale ist ein Sechserkreis. Die Schrittmatrix P auf den Schalen hat die Raten 1, -2/3, 1/3, und P^6 hat genau
1, 64/729, 1/729. Die Ladung der Paare ist der Eigenvektor zu -2/3. Den Exponenten setzt das Original aber fest:
`v124` schreibt p_2 = 6 = Zahl der positiven Wurzeln von A3 = "Größe des Sechsecks", `v487` setzt P2 = 6. v808 nennt
das Ergebnis selbst "nur das Spektrum" und begründet die Sechs nicht. Das ist also nichts Neues von mir.

## Was die Nachrechnung dazu bringt

1. **Die gleichmäßige Leiter 1, 2/3, 1/3 gibt es nur bei fünf Plätzen.** Für Paare aus n Plätzen sind die Raten
   1, (n-3)/C(n-2,2), 1/C(n-2,2). Gleiche Abstände verlangen n^2 - 9n + 20 = 0, also n = 4 oder 5; n = 4 zerfällt.
   Die 3 im Nenner ist dort die Zahl der Nachbarn eines Paars. Das Original schreibt (1 - n/N_fam)^6 mit der
   Familienzahl 3. Beide Dreien haben denselben Wert; dass sie dasselbe sind, ist nicht gezeigt.
2. **Die Sechs ist die Ordnung der Vertauschung "Farbe im Kreis, schwach getauscht".** (f1 f2 f3)(w1 w2) hat Ordnung
   kgV(3, 2) = 6, läuft auf den Paaren in Bahnen 1 + 3 + 6, und ihre Sechserbahn ist genau der Sechserkreis. In S5
   haben nur die 20 Elemente dieses Typs die Ordnung 6. Das Original liest dieselben sechsten Einheitswurzeln als
   Familie mal Blatt (`v118`); hier wären es Farbe mal schwach. Zwei Lesarten, nicht entschieden.
3. **Positivität verlangt nur eine gerade Schrittzahl.** P ist symmetrisch bezüglich der Gewichte 1 : 3 : 6; ungerade
   Potenzen haben einen negativen Eigenwert. Aber 2 und 4 Schritte wären auch positiv.
4. **Sechs ist der kürzeste gerade Rundweg ohne Umkehr, und nur bei fünf Plätzen.** Der Petersengraph hat keine
   Vierecke; kürzester Kreis 5, kürzester gerader Kreis 6. Bei sechs und sieben Plätzen ist der kürzeste gerade Kreis 4.
5. **Fünf Lesarten der Sechs fallen nur bei 3 + 2 zusammen** (Ringlänge, Ordnung der Drehung, Zahl gemischter Paare,
   zweimal Nachbarzahl, kürzester gerader Kreis). Bei 4+2, 5+2, 4+1, 3+3, 2+2 gehen sie auseinander.

## Urteil

Kandidat für eine Regel: "Ein Transferschritt ist der kürzeste Rundweg, der positiv ist." Zusammen mit Punkt 3 und 4
gäbe das die Sechs ohne Setzung, und nur für fünf Plätze. **Diese Regel steht nicht im Original und ist nicht
hergeleitet.** Offen bleibt die konkrete Abbildung: Der Transfer des Originals wirkt auf drei Spitzen einer Modulkurve,
P wirkt auf drei Schalen von Paaren. v808 zeigt gleiches Spektrum, aber keine gleiche Basis. Solange niemand die drei
Spitzen mit den drei Schalen identifiziert, ist die Sechs erklärt als Zahl des Objekts, nicht als Folge des Originals.

Grenzen: endliche Graphenrechnung; die Lesarten wurden nach Ansicht des Ergebnisses zusammengestellt; fünf
Vergleichsteilungen, keine vollständige Liste. Nichts befördert, kein Ledgerwechsel, kein Atlaseintrag, kein Commit.

---

# SECHSB20261119: Sind die drei Zustände des Originals die drei Schalen?

Aufruf: `../../../.venv/bin/python -B check_cusps.py` (20 von 20 Prüfungen, exakt, `results_cusps.json`).

Das Original (`v486`, `v487`) hat einen festen Zustand und ein tauschbares Paar: bleiben 1/2, wechseln 1/6, Verlust
1/3. Gerade Paarform Rate 2/3, ungerade 1/3. Dagegen gestellt: der Schritt auf den Schalen 1, 3, 6.

1. **Schalen als Zustände: ausgeschlossen.** Kein Eigenvektor des Schalenschritts hat einen Nulleintrag, also ist in
   keiner Funktion des Schrittes eine Schale fest oder abgekoppelt. Keine Vertauschung von Schalen lässt den Schritt
   unverändert, also gibt es kein tauschbares Schalenpaar. Der positive Schritt mit Raten 1, 2/3, 1/3 hat auf den
   Schalen negative Einträge, ist also keine Sprungregel. Gilt für jede Reihenfolge der Schalen.
2. **Was passt: Zustände des Originals als Schwingungsformen.** Die drei Formen sind 1, die Ladung Y, und ein Ausdruck
   zweiten Grades in Y (25 Z = 108 Y^2 - 6 Y - 27). Form n hat die Rate 1 - n/3 im Betrag. Die Stufe n der Leiter des
   Originals wäre damit der Grad in der Ladung. Bei nur zwei möglichen Zuordnungen ist das schwache Stütze.
3. **Was nicht passt: die Vorzeichen.** Auf den Schalen wechselt Form n je Schritt das Vorzeichen (-1)^n. Im Original
   sind Stufe 0 und 1 gerade, Stufe 2 ungerade. Nur gerade Schrittzahlen stimmen überein.
4. **Der Blatttausch hat auf dem ganzen Graphen kein Gegenstück.** Keine der 25 Spiegelungen in S5 ist gerade auf
   allen 2/3 Formen und ungerade auf allen 1/3 Formen.
5. **Auf dem Sechserkreis allein passt die Parität.** Der Schritt dort verliert je Zug genau 1/3 (wie im Original).
   Der Tausch der zwei schwachen Plätze ist der Gegenpunkt auf dem Kreis: gerade auf 2/3, ungerade auf 1/3, wie `v487`.
   Aber: es gibt dort auch -1/3 und -2/3, die Vielfachheiten sind 2 und 4 statt 1 und 1, einen festen Zustand gibt es
   nicht, und zusammengefasst nach schwachem Platz kommt "bleiben 0, wechseln 2/3" heraus, nicht 1/2 und 1/6. Die
   Rate 1/3 entsteht dort als (2/3) cos 60°, auf den Schalen als 1 - 2/3.

## Urteil

Die Gleichsetzung Spitzen = Schalen ist ausgeschlossen. Es bleiben zwei Teillesarten, die sich gegenseitig
ausschließen: Formen auf den Schalen (Raten und Leiter passen, Parität nicht) und Formen auf dem Sechserkreis (Verlust
und Parität passen, fester Zustand, Vielfachheit und Sprungregel nicht). Keine gibt den Transfer des Originals als
Matrix. Die Sechs bleibt im Original gesetzt; die Kandidatenregel aus SECHS20261119 hat damit keine Stütze gewonnen.

Grenzen: geprüft sind Funktionen des einen Schrittes auf Schalen und der Schritt auf dem Sechserkreis; andere
Operatoren auf den zehn Paaren (etwa mit Gewichten je Paartyp) sind nicht durchsucht. Nichts befördert, kein
Ledgerwechsel, kein Atlaseintrag, kein Commit.

---

# SECHSC20261119: Suche über alle Operatoren, die die Teilung 3 + 2 respektieren

Aufruf: `../../../.venv/bin/python -B search_transfer.py` (rund zehn Sekunden, `results_transfer_search.json`), danach
`check_search.py` (13 von 13, exakt).

Suchraum: alle 15 Operatoren auf den zehn Paaren, die mit den Vertauschungen der drei Farbplätze und der zwei schwachen
Plätze vertauschen. Geprüft: alle 2719 Dreiteilungen der zehn Paare (bis auf Symmetrie), jeweils ob es einen Operator
gibt, der zusammengefasst genau den Originaltransfer ergibt (fester Zustand, Paar mit bleiben 1/2, wechseln 1/6). Zwei
Lesarten des Verlusts von 1/3: er verschwindet, oder er fließt in den festen Zustand.

| Einschränkung | Verlust verschwindet | Verlust wird aufgesogen |
|---|---|---|
| lösbar überhaupt | 48 | 183 |
| mit erlaubter Sprungregel (nichts negativ) | 36 | 57 |
| nur aus dem Graphen (je ein Gewicht für gleich, Nachbar, Nichtnachbar) | 0 | 0 |
| Teilung symmetriefest | 4 | 4 |
| symmetriefest und Regel symmetrisch | 1 | 0 |
| Regel symmetrisch und nur über Nachbarn | 1 | 0 |

1. **Aus dem Graphen allein kommt der Originaltransfer nie heraus.** Kein Treffer in 2719 Teilungen.
2. **Mit freien Gewichten je Paartyp gibt es Dutzende Treffer.** Die Regel ist dann hineingeschrieben, nicht erzwungen.
   Auch die Schalen gehen so (das widerspricht SECHSB nicht: dort waren nur Funktionen des einen Schrittes geprüft).
3. **Eine Teilung sticht heraus:** fest = die drei Farbpaare und das schwache Paar; das tauschbare Paar = die zwei
   Seiten des Sechserkreises (Paare mit w1 gegen Paare mit w2). Nur sie ist symmetriefest mit symmetrischer Regel.
4. **Auf dem Sechserkreis, nur mit den Raten 2/3 und 1/3, bleiben genau zwei Regeln.** Regel A: bleiben 1/2, mit 1/6
   den schwachen Platz wechseln bei gleicher Farbe. Das ist die Paarregel des Originals, dreimal nebeneinander (einmal
   je Farbe). Regel B: bleiben 7/18, zu jedem der fünf anderen mit 1/18.
5. **Regel A ist kein Schritt des Petersengraphen.** Die beiden Paare teilen einen Platz. Die einzige Regel über
   Nachbarn (bleiben 1/2, je Kreisnachbar 1/12) bringt die fremden Raten 7/12 und 5/12 mit.
6. **Die Sechs legt keine dieser Regeln fest.** Die Schrittzahl bleibt von außen gesetzt.

## Urteil

Der Originaltransfer passt auf das Objekt, aber als "schwacher Platzwechsel bei fester Farbe", nicht als Lauf auf dem
Petersengraphen. Die Dreierzahl wäre dann die Zahl der Farben, nicht der Familien; das Original nennt den festen
Zustand Familienkanal. Damit trennen sich die zwei Befunde: Die Zahlen 1, 2/3, 1/3 hoch sechs kommen aus dem Lauf auf
dem Graphen (v808), die Sprungregel des Originals aus dem Platzwechsel. Dass beide dieselben Raten haben, ist exakt;
ein gemeinsamer Operator ist in diesem Suchraum nicht vorhanden. Sechs Schritte von Regel A geben auf den zehn Paaren
1 viermal, (2/3)^6 dreimal, (1/3)^6 dreimal; v808 sagt für den Lauf 1, (1/3)^6 fünfmal, (2/3)^6 viermal voraus. Die
Vielfachheiten unterscheiden die beiden Lesarten, falls das Original je einen Operator auf zehn Zuständen festlegt.

Grenzen: Suchraum sind lineare Operatoren auf den zehn Paaren mit der Symmetrie der Teilung und Zusammenfassen zu
Klassen von Paaren; Zusammenfassen zu Mischzuständen und Operatoren auf größeren Räumen (16 oder 32 Zustände) sind
nicht durchsucht. Die Auswahlkriterien (symmetriefest, symmetrisch, reine Raten) habe ich nach Ansicht der Trefferliste
gewählt. Nichts befördert, kein Ledgerwechsel, kein Atlaseintrag, kein Commit.

---

# SECHSD20261119: Eine Regel für den ganzen Transfer: drei Familienbits, sechs Züge

Aufruf: `../../../.venv/bin/python -B check_cube.py` (23 von 23, exakt, `results_cube.json`).

Anlass: Die Suche auf den zehn Paaren (SECHSC) fand keinen Operator. Der Blick zurück in die Originalstellen zeigt,
dass der Dreierraum des Transfers nicht auf den fünf Plätzen liegt, sondern beim Familienteil: `v69`, `v137`, `v95`
legen ihn als Raum mit Zähler 0, 1, 2 fest (Familienzahl mal Spitzengewicht), mit Wand bei 3 (`v124`).

## Die Regel

Drei Bits (die drei Familienmoden der Teilung (3+2)+3). Ein Zug setzt ein Bit auf einen Wert. Es gibt 3 mal 2 = 6
Züge; je Schritt wird einer gleich wahrscheinlich gezogen. Eine Hand sind so viele Schritte, wie es Züge gibt.

## Was daraus exakt folgt, ohne weitere Wahl

| Angabe im Original | Quelle | Aus der Regel |
|---|---|---|
| bleiben 1/2, wechseln 1/6 (dort aus dem Spektrum gewählt, Herleitung offen) | `v486` | von jedem Zustand: bleiben 1/2, je Nachbar 1/6 |
| Verlust 1/3 | `v486` | die zwei anderen Nachbarn; gleich 1/N nur bei N = 3 |
| Leiter 1 - n/3 | `v124` | Raten 1, 2/3, 1/3 mit Vielfachheit 1, 3, 3 |
| Wand bei n = 3, nur drei Stufen | `v124` | vierte Rate ist genau 0, einmal |
| gerade Paarform 2/3, ungerade 1/3; jede Stufe wechselt die Parität | `v487` | Umklappen aller Bits hat auf Stufe n das Vorzeichen (-1)^n |
| fester Zustand = demokratisches Gesetz | `v327` | Gleichverteilung |
| vier Marken = Familienzahl + 1 | `v69`, `v503` | nach Zahl gesetzter Bits: vier Zustände, jede Rate einmal |
| Hand der Ordnung 6 = 2 N = Zahl positiver Wurzeln von A3 | `v124`, `v327` | sechs Züge; als Operatoren sechs Majoranagrößen, Sektoren 1, 15, 15, 1 |
| lineare Uhr mit Steigung 2 je Stufe | `v124` | Schritt = 1 + G/6 mit G = -2n |
| -M0 hat Ordnung 6 | `v118` | Bits weiterrücken und Wert umklappen durchläuft alle sechs Züge |
| Spektrum 1, (2/3)^6, (1/3)^6 | `v54` | eine Hand: 1, (2/3)^6 dreimal, (1/3)^6 dreimal, 0 einmal |

Eindeutigkeit: Unter den 24 Regeln mit 1 bis 6 Stellen und 2 bis 5 Werten trifft nur (3, 2) die fünf Angaben bleiben,
wechseln, Verlust, zweite Rate und Hand; die zweitbeste trifft zwei. Die Raten legen die drei Stellen fest, bleiben
1/2 die zwei Werte.

## Was nicht hergeleitet ist

1. **Eine Hand = sechs Schritte.** Die Regel sagt "so viele Schritte wie Züge". Das ist eine Setzung, jetzt mit
   Bedeutung (jeder Zug im Mittel einmal), aber nicht erzwungen. Im Original ist dieselbe Stelle offen (`v124` Punkt 5).
2. **Dass die drei Bits die Familienmoden sind.** Die Zuordnung zu P1 und P2 ist nicht gezeigt.
3. **Verhältnis zur Farbseite.** Der Lauf auf dem Petersengraphen (v808, Atlasbefund zum Petersen Quantenkanal) hat
   dieselben drei Raten aus anderem Grund. Ein gemeinsamer Ursprung ist nicht gezeigt.

## Vorhersagen, die die Lesarten trennen

Nach einer Hand auf allen Zuständen: Würfel 1, 3, 3 und eine Rate 0 (acht Zustände); Petersenlauf 1, 4, 5 (zehn);
schwacher Platzwechsel 4, 3, 3 (zehn). Auf Operatoren: Würfel 1, 15, 15, 1; Petersenkanal 1, 24, 75. Nur die
Würfelregel hat eine Stufe, die in einem Schritt ganz verschwindet.

## Suche

Die Regel steht als allgemeine Objektfamilie `moves:{d}x{q}` in `research/originlab/structure_objects.py` (alle d von
2 bis 6, q von 2 bis 5), und jede Gruppe von n Fermionmoden trägt ihren Würfel als Teil. Lauf auf Atlasrevision
`3d4edf28f96138a983aa8968`, 88 Zahlen, 1000 Zufallsvergleiche, Saat 1:

| | ohne die Familie | mit der Familie |
|---|---|---|
| bestes Objekt | (3+2)+3 mit Teilen | (3+2)+3 mit Teilen |
| Wert / Zufallsmedian | 200,15 / 85,17 | 213,90 / 88,46 |
| Anteil am Fingerabdruck | 46,4 % | 51,6 % |
| zurückgehaltener Test p | 0,0008 | 0,0005 |
| unerklärt | 35 | 31 |

Neu getragen: 1/729, 64/729, 1/10, 1/16. Die Familie wurde nach Ansicht der Restliste gebaut; die p Werte
berücksichtigen das nicht.

## Urteil

Keine vollständige Lösung. Gewonnen ist: Zehn getrennt festgehaltene Angaben des Originals zum Transfer folgen aus
einer Regel mit zwei ganzen Zahlen (3 Stellen, 2 Werte). Die Sechs ist darin die Zahl der Züge. Offen bleibt genau eine
Setzung (Hand = Zahl der Züge) und die Bindung der drei Bits an den Ursprung. F1 bleibt offen. Nichts befördert, kein
Ledgerwechsel, kein Atlaseintrag, kein Commit.
