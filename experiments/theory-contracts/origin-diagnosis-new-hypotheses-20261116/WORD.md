# Das geordnete E8 Wort auf dem Materieraum

16. November 2026 · Fortsetzung von ORIGINDIAG20261116 · Kennung E8WORD20261116

Geltungsbereich: E8 kovarianter Stromvertrag im festen konformen Vakuum des
Niveau eins Netzes, D8 Chart. Marker U ist ein neutraler Cartanstrom, Y ein
(10,6) Wurzelstrom. Die Zuordnung zum ursprünglichen skalaren Higgsmarker und
zur physischen Zeitkarte ist nicht hergeleitet. Kein Ledgerwechsel.

## Die vier Marken

Für einen Projektor P auf Materiewurzeln und einen Mittelstrom γ aus (10,6):

    Z0  = Tr P
    ZU  = Tr P h h'
    ZY  = Tr P X_{−γ} P' X_γ
    ZUY = Tr P h X_{−γ} P' h' X_γ

P' = 1 lässt den Rückweg durch den ganzen Raum zu (bedingte Erwartung).
P' = P hält den Zwischenzustand im behaltenen Raum (Kompression). Alle vier
Größen sind quadratisch in den Strukturkonstanten und hängen nicht von der
Wurzelsektion ab. Der ortsabhängige Faktor ist getrennt: an den Marken 1, i, −1
ergibt der Zweipunktkern −1/4 und der Dreipunktkern −i/4.

## Ergebnis auf dem vollen Raum (128)

- ZY/Z0 = 1/4 für jeden der 60 Mittelströme.
- ZU/Z0 = |h|²/4 für jeden Cartanmarker. Für die Hyperladung 5/24, für T3 1/8.
- ZUY = ZY · (Kovarianz von h quer zu γ − (h·γ)²/4), exakt für alle 60 Ströme
  und vier geprüfte Marker.
- Hyperladungsmarker: ZUY/Z0 = 1/48 für schwache Mittelströme mit (Y·γ)² = 1/4
  und 11/288 für Farbströme mit (Y·γ)² = 1/9.
- Die verbundene gemischte Antwort ZUY − ZU·ZY ist für jeden Mittelstrom strikt
  negativ: −1/32 für schwache und −1/72 für Farbströme.

Im E8 kovarianten Vertrag sind die vier Marken damit ohne freien Parameter
festgelegt. Das Vorzeichen der gemischten Antwort folgt aus der Stromalgebra.

## Drei Familien

Für jede der vier Wahlen der ausgelassenen Familie und jeden Mittelstrom:

| Größe | bedingte Erwartung | Kompression | voller Raum |
| --- | --- | --- | --- |
| ZY/Z0 | 1/4 | 1/6 | 1/4 |
| ZU/Z0 (Hyperladung) | 5/24 | 5/24 | 5/24 |
| ZUY/Z0 (Hyperladung, schwacher Strom) | 1/48 | 1/72 | 1/48 |
| ZUY/Z0 (Hyperladung, Farbstrom) | 11/288 | 11/432 | 11/288 |

Die bedingte Erwartung reproduziert für alle Trägermarker die Werte des vollen
Raums. Die Kompression verliert ein Drittel der Stromantwort. Das fehlende
Drittel sind die Rückwege durch die vierte Familie. Familienmarker sehen die
ausgelassene Familie: ZU für den Familiencartan F1 ist 1/3 oder 2/3 statt 1/2.

## Einordnung

Der Befund vom selben Tag (R113PRIMITIVEHISTORY20261116) verlangt, erst nach
der ganzen Wortkomposition zu projizieren. Hier ist das am ursprünglichen
Gitter ausgerechnet: Die unbehaltene Richtung ist die vierte Familie, und mit
ihr als Rückweg entsteht der starre E8 Wert. Der native Aufruf auf Full96
entspricht der Kompression und kann diesen Wert nicht liefern.

Offen bleibt die Zuordnung. Der ursprüngliche Higgsmarker ist skalar und
nichtnegativ, ein Cartanstrom ist es nicht. Erst ein hergeleitetes Wörterbuch
zwischen beiden macht aus den Zahlen eine Aussage über R113. Ebenso offen sind
Feldmaß, Zeitkarte und die Auswahl, welche Wortpaare der Ursprung verkettet.

Reproduktion: `python3 -B check_word.py` (Exit 0, Ausgabe in `run_word.log`).
