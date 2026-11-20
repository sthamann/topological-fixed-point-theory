# Neuansatz, 20. November 2026: P1 am singulären Punkt

Anlass: Auftrag des Nutzers, von vorn zu beginnen und Ansätze zu suchen, die mit den bisherigen nichts zu tun haben
(Singularität, spätere Prozesse, Doppelüberlagerung mit iterativer Dynamik, Automaten, Holografie).
Rechnung: `check_singular_point.py`, 20 von 20 Kontrollen, Ausgabe `check_singular_point.json`.
Kein Ledgerwechsel, kein Atlaseintrag, nichts befördert. Die Mathematik ist klassisch (komplexe Multiplikation,
Lattès Abbildungen, Isogeniegraphen); neu ist nur die Anwendung auf das Viermarkenobjekt.

## Ausgangsbeobachtung

Die Sphäre mit vier Marken hat als Doppelüberlagerung einen Torus. Der Bestand hält fest, dass der quadratische
Torus (τ = i, j = 1728) gewählt ist. Dieser Punkt ist ein singulärer Punkt des Modulraums. Der Bestand nutzt dort
nur die umkehrbare Symmetrie der Ordnung vier. Am singulären Punkt existiert aber mehr: Der Torus hat dort
Selbstabbildungen, die nicht umkehrbar sind (Multiplikation mit einer Gaußschen Zahl). Im Atlas nicht gefunden.

## Was exakt gilt (alles gerechnet)

1. **Eine Regel, die iteriert.** Multiplikation mit 1 + i steigt auf die Sphäre ab als f(w) = (w² − 1)/(2iw),
   Grad 2. Die vier Marken sind genau die kritischen Werte und ihre Bilder: 1 → 0, −1 → 0, 0 → ∞, ∞ fest.
2. **Die Marken sind nicht gleichberechtigt.** Der Schritt ist auf ihnen nilpotent und liefert eine Fahne
   {Ursprung} ⊂ {Ursprung, c} ⊂ alle vier. Eine Marke ist fest, eine ist ausgezeichnet, zwei werden von der Symmetrie
   vertauscht. Das ist die Gestalt von Z0, ZUY und dem Paar ZU, ZY.
3. **Die Marken sind Durchgangszustände.** Nach zwei Schritten liegt jede Marke im Ursprung. Wiederkehrende Bahnen
   gibt es nur auf ungerader Torsion: Fixpunkte von n Schritten N((1+i)ⁿ − 1) = 1, 5, 13, 25, 41, 65, 113, 225.
   Periode 4 ist die 5 Torsion, Periode 8 die 15 Torsion. Der Schritt hat auf der 3 Torsion Ordnung 8, auf der
   5 Torsion Ordnung 4.
4. **Eine Dilatation, die die Marken erhält.** Der Bestand schließt eine Dilatation durch umkehrbare
   markenerhaltende Abbildungen aus. Der Schritt ist nicht umkehrbar, streckt um √2 und bildet Marken auf Marken ab.
   Entropie log 2 je Schritt, Ljapunowexponent (log 2)/2. Er ist eine Skalenstufe, keine Zeit mit Lichtkegel.
5. **1/(8π) am singulären Punkt.** Periode mal Quasiperiode ist dort exakt π (Legendrerelation plus Symmetrie der
   Ordnung vier), gleichwertig E₂(i) = 3/π. Also c₃ = 1/(8π) = E₂(i)/24. Wert und Begründung (die ergänzte Form verschwindet am Fixpunkt der
   Vierteldrehung) stehen schon im Bestand (P1KERN20261119, Brücke Kugel zu Torus). Nicht neu. Am zweiten
   singulären Punkt gilt E₂ = 2√3/π.
6. **Zwei Füllungen sind entartet.** Die Höhe der Füllung mit Steigung (c, d) ist bei τ = i gleich 1/(c² + d²).
   Die beiden obersten sind gleich hoch und werden von τ → −1/τ vertauscht; dessen Quadrat ist die Deckabbildung
   der Doppelüberlagerung. An einem allgemeinen Punkt ist die oberste Füllung eindeutig. Das ist der Punkt des
   Hawking Page Übergangs.
7. **Ein Baum als Inneres.** Die drei Marken neben dem Ursprung sind die drei Schritte vom Grad 2 aus einem Torus
   heraus. Die so erreichbaren Tori bilden einen Baum der Valenz 3. Bei j = 1728 ist die Ecke gefaltet: Ein Schritt
   führt zurück auf denselben Torus, die beiden anderen auf ein und denselben (j = 287496). Gleichförmiger Gang im
   Baum: nach außen 2/3, zurück 1/3. Die Kugel vom Radius 2 hat 10 Ecken; ohne kurze Kreise geschlossen ist das der
   Graph der Paare aus fünf mit Schrittspektrum 1, 1/3, −2/3.

## Antworten auf die gestellten Fragen

| Frage | Befund |
|---|---|
| Liegt P1 auf einer Singularität? | Ja, im Modulraum (τ = i). Der Bestand kennt die Wahl, nicht die Folgen 1 bis 7. |
| Origin Theorie geprüft? | Gelesen (Kern, Attraktor, Zyklus). Sie setzt einen Attraktor mit Raten (2/3)⁶, (1/3)⁶ und deutet einen kosmischen Zyklus; die Deutung ist dort selbst als nicht hergeleitet markiert. Ihre Skripte v53 bis v56 habe ich nicht ausgeführt. |
| Geschehen die Prozesse erst später? | In diesem Bild ja: Die Marken sterben nach zwei Schritten aus, Wiederkehr gibt es erst auf feinerer, ungerader Torsion. Kennzeichnung: Deutung. |
| Doppelüberlagerung mit iterativer Dynamik? | Ja, exakt: f ist die Regel, und sie existiert nur am singulären Punkt. |
| Automat oder Hypergraphregel? | Der Schritt ist auf jedem Torsionsgitter ein linearer Automat (x, y) → (x − y, x + y). Er erzeugt Skala, keine Zeit. Eine Hypergraphregel mit eigenem Gehalt habe ich nicht gefunden. |
| Holografie? | Der Baum der Valenz 3 ist ein diskretes hyperbolisches Inneres, die Marken sind seine Richtungen, zwei Füllungen sind entartet. Kennzeichnung: Analogie mit exakten Bausteinen. |

## Was das nicht zeigt

- Keine Zahl der Physik ist neu hergeleitet. c₃ ist neu begründet, nicht neu gewonnen.
- 2/3 und 1/3 treten im Baum als Übergangswahrscheinlichkeiten auf, im Bestand als Eigenwerte. Gleiche Zahlen,
  verschiedene Rollen; das ist schwach. Die Sechs der sechs Schritte folgt nicht.
- 5, 15 und die Ordnungen 4 und 8 sind kleine Zahlen. Ohne Gegenprobe gegen andere Schritte ist ihre Nähe zu den
  Uhren des Bestands kein Beleg.
- Dass die Dynamik von P1 dieser Schritt ist, ist eine Hypothese. Geprüft ist nur, dass der Schritt existiert.

## Der Test, der entscheidet: ausgeführt (`check_self_similarity.py`, 13 von 13)

Selbstähnlich heißt: Der mit dem Schritt zurückgezogene Kern ist die Summe des Kerns über die Blätter,
α·k(αz) = k(z) + k(z + c).

| Kern | Herkunft | Ergebnis |
|---|---|---|
| Kern der runden Sphäre | Calderón Seite des Bestands (P1KERN20261119) | nicht selbstähnlich. Zurückgezogen entsteht der Kern mal √((1+x²)(1+y²))/(1+xy): ein zweiter Pol auf dem anderen Blatt und Verzweigung an den kritischen Punkten ±i. |
| Drei Kerne vom Gewicht 1/2 auf dem Torus, gerade Spinstruktur | im Bestand bei τ = i geprüft | nicht selbstähnlich. Der Schritt schiebt sie weiter wie die Marken: die beiden vertauschten auf die dritte, die dritte in die ungerade Struktur. |
| Kern der ungeraden Spinstruktur ohne Nullmode, Z(z) = ζ(z) − π·z̄ | nicht im Bestand gefunden | selbstähnlich, unter 1 + i, unter 2 + i und unter 2. Der rein holomorphe Teil allein ist es nicht; der Term −π·z̄ ist nötig. |

**Urteil:** Für die aufgezeichneten Kerne ist der Ansatz widerlegt. Keiner von ihnen ist unter dem Schritt
selbstähnlich. Es gibt genau einen Fixpunkt der Dynamik auf den Kernen: den Kern des Sektors mit Nullmode. Jede
gerade Spinstruktur läuft in höchstens zwei Schritten dorthin, so wie jede Marke in den Ursprung läuft.

**Was am Fixpunkt steht:** Z(z) = ζ(z) − π·z̄. Der Koeffizient π ist 1/(8·c₃). Er ist die Hintergrundladung der
entfernten Nullmode (∂̄Z = π·(δ − 1/Fläche)). Der Bestand hält fest, dass die Nullmode +c₃ trägt und die Schwingungen
−c₃. Das passt, ist aber dieselbe Zahl E₂(i) = 3/π und kein zweiter Beleg.

**Kontrollen:** An einem allgemeinen Torus gilt dasselbe Gesetz für die ganze Zahl 2. Die Selbstähnlichkeit unter
ganzen Zahlen ist also keine Besonderheit von τ = i; besonders ist dort nur der Schritt vom Grad 2 und der Wegfall
des linearen Terms im holomorphen Teil.

**Offen und prüfbar:** Trägt P1 eine Nullmode? Die Quelle nennt eine Naht mit Einheitswindung. Ist der zugehörige
Sektor der ungerade, dann ist Z der Kandidat für den Kern, und die Dynamik ist der Schritt. Ist er gerade, ist der
Ansatz für P1 erledigt. Zu prüfen an der Kernquelle des Bestands (Randoperator B_Σ, Phasenterm (iπ/2)·Δη).

## Kernquelle gelesen: liegt P1 im Sektor mit Nullmode?

Quelle: `documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/001-…tex`.

- Die Quelle lässt beide Fälle zu. Die Randbedingung ist „modifiziert APS": ohne Kern des Randoperators die negative
  Hälfte, mit Kern zusätzlich ein gewählter Teilraum L des Kerns. Ein Kern ist also vorgesehen, aber nicht behauptet.
- Die Saatantwort ist 1/λ₁⁺(|B_Σ|), der Kehrwert des kleinsten positiven Eigenwerts. Das ist in beiden Fällen
  definiert und entscheidet nichts.
- Die im Bestand ausgeführte Realisierung (runde Sphäre, P1KERN20261119) hat ganzzahliges Spektrum ±1, ±2, … und
  keinen Kern. Die Vierteluhr exp(iπB/2) des Bestands braucht diese Ganzzahligkeit; auf dem quadratischen Torus ist
  das Spektrum 2π·|m + in| und nicht ganzzahlig.

**Urteil:** nicht entschieden. Die Quelle erlaubt die Nullmode, die ausgeführte Realisierung hat keine. In der
ausgeführten Realisierung gilt der Fixpunktkern Z also nicht. Offen bleibt die Realisierung auf der Sphäre mit vier
Marken selbst; dort wäre zu rechnen, ob der Randoperator einen Kern hat.

## Die drei Schritte und die drei Paarungen (`check_pairings.py`, 9 von 9)

Frage aus dem zweiten Fremdbericht: Erzeugt die Theorie selbst die zyklische Versetzung des ausgezeichneten Kanals?

- Jeder Schritt vom Grad 2 hat eine Marke als Kern. Er schickt sie in den Ursprung und verschmilzt die beiden
  anderen. Das ist genau die Gestalt „ein ausgezeichneter Kanal und ein gekoppeltes Paar". Die drei Schritte sind die
  drei Paarungen.
- Die Familienmatrix aus `v395` vertauscht mit dem Tausch ihres Paars. Ihre Eigenlinien (Kanal, Summe, Differenz)
  sind also schon durch diese Symmetrie festgelegt, nicht durch die Zahlen 1/2 und 1/6.
- **Am quadratischen Torus (τ = i)** vertauscht die Symmetrie zwei Marken und hält die dritte fest. Der geschlossene
  Vergleich zweier Lagen gibt 0, 1/4, 1/4 und trägt kein Vorzeichen.
- **Am sechseckigen Torus (τ = ω)** bewegt die Symmetrie die drei Marken zyklisch. Der geschlossene Vergleich gibt
  0, +1/8, −1/8. Die Symmetrie hat auf dem Torus Ordnung sechs, auf den Marken Ordnung drei.
- In der Modulgruppe geben S und T Vertauschungen, ST den Dreierzyklus. Der Zyklus hält die Lage der vierten Marke
  bei exp(±iπ/3) fest (j = 0), die Vertauschung bei −1, 2, 1/2 (j = 1728).

**Urteil:** Ja, aber nicht bei τ = i. Die zyklische Versetzung, die der Fremdbericht annimmt, liefert der Torus
selbst am zweiten singulären Punkt. Am ersten liefert er den nicht umkehrbaren Schritt. Die beiden singulären Punkte
tragen verschiedene Dinge: der quadratische die Skalenstufe, der sechseckige den Umlauf mit Vorzeichen. Der Bestand
ordnet die Familien schon dem Punkt τ = ω zu (Origin Theorie, SEAM.EQUIV.02); die Matrix aus `v395` ist die
Familienmatrix. Das passt zusammen.

**Nicht gezeigt:** dass die drei Kanäle der Familienmatrix die drei Marken sind (Zuordnung gesetzt); dass das
Vorzeichen eine spätere Antwort ändert (zweite Hälfte des Tests im Fremdbericht); wie ein Prozess von einem
singulären Punkt zum anderen kommt.


## Familienregel, Netz und eine gemeinsame Form (Nachtrag, 20. November)

Aufrufe: `check_family_rule.py` (18 von 19; die eine verfehlte Zeile ist eine bewusst festgehaltene falsche
Erwartung von mir), `check_network.py` (10 von 10), `check_unified.py` (12 von 12). Alles exakt.

**Familie.** Jeder der drei Schritte vom Grad 2 vertauscht die vier Marken paarweise (tau_t). Die gleichmäßige
Mischung der drei Operatoren (1 + tau_t)/2 ist (1/3)·I + (1/6)·J: Diagonale 1/2, Nebendiagonale 1/6. Auf jedem
Markenpaar ist das der Paarblock [[1/2, 1/6], [1/6, 1/2]] aus `v395`/`v486`/`v487` mit den Raten 2/3 und 1/3.
Die Zahlen 1/2 und 1/6 des Bestands folgen also aus „halb bleiben, halb einen der drei Schritte anwenden, jeden
gleich oft". Eine feste zyklische Reihenfolge der drei Schritte löscht dagegen alles Nichtkonstante in zwei
Schritten; Transfer mal zyklische Verschiebung hat das Polynom 18x³ − 3x² − 3x − 4 (reell 0,764672, Paar vom
Betrag 0,539084) und trifft 64/729 nicht. **Die Raten des Bestands gehören zur gleichmäßigen Mischung, nicht zur
zyklischen Ordnung.** Der Gruppencharakter von `v117` (3, −1, 0, 1) ist (Funktionen auf vier Marken minus Konstante)
mal Vorzeichen. Auf der Doppelüberlagerung hat die Vierteldrehung Ordnung 8, die Dreierdrehung Ordnung 6.

**Netz.** Die Punkte der Periode vier des Schritts sind die 5-Torsion. Ihre Symmetriegruppe SL2(F5) hat 120
Elemente in neun Klassen; die Spuren der zweidimensionalen Darstellung sind das Spektrum des erweiterten E8
Diagramms, die Dimensionen sind die Marken 1, 2, 3, 4, 5, 6, 4, 2, 3 (Quadratsumme 120). Das ist das Netz mit
Kac Marken aus `v312`/`v395`. Der Schritt selbst hat Determinante 2, kein Quadrat modulo 5, wirkt also von außen
und vertauscht die beiden Klassen der Ordnung 5 (√5 ↔ −√5). Stufe 3 gibt entsprechend das erweiterte E6.
Kontrollen: Stufe 7 (336) und Stufe 4 (48, nicht binär oktaedrisch) geben kein solches Diagramm.

**Gemeinsame Form.** Der Netzschritt des Bestands T_net = (A + 2)/4 hat genau die Eigenwerte cos²(θ/2) über die
neun Klassen, und cos²(θ/2) ist (1 + g)/2 mal seine Adjungierte für ein g mit Drehwinkel θ. Netz und Familie haben
damit dieselbe Regel: **„halb bleiben, halb die Symmetrie anwenden"**, einmal mit der Symmetrie der 5-Torsion,
einmal mit den drei Schritten auf der 2-Torsion. Umgeschrieben mit den Marken ist der Netzschritt stochastisch mit
stationärer Verteilung Marke²/120.

**Was dabei nicht herauskommt.**
- Die Sechs. Jede gerade Potenz des Paarblocks ist positiv; Kandidaten für sechs (Ordnung der Dreierdrehung auf der
  Doppelüberlagerung, sechs Elemente der Symmetrie der 2-Torsion, 2 mal Familienzahl im Original) zeichnen die
  Potenz nicht aus. Stand wie SECHS20261119.
- c₃. Der Anteil 1/6 in φ_base = 1/(6π) ist der Anteil der Elemente der Ordnung 3 in SL2(F5) (wie `v396`:
  20/120). Der Anteil 1/8 in c₃ = 1/(8π) ist keine Vereinigung von Klassen. π selbst kommt im Ansatz nur über den
  selbstähnlichen Kern ζ(z) − π·z̄ vor. Keine Herleitung von c₃, φ₀ oder α.
- Dass die 5-Torsion des Torus das Netz des Bestands *ist* (und nicht nur dieselbe Gruppe hat), ist eine Zuordnung.
- Die Kopplung T_net ⊗ M⁶ ist im Bestand ein Produkt zweier vertauschender Teile; der Ansatz liefert beide Teile,
  aber keinen Grund, warum auf einen Netzschritt sechs Familienschritte kommen.
- Raum, Zeit, Gravitation, die Verbindung der beiden singulären Punkte, die Nullmode an der Kernquelle: unverändert offen.


## Nullmoden und Chiralität an beiden singulären Punkten (Nachtrag, 20. November)

Aufrufe: `check_orbifold_zero_modes.py` (11 von 11), `check_zero_mode_reps.py`, `check_hex_point.py`,
`check_hex_point2.py`, `check_both_points.py`, `check_family_channel.py`. Exakte Gitterrechnung auf den 240 Wurzeln.

**Eingabe aus dem Bestand.** Die Viertelwirkung z -> iz auf den vier Marken hat im Bestand
(`TFPT-DREI-CHATS-QUARTER-ORDNUNG-20261111`) die Winkel d = (0,1,2,-1,0)/4 auf C^5 und a = (3,1,-1,-3)/8 auf der 4
von SU(4); der SU(4) Lift ist bis auf das Zentrum eindeutig. Die Graddimensionen (30,32,...) sind reproduziert.

**Geometrie.** Vier Marken auf mu_4 haben harmonisches Doppelverhältnis, das ist der Punkt tau = i. Die Drehung hebt
auf den Torus als u -> i·u + t; ein Spinor bekommt die Phase exp(±i·pi/4). Eine konstante Nullmode der periodischen
Struktur überlebt genau dann, wenn interne Phase mal Spinorphase eins ist.

**Ergebnis am quadratischen Punkt.**
1. Nullmoden gibt es nur in den 128 Verklebungsrichtungen, keine in den 120 übrigen. 32 je Händigkeit, die
   rechtshändigen sind die Konjugierten der linkshändigen.
2. Die von der Viertelwirkung festgelassene Algebra ist SU(4) × SU(3) × SU(2) × SU(2) × U(1) (30 Dimensionen).
3. Die 32 linkshändigen Nullmoden sind (4, 1, 2, 1) + (4-quer, 3, 1, 2): **eine Familie der einen Händigkeit und
   drei Familien der anderen, in der Form von Pati und Salam.** Das Spektrum ist unter dem halbeinfachen Teil nicht
   reell, also chiral. Mit dem anderen Vorzeichen der Spinhebung tauschen die beiden SU(2) die Rollen.
4. Unter der Farbe und Hyperladung des Bestands (Koordinaten 3 + 2 in D5) sind dieselben 32 Zustände vektorartig:
   16 aus der 16 und 16 aus der 16-quer. Die SU(4) und SU(3) der festen Algebra mischen Träger und Familie
   (Wurzeln aus (45) und aus (10,6)); die Farbe des Bestands wird von der Viertelwirkung gebrochen.
5. Die kubische Anomalie unter SU(4) und SU(3) ist nicht null. Die Nullmoden allein sind also kein vollständiges
   Spektrum; es fehlen mindestens die Beiträge der drei Kegelpunkte (4, 4, 2) des Quotienten. Nicht gerechnet.

**Am hexagonalen Punkt.** Die Dreierwirkung auf D5 steht nicht im Bestand; angesetzt ist dieselbe Regel (Winkel m/3
der Funktionen z^m/(z³−1)) und der kanonische Lift des Dreierzyklus, der anders als beim Viererzyklus nicht
eindeutig ist. Unter dieser Annahme: Hebung der Ordnung 6 lässt SU(5) × SU(4) × U(1) fest, Nullmoden
(10, 4-quer) + (1, 4): vier Familien der 10, keine 5-quer. Hebung der Ordnung 3 lässt SU(9) fest, Nullmoden 84.

**Beide Punkte zusammen.** Die gemeinsam feste Algebra hängt von der relativen Lage der beiden Wirkungen ab, die der
Bestand nicht festlegt. Über alle 2880 Lagen kommt A2 + A1 (die Form der Standardmodellalgebra) in 288 Fällen vor,
in der naheliegenden Lage nur A1 + A1. **Das Standardmodell folgt nicht als gemeinsamer Rest.**

**Familienkanal.** Der im Fremdbericht verlangte Kanal (Sektorwerte 1, 1/3, eta, 2/3) ist keine Mischung der S4
Wirkungen selbst. Die Regel des Neuansatzes als Quantenkanal gibt (1, 1, 1/3, 1/3). Der Kanal des vierten
Fremdberichts (1/2 Identität + 1/3 Transponieren + 1/6 Mitteln) gibt B exakt, aber nur in der komplexen Clockbasis.

**Einordnung.** Neu und exakt: Die Viertelwirkung des Bestands erzeugt am singulären Punkt ein chirales Spektrum mit
der Teilung eins zu drei. Das beantwortet die Nullmodenfrage für diese Wirkung: der ungerade Sektor ist da, und er
sitzt genau in der Verklebung. Nicht gezeigt: dass dieses Spektrum das der Natur ist (falsche Eichgruppe, Anomalie
offen, je Händigkeit nur eine statt drei Familien), dass der Bestand diese Struktur als Felder trägt, alle reellen
Größen, Raum, Zeit, Gravitation.


## Drei neue Ansätze: Gewichtsgesetz, Körper mit zwei Elementen, E8 aus Einheitsquaternionen (Nachtrag, 20. November)

Aufrufe: `check_new_ansaetze.py` (13 von 13), `check_icosians.py` (3 von 3).

**A. Gewichtsgesetz für c₃.** Für jede Form f vom Gewicht k, deren Betrag unter tau -> −1/tau erhalten bleibt, gilt
am Punkt tau = i exakt (1/2πi)·d/dtau log f = k/(4π). Geprüft an eta, theta_3 und E4 auf 15 Stellen. Damit ist
c₃ = 1/(8π) der Fall k = 1/2: die Antwort einer einzelnen Halbdichte (eines freien Bosons) am selbstdualen Punkt. Die
Gittersumme von E8 antwortet mit 8·c₃ = 1/π, also Rang mal c₃. Am hexagonalen Punkt wird daraus 1/(8π·Im tau).
Das gibt P1 einen Inhalt („ein Boson am selbstdualen Torus"), leitet aber nicht her, dass der Rand genau das ist.
phi_base = 1/(6π) wäre Gewicht 2/3; ein Objekt des Bestands mit diesem Gewicht kenne ich nicht. d_top = 48·c₃⁴ ist
kein Wert dieser Art. Die Kubik für alpha folgt daraus nicht.

**B. Der Schritt ist ein Frobenius.** Die Kurve y² + y = x³ + x + 1 hat über den Körpern mit 2, 4, …, 256 Elementen
genau 1, 5, 13, 25, 41, 65, 113, 225 Punkte (ausgezählt). Das sind die Zahlen N((1+i)^n − 1) der Wiederkehrpunkte des
Schritts. Der nicht umkehrbare Schritt am singulären Punkt ist also der Frobenius einer elliptischen Kurve über dem
Körper mit zwei Elementen, die dort nur einen einzigen Punkt hat. Bekannte Sätze dazu, hier nicht nachgerechnet: diese
Kurve ist supersingulär, ihre 2-Torsion fällt in einen Punkt zusammen (die vier Marken werden eine), ihre
Symmetriegruppe hat 24 Elemente (SL2(F3)), ihr Endomorphismenring sind die Hurwitz Quaternionen.

**C. Auswahl der Trägerzahl.** Die Zahl der Punkte der Periode g ist gleich g·2^(g−2) + 1 (das ist 10·b₁ des
Bestands als Funktion von g) genau bei g = 5, geprüft exakt für g < 400. Bei g = 5 sind es 41 Punkte über dem Körper
mit 32 Elementen. Art: Zahlenbeobachtung nach Ansicht des Ergebnisses, keine Herleitung; die Formel für b₁ ist die
des Bestands, nicht meine.

**D. E8 direkt aus der Symmetrie der Periode vier.** Die Symmetriegruppe der 5-Torsion, als 120 Einheitsquaternionen
über Q(√5) geschrieben, spannt ein Gitter vom Rang 8 auf. Mit der Form 2(a ± b) für |x|² = a + b√5 ist es gerade,
ganzzahlig und hat Determinante 1, ist also E8 (bekannte Konstruktion der Ikosianer, hier exakt nachgerechnet). Die
beiden Formen gehen durch √5 -> −√5 ineinander über; genau diesen Tausch macht der Schritt (siehe Netz). E8 ist in
diesem Ansatz kein Eingang mehr, sondern Folge: Kurve über F₂, Punkte der Periode vier, ihre Symmetrie, deren Gitter.

**Was das zusammen ist.** Ein Kandidat für einen einzigen Ausgangspunkt: die elliptische Kurve mit einem Punkt über
dem Körper mit zwei Elementen. Aus ihr folgen exakt der Schritt, die Zahlen 5, 25, 41, das Netz mit Kac Marken, E8.
**Nicht gezeigt:** die Teilung D5 + A3 und die Familienzahl aus dieser Kurve, phi_0, d_top, alpha, die Massen, Raum,
Zeit, Gravitation; dass die Natur diesen Ausgangspunkt hat. Keine fundamentale Lösung.
