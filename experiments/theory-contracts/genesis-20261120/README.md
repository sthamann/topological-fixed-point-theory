# GENESIS20261120 · Suche nach einer Ursprungsregel über Regelklassen hinweg

Stand des Repositorys: `da29a4fcc` (lokal gleich `origin/main`). Atlas Revision `3b1bd8ffcddefa5436f0227c`.
Kein Ledgerwechsel, kein Atlaseintrag, nichts befördert, kein Commit.

**Urteil vorweg:** Es gibt keine neue Gesamttheorie. Zwei kurze Regeln liefern je für sich den festen diskreten Kern des
Bestands (fünf Vorzeichen und 16 Zustände; Partner vom Rang 3, E₈, Klebegruppe Z₄), ohne dass eine Zahl des Bestands
eingesetzt wird. Als eine gemeinsame Regel gelesen widersprechen sie sich: Das ergänzte Objekt trägt keine einzige
spiegelungerade Markierung mehr, die für Proben in allen acht Richtungen unsichtbar ist (exakt, 0 von 64). Sie bleibt
nur bestehen, solange mindestens eine der drei Partnerrichtungen nie geprüft wird. Fast alle Einzelaussagen stehen
bereits im Bestand oder im Lehrbuch. Raum, Zeit, Dynamik, Wahrscheinlichkeiten, drei Familien und jede reelle Konstante
folgen aus nichts davon.

## Was gebaut und ausgeführt wurde

| Teil | Ort | Ergebnis |
|---|---|---|
| Regelsuche (neue Spur in OriginLab) | `research/originlab/genesis_rules.py`, eingebunden in `main.py` und `autonomy.py` | 177 965 Strukturen in sechs Universen, 162 erzeugte Regeln, Laufzeit etwa 10 s |
| Tests | `research/originlab/test_genesis_rules.py` | 14 von 14 bestehen (`test_genesis_rules.log`) |
| Rohausgabe | `genesis_rules_run.json` | alle Regeln, Urteile, Ablationen |
| Sechs blinde Agenten | `blind/b1` bis `blind/b6` | 25 Hypothesen, 12 ausgeführte Aufzählungen mit Gegenproben |
| Red Team | `redteam/report.md`, `redteam/controls.py`, `redteam/controls.out` | acht Gegenproben; Verbindung der beiden Sätze widerlegt |
| Katalog | `HYPOTHESES.md` | 32 Hypothesen mit Art des Belegs |

Wiederholen:

```bash
.venv/bin/python -B research/originlab/genesis_rules.py > /tmp/genesis.json
```

```bash
cd research/originlab && ../../.venv/bin/python -B -m unittest test_genesis_rules
```

Nicht ausgeführt: Lean, Wolfram, die volle Verifikationssammlung `run_all.py`, der autonome Supervisor als Ganzes.
Die blinden Berichte von B1 und B6 konnten nicht als Datei gespeichert werden; ihre Skripte und Ausgaben liegen vor,
der Berichtstext steht nur im Sitzungsprotokoll und ist in `HYPOTHESES.md` zusammengefasst.

## Wie die Regelsuche arbeitet

Ein Universum ist eine vollständig aufgezählte Klasse kleiner Strukturen mit einem Größenmaß und einem Vorrat an
Eigenschaften. Eine Regel ist eine Teilmenge dieses Vorrats. Die Suche erzeugt jede Regel, nimmt die kleinsten
Strukturen, die sie erfüllen, und urteilt: leer, eindeutig oder mehrdeutig. Für jede eindeutige Auswahl meldet sie die
kürzesten Regeln und was geschieht, wenn man einen Satz weglässt. Der Bestand wird erst danach genannt.

| Universum | Strukturen | Kürzeste Regel mit eindeutiger Auswahl, die den Bestand trifft | Ohne einen der Sätze |
|---|---|---|---|
| Binäre Codes bis Länge 8 | 161 640 | selbstdual + doppelt gerade → Länge 8, Hamming Code, 30 Codes in einer Bahn | Länge 2 oder Länge 4 |
| Zweiteilige Schnitte eines geraden selbstdualen Gitters | 11 (Rang 8; darunter keine) | Klebegruppe mit konjugiertem Paar + neutraler Klasse → A₃ \| D₅ | A₂ \| E₆ oder A₁ \| E₇ |
| Ein Vektor und eine Uhr | 3 850 | selbstdual + gerade → Rang 8, 240 Vektoren der Norm 2 | Rang 1 |
| Faltungen auf die Sphäre | 12 298 | mit Formparameter → (2,2,2,2); teilerfremde Drehordnungen → (2,3,5) | (2,2,2) |
| Vorzeichenwörter bis 7 Vorzeichen | 127 | blind bis dritte Ordnung → fünf Vorzeichen, 16 Wörter | drei Vorzeichen, 4 Wörter |
| Umbenennungen von 2 bis 40 Dingen | 39 | es gibt eine zweite Beschreibung → sechs | zwei |

## A. Die überraschendsten Befunde

1. **Drei blinde Agenten kamen unabhängig auf dieselbe Regel und dasselbe Objekt.** B1 (Darstellungen), B2 (Prüfcodes)
   und B6 (Punktmengen im Gitter) fragten nach der kleinsten Menge, die nicht ihr eigener Spiegel ist und für jede
   Probe bis zur dritten Ordnung unsichtbar bleibt. Alle drei erhielten die 16 Gewichte des Halbspinors von D₅. Die
   Regelsuche bestätigt das als Satz über Charaktere des Würfels: Die Lösungen werden von den Charakteren ungeraden
   Grades über 3 aufgespannt, der erste hat Grad 5. Exakt. Folge: Jede lineare Ladung auf den 16 Wörtern ist von selbst
   frei von linearer und kubischer Anomalie. Das ist die Lehrbuchaussage, dass die 16 von SO(10) der kleinste komplexe
   anomaliefreie Spinor ist. Neu ist nur, dass drei verschiedene Ontologien ohne Kenntnis voneinander dort landen.
2. **Der Partner ist erzwungen.** Die fünf Vorzeichen spannen ein Gitter der Determinante 4 auf. Die Gaußsumme seiner
   Klebegruppe hat Phase 5/8; ein Partner in einem geraden selbstdualen Gitter hat deshalb Rang 3 modulo 8. Die
   Aufzählung aller geraden Grammatrizen bis Rang 3 findet genau einen Partner: A₃ mit 12 Wurzeln. Exakt in den
   Schranken.
3. **Der Schnitt A₃ | D₅ braucht genau zwei Eigenschaften der Klebegruppe.** Ein konjugiertes Paar allein wählt
   A₂ | E₆, also die bekannte Alternative mit Familiengruppe SU(3). Erst die zusätzliche neutrale Klasse wählt
   A₃ | D₅. Einschränkung des Red Teams: A₂ + A₁ | A₅ und A₃ + A₁ | A₃ + A₁ tragen beide Eigenschaften ebenfalls;
   A₃ | D₅ gewinnt nur über die kleinste Klebegruppe. Der Satz lautet also: Z₄ ist die kleinste Gruppe mit beiden
   Klassenarten. Der Bestand hat die Eindeutigkeit in v92, v993 und BOOT.CLASS.01.
4. **Vier Marken sind die einzige Faltung mit Formparameter**, und nur bei vier Marken gelten zugleich: Überlagerung
   flach, Moden = zweimal Marken, ungerade Moden = Schleifen der punktierten Sphäre. Im Bestand vorhanden (v189, v214,
   v216, GESAMT20261120); hier nur gemeinsam aufgezählt.
5. **Die Zahl Sechs hat als einzige eine zweite Beschreibung ihrer Umbenennungen.** Diese tauscht den Typ 3+2+1 mit dem
   Sechserzyklus. Im Bestand ist die Sechs die Ordnung von (f₁f₂f₃)(w₁w₂) und zugleich die Länge des Sechserkreises
   (SECHS20261119); dort stehen beide Lesarten nebeneinander. Hier ist die Abbildung zwischen ihnen explizit, auf allen
   720 Elementen gerechnet. Die Sonderstellung der Sechs steht laut Red Team bereits im Bestand (v796). Kandidat, keine
   physikalische Aussage.
6. **Uhren.** Ein Vektor und eine Uhr erzeugen ein gerades selbstduales Gitter erst im Rang 8. Uhren ohne Teiluhr haben
   dann die Ordnungen 15, 20, 24, 30; insgesamt kommen 8, 12, 14, 15, 18, 20, 24, 30 vor. Der Bestand kennt die
   Krylovbasis des Coxeterelements (v223). Die Uhren der Ordnungen 15, 20 und 24 mit zyklischem Vektor fanden weder
   ich noch das Red Team im Bestand.
7. **Das ergänzte Objekt verliert die blinde Chiralität vollständig.** Auf den 128 Klebevektoren (gerade Wörter aus acht
   Vorzeichen, 64 Spiegelpaare) gibt es 49 unabhängige spiegelungerade Markierungen, die für Proben erster und dritter
   Ordnung in den fünf Trägerrichtungen unsichtbar sind. Prüft man zusätzlich eine, zwei oder drei Partnerrichtungen,
   bleiben 38, 22 und 0. Exakt (Rang der Probenmatrix). Folge: Mindestens eine Partnerrichtung darf keine geprüfte
   Ladung sein. Das ist die Gitterfassung des Spiegelfamilienproblems. Für die beiden naheliegenden Markierungen:
   Die Markierung „16 gegen Spiegel" mit vier Familiengewichten (64 Wörter) bleibt blind, solange höchstens zwei
   Partnerrichtungen geprüft werden. Die Markierung mit drei Familiengewichten (48 Wörter, das vierte ausgelassen)
   bleibt nur blind, wenn keine einzige Partnerrichtung geprüft wird.
8. **Aus den blinden Läufen, ohne Bezug zum Bestand:** Unter den einfachen reellen assoziativen Algebren hat die
   Spurform tr(x²) genau eine Richtung abweichenden Vorzeichens nur bei M₂(R) mit (3,1), bei H mit (1,3) und bei C.
   Die Dimension 4 ist die Aussage; welche der beiden die Zeit trägt, hängt an der Vorzeichenwahl. Die Bedingung
   „Zahl der Ausgangsstatistiken = Zahl der wirksamen umkehrbaren Bewegungen" erfüllen nur komplexe Amplituden mit der
   2 Norm (B2, numerisch, 45 Fälle).

## B. Die zehn ungewöhnlichsten Konzepte

| Konzept | Worin es sich unterscheidet |
|---|---|
| Unsichtbarkeit statt Symmetrie (B2.2, G6) | Materie ist, was keine Probe niedriger Ordnung von seinem Spiegel unterscheiden kann. Es wird keine Gruppe vorausgesetzt. |
| Wörter und Umschaltungen sind dasselbe (Variante von G6) | Zustände und Übergänge haben nur bei acht Vorzeichen dieselbe Norm. Materie und Kraft sind dann eine Sorte. |
| Ein Ding ist die Menge seiner eigenen Änderungen (B5.1) | Lie(Aut(V,B)) ≅ V. Wählt Dimension 3, nicht die Signatur. |
| Wirklichkeit als Unverträglichkeit von Projektionen (B5.2) | Je zwei Zerlegungen sind gemeinsam vollständig. Erzeugt Körper für n ≤ 5 und verwirft 6. |
| Jedes Syndrom ist ein Erzeuger (B2.1) | Bedingung an ein einzelnes System, keine Zusammensetzung nötig. Liefert C und die 2 Norm. |
| Zeit als Eliminationsordnung (B3.3) | Es gibt nur eine ungeordnete Bedingung; jede Ordnung, die sie löst, ist eine Zeit. Ungetestet. |
| Zweimal dasselbe tun (B3.1) | Geometrie ist das Vorzeichenmuster von tr(x²). Liefert (3,1) für die kleinste nichtkommutative einfache Algebra. |
| Selbstdeuter (B4.1) | Die Verknüpfung ist ihr eigener Interpreter. Sehr selektiv, aber nicht eindeutig. |
| Starres Ladungsbuch (B6.1) | Ladungsquantelung als isolierte Lösung. Trifft das Standardmodell erst mit der Zusatzklausel „alles geladen". |
| Zweite Beschreibung (G5) | Eine Struktur, deren Umbenennungen sich auf zwei nicht ineinander überführbare Arten lesen lassen. Nur bei sechs. |

## C. Die fünf aussichtsreichsten Modelle

### C1. Blinde Vorzeichen mit selbstdualem Abschluss (G6)

- **Objekte:** endlich viele Vorzeichen; Wörter aus ihnen.
- **Regel:** (i) Die Wortmenge ist nicht ihr eigener Spiegel und für jede Probe erster und dritter Ordnung unsichtbar;
  nimm die kleinste. (ii) Ergänze zum kleinsten positiven geraden Gitter, das sein eigenes Dual ist.
- **Folgen, exakt:** fünf Vorzeichen; 16 Wörter; Gitter der Determinante 4 mit 40 Wurzeln (D₅); Wortnorm 5/4; Partner
  vom Rang 3 modulo 8, eindeutig A₃; zusammen 240 Wurzeln (E₈); Klebegruppe Z₄ mit konjugiertem Paar und neutraler
  Klasse. Mit einer Teilung der fünf Vorzeichen in 3 + 2 und Y = −⅓ auf drei, ½ auf zwei ergeben sich die Ladungen
  {−1, −⅓, −⅙, 0, ½, ⅔} mit Vielfachheiten {1, 3, 6, 1, 2, 3}, eine Familie mit rechtshändigem Neutrino bis auf
  Konjugation.
- **Verbindung zum Bestand:** D₅ Träger, A₃ Familienseite, Z₄ Klebung, E₈. Der Bestand hält fest, dass Gitterdaten den
  Trägerrang nicht erzwingen (D₁₃ + A₃ gibt D₁₆⁺). Satz (i) liefert die Fünf ohne Marken.
- **Gesetzt, nicht hergeleitet:** „dritte Ordnung" (entspricht vier Raumzeitdimensionen), zweiwertiges Alphabet, die
  Teilung 3 + 2 und die Koeffizienten ⅓ und ½, dass Satz (ii) auf Satz (i) folgt. Gegenproben des Red Teams: dreiwertiges
  Alphabet gibt 12 Zustände bei drei Stellen; eine einzelne Probe mit ganzen Ladungen gibt 5 Zustände; ΣY = ΣY³ = 0
  gilt für jede lineare Ladung und ist kein Beleg; die Teilung 5 + 0 gibt Vielfachheiten (1, 5, 10), 4 + 1 gibt
  (1, 1, 4, 4, 6).
- **Widerlegt in wörtlicher Lesart:** Satz (i) gilt für das Ergebnis von Satz (ii) nicht mehr. Die 128 Klebevektoren
  sind spiegelsymmetrisch (64 + 64); mit der D₅ Klasse als Chiralität ist die Probe dritter Ordnung in den drei
  Partnerrichtungen 128 statt 0, nach Auslassen des vierten Familiengewichts sind auch Proben erster Ordnung −32.
  Befund 7 oben zeigt, dass keine andere Markierung das rettet.
- **Offen:** drei gegen vier Familien (A₃ trägt vier Gewichte der Norm ¾), Spiegelfamilien (im Bestand
  CHIRAL4D.MIRROR.DET16.01 als Kandidat), alles Kontinuierliche.
- **Nächster Beweis:** siehe G.

### C2. Wörter und Umschaltungen sind dieselbe Sorte (Einsatzvariante)

- **Regel:** Bei n Vorzeichen haben Wörter die Norm n/4 und Umschaltungen zweier Vorzeichen die Norm 2. Verlange
  Gleichheit.
- **Folgen, exakt:** n = 8; 128 Wörter und 112 Umschaltungen, zusammen 240. Ganzzahlig und gerade ist die Vereinigung
  bei n = 8, 16, 24. Der kleinste blinde chirale Teil hat fünf Vorzeichen, der Rest drei. Die geraden Wörter aus drei
  Vorzeichen bilden die Kleinsche Vierergruppe: ein Einselement und drei weitere.
- **Verbindung:** acht Moden, (3 + 2) + 3, Vierergruppe der Marken (GESAMT2-20261120).
- **Offen:** Die Gleichsetzung der vier geraden Wörter mit den vier Marken ist eine Gruppenisomorphie, keine
  konstruierte Abbildung auf die Marken des Originals. Die Teilung 5 + 3 folgt hier nur aus C1 (i).

### C3. Faltung mit Formparameter (G4)

- **Regel:** eine flache Fläche, auf eine Sphäre gefaltet, die noch eine Form haben kann.
- **Folgen, exakt:** vier Kegelpunkte der Ordnung 2; Überlagerung ein Torus; 8 Moden = 5 gerade + 3 ungerade; nur hier
  sind die ungeraden Moden zugleich die Schleifen.
- **Verbindung:** gemeinsames Viermarkenobjekt. Bekannt; die Regel „hat einen Formparameter" als Auswahl ist die
  kürzeste, die die Suche fand.
- **Offen:** Bündel und Grad sind gesetzt (Rückzug des Markendivisors); τ = i ist im Bestand als kritischer Punkt
  belegt, nicht als Auswahl.

### C4. Drei teilerfremde Uhren (B3.2, G4)

- **Regel:** Uhren teilerfremder Ordnung, deren Produkt die Identität ist, erzeugen eine endliche Gruppe.
- **Folgen, exakt:** nur (2,3,5), Ordnung 60; Doppelüberlagerung mit McKay Graph Ẽ₈; Produkt 30.
- **Verbindung:** Coxeterzahl 30 = 2·3·5 und Z₃₀ ≅ Z₅ × Z₆ im Clocknetz.
- **Offen:** Endlichkeit ist ein eigenes Axiom. Eine Abbildung von den drei Uhren auf C₅, C₆ und μ₄ des Bestands fehlt;
  μ₄ liegt dort auf der Galoisseite.

### C5. Zwei Spuren und eine Statistikbedingung für Raumzeit und Quantenmechanik (B3.1, B2.1)

- **Regeln:** tr(x²) mit genau einer negativen Richtung auf einer einfachen nichtkommutativen Algebra; Ausgangsstatistiken
  und wirksame Bewegungen gleich zahlreich.
- **Folgen:** M₂(R) mit (3,1) exakt in n ≤ 8; komplexe Amplituden mit 2 Norm numerisch.
- **Verbindung zum Bestand:** keine gefunden. Das ist der einzige Kandidat der Sitzung für Signatur und Amplituden, und
  er steht unverbunden neben C1 bis C4.
- **Offen:** (3,1) gegen (1,3); „einfach" trägt die Auswahl; die Form |·|^p ist gesetzt.

## D. Verbindungen

| Verbindung | Stand |
|---|---|
| Blinde Vorzeichen ↔ Halbspinor von D₅ ↔ kleinste anomaliefreie chirale Bahn | bewiesen (drei unabhängige Rechnungen, Charaktersatz) |
| Fünf Vorzeichen → Partner A₃ → E₈ | bewiesen in den Schranken (Gaußsumme, Aufzählung bis Rang 3); im Bestand als q(D₅) + q(A₃) = 2 |
| Ergänzung ↔ Verlust der blinden Chiralität (49, 38, 22, 0) | bewiesen; im Atlas nicht gefunden |
| Hamming Code → Konstruktion A → 240 Vektoren; Uhrgitter → 240 Vektoren | bewiesen; im Bestand bekannt |
| Sechs: Typ 3+2+1 ↔ Sechserzyklus | bewiesen als Abbildung; Deutung für den Transfer Kandidat |
| Vier Marken ↔ gerade Wörter aus drei Vorzeichen; sechs Marken ↔ gerade Wörter aus fünf | Gruppenisomorphie bewiesen; geometrische Deutung Kandidat |
| (2,3,5) ↔ Clocknetz | gleiche Zahlen, keine Abbildung |
| M₂(R) (3,1) ↔ Bestand | keine |

## E. Gescheitert oder eingeschränkt

- Kleinste anomaliefreie Ladungsmenge ohne weitere Bedingung: nicht eindeutig (B1.4, B6.1).
- Starrheit der Ladungen: wählt 14 Zustände, nicht das Standardmodell (B6.1).
- Zyklischer Mischungsansatz mit drei Familien: θ₁₂ = θ₁₃ = 0, widerlegt (B6.4).
- Selbstdeuter: zwei Klassen bei n = 4, keine bei 5; keine Auswahl (B4.1).
- Zulässige Geschichten mit Kreuzungssymmetrie: rein topologisch, keine Ausbreitung (B4.4).
- S² = (n−1)·1: erfüllt auch n = 14, 18, 30; Signatur (3,3) statt (3,1) (B1.2).
- Autologie: Signatur nicht gewählt; entartete Lösung bei n = 2 (B5.1).
- Chiralität allein wählt im Gitterschnitt E₆ × SU(3), nicht den Bestand (G2).
- Vorzeichenregel ohne Zweiwertigkeit: 15 statt 16 Zustände (B2, B6).
- Blinde Vorzeichen und selbstdualer Abschluss als eine Regel: widerlegt (Befund 7).
- „ΣY = ΣY³ = 0 bestätigt die Hyperladung": kein Beleg, gilt für jede lineare Ladung.
- Keine Regel der Sitzung liefert drei Familien, eine Dynamik, Gravitation oder eine reelle Konstante.

## F. Beste verbleibende Gesamtfassung

Die Zweisatzregel aus C1 ist als Gesamtregel gefallen. Was bleibt, ist eine Fassung mit einer benannten Trennung:

> Es gibt acht Vorzeichen; Wörter und Umschaltungen sind dieselbe Sorte Vektor. Fünf der acht Richtungen werden
> geprüft, von den übrigen drei höchstens zwei. Wirklich ist eine Wortmenge, die sich von ihrem Spiegel unterscheidet
> und die keine Probe bis zur dritten Ordnung in den geprüften Richtungen davon unterscheiden kann.

Hergeleitet, exakt: 8 als einzige Stellenzahl mit gleicher Norm; 240 Vektoren; die kleinste blinde chirale Teilmenge
braucht fünf Richtungen und hat 16 Wörter; der Rest hat drei Richtungen, seine geraden Wörter bilden die Vierergruppe;
blinde Chiralität überlebt genau dann, wenn nicht alle drei Restrichtungen geprüft werden.

Gesetzt: dass es Proben und geprüfte Richtungen gibt; die dritte Ordnung; die Auswahl der fünf; welche der 49
Markierungen wirklich ist; die Teilung 3 + 2.

Nicht erfüllt: Quantenformalismus und Wahrscheinlichkeiten; Raum, Zeit und Lorentzsymmetrie; Eichdynamik; Gravitation;
drei Familien; Massen und Mischungen; c₃, φ₀, α. Die Fassung ist eine Auswahlregel für den diskreten Kern und kein
physikalisches Gesetz. Gegenüber dem Bestand sagt sie nichts Neues über die Natur; sie gibt dem Unterschied zwischen
Trägerseite (geprüft) und Familienseite (nicht vollständig geprüft) eine rechenbare Form.

## G. Das entscheidende Experiment: ausgeführt

Frage: Koppelt ein Strom des ergänzten Objekts an eine Familienmarkierung? Gerechnet mit `current_probes()` in
`genesis_rules.py` für alle 112 Wurzelströme außerhalb der 128 Klebevektoren (40 des Trägers, 12 des Partners, 60
gemischte vom Typ (10,6)). Probe eines Stroms r auf einer Markierung f: Σ f(w)·⟨r,w⟩ und Σ f(w)·⟨r,w⟩³, dazu die
gemischte dritte Ordnung für mehrere Ströme.

| | drei Familien (48 Wörter) | vier Familien (64 Wörter) |
|---|---|---|
| 40 Ströme des Trägers | blind | blind |
| 12 Ströme des Partners | 6 sehen sie in erster Ordnung (±32) | einzeln blind; zusammen sichtbar (gemischte dritte Ordnung) |
| 60 gemischte Ströme (10,6) | alle 60 sehen sie in erster und dritter Ordnung (±16) | einzeln blind; drei Ströme in drei Partnerrichtungen: 128 |
| Wirkung der gemischten Ströme als Übergang | 16 Wörter in den Spiegel, 8 aus der Markierung heraus | alle 32 in den Spiegel |
| größte blinde Strommenge unter den geprüften | nur die 40 des Trägers (Rang 5) | Träger und gemischte Ströme in zwei Partnerrichtungen: 84 Ströme, Rang 7 |

**Ausgang:** Der ungünstige Fall ist eingetreten. Jeder der 60 Ströme vom Typ (10,6) sieht die Dreifamilienmarkierung
schon in erster Ordnung. Auch die sechs Partnerströme, die in erster Ordnung blind sind (die Familiengruppe mit Rang 2),
sehen sie in gemischter dritter Ordnung. Drei chirale, blinde Familien vertragen im ergänzten Objekt also nur die
Ströme des Trägers als Proben. Das stimmt mit dem Satz des Bestands überein, dass kein (10,6) Strom den Raum der drei
Familien erhält, und beziffert ihn. Vier Familien bleiben für eine Algebra vom Rang 7 mit 84 Strömen blind.

**Was daraus folgt:** Der E₈ Rahmen kann drei blinde chirale Familien nur tragen, wenn die 60 gemischten und die 12
Partnerströme keine Ladungen messen, sondern ausschließlich Übergänge sind. Ob der Bestand das so festlegt, ist die
nächste Frage; sie steht in den Stromverträgen (`seam-collar-multifield-20261016/e8_boundary_response.py`,
`kernel-followups-20261120`), nicht in dieser Aufzählung. Standardwissen dahinter: E₈ hat nur reelle Darstellungen,
die 3 und die 4 der Familiengruppen haben eine kubische Anomalie.

### Lesart der Ströme im Bestand (`check_current_reading.py`, 9 von 9)

Gerechnet mit den Antwortfunktionen des Bestands selbst (`seam-collar-multifield-20261016/e8_boundary_response.py`,
nur gelesen) und den Markenformeln des E₈ Wortvertrags.

- **Im Markenvertrag** (Atlas `premise:e8-covariant-current-word-cartan-marker-20261116`): Der gemessene Marker ist ein
  Cartanstrom des Trägers (6Y = (−2,−2,−2,3,3) oder T₃). Der (10,6) Strom tritt nur quadratisch auf, als
  Tr P X₋γ X_γ über beide Chiralitäten. Der spiegelungerade Anteil von Z0 und ZU ist null. Der von ZY wird im Bestand
  nicht gebildet; getrennt gebildet ist er ±8 für einen einzelnen Strom und hebt sich zwischen γ und −γ auf. Keine Größe
  erster oder dritter Ordnung in einer Partnerrichtung kommt vor. In diesem Vertrag bleibt die Dreifamilienmarkierung
  blind.
- **Im vollen lokalen E₈ Netz** desselben Bestands: Alle acht Cartanströme sind lokale Felder mit exakter Wardregel
  (Residuum d·w am Materiepaar, nachgerechnet). Die drei Partnerströme geben auf der Markierung je −16. Das Paar
  V_α V₋α eines (10,6) Stroms misst in erster Ordnung die Ladung α·w, auf der Markierung summiert −16. Ein
  Familienmarker (1,1,1) gibt −48.

**Urteil:** Die Antwort hängt davon ab, welche der beiden Algebren die physische ist. Im Markenvertrag sind die
Partnerrichtungen ungeprüft, und drei chirale Familien sind blind. Im vollen Netz sind sie geprüft, und die
Blindheit fällt. Der Bestand führt beide als verschiedene Algebren (Kompression auf 96 Zustände ist kein
Homomorphismus, Atlas `finding:no-bifundamental-e8-current-preserves-three-family-space-20261116`). Neu ist die
Schärfung: Der Unterschied ist genau die Frage, ob die drei Cartanströme des Partners lokale Observablen sind. Der
Bestand kennt dazu schon einen Hinweis in dieselbe Richtung (Sitzungsdokumentation vom 20. September: die
Quellentransformation J hat auf den acht Cartanrichtungen keinen Fixvektor); er ist nicht geprüft und hier nicht
verwendet.

**Was eine Ursprungsregel also leisten muss:** die fünf Cartanströme des Trägers als Observablen liefern und die drei
des Partners nicht. Das ist eine Bedingung an die Observablenalgebra, die sich prüfen lässt.

**Nächster Test:** die Abbildung der fünf Vorzeichen auf die fünf geraden Moden des Viermarkenobjekts, unter der V₄ die
16 Wörter erhält; und im Stromvertrag des Bestands nachlesen und nachrechnen, ob die (10,6) Ströme dort eine Ladung
tragen, die gemessen wird.

## Befunde des Bestands nach Art (Kurzfassung des interpretationsfreien Katalogs)

- **Mathematische Tatsachen:** D₅ ⊕ A₃ mit Z₄ Klebung ist E₈; ein Divisor vom Grad 4 gibt Rang 3 und Dimension 5;
  Clockidentitäten; G₃₁ als Stabilisator; Lean Dateien prüfen ganzzahlige Identitäten, die Brücken zur Naht stehen dort
  als Axiome (`SeamIsE8`, `mmst_scaling_limit`, `rigidity_from_clock_invariance` und weitere).
- **Modellbefunde:** Werte der Gruppen A1 bis A11 in `THEORY-BRIEF.md`, je in ihrem Vertrag; Struktursuche: (3+2)+3
  trägt rund die Hälfte der wiederkehrenden Zahlen des Bestands.
- **Deutungen:** c₃ = 1/(8π) als Nahtnormierung, φ₀, die α Gleichung, Gravitation als Ausgabekanal; 49 Zeilen
  `predictions.txt`, davon die meisten bedingt.
- **Offen und Gegenbeispiele:** eine gemeinsame ausgeführte Quelle; die vier Werte des R113 Endpunkts; Zeit mit Länge;
  drei gegen vier Familien; Spiegelfamilien; isospektrale Verklebungen mit gleichen Einzelmarken.

## Red Team

Vollständiger Bericht mit Rohausgaben unter `redteam/`. Urteile:

| Aussage | Urteil |
|---|---|
| Blinde Vorzeichen → 16 | hält; bekannt (Lehrbuch). Binäres Alphabet und dritte Ordnung tragen |
| ΣY = ΣY³ = 0 als Beleg | widerlegt als Beleg |
| Vielfachheiten einer Familie | geschwächt: folgen aus der gewählten Teilung 3 + 2, Werte aus gesetzten Koeffizienten |
| Selbstdualer Abschluss → A₃, E₈ | hält; bekannt (v92, v993, BOOT.CLASS.01, E8.GLU.01) |
| Beide Sätze als eine Regel | widerlegt; hier exakt nachgerechnet (49, 38, 22, 0; drei Familien nur bei 0 geprüften Partnerrichtungen) |
| Uhren, Faltungen, (2,3,5), Sechs | bekannt (v223, v214, v216, v189, v219, v1004, v796); Uhren 15, 20, 24 nicht gefunden |
| M₂(R) mit (3,1) | geschwächt: H hat (1,3); Dimension 4 hält |
| Sätze von Distler und Garibaldi, Spiegelfamilien | treffen die Verbindung der Sätze |
| Nielsen und Ninomiya, Weinberg und Witten, Coleman und Mandula | nicht berührt, weil nichts Dynamisches behauptet wird |

Schätzung des Red Teams zum Auswahlrisiko: Acht von neun Übereinstimmungen mit dem Bestand sind durch Minimalität in
Universen erzwungen, die in Kenntnis des Bestands gewählt wurden. Die Universen der blinden Agenten sind davon
ausgenommen; ihre Übereinstimmung betrifft nur die 16.

## Quellen und Suchanfragen

- Atlas: `graph_status`, `theory_overview`; Suchen nach Springer und zyklotomischen Idealgittern, nach der
  Galoiswirkung auf Z/5, nach der Eindeutigkeit der Z₄ Klebung, nach Kummer und S₆, nach kubischer Blindheit.
  Gefunden und berücksichtigt: SEAM.GALOIS.01, E8.ST31.01, BOOT.CLASS.01, v92, v993,
  `finding:coxeter-seventh-power-has-complete-thirty-e8-lattice-lifts`,
  `finding:tfptfuture20261109-paired-z4-glue-without-rank-selection-has-d16plus-counterexample`,
  `finding:source-inversion-forces-six-line-theta-completion`.
- „Nicht gefunden" ist kein Beleg für Neuheit.
- Gelesen: `AGENTS.md`, Skill `tfpt-knowledge`, `research/originlab/THEORY-BRIEF.md`,
  `STRUCTURE-SCENARIOS-20261119.md`, `RESULTS.md` (Einträge vom 19. und 20. November), `predictions.txt`,
  `verification/status_ledger.csv` (788 Zeilen, Statusspalten), Axiomzeilen der Lean Dateien.
