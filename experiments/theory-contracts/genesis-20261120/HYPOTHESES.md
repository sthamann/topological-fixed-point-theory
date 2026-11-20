# Genesis, 20. November 2026: Hypothesenkatalog

32 Hypothesen. B1 bis B6 stammen von sechs Agenten, die ohne Zugriff auf Repository, Atlas und Netz gearbeitet haben
(Rohberichte und Skripte unter `blind/`). G1 bis G7 stammen aus der Regelsuche `research/originlab/genesis_rules.py`.
Kennzeichnung: **E** exakt hergeleitet oder vollständig aufgezählt, **A** unter benannter Zusatzannahme, **N** numerische
Beobachtung, **H** Analogie, **U** ungetestet. „Gerechnet" heißt: in dieser Sitzung ausgeführt.

| Nr. | Was existiert | Regel | Ergebnis | Art | Gerechnet |
|---|---|---|---|---|---|
| B1.1 | eine irreduzible unitäre Wirkung | kleinste Darstellung, die nicht ihr Spiegel ist und deren kubische Spur verschwindet | 16 von Spin(10), eindeutig; dann 27, 64 | E in der Klasse | ja, alle einfachen Typen bis Dimension 1000 |
| B1.2 | eine symmetrische ±1 Relation | S² = (n−1)·1, kleinstes n | n = 6, eine Klasse, sechs gleichwinklige Geraden im R³, Gruppe der Ordnung 60 | E | ja, n ≤ 8 |
| B1.3 | eine Verknüpfung, deren Elemente ihre eigenen Endomorphismen sind | x∗(y∗z) = (x∗y)∗(x∗z), ein Erzeuger | Laver Tafeln der Größe 2^k | U | nein |
| B1.4 | endliche Ganzzahlmengen | Σq = Σq³ = 0, keine Teillösung | Größe 5, nicht eindeutig | E, wählt nichts | ja |
| B2.1 | Unterscheidungen mit Größenfunktion | Zahl der Ausgangsstatistiken = Zahl der wirksamen umkehrbaren Bewegungen | nur komplexe Amplituden mit 2 Norm | N, Form |·|^p gesetzt | ja, 45 Fälle |
| B2.2 | Syndrome von n Vorzeichenprüfungen | kürzester Vorzeichencode, nicht spiegelsymmetrisch, ungerade Momente bis Ordnung 3 null | n = 5, 16 Zustände, eindeutig bis Spiegel | E, binäres Alphabet trägt | ja, vollständig n ≤ 5 |
| B2.3 | eine Wiederherstellbarkeitsrelation | genau einer von Teil und Komplement ist wiederherstellbar, lokal sphärisch | Dimension in {0, 2, 4, 8, 16} (erinnerter Satz) | U | nein |
| B2.4 | Menge mit fixpunktfreien Involutionen | fünf Paarungen, 5 Würfel als 4 Sphäre | D/2 + 3 = D + 1 nur bei D = 4 | H, schwach | nein |
| B3.1 | Akte als einfache reelle assoziative Algebra | Q(x) = tr(x∘x) ist die Geometrie | genau eine negative Richtung nur M₂(R) (3,1), C (1,1), H (1,3) | E; „einfach" trägt | ja, n ≤ 8 |
| B3.2 | k zyklische Uhren teilerfremder Ordnung | x₁⋯x_k = 1, Gruppe endlich und nicht trivial | nur (2,3,5), Ordnung 60; Doppelüberlagerung mit McKay Graph Ẽ₈ | E; Endlichkeit gesetzt | ja, 120 Tripel |
| B3.3 | eine ungeordnete dünne Bedingung Aφ = 0 | Zeit = jede Eliminationsordnung, die A dreieckig macht | Signatur mit einer Zeitrichtung | U | nein |
| B3.4 | vollständige Paarungen | ⟨a\|b⟩ = d^Schleifen, positiv | d = 2cos(π/k) oder d ≥ 2 | bekannte Sätze | nein |
| B4.1 | Menge mit einer Verknüpfung | x → (y → x∗y) ist Bijektion auf die eigenen Endomorphismen | 1, 1, 1, 2 Klassen bei n = 1..4, keine bei 5 | E, nicht eindeutig | ja |
| B4.2 | F₂ Raum von Unterscheidungen | V ≅ V* durch eine alternierende Form | gerade Dimension, Pauli Darstellung | E, nur Stabilisatortheorie | nein |
| B4.3 | ein unendliches Wort | Wort = eigene Lauflängenbeschreibung | Kolakoski Wort | U | nein |
| B4.4 | zeitlose zulässige Geschichten | je zwei Beine bestimmen das dritte, Bein und Kreuzungssymmetrie | Alphabet F₂^k; topologisch, keine Ausbreitung | E, spricht gegen Gravitation | ja, n ≤ 9 |
| B5.1 | Raum V mit einer Bilinearform | Lie(Aut(V,B)) ≅ V | n = 3: so(3) und so(2,1); n = 4, 5 keine gefunden | N für n ≥ 4 | ja |
| B5.2 | Familie von Zerlegungen einer Menge | je zwei gemeinsam vollständig und nicht redundant | n ≤ 5: Körperstruktur; n = 6: kein Partner | E, bekannte Sätze | ja, n ≤ 6 |
| B5.3 | eine Relation M mit zwei Randverteilungen | M Mᵀ = Mᵀ M | 2×2 Blöcke a + ib | E, vermutlich klassisch | nein |
| B5.4 | Relation R | K = X ∖ R[K] | Auswertungsordnung; ungerade Zyklen erzwingen 1/2 | U, CHSH Test offen | nein |
| B5.5 | mögliche Wörter | Kugeln um die wirkliche Menge pflastern den Raum | perfekte Codes, kein Kontinuum | H | nein |
| B6.1 | Buch ganzzahliger Ladungen | kleinste chirale, anomaliefreie, starre Ladungsmenge | 14 Zustände; Standardmodell (15) erst mit „alles geladen" | E; Zusatzklausel trägt | ja |
| B6.2 | eine endliche Punktmenge im Gitter | kleinste einzelne Weylbahn, chiral, erste und dritte Momente null | 16 von D₅; zwei Bahnen: 15 von A₄ | E in der Klasse | ja, Rang ≤ 8 |
| B6.3 | geschlossene 4D Geschichten | Amplitude auf jeder geschlossenen Mannigfaltigkeit eindeutig | Zustandszahl 0 mod 16 | A (bekannter Z₁₆ Satz) | nein |
| B6.4 | eine zyklische Ordnung | kleinste Uhr mit nichtreellem Charakter | n = 3; Mischungsansatz in exakter Form widerlegt | A, teils widerlegt | ja |
| G1 | binäre Prüfungen, die sich gegenseitig prüfen | selbstdual und doppelt gerade | Länge 8, 30 Codes in einer Bahn (Hamming), Konstruktion A: 240 Vektoren der Norm 2 | E | ja, 161 640 Codes |
| G2 | Schnitt eines selbstdualen geraden Gitters in zwei Teile | Klebegruppe hat ein konjugiertes Paar und eine neutrale Klasse | A₃ \| D₅; ohne neutrale Klasse A₂ \| E₆ | E | ja, 11 Schnitte, Rang < 8 leer |
| G3 | ein Vektor und eine Uhr | sie erzeugen ein selbstduales gerades Gitter | erst im Rang 8; Uhren ohne Teiluhr: 15, 20, 24, 30 | E in den Schranken | ja, 3 850 Gitter |
| G4 | eine auf eine Sphäre gefaltete Fläche | flach und mit Formparameter | (2,2,2,2), vier Marken; Moden = 2·Marken nur bei vier | E | ja, 12 298 Faltungen |
| G5 | Umbenennungen von n Dingen | es gibt eine zweite Beschreibung | n = 6; tauscht Typ 3+2+1 mit dem Sechserzyklus | E | ja, n ≤ 40 und alle 720 |
| G6 | Vorzeichenwörter, dann Gitter | (i) nicht eigener Spiegel, blind bis Ordnung 3; (ii) selbstdual und gerade | 5 Vorzeichen, 16 Wörter, Partner eindeutig A₃, zusammen E₈ | E je Schritt; Verbindung der Schritte gesetzt | ja |
| G7 | Marken einer verzweigten Überlagerung | Zustände = gerade Teilmengen der Marken modulo Komplement | 4 Marken ↔ 3 Vorzeichen (V₄), 6 Marken ↔ 5 Vorzeichen (16) | Gruppenisomorphie E; Deutung H | nur die Isomorphie |

## Ontologien (mindestens zwölf verschiedene)

Wirkung einer Gruppe (B1.1), Relation (B1.2, B5.3, B5.4), selbstbezügliche Verknüpfung (B1.3, B4.1), Statistik mit
Größenfunktion (B2.1), Prüfcode (B2.2, G1), Wiederherstellbarkeit (B2.3), Algebra der Akte (B3.1), Uhren (B3.2, G3, B6.4),
ungeordnete Bedingung (B3.3), Schleifenpaarung (B3.4), Unterscheidungsraum (B4.2), selbstbeschreibendes Wort (B4.3),
zulässige Geschichten (B4.4, B6.3), Projektionen (B5.2), Ladungsbuch (B6.1), Punktmenge im Gitter (B6.2), Schnitt eines
Gitters (G2), Faltung (G4, G7), Umbenennung (G5).
