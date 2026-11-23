Gebiete: endliche Quantenmechanik und Weil-Darstellungen (diskreter Phasenraum) · holografische Tensornetze und hyperbolische Pflasterungen · Topologie von 3-Mannigfaltigkeiten · nichtassoziative Algebra (Jordanalgebren) · Renormierung als dynamisches System · Thermodynamik der Berechnung

# Exploration 23. November 2026 (zweite Runde): Phasenraum über 𝔽₅ und minimaler holografischer Code

Ausgangslage: `docs/research/TFPT-STAND.md` mit den konsolidierten Befunden vom 20. bis
23. November. Fast jede verbleibende Eingabe läuft auf zwei Fragen hinaus: woher die Vier
kommt (μ₄, vier Marken, d = 4, Viertelzählung in c₃) und welches Prinzip Minimalität
verlangt. Gesucht sind Prinzipien, die beides zugleich liefern, ohne ein TFPT Objekt
vorauszusetzen.

Vielfaltsprüfung: Die letzten Ordner nutzten Operatoralgebren, Quanteninformation und endliche
Geometrie, Fehlerkorrektur (perfekter Code), Kobordismus und Anomalien, Arithmetik über 𝔽₂,
Extremalprinzipien, Quantenstatistik, Zufallsmatrizen, ART in Spinorform, Kristallographie,
AQFT. Neu in dieser Runde: holografische Tensornetze, Weil-Darstellungen endlicher
Phasenräume, 3-Mannigfaltigkeiten, Jordanalgebren, Renormierung, Thermodynamik der Berechnung.

## Ansätze

**W. Kleinster endlicher Phasenraum mit Spinor (Weil-Darstellung).**
Grundobjekt: ein Quantensystem mit endlichem Phasenraum 𝔽_p × 𝔽_p, Heisenberg–Weyl
Verschiebungen und der vollen symplektischen Symmetrie SL(2, 𝔽_p), die auf ℂ^p durch die
Weil-(Oszillator-)Darstellung wirkt. Die Vierteldrehung des Phasenraums ist die diskrete
Fouriertransformation S, S² ist die Parität x ↦ −x, S⁴ = I.
Prinzip: Der „fermionische“ Teil der Darstellung (auf dem das zentrale −I als −1 wirkt) soll
ein einzelner Spin-½ Dublett sein. Das legt p fest.
Soll liefern: 5 Plätze, 3 + 2, die Hyperladung als einzige spurfreie mit der ganzen Symmetrie
verträgliche Größe, μ₄ als Fourier-Vierteldrehung mit Spektrum {1,1,−1,i,−i} (TFPT R₅),
die binäre Ikosaedergruppe 2I ≅ SL(2,5) und damit E₈ (Ikosiane, McKay).
Gemessene Größe: keine Zahl, nur Struktur (5, 3+2, Y Verhältnis −2/3, E₈).
Freie Parameter: keiner. Erste freie Wahl: das Auswahlprinzip „Spinorteil = ein Dublett“.

**H. Minimaler holografischer Quantencode.**
Grundobjekt: ein Tensornetz aus perfekten Tensoren auf einer hyperbolischen Pflasterung
{p, q} (wie HaPPY), das einen gekrümmten Innenraum auf einen Rand abbildet.
Prinzip: kleinster Code, der eine logische Einheit gegen jeden Einzelfehler schützt
(Quanten-Singleton: n ≥ 5 Plätze), auf der kleinsten hyperbolischen Pflasterung mit diesem
Fünfeck ((p−2)(q−2) > 4 gibt q = 4).
Soll liefern: 5 und 4 aus Minimalität, die Dreiecksgruppe (2, 4, 5), ihr kleinstes endliches
Bild (zu prüfen: S₅ ≅ PGL(2, 5), die schon bekannte TFPT Uhrengruppe), eine Flächen-
Entropie-Beziehung (Ryu–Takayanagi) als Gravitationsanschluss.
Gemessene Größe: keine. Erste freie Wahl: Holografie (negativ gekrümmter Innenraum) und
„eine logische Einheit, Abstand 3“.

**P. Poincaré-Sphäre als einzige ladungsfreie nichttriviale 3-Mannigfaltigkeit.**
Prinzip: Der Rand eines Raumgebiets hat keine Homologie (keine Ladung sichtbar), aber
nichttriviale endliche Fundamentalgruppe. Einzige Lösung Σ(2,3,5) = S³/2I; Spin-Füllung mit
E₈ Schnittform, Rokhlin 16.
Abgleich vorab: Das steht schon in TFPT (v219, McKay, „2I perfekt ⇒ Poincaré-Sphäre“).
Daher nicht gerechnet; als bereits erledigt geführt.

**J. Ausnahme-Jordanalgebra J₃(𝕆).**
Prinzip: Die Observablen eines elementaren Systems bilden eine formal reelle Jordanalgebra;
die einzige nicht darstellbare ist J₃(𝕆) (Jordan–von Neumann–Wigner). Literatur: das
Standardmodell als Schnitt zweier maximaler Untergruppen von F₄ (Dubois-Violette und Todorov),
drei Generationen aus 3 × 3.
Soll liefern: drei Familien, (SU(3)×SU(2)×U(1))/ℤ₆. Fraglich: 5 = 3 + 2 (𝕆 = ℂ ⊕ ℂ³ gibt 1 + 3).
Kleinster Test: Wirkt Spin(9) ∩ SU(3)-Stabilisator auf einem natürlichen Fünferraum mit 3 + 2?
Nicht gerechnet in dieser Runde.

**R. Renormierung als dynamisches System.**
Prinzip: Die Physik ist der Fixpunkt einer Vergröberung von Codes; die Schichten einer
hyperbolischen Pflasterung sind eine diskrete Renormierungsgruppe (MERA-artig). Die
Wachstumsrate der Schichten von {5,4} ist ein RG Eigenwert.
Soll liefern: Skalendimensionen, eventuell die Transferzahlen (2/3)⁶, (1/3)⁶.
Kleinster Test: Wachstumsrate der {5,4} Schichten berechnen und mit den Atlaszahlen
vergleichen. Nur als Nebenrechnung in H.

**T. Thermodynamik der Berechnung.**
Prinzip: c₃ = 1/(8π) = (1/4)·(1/2π) als Landauer/Bekenstein-Viertel mal Unruh-Periode.
Abgleich vorab: Die Lesart „Einstein 8π = 1/c₃, Bekenstein 1/4 = 1/|μ₄|“ ist in TFPT bereits
erfasst. Nicht gerechnet.

## Begründung der Auswahl zum Rechnen

W und H setzen kein TFPT Objekt voraus und könnten 5 und 4 zugleich erzeugen. P und T sind
bereits in TFPT vorhanden und zählen nicht als neuer Weg. J braucht große Rechnungen in F₄ und
hat ein frühes Gegenargument (1 + 3 statt 3 + 2); R ist eine Nebenfrage zu H.

## Kleinste entscheidende Tests

- **W:** Weil-Darstellung von SL(2, p) für p = 3, 5, 7, 11, 13 exakt aufbauen; Gruppenordnung,
  Zerlegung in gerade/ungerade, Treue und SU(2)-Wertigkeit des ungeraden Teils, reelle
  3D-Drehung im geraden Teil, Kommutant, Spektrum der Fouriervierteldrehung. Abbruch, wenn der
  ungerade Teil bei p = 5 kein treues SU(2) Dublett ist.
- **H:** [[5,1,3]] als perfekten Tensor prüfen; Tabelle der kleinsten hyperbolischen q; alle
  erzeugenden Tripel von (2,4,5) in S₅; kleinstes glattes endliches Bild; Quadrat des
  Ordnung-4-Erzeugers auf ℂ⁵. Basislinie: dieselbe Rechnung für (2,3,7), (2,4,6), (2,5,5),
  (2,3,8).
