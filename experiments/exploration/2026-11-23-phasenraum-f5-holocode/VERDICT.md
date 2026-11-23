# Urteil

## Aktueller Stand: H6h13 — die vier Entscheidungen und die TOE-Lückenkarte (`z21`–`z24`)

Bilanz unverändert bei 15 P und 5 L von 20. Drei Fragen liefen parallel und wurden danach selbst
nachgerechnet: z21 (Lückenkarte), z22 (Prinzip U(1)_J) und z23 (die zwei 16er). Alle drei liefern
identische Ausgaben; die JSON von z22 weicht nur in der 15. Stelle ab (BLAS). z24 ist eine eigene
Rechnung, die erst beim Zusammenlegen von z22 und z23 nötig wurde.

**In einem Satz:**
- Entscheidung 3 löst sich auf.
- Entscheidung 4 ist das Gitter-Chiralitätsproblem selbst; kein Prinzip löst sie ohne Abstimmung.
- Entscheidung 1 bleibt eine Lesart mit genauer Eingabeliste.
- Entscheidung 2 bleibt offen.

**1. Die zwei 16er (`z23`, 27 Prüfungen, 11 Befunde).**
- **Keine Vereinigung.** Die Naht-16 (Vektor von SO(16)) liegt in keinem E₈-Multiplett:
  exp(2πi·(½)⁸) ist in E₈ trivial, auf 16_v aber −1. Kein Automorphismus, keine Trialität und keine
  Wahl der Untergruppe führt sie in eine Generations-16 über. Grund: Reell bleibt reell, und die
  Generation ist nicht selbstkonjugiert.
- **J-Satz.** Ein reelles Multiplett mit Paarstruktur J gibt ξ₋ = konj(ξ₊). Die Familien aus Z18
  wären also 3 × [(10,1)+(1,6)], vektorartig. Umgekehrt erzwingen Familien 3 × 16 eine Naht in 16b,
  und die gibt so(32)₁ statt E₈.
- **Wörterbuch Ströme → Fermionen.** Die 128 ist der gerade Fockraum der 16 Nullmoden und enthält
  (16,3). Die Quantenzahlen stimmen also, das Wörterbuch scheitert aber zweifach:
  - an der Statistik: h ist ganzzahlig, und die E₈-Klammer ist antisymmetrisch;
  - an der Chiralität: Unter Fluss hebt die vierte (16,1) die drei (16,3) auf, und das ℤ₄-Orbifold
    gibt netto stets eine gerade Zahl, nie 3.
- **Konsistent ist nur die direkte Summe.**
  - Familien F_k in der 16: leicht, mit den Charakteren {−1, i, −i}.
  - Eltern P_a in (10,1)+(1,6): schwer, über Φ P P.
  - Das Ganze ist anomaliefrei und trägt ein gemeinsames Spin-ℤ₄ J = z_X⁻¹ ⊗ R_Paar mit ν ≡ 0 mod 16.

**2. Prinzip U(1)_J (`z22`, 22 Prüfungen, 16 Befunde).** Hypothese: Das interne Deck ist ℤ₄ ⊂ U(1)_J.
- **Paarung.** ψχ = (i/2)(ξ₊ξ₊ − ξ₋ξ₋) ist nicht kovariant, Φ koppelt nur an ξ₊ξ₊. U(1)_J erlaubt
  aber y′ Φ\* ξ₋ξ₋. Satz: Jede Symmetrie, die das Wilson-Glied ξ₊ξ₋ und Φξ₊ξ₊ erhält, erhält auch
  Φ\*ξ₋ξ₋.
- **Gitter.** Das Phasenpotential ist exakt flach, und Wilson r ≤ 1 ist für alle Phasen linkpositiv.
- **Index.** Mit y′ = 0 zählt der Kern netto eine Mode je Paar (L = 14 … 38). Geschützt ist sie nicht:
  - det M am Knoten ist 0, ξ₋ bleibt lückenlos;
  - die Nahtmode hybridisiert, das Kerngewicht fällt von 0,84 auf 0,73;
  - jedes y′ > 0 gibt netto 0.
- **Eine Naht.**
  - Mit ℤ₄-wertigen Holonomien ist genau ein Schlauch minimal, die Zahl der
    Nullstellen-Schläuche ist ungerade, und die Z11-Parität gilt automatisch.
  - Neu: U(1)_J quantisiert die Holonomie nicht.
  - Unter U bevorzugt die Fermionenergie die Orientierung ohne Mode (+0,0037 je Paar). Unter G
    tritt das Problem nicht auf: Beide Orientierungen tragen eine Mode (+1 bzw. −1), die
    Orientierung legt nur das Vorzeichen der Chiralität fest.
- **Darstellungen.** Ein Singulett-Φ verlangt ein reelles R. Für die Spin(10)-16 hat 16 ⊗ 16 keine
  Invariante.
- **Ein Block.** Weder die U(1)_J-Anomalien noch Spin-ℤ₄ noch die Klebung legen N fest.

**3. Die Wand unter G (`z24`, 11 Prüfungen, eigene Rechnung).** Nach z23 braucht der
Eltern-Sektor allein eine voll gegappte Paarung mit geschütztem Index, also G. Damit kehrt die Frage
aus Z19 zurück: Ist die Wand ein Gitterartefakt?
- **Kein Artefakt.** Hamiltonsch gilt ΔV = 0,254 m², konvergiert auf 4·10⁻⁷. Euklidisch ist der
  Exponent 2,00, für die lokale GJ-Form ebenfalls 2,00. Die Fourier-Harmonischen gehen wie m², m⁴
  und m⁶. Das ist der relevante Operator Λ² Re(Φ²): J-gerade, von G nicht verboten. In physischen
  Einheiten wächst ΔV/M⁴ wie (Ma)⁻².
- **Erklärung.** Die Kontinuumstheorie von G hat eine axiale U(1)_A: ξ₊ und ξ₋ tragen beide die
  Ladung +1, die Symmetrie ist gravitativ anomal, und daher ist der Index geschützt. Wilson bricht
  sie so, wie bei Wilson-QCD die chirale Symmetrie gebrochen wird. Zum Wiederherstellen braucht es
  zwei abgestimmte Gegenterme, Re Φ² (relevant) und Re Φ⁴ (marginal), analog zur kritischen Masse.
- **Satz.** Eine exakte, ortslokale U(1), unter der Φ geladen ist, lässt einem gegappten Wirbel netto
  den Index 0. Der Grund: det M(k) windet an allen Ecken gleich, und Σχ_c = 0 (Nielsen–Ninomiya).
  Geprüft an drei Modellen:
  - G: Index 2 je Umdrehung, also 1 je Halbwindung;
  - U mit y′: Index 0;
  - ein anomal geladenes Modell (+1, +1, −1): Index 0, mit k-unabhängiger Windung.
- **Korrektur zu z22.** „U löst Entscheidung 4 vollständig“ gilt nur im Gitterteil. Die flache Phase
  kommt dort von einer vektorartigen U(1)_J, und eine solche kann den Index nach dem Satz nicht
  schützen.

**Was das für die vier Entscheidungen heißt.**

| Entscheidung | Stand H6h13 | Was bleibt |
|---|---|---|
| 1: genau eine Naht | Lesart: Z20 und Z22 geben genau einen Schlauch als Minimum, die Parität stimmt | Deck intern, Holonomie physisch und ℤ₄-wertig, das Reziproke spontan gebrochen; empirisch N_ν (LEP) |
| 2: genau ein Block | offen | die Forderung „E₈ auf der Naht“ (N = 1); keine Rechnung legt N fest |
| 3: Paarung nur auf der Nahtkomponente | aufgelöst | Die Frage entstand nur im gemeinsamen Multiplett (Z18), das Z23 ausschließt. Getrennt kann ein Singulett-Φ die chiralen Familien gar nicht paaren (16 ⊗ 16 ohne Invariante); die Familienmassen 0, ½, ½ entfallen. Preis: kein gemeinsamer Operator. Gemeinsam bleiben J, die Abbildung T und die Geometrie (Marken, Deck, μ₄). |
| 4: Wand oder Positivität | als Gitter-Chiralitätsproblem erkannt | Geschützte Naht, lokale positive Transfermatrix und flache Phase ohne Abstimmung schließen sich aus. Zur Wahl stehen: G mit zwei abgestimmten Gegentermen (der Standardweg), U mit ungeschütztem Index oder nichtlokale bzw. nicht positive Formen (Z19). Bei Ma ~ 1 ist die Abstimmung eine physische Wahl. |

**4. TOE-Lückenkarte (`z21_toe_luecken.md`, Zensus `z21_ledger_zensus.py`, 14 Prüfungen).**
- **Zensus.** Das Ledger hat 784 aktive Zeilen (E 192, C 216, O 325). Nur 38 der offenen sind ganz
  offen.
- **Der Dachvertrag.** TFPT.TOE.COMPLETE.01 [O] zerlegt die Vollständigkeit in acht Tore T1–T8. Alle
  acht sind [O]. H6 trägt fast nur zu T1 (Ursprung von P1) und T2 (Naht = (E₈)₁) bei. Die Bilanz
  15 P misst also T1/T2-Struktur, nicht TOE-Vollständigkeit.
- **Was am meisten blockiert, in dieser Reihenfolge:**
  1. S2: von der 2D-Naht zur lokalen, wechselwirkenden 3+1D-Theorie (T3/T5);
  2. S3: chirale Materie ohne Spiegelgeneration (T4);
  3. S0: Ursprung (T1);
  4. S5: Gravitation über Jacobson hinaus (T7);
  5. S1: E₈ auf Netzebene (T2);
  6. S4: Kopplungen, Higgs, Flavour (T6);
  7. S6: Anfangszustand (T8);
  8. S7: Statistik.
- **Durch Rechnung allein nicht schließbar:**
  - die Einheit v_geo;
  - die Rahmenprämissen;
  - die Forderung „E₈ auf der Naht“;
  - die Wahl zwischen D und einem mikroskopischen Spin-2;
  - eine nichtentropische Anfangszustandswahl;
  - die absolute Neutrinoskala;
  - δ_PMNS und θ₂₃;
  - die Kosmologiegabel;
  - α_s.
- **Querbefunde:**
  - Z20s Einzelnaht steht in Spannung zur räumlichen Linse S³/ℤ₄ von CCC.SEAM.CROSSOVER.01 (D1).
  - FLAV.STRONGCP.02 braucht Reflexionspositivität, gilt also im G-Zweig mit abgestimmten
    Gegentermen.
  - sin²θ_W = 3/8 an der Gitterskala trifft auf den [X]-Befund, dass reine SM-Läufe nicht vereinigen.
- **Bester Ertrag je Aufwand:**
  - der H6-RG-Test (zwei Schleifen von M_Z bis zur Gitterskala, SM gegen Pati–Salam);
  - Lüschers projizierte Weylfelder auf Reflexionspositivität;
  - Sakharov gegen Jacobson (macht a/ℓ_P zur Vorhersage oder widerlegt D);
  - der explizite H6-Hamiltonian für T3.

  Die D3-Suche, die z21 noch empfahl, ist durch z23 und z22 erledigt.

## Davor: H6h12 — der Deck-interne Zweig von v957 (`z20_deck_intern.py`)

v957 (CCC.SEAM.KINEMATICS.01) lässt eine Gabel offen: Das Deck D = diag(i, −i) wirkt entweder als
Raumquotient (Linse X = S³/ℤ₄) oder als Strukturgruppe eines internen Bündels über S³. Z17 zeigte,
dass im räumlichen Zweig eine einzelne Naht unmöglich ist. Z20 rechnet den internen Zweig; alle zehn
Prüfungen bestehen.

**Nachgebaut aus v957.**
- Das Deck ist die Halbdrehung z → −z der Nahtsphäre mit D² = −1, die Uhr die Vierteldrehung, das
  Reziproke z → −1/z.
- D wirkt frei auf S³.
- Jede der vier Markenfasern verschlingt die Polfaser genau einmal, und die beiden Polfasern bilden
  eine Hopf-Verschlingung (Gauß-Integral ±1,000).

**Räumlicher Zweig.** Monome z₀ᵃ z₁ᵇ haben Hopf-Ladung k = a + b und Deck-Charakter i^{a−b}. Bei k = 1
gibt es kein invariantes Monom, das kleinste invariante k > 0 ist 2. Generisch entstehen also zwei
Nullfasern: Nähte kommen paarweise.

**Interner Zweig.**
- **Flussschlauch erzwungen:** π₁(S³) = 0, also ist ein flaches internes ℤ₄-Bündel ohne Defekte
  trivial. Soll die Holonomie χ(D) = i aus v957 physisch sein, muss sie um eine Defektlinie liegen,
  einen Flussschlauch. Das Komplement eines Kreises erlaubt jede Meridian-Holonomie, auch J mit
  J² = −1 = (−1)^F.
- **G folgt:** Φ ψχ ist unter J ungerade (Jᵀ M J = −M), also trägt Φ die J-Ladung 2 und wechselt um
  den Schlauch das Vorzeichen. Eine stetige Funktion mit dieser Monodromie hat halbzahlige Randwindung
  (in 200 von 200 Proben), muss also im Inneren verschwinden. Φ verschwindet auf dem Schlauch und
  windet um π: **Das ist genau der verdrillte Defekt aus G, jetzt hergeleitet statt gesetzt.**
- **Lage:** Unter der Uhr sind nur die Pole fest. Ein einzelner uhreninvarianter Schlauch liegt also
  auf einer Polfaser, wo v957 genau die volle μ₄-Holonomie χ(D) = i führt.
- **Zahl:** Verlangt man das Reziproke von der Konfiguration, braucht es beide Polfasern (zwei Nähte,
  Hopf-verschlungen). Die Zustände mit einer Naht bilden aber ein Dublett unter dem Reziproken; der
  Grundzustand wählt einen Pol und bricht das Reziproke spontan. Minimal ist genau eine Naht.

| Zweig | Nähte | G hergeleitet | c | N_fam |
|---|---|---|---|---|
| räumlich (Linse) | 2 | nein | 16 | 6 |
| intern, Reziprokes verlangt | 2 | ja | 16 | 6 |
| intern, Reziprokes spontan gebrochen | 1 | ja | 8 | 3 |

(c und N_fam unter der Lesart „je Naht ein 16er-Satz und drei Familien nach C″“.)

**Urteil: Lesart, deutlich gestärkt.** Der interne Zweig leitet den verdrillten Defekt her und
erlaubt genau eine Naht. Eingaben sind:
- das Deck wirkt intern;
- seine Holonomie ist physisch (der v522-Umlauf);
- die Konfiguration muss das Reziproke nicht selbst tragen.

Die spontane Wahl eines Pols ist eine Wahl zwischen den beiden Äonen-Enden, also eine Zeitrichtung.

## Davor: H6h11 — die drei offenen Verbindungen (`z17`, `z18`, `z19`)

Bilanz unverändert bei 15 P und 5 L von 20. Die drei Verbindungen, die H6h10 als offen nannte, wurden
parallel gerechnet und danach selbst nachgeprüft. Alle drei Skripte laufen reproduzierbar mit
identischer Ausgabe (35, 23 und 20 Prüfungen); negative Ergebnisse stehen in `befunde`.

**1. Wahl des Sektors 1 (`z17`).** Kein Kandidat wählt den primitiven Sektor und genau einen
windenden Block als Eigenschaft des Grundzustands. Beide Fragen verengen sich auf je eine Lesart.

| Kandidat | Ergebnis | Urteil |
|---|---|---|
| Ordnungsfeld als Datum der Doppelüberlagerung w² = q | Torus mit nirgends verschwindender Form nur für k ∈ {3, 4}; bei k = 1 gibt es keine holomorphe Form | gescheitert |
| Spin und Statistik | die Monodromie von √q ist (−1)^k, also zirkulär; nach Spin^c-Buchhaltung unter G gibt es nie ungerade viele Nähte | gescheitert |
| Energie und Entropie der Hopf-Sektoren | E(k) steigt streng, k = 0 ist Grundzustand; die Entropie log 12 ist vernachlässigbar (βE/log 12 ≥ 2,8·10³); ein Fenster für genau eine Naht gibt es nur mit Feinabstimmung | gescheitert |
| Ordnungsfeld ohne Mittelwert | jede generische l = 1 Mode verschwindet auf genau einem Großkreis mit Windung 1; nahe der Schwelle wählt der Quartterm die Kreismode; trägt bis R/ξ ≈ 82 | Lesart (Mittelwert null aus A\* übertragen) |
| Primitive Kopiendeterminante | primitive Determinante gibt genau einen Rangverlust; unter G winden aber alle gekoppelten Kopien am selben Ort, und genau ein Block braucht Rang Y = 1 (Kodimension 3); mit Familiencharakteren hat Y Rang 2 oder 3 | Lesart |
| TFPT-Quellen | die Einheit ist ein Postulat; tfpt-42 029:203 setzt N_fam = 3n mit kanonischem n = 1. Familien je Naht geben N_fam = 3k; mit LEP N_ν = 2,984 ± 0,008 (extern) liegt k = 1 bei +2σ, k = 2 bei +368σ | Lesart (empirisch) |
| Gemeinsame GSO-Summe unter A\* | n Blöcke geben D₈ₙ⁺; n = 3 ist so(48)₁ (D₂₄⁺), ohne E₈, mit sichtbarem Niveau 1 und c = 24. α schließt es nicht aus, nur die Forderung „E₈ auf der Naht“ | Lesart |

**Neu und wichtig:** In TFPTs eigenem Hopf-Lift (v957, CCC.SEAM.KINEMATICS.01) sind die vier Marken
vier verschlungene Hopf-Fasern. Jedes Feld mit Hopf-Ladung 1 trägt dort die Deck-Holonomie ±i; das
kleinste Deck-invariante Feld hat Ladung 2 und gibt zwei Nahtstränge. Wirkt das Deck auf den Raum,
ist eine einzelne Naht dort also unmöglich, im Einklang mit Z11. Eine einzelne Naht verlangt den
Zweig von v957, in dem das Deck intern wirkt.

**2. Gemeinsamer Operator mit benannten Kopplungen (`z18`).** Gerechnet ist ein Multiplett auf dem
Nahttorus E. Der Dirac koppelt an K_E^{1/2} ⊗ π\*K(D) = K_E(R); R sind die vier Marken auf E, die
Spinstruktur ist (P,P). Physisch sind nur Zustände, die unter σ̃ ⊗ J⁻¹ invariant sind (geeichtes
Deck).
- **Hergeleitet:**
  - J² = −1 ist erzwungen, weil σ an den Marken um π dreht.
  - Ein Bündel O(R) ersetzt die zwei Bündel aus Z12.
  - Die Nullmoden teilen sich in 3 deckgerade (die Familien mit {−1, i, −i}) und 1 deckungerade
    (ω = dz/w = √q). Laut Lefschetz geht 3 + 1 nur bei m = 1.
  - Die Eichbedingung ordnet ξ₊ der Naht und ξ₋ den Familien zu; 16 Paare geben 16 Nahtmoden und
    48 Familienkomponenten.
  - E₈ entsteht nur bei richtiger Orientierung; die falsche gibt so(48).
  - Ein Verflechter ι (z → 1/z) erfüllt ι² = σ.
- **Exakte Obstruktion:** Der Nahtkanal liegt in den Familienkanälen (q = f₂² − f₀²), und das
  Paarungsprofil ist eindeutig s ∝ q. Mit Gs Paarung ψχ bekommen die Familien Massen 0, ½ und ½ mal
  die Nahtmasse. Keine Symmetrie trennt die Kanäle.
- **Weiter Eingabe:**
  - die Paarung Φ ξ₊ξ₊ (Projektor auf J = +i), die ψχ aus G widerspricht;
  - die Eichdarstellungen (chirale 16 von SO(10) gegen reelle 16 von SO(16));
  - die Spinstruktur (P,P);
  - m = 1.

  Urteil: Lesart mit klarer Obstruktion.

**3. Ginsparg–Wilson-Paarung und Positivität (`z19`).**
- Für Majoranafelder gibt es nur die symmetrische chirale Drehung δλ = iαXλ mit X = γ₅(1 − D/2),
  deren Spektrum nicht quantisiert ist. Die exakte Drehidentität aus Z16 (Dirac) überträgt sich nicht.
- **Phasenpotential:**
  - Wilson: 0,083 (euklidisch) bzw. 0,088 (wie Z9).
  - Lokale, eindeutige Majorana-Form (GJ): 0,0073 mit Periode π.
  - Formen X und Xq: flach. X ist nicht 2π-periodisch und kostet am Schnitt; Xq ist nichtlokal.

  **Korrektur zu H6h10:** Die Barriere ist nicht bloß ein Wilson-Artefakt; ohne Wand geht es nur
  nichtlokal oder mit Kosten am Schnitt.
- **Index:** eine chirale Mode je Paar mit allen Formen; die Overlap-Verdoppler tragen nichts bei.
- **Reflexionspositivität:**
  - Naiv und Wilson r ≤ 1 sind linkpositiv, r = 2 nicht.
  - Overlap bei θ = 0 ist positiv (Pol plus Schnitt mit positivem Gewicht). Die Behauptung von Lu–Shao
    ist für das freie vektorartige bzw. Majorana-Overlap nicht reproduziert; Lüschers projizierte
    Weylfelder sind nicht geprüft.
  - Die flachen Formen verletzen die Positivität bei gedrehter Phase, also dort, wo der Defekt
    durchläuft (bis −0,23).

  **Abwägung für H6:** entweder Wand (Wilson, positiv) oder Positivitätsverletzung auf Gitterskala
  (flache Formen).
- **Gravitation:** Kein chiraler Regulator macht a/ℓ_P bestimmter. SLAC hat kein Flächengesetz,
  Overlap hat keinen endlichen lokalen Hilbertraum.

**Was das für H6 heißt.** Die drei Verbindungen sind jetzt genau benannt:
- Für genau eine Naht braucht es ein Ordnungsfeld ohne Mittelwert oder den Deck-internen Zweig von
  v957.
- Für genau einen Block braucht es die Forderung „E₈ auf der Naht“.
- Für den gemeinsamen Operator braucht es eine Paarung, die nur ξ₊ξ₊ trifft, statt ψχ.

Für die Zeitrekonstruktion bei endlichem Gitterabstand muss H6 wählen zwischen einer Wand mit
Positivität und einer Positivitätsverletzung ohne Wand. Keine dieser Entscheidungen folgt bisher aus
einem Prinzip.

## Davor: H6h10 — dritte externe Prüfrunde (`z16_kritik_drei.py`)

Bilanz unverändert bei 15 P und 5 L von 20. Drei externe Berichte vom selben Tag wurden geprüft. Ihre
Pakete (PDF, ZIP, ein Lean-Zertifikat mit 15 Lemmas) liegen nicht im Repo; nachgerechnet ist alles,
was sich aus dem Text rechnen lässt. Alle 21 Prüfungen in Z16 bestehen.

**Bestätigt und übernommen.**
- *Einheitswindung als relative Chernzahl.* Für die explizite Masse
  M(z) = ½[[1+z, i(z−1)], [i(z−1), −1−z]] gilt exakt:
  - det M = −z, M² = zI, Singulärwerte {1, |z|};
  - Pf(M ⊗ iσ₂) = det M;
  - Windung 1 auf dem Rand und genau ein Rangverlust bei z = 0.

  Im Sektor 1 ist der Rangverlust also erzwungen. Die Wahl des Sektors bleibt offen.
- *Äquivariante Abbildung vom Träger auf die Familien.* T(f) = (f − P₀f) dz/z ist surjektiv von
  H⁰(O(D)) auf H⁰(K(D)), hat als Kern genau die zwei invarianten Trägermoden und vertauscht mit der
  Vierteldrehung. Die Familiencharaktere (i, −1, −i) sind damit der nichttriviale Charakteranteil
  des Fünferträgers (1, i, −1, −i, 1).
- *Geometrischer Ursprung der Marken.* Eine verzweigte Doppelüberlagerung w² = s(z) der Sphäre mit
  Torus braucht nach Riemann–Hurwitz genau r = 4 Verzweigungspunkte. Unter z → iz invariant ist nur
  s = a z⁴ + e; das gibt w² = z⁴ − 1 mit Doppelverhältnis −1 und j = 1728. Vier Marken, zwei Blätter,
  Torus und τ = i stammen so aus einem Objekt (Ordnung-4-Symmetrie vorausgesetzt).
- *Nahttorus.* Zurückgezogen hat O(2) dort 3 gerade und 1 ungerade Schnitte, O(4) hat 5 gerade und 3
  ungerade (wie im Atlasbefund vom 19.11.). Ein Elternoperator muss die leichten Fermionen an den
  geraden Teil koppeln.
- *Hopf-Modell.* Ein Ordnungsfeld mit Ladung 1 unter der Hopf-U(1) auf S³ hat als tiefste Moden genau
  A z₀ + B z₁. Ihre Nullmenge ist ein Großkreis mit Windung 1: **genau eine Naht ist erzwungen**, ohne
  Lage, Form oder Anzahl vorzugeben. Das ist verträglich mit Z11, weil H₂(S³) = 0 ist. Offen bleibt,
  warum das Feld primitive Ladung trägt.
- *Gibbs-Variationsprinzip.* Unter allen Zuständen mit mittlerer Energie null hat der KMS-Zustand bei
  β = L eindeutig die größte Entropie, S(ρ\*) − S(ρ) = D(ρ‖ρ\*) > 0 (gerechnet). K6 wird damit ein
  Variationsprinzip.
- *Gravitativer Zufluss.* Der Bericht findet thermisch c = 0,49997 je Paar und 7,99995 für 16 Paare;
  das bestätigt Z13 unabhängig.

**Korrekturen an früheren Aussagen.**
1. *Z9, Wand am Schnitt:* Die Phasenbarriere ist ein Artefakt des Wilson-Glieds. Mit einem
   Ginsparg–Wilson-Operator gilt die Drehidentität exakt (Rest 1,3·10⁻¹⁵), und das Phasenpotential ist
   flach (Spannweite 0 gegen 0,56 beim Wilson-Glied). Nielsen–Ninomiya erzwingt die Wand also nicht.
   Für die Majoranapaarung von H6 ist das nicht gerechnet. Laut Bericht haben Overlap-Weylfermionen bei
   endlichem Gitterabstand ein Problem mit der Reflexionspositivität (Lu–Shao); hier nicht geprüft.
2. *Z8, die 16 aus ℤ₁₆:* Im Paarmodell tragen die Partner die Ladung −1, also ν = n₁ − n₃ = 0 für jede
   Paarzahl. Die Elternanomalie erzwingt dort keine 16. Nur in der partnerfreien Fassung (16 Weyl mit
   Ladung 1, Masse aus Φ mit Ladung 2) ist ν = 16 die ℤ₁₆-Bedingung. Sonst wählt erst die
   E₈-Vervollständigung (Twistgewicht 1) die 16.
3. *Z15, Untergrenze a/ℓ_P ≥ 3,3:* Sie gilt nur für die dort gerechneten Regulatoren. Wilson-Fermionen
   geben je Feld κ = 0,221 (r = 1) bis 0,017 (r = 8). Robust bleibt nur der bosonische Sockel ohne
   Randterm, a/ℓ_P ≥ 1,65. Das Gitter liegt also oberhalb der Planckskala, aber ohne feste Zahl.
4. *Daten:* Für θ₁₃ sind Reaktorexperimente maßgeblich, nicht JUNO. Mit den Werten aus den Berichten
   (hier nicht geprüft) ergeben sich folgende Abstände:
   - sin²θ₁₂ gegen JUNO: 0,49σ;
   - sin²θ₁₃ gegen NuFIT 6.1: 1,08σ, gegen Bari: 1,69σ;
   - sin²2θ₁₃ gegen Daya Bay: 2,17σ.

   δ = 240° ist im Register ein zugewiesener Texturwert, nicht hergeleitet.

**Nicht prüfbar:** Die Lean-Zertifikate der Berichte und die dort behauptete Lean-Aussage, alle
Nullstellen von E₆ lägen im Orbit von i. Im Repo ist nur E₆(i) = 0 bewiesen.

**Was die Berichte gemeinsam schärfen.** Offen sind jetzt drei Verbindungen:
- eine Regel, die den primitiven Sektor (Ladung oder Windung 1) und genau einen windenden Block wählt;
- ein gemeinsamer Operator für schwere Nahtfelder, E₈ und drei leichte Familien, mit benannten
  Kopplungen (gerader Teil auf dem Nahttorus);
- Positivität und Gravitation bei endlichem Gitterabstand, für einen chiralen Regulator.

## Davor: H6h9 — sechs Folgeprüfungen (`z11`–`z15`, Lean-Vertrag H6EXAKT20261123)

```
        K1 … K12 K13 K14 K15 K16 K17 K18 K19 K20    P  L  F  O  E
   H6h8  P …  P   L   P   P   L   L   L   L   O    14  5  0  1  0
   H6h9  P …  P   L   P   P   L   P   L   L   L    15  5  0  0  0
```
K20 (gemeinsamer Elternoperator, der zentrale Punkt der zweiten Kritik) ist neu im Katalog. Die
sechs Prüfungen liefen parallel und wurden danach selbst nachgerechnet. Alle Skripte laufen
reproduzierbar; negative Ergebnisse stehen jeweils in `befunde` und sind nicht als Prüfung erzwungen.

**1. Einheitswindung als c₁ (`z11`, 15 Prüfungen).** TFPTs Einheitswindung ist im Original
[u_Σ] = 1 ∈ K¹(S¹) (tfpt-42, Abschnitt 006, Z. 253) bzw. c₁(L₂) = deg g₂ = 1 auf der kompaktifizierten
Normalensphäre (Abschnitt 012, Z. 27).
- **Lokal: Lesart.** Diese Klasse ist genau die Windung von det M um die Naht, also der Index. Sie
  wählt den verdrillten Defekt (Windung 1) und nicht den gewöhnlichen Wirbel (2).
- **Als Zwang zu genau einer Naht: gescheitert.** Erstens liefert TFPT keinen geschlossenen
  transversalen 2-Zyklus. Zweitens bricht das Wilson-Glied die Paarungsphase (Z9). Drittens ist
  unter G der Grad von det M gerade (M ∈ Sym²E\*), also kommen verdrillte Nähte durch eine
  geschlossene Fläche nur paarweise vor. c₁(Φ) = 1 erzwingt zwei, die Viereckenbilanz vier.
  Gitter-Ginzburg–Landau auf Torus und Kugel bestätigt das.
- Erzwungene Defekte liegen nicht quadratisch (auf der Kugel tetraedrisch). Die quadratische Lage
  kommt in H6 aus A\*, nicht aus der Abstoßung.

**2. Elternoperator (`z12`, 16 Prüfungen, 10 negative Befunde).** Gebaut ist ein BdG-Operator auf
(Querschnitt mit Defekt) × k_z × Markensphäre. Der innere Dirac ist an K^{1/2} ⊗ L gekoppelt, und der
Lift der Vierteldrehung wird aus der Invarianz von u⁴ − v⁴ berechnet.
- **Hergeleitet:** Index 3 bei Grad 3; Charaktere {−1, i, −i}, die nur bei einer freien μ₄-Bahn der
  Marken entstehen. Nahtmoden = Paare × Index des inneren Operators, alle mit Geschwindigkeit 1.
- **Mit einem einzigen Multiplett gescheitert:**
  - keine E₈, denn 48 Majorana-Weyl geben so(48);
  - die Familien würden schwer;
  - keine Vierteldrehung erfüllt beides zugleich, Familiencharaktere und J² = −1 auf der Naht.
- **Mit zwei Multipletts Lesart.** Die Naht koppelt an K^{1/2} ⊗ M mit M² = O(D) (Index 1, Nullmode
  √q, q → −q ist der Blattwechsel), die Materie an K^{1/2} ⊗ O(D) (Index 3, ungepaart). Dann gelten
  Naht, E₈ und Familien zugleich. Die Aufteilung ist eingebaut.
- Die Warnung des Atlas ist bestätigt: Auf den 248 Strömen ist J eine Involution (120/128 oder
  136/112), nie i. Geometrische Drehung, Monodromie und Markendrehung sind drei verschiedene
  ℤ₄-Wirkungen; ihre Identifikation braucht einen Verflechter.

**3. Gravitativer Zufluss (`z13`, 16 Prüfungen).** **Hergeleitet** im freien Gitter.
- Aus Casimirimpuls, Impulspolarisation und Energiestrom übereinstimmend: c = 0,4994 ± 0,0024 und
  h_R = 0,0624 je verdrilltem Paar, also c = 7,99 und h_R = 0,999 für 16 Paare.
- Eine Phasenwand im Volumen hat Chernzahl 1 je Windung, also c_Wand = 1/2 = c_Naht; der Viertelfluss
  ändert sie nicht.
- Die Wand trägt damit genau den Ausgleich von I₄ = (8/24) p₁.
- Nicht gerechnet ist der ortsaufgelöste Fluss in gekrümmtem Hintergrund.

**4. Satzzahl (`z14`, 26 Prüfungen).** Eine gemeinsame Sonde misst auf dem Gitter das Niveau n, auch
bei gekoppelten Sätzen; die SO(2) hat in E₈ den Einbettungsindex 1.
- *Gleichförmig drei Sätze* widersprechen α: Unter der Lesart EM1 (α³-Koeffizient = Niveau) gäbe
  Niveau 3 α⁻¹ = 198,24.
- *Drei Sätze, zwei verborgen* (sichtbares Niveau 1, c = 24) bestehen alle Zahlen; dagegen stehen nur
  Lesarten.
- Die Entscheidung bringt c₋/k der Naht (8 gegen 24), gleichwertig die eigene T-Phase der Naht.

**5. a/ℓ_P mit vollem Feldinhalt (`z15`, 16 Prüfungen).** **Nicht festgelegt.**
- Im H6-Schema (Wilson r = 1, schwerer Sektor, erweiterter Hilbertraum, g = 0,55) gilt a/ℓ_P ≈ 10,8.
  Über alle Schemata reicht es von 3,3 bis 14,3.
- Dominant ist der Shannon-Randterm des Eichflusses (etwa 1,8 nat je Kante und Eichboson), gefolgt von
  der Behandlung der Verdoppler.
- Robust ist nur a/ℓ_P ≥ 3,3: Das Gitter liegt nicht bei der Planckskala.

**6. Lean-Vertrag (`experiments/theory-contracts/h6-exakte-teile-20261123/`).** 87 Theoreme
kompilieren mit Mathlib, ohne `sorry`, ohne neue Axiome und ohne `native_decide`.
- Viertelsektoren: g̃⁴ = −1, char. Polynom (X⁴+1)², jeder Sektor genau ein Viertel jedes Eigenraums.
- Familiencharaktere und die Teilung 3+2 des Trägers.
- Arithmetik: c ∈ 8ℤ ⇔ 16 | ν, Twistgewicht 1 ⇔ N = 16, 48 | N.
- Die X-Ladungen des Standardmodells.
- E₆(i) = 0.

Offen bleibt die Eindeutigkeit der Nullstelle, weil Mathlib die Valenzformel fehlt. Das Ledger ist
unverändert.

**Prämissen von H6h9:** A\*, C″, D, E, G.

**Was sich verschoben hat:**
- Gravitativer Zufluss und Twistgewicht sind jetzt gerechnet.
- Die Naht als genau eine erzwungene Naht ist widerlegt; Nähte kommen paarweise.
- Der Elternoperator gelingt nur mit zwei Multipletts.
- Die Satzzahl hängt an einer einzigen Größe, c₋/k.

**Offen:**
- Warum die Aufteilung in Naht- und Materiemultiplett.
- Der Grundzustand.
- c₋/k.
- Der Fluss in Krümmung.
- Die Valenzformel in Lean.

## Davor: H6h8 — Zufluss, Grundzustand, zweite Kritik (`z9_zufluss_grundzustand.py`, `z10_kritik_zwei.py`)

Bilanz unverändert: 14 P und 5 L von 19, nichts gescheitert, offen oder gesetzt. Geändert haben sich
die Grundlagen von K9, K10, K17 und K19, und drei frühere Aussagen werden zurückgenommen. Alle
Prüfungen in Z9 (6) und Z10 (8) bestehen.

**Zufluss (Z9, Teil A).** Zwei Komponenten des reellen 16-Pletts bilden eine exakte SO(2)-Ladung
([H, Q] = 0). Im Sektor Q = +1 ist eine Verschiebung von k_z ein SO(2)-Eichfeld entlang der Naht. Die
Zahl der besetzten Zustände bleibt konstant (576).
- Beim Nulldurchgang verliert der Kern genau eine geladene Mode (Sprung −0,900 bei Kerngewicht
  0,922, Verhältnis 0,976).
- Der äußere Rand gewinnt sie (+0,930).
- Über die übrige Periode gleicht das lückenhafte Volumen stetig aus (+0,900).

Das ist der Zufluss für die SO(2)-Ladung. Der gravitative Zufluss (c = 8) ist nicht gerechnet.

**Explizite Brechung und Wand (Z9, Teil B, neu).** Das Wilson-Glied, das Nielsen–Ninomiya erzwingt,
bricht die U(1) der Paarungsphase. Die Energiedichte hängt von der relativen Phase ab:
- Periode π in der Massenphase (J-Symmetrie).
- Minimum bei φ₀ = π, Barriere 0,088 je Gitterzelle.

Der Defekt dreht die Phase einmal herum und muss über die Barriere, also **berandet er eine
Domänenwand. Die Wand ist der Verzweigungsschnitt.** Korrektur zu H6h5: Der Schnitt ist nur
topologisch unsichtbar (Anomaliephase, c ∈ 8ℤ), energetisch ist er eine Wand mit Spannung auf
Gitterskala. Wie bei Axionen mit Domänenwandzahl 1 enden Wände auf Nähten, und Nähte beranden Wände.
Der frühere Vergleich mit der Phase 0 (negative Defektenergie, mit der Fläche wachsend) war ein
Vergleich mit dem Maximum. Die Kritik hat dieselben Zahlen unabhängig erhalten und gegen die Phase π
positive Werte gefunden.

**Grundzustand.** Im freien Modell erzwingt nichts einen Defekt: Er kostet Kern- und Wandenergie, und
der homogene Zustand ist erlaubt. Es gibt zwei Wege.
- Topologisch: Ein nichttriviales Bündel des Massenfelds (c₁ ≠ 0) erzwingt eine Nullstellenfläche
  [Σ] = PD(c₁). Was c₁ ≠ 0 erzwingt, ist offen.
- Als Fluktuation: Die euklidische Summe über Defektflächen gewichtet mit
  exp(−σ_W V − σ_N A)·χ_E₈(i). Flächen mit Wirkung unter log 12 sind nicht unterdrückt, also
  Nähte auf Gitterskala.

Die Kritik meldet ein Zusatzmodell auf S³ mit erzwungener Nulllinie oberhalb von δ_c = 2 − √(8/3) ≈
0,367. Die Arithmetik ist geprüft; das Modell liegt nicht vor.

**Zweite Kritik, Punkt für Punkt (Z10).**
1. *J folgt nicht aus der Modenzahl.* **Richtig.** Auch ohne J trägt der Defekt (1, 0) genau eine Mode
   mit v = 1. Der Index ist die Windung von det M an der physischen Ecke, n₊ + n₋. Die Verdoppler
   tragen wegen W ≥ 2r_w > m nichts bei, und ohne Wilson-Glied heben sie sich auf. J leistet etwas
   anderes: Mit J genügt eine einzige komplexe Masse, deren Phase um π windet. Die Masse ist dann
   zweiwertig (Φ → −Φ); das sind die zwei Blätter auf Bosonischem, und auf Fermionen wirkt μ₄. J wird
   durch die Blattstruktur von P1 gewählt, nicht durch die Zählung. Die Klassifikation über ein
   Energiefenster ist durch die Windungsinvariante ersetzt.
2. *Die Familiencharaktere entscheiden die Satzzahl nicht.* **Richtig.** Bei drei identischen Sätzen
   ist der Unterraum „alle tragen dasselbe Label“ invariant und hat dieselben Eigenwerte {−1, i, −i}.
   Die Behauptung aus Z8 ist zurückgenommen. n = 1 stützt sich nur noch auf Minimalität mit Zufluss und
   auf TFPTs ein E₈. Unterscheiden könnte das Niveau eines gemeinsamen E₈ Stroms (1 gegen 3).
3. *Familien als chiraler Index der vier Marken.* **Bestanden und übernommen.** Für D = μ₄ auf CP¹
   gilt H⁰(K(D)) = H⁰(O(2)) mit Dimension 3 und h¹ = 0. Spin-Dirac an K^{1/2} ⊗ O(D) = O(3) hat
   ker D₊ = 3 und ker D₋ = 0. Die Vierteldrehung z → iz wirkt auf ω = P(z) dz/(z⁴ − 1) als
   P(z) ↦ i P(iz), mit Eigenwerten **genau {−1, i, −i}**. Dieselben Marken geben H⁰(O(D)) mit
   Dimension 5 (TFPTs Fünferträger), und die Halbdrehung teilt ihn exakt **3 + 2**.
   Λ^gerade ℂ⁵ = 1 + 10 + 5 = 16, also Materie = H⁰(K(D)) ⊗ Λ^gerade H⁰(O(D)) mit 48 Komponenten. Das
   ersetzt C′ (Spinstrukturen) durch **C″: Familien sind der chirale Index der vier Marken**. Er ist
   topologisch geschützt und kommt ohne Wahl des Drehzentrums aus. Offen ist, warum die physischen
   Fermionen an dieses Bündel koppeln.
4. *Zwei verschiedene 16.* **Richtig und jetzt getrennt.** Die 16 der Naht sind ein schweres reelles
   16-Plett (Vektor von SO(16)); die 16 einer Generation sind Λ^gerade des Fünferträgers (komplexer
   Spin(10)-Spinor).
5. *Twistgewicht h_R = N/16.* **Bestanden.** Nur N = 16 gibt Ramond-Grundzustände vom Gewicht 1, also
   128 Spinorströme und 120 + 128 = 248. Das wählt Blöcke aus 16 und schließt eine gemeinsame 48er-Summe
   aus (h_R = 3). Zwischen einem und drei Blöcken entscheidet es nicht.
6. *Stabilisierung einer Schleife durch (h − c/24)/R.* Im Ramond-Grundzustand ist der Koeffizient +2/3.
   Im KMS-Zustand von H6 ist er null (Z4: ⟨L₀⟩ = c/24). Die Stabilisierung gilt also nur für eine
   andere Zustandswahl.
7. *Raum mal periodische Zeit.* Der Torus setzt euklidisch periodische Zeit voraus. In H6 ist das der
   KMS-Zustand (K6); das ist eine Prämisse, keine versteckte Annahme.

**Prämissen von H6h8:** A\*, C″, D, E, G.

**Offen.**
- Was ein nichttriviales Massenbündel erzwingt, oder eine Rechnung zur Defektsumme.
- Der gravitative Zufluss aus der Fermiondeterminante.
- Die Kopplung der physischen Fermionen an K^{1/2} ⊗ O(D).
- Die Satzzahl über das Niveau eines gemeinsamen E₈ Stroms.
- Ein gemeinsamer Elternoperator für Naht, E₈ und Materie.

## Davor: H6h7 — die 16 und die Zahl der Sätze (`z8_sechzehn_und_saetze.py`)

Die Bilanz bleibt bei 14 P und 5 L von 19. Es gibt aber eine Prämisse weniger: F (E₈ Phase auf dem
Schnitt) ist jetzt eine Lesart innerhalb von G. Die Prämissen sind **A\*, C′, D, E, G**. Alle sechs
Prüfungen bestehen; zwei davon sind reine Arithmetik (c = ν/2 und die Minimalität).

**Was gezählt wird.** Im Paarmodell aus Z7 windet nur ξ₊ = (ψ − iχ)/√2. ξ₊ trägt unter J die Phase i,
ist also ein Weyl-Fermion mit Spin-ℤ₄-Ladung 1 (J² = (−1)^F). Der Partner ξ₋ (Ladung −1) windet
nicht; das Gitter braucht ihn nur, damit ein J-invariantes Wilson-Glied die Verdoppler entfernt. ν
Fermionen mit Ladung 1 geben c = ν/2. Der Schnitt ist genau bei c ∈ 8ℤ unsichtbar (Z5), also bei
ν ∈ 16ℤ. **Das ist dieselbe Bedingung wie die ℤ₁₆-Anomaliefreiheit von Spin-ℤ₄** (Garcia-Etxebarria und
Montero 2019).

**Welche 16?**
- *Die Standardmodell-Generation hat die richtige Zahl.* Mit X = 5(B−L) − 4Y gilt X ≡ 1 mod 4 für alle 16
  Weyl-Fermionen einschließlich ν^c. Ohne ν^c wären es 15, also anomal.
- *Sie bildet die Naht trotzdem nicht.* In einer Generation ist der einzige eichinvariante
  Majorana-Bilinear ν^c ν^c. Nur ν^c kann am Defekt eine windende Masse tragen: eine Mode je
  Generation, drei Generationen geben c = 3/2. Die leichten SM-Fermionen scheiden auf der Ebene von
  Bilinearen aus (Wechselwirkungen, etwa symmetrische Massenerzeugung, sind nicht geprüft).
- *Ein reelles 16-Plett funktioniert.* Für den Vektor von SO(16) ist die Masse Φ ψ_a ψ_a mit
  komplexem Singulett Φ eichinvariant. Φ trägt ℤ₄-Ladung 2, und **sein Vorzeichen ist das Blatt**.
  Alle 16 winden, also c = 8 und ν = 16, anomaliefrei. Unter SO(16) ⊃ SO(10) × SO(6) ist 16 =
  (10,1) + (1,6). Das sind genau die Quantenzahlen der 16 Nahtmoden in E₈ ⊃ SO(16): Die 4D Eltern
  tragen dieselben Quantenzahlen wie die Nahtfermionen. Sie sind schwer und vektorartig, ohne
  Eichanomalie.
- *Wo die Generationen sitzen.* E₈ = 120 + 128 unter SO(16), und 128 = (16,4) + (16̄,4̄). Mit SU(4) ⊃
  SU(3) × U(1) wird (16,4) zu (16,3) + (16,1). Die SM-Generationen sind Zustände im Spinorteil des
  Nahtkerns, nicht seine 16 Fermionen. Das passt zu TFPTs Familien in einem E₈.

**Zahl der Sätze: n = 1.**

| n | c | ν | Naht allein T-invariant | mit Zufluss | Familiencharaktere aus der Naht |
|---|---|---|---|---|---|
| 1 | 8 | 16 | nein (e^{−2πi/3}) | ja | 3 Sektoren, {−1, i, −i} |
| 2 | 16 | 32 | nein | ja | 6 Sektoren |
| 3 | 24 | 48 | ja | ja | 9 Sektoren |

Unter G ist die Naht nie allein. Die schweren 4D Fermionen, deren Masse am Defekt windet, liefern den
Zufluss (Callan–Harvey), deshalb wird strenge T-Invarianz der Naht für sich nicht verlangt. Das
kleinste n mit unsichtbarem Schnitt ist dann 1. Nur bei n = 1 trägt die Naht genau drei
Spinstruktursektoren mit {−1, i, −i}. n = 3 wäre die kleinste Zahl ohne Zufluss, verlöre aber die
Familiencharaktere und passt nicht zu TFPTs einem E₈.

**Nicht gezeigt.**
- Der Callan–Harvey-Zufluss ist hier nicht gerechnet.
- Warum es genau ein schweres reelles 16-Plett gibt.
- Warum der Grundzustand einen geschlossenen Defekt enthält.
- Ob Wechselwirkungen die SM-Generation selbst zur Quelle machen könnten.

## Davor: H6h6 — Ursprung der Naht (`z7_ursprung_naht.py`)

```
        K1 … K15 K16 K17 K18 K19    P  L  F  O  E
   H6h5  P … P   L   L   L   E    14  4  0  0  1
   H6h6  P … P   L   L   L   L    14  5  0  0  0
```
K19 (Existenz der Naht) ist neu im Katalog; bei allen früheren Hypothesen war die Naht Eingabe. In
H6h6 ist nichts mehr als Eingabe gesetzt; alles steht unter benannten Prämissen.

**Hypothese G.** Die Naht ist der verdrillte Defekt eines komplexen Paarungsmassenfelds
m e^{iθ} ψχ zwischen Weyl-Fermionen ψ und Partnern χ. Um den Defekt werden ψ und χ mit der
Vierteldrehung J: ψ → χ → −ψ vertauscht, und die Massenphase windet um π. Eindeutig ist das, weil J die
Masse ψχ umpolt.

**Gerechnet (alle sieben Prüfungen bestehen).** 3D Gitter, Defektlinie entlang z, Bogoliubov-de-Gennes
mit Weyl-Kinetik und unter J invarianter Wilson-Majoranamasse gegen die Verdoppler. Kern- und
Randmoden sind über den Kernprojektor getrennt, die Chiralität kommt aus ⟨∂H/∂k_z⟩.

| Defekt | chirale Majorana-Weyl Moden am Kern |
|---|---|
| verdrillt (J, Phase windet um π), ein Paar | **1**, Geschwindigkeit +1,00 |
| verdrillt, zwei Paare | **2** |
| gewöhnlicher Wirbel (ohne J, Windung 1), ein Paar | 2 |
| nur J, ohne Massenwindung | 0 |

Der äußere Rand trägt jeweils dieselbe Zahl mit umgekehrter Richtung. Das Ergebnis ist robust bei
L = 14 und 18 und für Massen von 0,4 bis 0,9.

**Was daraus für die Naht folgt.**
- Kodimension 2: Eine Defektlinie im Raum ist eine Fläche in der Raumzeit.
- Geschlossen und flach heißt Torus (Gauß-Bonnet, χ = 0).
- J² = −1 auf Fermionen, also wirkt J auf allem Bosonischen mit Ordnung 2. Das sind die **zwei
  Blätter**; auf Fermionen ist es **μ₄**, wie die gelifte Vierteldrehung (g̃⁴ = −1, Z1).
- Der Verzweigungsschnitt ist der Ort, an dem die halbe Phase springt. Seine Lage ist eine Eichwahl;
  das ist der unbeobachtbare Schnitt aus F (Z5).
- Die chiralen Moden laufen mit Lichtgeschwindigkeit. Ihre Zahl ist die Zahl der Paare, deren Masse am
  Defekt windet; **16 Paare geben 16 Majorana-Weyl Moden, also c = 8**.
- Gegenprobe: Die 16 Verdoppler eines naiven Feldes liefern die Moden nicht. Ohne Wilson-Glied haben
  sie abwechselnde Chiralität, und der verdrillte Defekt trägt dann keine Mode. Gebraucht werden 16
  echte Paare.

**Korrektur zur Gabel (Z5).** Z5 hat nur die Operatorlesart (H|ψ⟩ = 0) ausgeschlossen. c = 24 mit der
Energiebedingung im Mittel ist damit nicht ausgeschlossen. Z3 und Z4 gelten für E₈³ genauso: Dort ist
nur i stationär, und die Energie ist null. Die Wahl c = 8 stützt sich auf TFPTs Familien in einem E₈
und auf C′ (Familien als Spinstruktursektoren). Unter G ist c = 8 × (Zahl der 16er-Sätze, deren Masse
windet).
- Windet ein Satz, ist c = 8 und F nötig.
- Winden drei, ist c = 24. Dann verschwindet die Rahmungsphase ohne F, und drei ist die kleinste Zahl
  mit 8N ∈ 24ℤ. Die Familiencharaktere {−1, i, −i} kämen dann aber nicht mehr aus der Naht, denn J
  wirkt innerhalb der Paare und vertauscht keine Sätze.

Die Gabel ist damit eine Frage nach der Zahl der windenden 16er-Sätze.

**Offen.**
- Die Identität der 16 Paare. Sie sind im Volumen schwer (Masse m) und daher nicht die leichten
  Fermionen des Standardmodells.
- Warum genau ein (oder drei) Satz windet.
- Warum der Grundzustand einen geschlossenen Defekt enthält.

## Davor: H6h5 — Gabel entschieden, Kontinuum (`z5_gabel.py`, `z6_kontinuum.py`)

```
        K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 K11 K12 K13 K14 K15 K16 K17 K18    P  L  F  O  E
   H6h4  P  P  P  P  P  P  P  P  P  P   P   P   L   P   P   O   L   L    14  3  0  1  0
   H6h5  P  P  P  P  P  P  P  P  P  P   P   P   L   P   P   L   L   L    14  4  0  0  0
```
Erstmals ist kein Kriterium mehr offen. Vier stehen unter benannter Lesart: Gravitation (D),
Kontinuum (physisches Gitter, perturbativ), Anomalie (F) und Neutrinomischung (TFPT-Kanäle).

**Gabel c = 8 gegen c = 24: entschieden für c = 8 mit F.** Alle fünf Prüfungen bestehen.
1. Die Operatorlesart H|ψ⟩ = 0 widerspricht H6 selbst. Der KMS-Zustand bei β = L hat Var(L₀) = 0,2426,
   ist also kein Eigenzustand. Unter H|ψ⟩ = 0 wäre e^{−βH} auf allen physischen Zuständen eins, und β,
   die Kragenbreite L/2 und die Zeitneutralität aus Z3 und Z4 wären unbestimmt (K6 fiele).
2. Ein E₈ trägt genau drei Spinstruktursektoren, die Terme θ₂⁸, θ₃⁸, θ₄⁸ der GSO-Summe. S hält (A,A)
   fest und tauscht AP mit PA; T tauscht (A,A) mit PA, mit der gemeinsamen Phase e^{−2πi/3}. Bis auf
   diese Phase wirkt SL(2,ℤ) als S₃ = SL(2, ℤ/2) auf genau drei Objekten, passend zu C′. E₈³ hätte neun
   Sektoren. TFPT bettet die Familien als SU(3)-Triplett in ein E₈ ein: 248 = 78 + 8 + 2·27·3.
3. **Neu: Die 3D Wand ist genau bei c ∈ 8ℤ unsichtbar.** Verschiebt man sie über ein 4D Gebiet, ändert
   sich die Phase um exp(2πi c σ/8) (Anomaliepolynom (c/24)p₁, ∫p₁ = 3σ). Weil die GSO-Summe die Naht
   spinstrukturfrei macht, zählen alle orientierten 4-Mannigfaltigkeiten, und c = 8 ist der kleinste
   Wert mit trivialer Phase. Das ist dieselbe 16 wie in Z1. Lesart: Die zwei Blätter von P1 sind
   entlang eines 3D Schnitts verklebt, der an der Naht endet. Trägt der Schnitt die invertierbare E₈
   Phase, ist seine Lage unbeobachtbar wie ein Dirac-String. F ist dann keine zusätzliche Struktur,
   sondern der Verzweigungsschnitt der zwei Blätter.

**Kontinuum (K16).** Das Gitter ist physisch (Z2: a ≥ 3 ℓ_P). Verlangt ist daher kein Limes a → 0,
sondern Universalität im Infraroten. Alle acht Prüfungen bestehen.
- **Lorentz ohne Feinabstimmung.** Auf dem isotropen Gitter ist die Lichtgeschwindigkeit für jede Art
  und jede Richtung exakt eins. Die Abweichung ist ~ p² (Dimension 6) und richtungsabhängig, zum
  Beispiel −1/3 entlang der Achse und −2/9 entlang der Raumdiagonale für das naive Fermion. Die
  Ein-Schleifen-Vakuumpolarisation ist entlang Zeit- und Raumachse gleich (−0,252432). Entlang der
  Diagonale weicht sie um 5·10⁻⁴ ab und skaliert mit q² (Verhältnis 0,260 gegen 0,25).
- **Gegenprobe mit ausgezeichneter Zeit.** Mit einem eigenen Zeitsprung ν je Art bekommt jede Art ihre
  eigene Geschwindigkeit, und das induzierte Photon liegt dazwischen: Zwei Arten mit ν = 1 und 0,9
  geben c_γ² = 1,122 gegen 1 und 1,235. Das ist Lorentzverletzung bei Dimension 4, nur durch
  Feinabstimmung zu beheben. Schon eine einzelne Art bekommt einen nicht-universellen Gitteranteil von
  etwa 2 %. **Die Gleichwertigkeit der vier Achsen, die die Kritik als Hindernis für Lorentz las, ist
  genau die Bedingung für Lorentz ohne Feinabstimmung.** Die Vierteldrehung zwischen Zeit- und
  Raumachse verbietet jedes ν ≠ 1 in jeder Ordnung.
- **Reflexionspositivität.** Naive Fermionen haben die Transfereigenwerte ±e^{−E}. T ist nicht positiv,
  T² schon (G(t+2)/G(t) = e^{−2E} exakt). Die OS-Rekonstruktion läuft über Doppelschritte; der
  Zeitverdoppler ist eine eigene Art.
- **Größenordnungen.** Quadratische Lorentzverletzung setzt bei etwa 2–7·10¹⁸ GeV ein; bei 100 TeV ist
  δv ~ 10⁻²⁷. Die Beobachtungsschranken liegen viele Größenordnungen darunter. Der Landau-Pol von
  U(1)_Y liegt bei etwa 10⁴¹ GeV, weit über der Gitterskala. Weil das Gitter physisch ist, ist die
  Trivialität von U(1) kein Hindernis; die nichtabelschen Faktoren sind asymptotisch frei.

**Nicht gezeigt.**
- Nichtperturbative Universalität (Confinement, Higgs-Sektor).
- Die chirale Gitterformulierung. Naive Fermionen haben 16 Geschmäcker in Paaren entgegengesetzter
  Chiralität (Nielsen–Ninomiya). Ein Weg aus der Literatur, hier nicht gerechnet: die Spiegelpartner
  symmetrisch massiv machen. Das geht nur, wenn die ℤ₁₆-Anomalie einer Generation verschwindet, also
  mit 16 Weyl-Fermionen je Generation einschließlich ν_R.

**Prämissen von H6h5:** A\*, C′, D, E und F. F ist jetzt an K1 gebunden: Es ist die E₈ Phase auf dem
Verzweigungsschnitt.

## Davor: Stabilität der Lage τ = i (`z4_stabilitaet.py`)

**Ergebnis: Die Lage ist stabil. Der Sattel aus Z3 ist keine Instabilität**, sondern die gewöhnliche
Struktur zweier konjugierter Paare. Alle neun Prüfungen bestehen.

In Rekonstruktion A (Raum L, Zeit β) ist τ = x + iy. Die beiden Richtungen des Sattels bedeuten
Verschiedenes.
- **Verdrillung x** koppelt an den Impuls: Z = Tr e^{−βH} e^{iθP}. Bei positivem Spektrum gilt
  |Z(x + iy)| ≤ Z(iy), mit Gleichheit nur bei ganzem x. Verdrillung null ist für jede unitäre Naht ein
  echtes Maximum von |Z|. Für E₈ gilt |Z(x+i)|/Z(i) = 0,953 bei x = 0,1 und 0,373 bei x = ½, also
  F_xx = +9,58: **stabil**. Die „eingefrorene Verdrillung“, die frühere Befunde als Zusatzannahme
  brauchten, ist damit automatisch.
- **Streckung y = β/L** ist die inverse Temperatur. log Z ist in β konvex, denn die zweite Ableitung
  ist die Energievarianz. Für E₈ gilt bei y = 1 exakt F_yy = −(2π)²·Var(L₀) = −9,58, mit
  Var(L₀) = 0,2426. Längs der Temperatur ist −log Z daher **für jede unitäre Naht** konkav. Das ist
  positive Wärmekapazität, keine Instabilität. Ein Minimum über die Temperatur gibt es grundsätzlich
  nicht; die Temperatur ist eine Legendre-Variable.
- **Mikrokanonisch:** Eine hintergrundfreie, geschlossene Naht hat kein äußeres Wärmebad. Ihre
  Temperatur wird durch die Energiebedingung festgelegt, nicht gewählt. S(E) = min_β[log Z + βE] ist
  konkav (stabil). Bei E = 0 liegt genau β = L, also τ = i, mit S(0) = log χ(i) = log 12. Dort gilt
  ⟨L₀⟩ = c/24 = 1/3 exakt: Casimir- und Anregungsenergie heben sich auf.

Die Lage τ = i ist damit der stabile Gleichgewichtszustand einer geschlossenen Naht bei Energie null
und Impuls null: stationär (Z3) und stabil (Z4).

**Neu sichtbar: die Gabel F gegen F′ als Frage nach der Energiebedingung.** Bei c = 8 liegen die
Niveaus bei n − 1/3; einen Zustand mit Energie genau null gibt es nicht. E = 0 gilt also nur im
Mittel des thermischen Zustands. Soll die Bedingung als Operatorgleichung gelten (H|ψ⟩ = 0,
Wheeler-DeWitt-artig), muss c durch 24 teilbar sein. Bei E₈³ (Z = j) hat das Nullniveau 744 = 3 × 248
Zustände. Das ist dieselbe Bedingung wie die T-Invarianz aus Z1. Die Lesart im Mittel passt zu TFPT,
das ein E₈ mit 248 = 78 + 8 + 2·27·3 nutzt (SEAM.EQUIV.02). Entscheiden muss das eine eigene Prüfung.

**Prüfstand.** Die Bilanz bleibt 14 P, 3 L, 0 F, 1 O. K2 trägt jetzt „stationär und stabil“.

**Offen.**
- K16: wechselwirkender Kontinuumslimes, Reflexionspositivität und SO(4).
- a/ℓ_P mit vollem Feldinhalt.
- Die Energiebedingung im Mittel (F, c = 8) oder als Operator (F′, c = 24).
- Warum es überhaupt eine geschlossene 2D Naht mit chiralen Fermionen gibt.

## Davor: H6h4 — Prämisse A hergeleitet (`z3_hintergrundfreiheit.py`)

```
        K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 K11 K12 K13 K14 K15 K16 K17 K18    P  L  F  O  E
   H6h3  P  P  P  P  P  P  P  P  P  P   P   P   L   P   P   O   L   L    14  3  0  1  0
   H6h4  P  P  P  P  P  P  P  P  P  P   P   P   L   P   P   O   L   L    14  3  0  1  0
```
Die Bilanz ist gleich; geändert hat sich die Grundlage. Die setzende Prämisse A („der orientierte
Tausch der Zyklen ist eine Symmetrie“) ist ersetzt durch:

**A\*: Die Naht hat keine Hintergrundgeometrie.** Ihre Basiswechsel (große Diffeomorphismen) sind
Symmetrien, und ihre Form ist nicht vorgegeben, sondern erfüllt ihre Bewegungsgleichung: Die Antwort
des Zustands auf jede Deformation der komplexen Struktur verschwindet, ⟨T⟩ = 0. Das ist die
Bedingung einer dynamischen Nahtmetrik, wie die Virasoro-Bedingung.

**Was daraus folgt (alle sieben Prüfungen bestanden).**
- Für die E₈-Naht gilt exakt d log χ/dτ = −(2πi/3)·E₆/E₄ (geprüft auf 10⁻²¹). ⟨T⟩ = 0 heißt also
  E₆(τ) = 0. Nach der Valenzformel hat E₆ im ganzen Fundamentalbereich genau eine Nullstelle, die
  halbe am elliptischen Punkt i. **Die Form der Naht ist damit eindeutig τ = i**, ohne
  Orientierungsvoraussetzung und ohne eingefrorene Verdrillung. Am Sechseck ist χ = 0; das ist keine
  zulässige Form.
- Dort ist die Gesamtenergie in beiden Rekonstruktionen exakt null: Die thermische Energie hebt die
  Casimirenergie auf. Bei β/L = 0,8 und 1,25 vertauschen sich die Energien (+2,40 und −1,54). Beide
  Rekonstruktionen sind am stationären Punkt dasselbe System (β = L). **Zeitneutralität ist damit
  eine Folge, keine Prämisse.** Der Stabilisator von i ist μ₄; daraus folgen wie bisher die
  zyklischen Marken, die Viertelzählung und c₃.
- **Verbot (neu):** Für eine chirale Naht ist log|Z| harmonisch (Hesse +9,58 / −9,58, Spur 0). Es gibt
  daher nirgends ein lokales Minimum. Energieminimierung kann die Form grundsätzlich nicht wählen; nur
  Stationarität kann es. Das erklärt, warum frühere Minimierungsansätze (KERNSIX, KERNSYM) bei τ = i
  einen Sattel fanden.
- **Gegenproben.** Bei der Dirac-Naht mit fester Spinstruktur (A,A) ist i im Modulraum der
  Thetagruppe ebenfalls die einzige stationäre Form. Bei E₈³ (Z = j) ist es nur i. Beim freien Boson
  und beim Leech-Gitter ist zusätzlich das Sechseck ρ stationär; dort entscheidet erst die
  Spinstruktur (Z1: am Sechseck bleibt nur PP mit Nullmoden). Die Auswahl ist also kein Artefakt
  jeder beliebigen Naht.

**Was davon schon bekannt war.** Dass die Energie des E₈-Charakters genau bei τ = i verschwindet und
dass E₆ dort seine einzige Nullstelle hat, steht im Atlas (Befund vom 18.11., „c₃ ist die
Gitterladungsenergie …“). Ebenso ist dort festgehalten, dass τ = i ein Sattel ist. Neu ist:
- diese Nullstelle als Bewegungsgleichung einer hintergrundfreien Naht zu lesen, die den gesetzten
  Fixpunkt ersetzt;
- das Verbot für Minimierung über die Harmonizität;
- die Gegenproben.

**Prämissen von H6h4.**
- **A\*** (hintergrundfreie Naht) ersetzt A.
- C′, D, E, F bleiben wie in H6h3.
- A\* trägt zweimal: Die großen Diffeomorphismen geben die GSO-Summe und damit E₈ mit T-Phase (daher
  F). Die Bewegungsgleichung gibt τ = i.

**Was A\* nicht leistet.**
- (Nachtrag: Z4 zeigt, dass der Sattel keine Instabilität ist; siehe oben.) Die Lage ist stationär, aber nicht stabil: Sattel, instabil gegen Streckung. Wenn die Form wirklich
  dynamisch ist, braucht die Stabilität einen weiteren Grund. Kandidaten sind die 3D E₈ Phase (F),
  die Ankopplung an die 4D Geometrie oder das Modulintegral, in dem ein Sattel der maßgebliche
  Punkt ist.
- Dass es überhaupt eine geschlossene 2D Naht mit chiralen Fermionen gibt, bleibt Ausgangspunkt.

## Davor: H6h3 — externe Kritik geprüft (`z1_kritik_pruefung.py`, `z2_verschraenkungsdichte.py`)

```
        K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 K11 K12 K13 K14 K15 K16 K17 K18    P  L  F  O  E
   H6h2  P  P  P  P  P  P  P  P  P  P   P   P   L   P   P   O   O   L    14  2  0  2  0
   H6h3  P  P  P  P  P  P  P  P  P  P   P   P   L   P   P   O   L   L    14  3  0  1  0
```
K17 (Rahmungsanomalie der Naht) und K18 (Neutrinomischung gegen Daten) sind neu im Katalog, damit
die Punkte der Kritik sichtbar bleiben. H6h2 wird rückwirkend ehrlich bewertet: K17 offen.

**Die Kritik Punkt für Punkt.** Alle 16 Rechenprüfungen in `z1` bestehen.
1. *Kovarianz wählt kein Quadrat.* Richtig für die Zustandssumme: Z_AA(−1/τ) = Z_AA(τ) gilt auch bei
   τ = 2i. Prinzip A verlangt aber mehr, nämlich dasselbe thermische System in beiden
   Rekonstruktionen. Bei τ = 2i sind es verschiedene Systeme (β/L = 2 gegen 1/2, Besetzung der
   tiefsten Mode 0,0019 gegen 0,172). A ist also die Fixpunktbedingung. Die von der Kritik
   vorgeschlagene Verschärfung (echte Drehung der Ordnung größer zwei, Spinstruktur ohne Nullmoden
   erhalten) ist dieselbe Bedingung. Drehungen der Ordnung > 2 gibt es nur bei i und ρ; bei ρ
   bleibt nur PP erhalten, das Nullmoden hat; also τ = i mit (A,A). Das deckt sich mit dem
   Atlasbefund vom 16.11. (nur j = 1728 erlaubt Ordnung vier auf den Marken). Offen bleibt, warum
   die Natur den Fixpunkt realisiert: Das ist Prämisse A selbst.
2. *Spinorwirkung der Vierteldrehung.* Gerechnet mit der gelifteten Drehung g̃ um den Viertelpunkt:
   g̃ vertauscht mit dem Dirac, g̃⁴ = −1, es gibt keinen invarianten Sektor, und der naive Projektor
   (1/4)Σg̃ᵏ ist keiner. Die Sektoren sind die vier Eigenwerte exp(iπ(2m+1)/4). Weil kein Impuls in
   (ℤ+½)² fix ist, trägt jeder Sektor **exakt ein Viertel jedes Eigenraums**, ohne
   Fixpunktkorrektur. Damit ist c₃ = 1/(8π) je Spinorsektor (Wärmespur 0,039788735773). Die
   Viertelzählung überlebt die fermionische Drehphase.
3. *Woher das Familienvorzeichen?* Neu und exakt: Die drei geraden Spinstrukturen der Naht sind die
   nichttrivialen Charaktere der Markengruppe E[2]. Die zyklische Markendrehung g(z) = iz + ½ ist
   affin. Ihr Verschiebungsteil wirkt auf diesen Charakteren über die Weilpaarung als Vorzeichen,
   und das Spektrum ist **genau {−1, i, −i}**: (A,A) mit −1, AP und PA mit ±i. Kein Vorzeichen ist
   eingesetzt. Eine Drehung um eine Marke statt um den Viertelpunkt gäbe {1, 1, −1}. Dasselbe g auf
   H¹ der vier Marken gibt dasselbe Spektrum. Für alle 24 Elemente von S₄ = ⟨R, C⟩ = Aut(Q₈) gilt
   Charakter auf Im H = Signum × Standard = Λ²(Standard). Die chiralen 2-Formen der Raumzeit sind
   also das Λ² des Familienraums, und das Signum ist dessen Determinante. Das ist genau der
   orientierte Viermarkenzustand, den die Kritik vermutet. Prämisse C wird damit zu
   **C′: Familien = gerade Spinstrukturen der Naht, chirale 2-Formen = ihr Λ²**.
4. *E₈ und die T-Phase.* Richtig: Die Vierteldrehung braucht nur S, und S ist für jedes N exakt.
   T-Kovarianz der GSO-Summe verlangt N ≡ 0 mod 16, strenge T-Invarianz N ≡ 0 mod 48. E₈ (N = 16)
   ist daher nur bis auf exp(−2πi/3) basisunabhängig. Es gibt zwei konsistente Ergänzungen:
   **(F)** eine invertierbare 3D E₈ Phase mit Rand c = 8. Ihr T²-Hilbertraum ist eindimensional,
   also bricht die Wahl der Füllung die Vierteldrehung nicht. Oder (F′) drei identische E₈ Blöcke,
   c = 24, mit Zustandssumme χ_E8³ = j (χ_E8(i) = 12, j(i) = 1728). H6h3 nimmt F als neue Prämisse;
   F′ würde TFPTs c = 8 in der Gravitation ändern. Nebenbefund: arg exp(−2πi/3) = 240° ist TFPTs
   δ_PMNS. Die Zahl ist gleich; ein Zusammenhang ist nicht hergeleitet.
5. *Lorentzsignatur.* Richtig: Die Gleichwertigkeit der vier Achsen ist eine euklidische Aussage.
   H6 nutzt sie nur so. Die lorentzsche Zeit entsteht durch OS-Rekonstruktion entlang der Achse, die
   der Zustand wählt. Dafür braucht es Reflexionspositivität der Gitterwirkung und die
   Wiederherstellung von SO(4) im Limes. Beides gehört zu K16 und ist nicht gezeigt.
6. *Eichdynamik.* Verbindung, Krümmung und induzierte Wirkung sind in `y1` gerechnet. Die gemeinsame
   Normierung kommt aus denselben Fermionschleifen bei vollständigen Multipletts. Schwellen durch
   ungleiche Massen sind nicht behandelt.
7. *Gravitation braucht η mikroskopisch (Gegentest D, `z2`).* Gerechnet: Ein naives Diracfeld auf
   ℤ³ trägt je Gitterfläche κ = 0,78657 nat (konvergiert in Ring und Querimpulsen; nicht π/4),
   also 0,0983 je Geschmack. η = 1/4 entspricht 2,54 Dirac-Geschmäckern. Drei Generationen mit ν_R
   geben schon ohne Eichfelder κ ≥ 2,36. Unter D ist daher a/ℓ_P = 2√κ ≥ 3,07 (naiv mit Verdopplern
   8,7). **Korrektur:** Das Lipschitzgitter sitzt nicht in Planckeinheiten. Unter D ist η = 1/4 die
   Definition der Planckeinheit, und a/ℓ_P wird eine Vorhersage, sobald der volle Feldinhalt samt
   Kontaktterm der Eichfelder feststeht (Sakharovs Artenproblem in Gitterform).

**θ₁₃ neu eingeordnet.** TFPT hat θ₁₃ bereits als eigenen Trägerkanal: sin²θ₁₃ = φ₀·exp(−tr Y²) mit
tr Y² = 5/6 über der Hyperladung der fünf Plätze (FLAV.TH13.01). In H6 ist Y der Generator, den der
lokale 3+2-Split (Prämisse E) liefert. H6 gibt die führende Ordnung (tribimaximal aus S₄), TFPT die
Korrekturen: sin²θ₁₂ = 0,3067, sin²θ₂₃ = 1/2, sin²θ₁₃ = 0,0231. Laut TFPTs eigener Prüfliste liegt das
etwa 2σ von den Daten; JUNO entscheidet. θ₁₃ ist damit keine Lücke von H6 mehr, sondern eine
geerbte Vorhersage (K18 unter Lesart).

**Der gemeinsame Nahtoperator (Vorschlag 5 der Kritik).** Ein einziger Operator, die geliftete
Vierteldrehung g̃ um den Viertelpunkt, leistet zwei der drei geforderten Eigenschaften exakt.
- Er permutiert die vier Marken zyklisch und gibt c₃ = 1/(8π) je Spinorsektor, mit Spinphasen.
- Auf den geraden Spinstrukturen gibt er die Familien mit {−1, i, −i}, ohne eingesetztes Vorzeichen.

Die dritte Eigenschaft, E₈ mit vollständig behandelter Anomalie, gilt für g̃ selbst (S ist exakt).
Für vollständige Basisunabhängigkeit braucht sie die Ergänzung F.

**Prämissen von H6h3.**
- **A:** Der orientierte Tausch der Zyklen ist eine Symmetrie (Fixpunkt statt Kovarianz).
- **C′:** Familien = gerade Spinstrukturen der Naht; die chiralen 2-Formen sind ihr Λ².
- **D:** Gravitation als Zustandsgleichung (Jacobson).
- **E:** Der 3+2-Split ist lokal.
- **F:** Die Rahmungsanomalie wird durch eine invertierbare 3D E₈ Phase kompensiert.

**Offen.**
- K16: wechselwirkender Kontinuumslimes, Reflexionspositivität und SO(4).
- a/ℓ_P mit vollem Feldinhalt.
- F gegen F′.
- Warum der Fixpunkt.
- Der Weg von der T-Phase zur CP-Phase, falls es einen gibt.

## Davor: H6h2 — Eichdynamik und Prämissen (`y1_eichdynamik_praemissen.py`)

```
        K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 K11 K12 K13 K14 K15 K16    P  L  F  O  E
   H6h   P  P  P  P  P  P  P  P  P  P   P   P   L   P   O   O    13  1  0  2  0
   H6h2  P  P  P  P  P  P  P  P  P  P   P   P   L   P   P   O    14  1  0  1  0
```

**Eichdynamik (K15) hergeleitet unter einer kleinen Prämisse E** (der 3+2-Split ist lokal):
- Der örtlich variierende Split definiert auf jeder Gitterkante eine Paralleltranslation
  U = polar(V_x* V_y); sie ist unitär und eichkovariant (U ↦ h_x* U h_y).
- Fermionschleifen auf dem Lipschitzgitter induzieren den Yang–Mills-Term. Mit Fermionmasse als
  Infrarotregler und im Schwachfeldbereich (F/M² ≤ 0,12) ist die Wirkung **quadratisch** in F
  (Fluss 1 und 2 gleich auf 2 %), **positiv**, **eichinvariant**, und die Kopplung **läuft
  logarithmisch** mit d(1/g²)/d log(1/M) = 0,287 gegen 0,270 im Kontinuum (16 Diracarten).
  Ohne Masse setzt das Feld selbst die Infrarotskala (F² log F), wie erwartet.
- Eine Generation hat für SU(3), SU(2) und U(1) (GUT-normiert) dieselbe Dynkin-Summe 2. Die
  induzierten Kopplungen sind daher gleich, **sin²θ_W = 3/8 an der Gitterskala**; TFPT hat 3/8 an
  der Trägerskala.

**Prämissen geprüft.**
- **A** (keine ausgezeichneten Zyklen) ist tragend. Ohne A wären τ, c₃ (1/(2π) oder 1/(4π)), die
  Kragenbreite, die Fermionenzahl der Naht, E₈ und die Spinstruktur frei.
- **B entfällt.** Drei chirale Mittelformen gibt es nur in d = 4. Mit der 2D Naht ist die
  Kodimension 2, und A auf allen vier Zyklen macht die Normalebene gleichwertig. B folgt also aus
  A, C und der Naht.
- **C gegen Daten.** Die Familiengruppe S₄ (Ordnung 24) mit Restsymmetrien Z₃ für die geladenen
  Leptonen und Z₂ × Z₂ für Majorana-Neutrinos gibt genau die **tribimaximale Mischung**:
  sin²θ₁₂ = 1/3 (TFPTs Leitterm; mit −φ₀/2 sind es 0,3067 gegen gemessen 0,307), sin²θ₂₃ = 1/2
  (gemessen 0,47 bis 0,57), aber **sin²θ₁₃ = 0 gegen gemessen 0,022**. C besteht zwei von drei
  Winkeln in führender Ordnung; θ₁₃ verlangt eine Brechung, die H6 noch nicht liefert. Die
  Majorana-Forderung passt zu TFPTs Kandidaten vom 18. November.
- **D** bleibt der Rahmen (Jacobson wie TFPT). Entscheidender Gegentest für später: die
  Flächendichte der Verschränkung des Gittervakuums ausrechnen; sie legt fest, ob die Gitterskala
  die Planckskala ist (2π/η = 8π = 1/c₃).
- **E** (neu, mild): der Split ist ein lokales Feld.

**Was jetzt noch offen ist.**
1. **K16:** Der freie Sektor hat die renormierbare Struktur (logarithmischer Lauf). Ein
   wechselwirkender Kontinuumslimes ist nicht gezeigt.
2. **θ₁₃:** Die S₄-Restsymmetrie gibt null; der gemessene Wert verlangt eine Korrektur. Das ist
   jetzt die schärfste Datenprüfung von H6.
3. **D:** Gravitation nur im Jacobson-Rahmen; Gegentest oben.

## Davor: H6 gehärtet (H6h)

Gleicher Prüfstand, Katalog um K15 (Eichdynamik in 4D) und K16 (Kontinuumslimes) erweitert,
damit nichts Fehlendes verdeckt wird.

```
        K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 K11 K12 K13 K14 K15 K16    P  L  F  O  E
   H6    P  P  L  P  P  P  L  P  L  L   P   L   O   P   O   O     8  5  0  3  0
   H6h   P  P  P  P  P  P  P  P  P  P   P   P   L   P   O   O    13  1  0  2  0
```

**Was die Härtung leistet.** Die fünf Lesarten von H6 sind durch Rechnungen ersetzt, die nur noch
auf vier offen genannten Prämissen stehen:
- **A** Keine ausgezeichneten Zyklen (zeitneutral und basisunabhängig).
- **B** Die Naht ist eine Kodimension-2-Fläche mit gleichwertiger Normalebene.
- **C** Familien = chirale Formen mittlerer Dimension der Raumzeit.
- **D** Gravitation als Zustandsgleichung (Jacobson, wie TFPT).

Die Ersetzungen im Einzelnen:
1. **GSO, 16 Fermionen und E₈ aus A.** Eine chirale Fermion-Naht, die von der Wahl der
   Basiszyklen nicht abhängt, ist nur für N ≡ 0 mod 16 modular konsistent (S exakt, T bis auf eine
   Phase; geprüft N = 8, 16, 24, 32, 48). Das kleinste N = 16 gibt c = 8 und das Gitter D₈⁺ = E₈.
   Der (P,P) Sektor trägt die 8 flachen Nullmoden des E₈ Worts. Die frühere Lesart „B ist über
   Spinstrukturen summiert“ ist damit hergeleitet.
2. **Familien, ihre Charaktere und das Drehzentrum aus C.** σ und die Vierteldrehung erzeugen auf
   den drei chiralen 2-Formen eine Gruppe der Ordnung 24 mit dem Charakter Standard × Signum von
   S₄. Das ist genau die Darstellung, die TFPT für die Familienmonodromie aufgezeichnet hat (v117,
   „signgetwistete Standarddarstellung von S₄“). H¹ der Marken ist die Standarddarstellung (Satz).
   Der früher „gewählte“ q-Charakter ist also das Signum von S₄. Die beiden Lesarten derselben
   Familien stimmen nur bei zyklischer Wirkung auf den Marken überein; das erzwingt das
   Drehzentrum auf dem Viertelpunkt ({−1, i, −i} gegen {−1, 1, 1}).
3. **3 + 1 aus B und C.** Drei chirale Mittelformen gibt es nur in d = 4 (d = 8: 35, d = 12: 462).
   Naht plus gleichwertige Normalebene ist ℤ[i] ⊕ ℤ[i]·j, also das Lipschitz-Gitter.
4. **Gravitation unter D mit Eingaben aus H6.** TFPTs Herleitung braucht c = 8, |μ₄| = 4 =
   |ℤ₂|·χ und 2π/η = 8π = 1/c₃. H6 liefert alle drei; c = 8 ist jetzt abgeleitet (N = 16), nicht
   gesetzt, und stimmt mit TFPTs g_car + N_fam = 8 überein.

**Neue Folgerung, zu prüfen.** Eine chirale Naht mit c = 8 hat die Rahmungsanomalie
exp(−2πi c/24) = exp(−2πi/3). Sie verlangt einen Zufluss aus einem 3D Volumen, also eine
E₈-artige Wand, die an der Naht endet.

**Was noch fehlt (ehrlich).**
- **K15 Eichdynamik:** Die Erhaltungsgruppe (SU(3) × SU(2) × U(1))/ℤ₆ folgt aus 3 + 2 (TFPT,
  exakt). Yang–Mills-Felder in 4D folgen aus H6 nicht.
- **K16 Kontinuumslimes:** Das Lipschitzgitter sitzt in Planckeinheiten; ein wechselwirkender
  Skalenlimes ist nicht gezeigt.
- **Prämissen B, C, D** sind physikalische Postulate, keine Herleitungen. Die quantitativen
  Ergebnisse (α, φ₀, Massen) kommen wie bisher aus dem TFPT Compiler mit dessen eigenen
  bedingten Schritten.

**Nächste Schritte am Ganzen (kein Nebenschauplatz).** K15 auf dem Lipschitzgitter angehen:
TFPTs Befund „Grassmannkrümmung spannt s(u(3) ⊕ u(2))“ als Gitter-Eichfeld auf den Verbindungen
von ℤ⁴ formulieren und prüfen, ob die Wirkung Yang–Mills ergibt. Danach K16 und für jede Prämisse
einen entscheidenden Gegentest. Nach jedem Schritt läuft der ganze Prüfstand neu.

## Vorher: vollständige Hypothesen auf einem festen Prüfstand (`v1_vollhypothesen.py`)

Arbeitsweise seit dem Hinweis des Nutzers: Jede Hypothese deckt die ganze Lösung ab und wird
gegen denselben Katalog K1…K14 geprüft. Ein Fehlschlag beendet die Prüfung nicht.

Katalog: K1 vier Marken mit zwei Blättern · K2 quadratische Lage τ = i · K3 μ₄ permutiert die
Marken zyklisch · K4 Spinstruktur (A,A) · K5 c₃ = 1/(8π) · K6 Zustand und Zeit mit festgelegter
Breite · K7 flacher Sektor für das E₈ Wort · K8 P2 über P1 · K9 drei Familien als gleiche Kopien ·
K10 Familiencharaktere {−1, i, −i} · K11 E₈ · K12 3+1 · K13 Gravitation · K14 keine freien
Parameter.

```
        K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 K11 K12 K13 K14    P  L  F  O  E
   H1    E  P  L  P  P  P  F  P  O  P   E   F   F   P     7  1  3  1  2
   H2    P  P  L  P  P  F  F  P  P  F   P   L   O   P     8  2  3  1  0
   H3    E  E  E  P  O  O  F  E  E  E   E   E   P   F     2  0  2  2  8
   H4    E  P  L  P  P  P  L  P  O  P   P   F   F   P     8  2  2  1  1
   H5    P  P  L  P  P  F  L  P  P  F   P   L   O   P     8  3  2  1  0
   H6    P  P  L  P  P  P  L  P  L  L   P   L   O   P     8  5  0  1  0
```
P hergeleitet, L unter benannter Lesart, F gescheitert, O offen, E Eingabe.

**H6 ist die erste vollständige Hypothese ohne Fehlschlag.** In Worten:
- **Raumzeit** ist das ganzzahlige Gitter der größten assoziativen Divisionsalgebra, die
  Lipschitz-Quaternionen ℤ⁴ = ℤ⟨1, i, j, k⟩. Ihre Einheiten ±1, ±i, ±j, ±k machen alle vier
  Achsen gleichwertig (keine vorab ausgezeichnete Zeit). Die Fermionen sind auf allen vier
  Zyklen antiperiodisch (invariant, ohne Nullmoden). 3 + 1 ist Re ⊕ Im.
- **Die Naht** ist eine Gaußsche Ebene ℤ[i] darin: quadratisch, τ = i, Spinstruktur (A,A). Die
  π-Drehung −1 gibt zwei Blätter und vier Marken. Die Vierteldrehung um einen Viertelpunkt
  permutiert die Marken zyklisch; c₃ = 2/(4π·4) = 1/(8π); die Kragenbreite ist ℓ = L/2.
- **P2** folgt aus der Naht wie in TFPT (P2.PARTITION.01).
- **E₈** ist dasselbe Gitter zweimal: aus ℤ⁸ = (ℤ⁴)² mit Halbzahlkleber (112 + 128 = 240) und als
  GSO-Summe von 16 Majorana-Weyl-Fermionen auf der Naht (120 + 128 = 248). Der (P,P) Sektor hat
  genau 8 komplexe Nullmoden, flach für jede Breite; dort lebt das E₈ Wort. Das löst die Gabel
  aus Abschnitt N.
- **Drei Familien** sind die drei selbstdualen 2-Formen (die drei Gaußschen Geraden ℤ[i], ℤ[j],
  ℤ[k]), von σ zyklisch permutiert, also gleiche Kopien. Die Vierteldrehung wirkt auf ihnen mit
  {1, i, −i}; zusammen mit dem Vorzeichencharakter von q, den TFPT auch für die Trägermoden
  nutzt, ergibt das genau {−1, i, −i}.
- **Gravitation**: 4D lokal; die Einsteingleichung mit Koeffizient 1/c₃ liefert nur TFPTs
  bedingte Herleitung (GRAV.NONLINEAR.01).

**Die fünf Lesarten von H6 (dort kann sie noch falsch sein).**
1. K3: Keine Marke ist ausgezeichnet, also liegt das Drehzentrum auf einem Viertelpunkt.
2. K7: Der Randoperator B ist der Naht-Dirac summiert über die Spinstrukturen (GSO).
3. K9: Familien = selbstduale 2-Formen.
4. K10: Die Familien tragen den Vorzeichencharakter von q. Diese Lesart ist nach Ansicht des Ziels
   gewählt; derselbe Charakter steht aber schon in TFPTs Trägermoden.
5. K12: Raumzeit = größte assoziative Divisionsalgebra (Frobenius).

**Ausschlüsse dieses Prüfstands (wiederverwendbar).**
- *Dreierdrehungen verbieten zeitneutrale invertierbare Fermionen.* Das hexagonale 2D-Gitter und
  das Hurwitz-4D-Gitter haben keine einheiteninvariante Spinstruktur ohne Nullmoden; das
  quadratische und das hyperkubische (Lipschitz) haben sie. Deshalb μ₄ und nicht μ₆, und deshalb
  Lipschitz statt Hurwitz.
- *Tensorprodukt-Spektraltripel* (H3) verlieren den flachen Sektor: mit invertierbarem
  geometrischem Dirac gibt es keine Nullmoden.
- *Naht allein* (H1, H4) liefert kein 3+1 und keine Gravitation.

**Nächster Schritt (für H6 als Ganzes).** Jede der fünf Lesarten durch eine Herleitung oder einen
entscheidenden Test ersetzen und danach den ganzen Prüfstand neu laufen lassen. Dazu die
Gravitation (K13) auf der Lipschitz-Raumzeit selbst prüfen, statt sie aus TFPT zu übernehmen.

---

Gerechnet am 23. November 2026 mit `.venv/bin/python` (numpy, scipy, sympy, mpmath). Zehn Skripte,
alle mit Rückgabewert 0; Ausgaben in `ausgabe_*.txt`. PLAN, progress, Ledger und Atlas sind
unberührt. Ausgangspunkt war der konsolidierte Stand in `docs/research/TFPT-STAND.md`
(Abschnitt „Befunde 20.–23. November“): Die offenen Eingaben laufen auf „woher die Vier“ und
„warum minimal“ hinaus.

## Übersicht

| Nr | Urteil | Stufe | Erste freie Wahl |
|---|---|---|---|
| **N Zeitneutrale Fermionnaht (Ursprung von P1)** | **Kandidat** | 1 für P1: Lage der Marken, Viertelzählung in c₃, Zustand; nicht 2, weil Naht und zwei Blätter Eingabe bleiben | Prinzip „keine vorab ausgezeichnete Zeitrichtung“ |
| **G Gaußscher Torus ℂ/ℤ[i]** | **Kandidat, korrigiert**: P2-Teil ist Doppelherleitung (P2 folgt aus P1) | 1 für P1 (Lage der Marken, c₃) | Startklasse CM-Torus; Auswahl über eine P2-Regel |
| W Weil-Darstellung von SL(2, 𝔽₅) | **Kandidat** | 2 für P2 unter benannter Lesart | Auswahlregel „fermionischer Teil = treues SU(2) Dublett“; Lesart „fünf Plätze = Ortsbasis eines Fünfersystems“ |
| W2 Clifford-Gruppe eines Fünfersystems | Uhrenrelation **bestätigt**, sechs Fasern **verworfen** | – | – |
| W3 die 16 unter 2I | **bestätigt**, exakt | – | – |
| H minimaler holografischer Code | **offen** | 1, an der Basislinie | Abstand 3; Holografie |
| P Poincaré-Sphäre | bereits in TFPT (v219), kein neuer Weg | – | – |
| T c₃ als Bekenstein-Viertel mal Unruh-Periode | bereits in TFPT, kein neuer Weg | – | – |
| J Jordanalgebra J₃(𝕆), R Renormierung | **offen**, nicht gerechnet | – | siehe IDEAS.md |
| F Familientest zu W (Nachtrag) | Kopienweg **verworfen**; Wörterbuch über 𝔽₅ **offen** | 1 für die Zahl 3 | siehe Abschnitt F |
| K Familien als étale Kohomologie (Nachtrag) | Test **bestanden**, aber nicht unabhängig; Frobeniusgleichheit **neu** | 1 | siehe Abschnitt K |

Keine vollständige Lösung (Stufe 4 fehlt: physische Zeit und Zustand, 3 + 1, Kontinuum,
Gravitation).

**Korrektur (Hinweis des Nutzers, mit Quellen geprüft).** P2 folgt in TFPT bereits aus P1: Aus der
markierten Viermarkennaht folgen Fünferträger und Teilung 3 + 2 exakt (P2.PARTITION.01,
AX.P2.01). W, F und der P2-Teil von G sind daher eine zweite Herleitung von schon Hergeleitetem,
kein neuer Ursprung. Was von G bleibt, ist die P1-Seite: die quadratische Lage der Marken und
c₃ aus demselben Torus. Dort wählt G den Torus aber über eine P2-Regel (Spinordublett im
Weil-Oszillator); eine Begründung von der P1-Seite selbst fehlt. Daran setzt Abschnitt N an.

## N: Ursprung von P1 — die zeitneutrale Fermionnaht (`n1_zeitneutrale_naht.py`)

**Prinzip.** Die Naht ist eine geschlossene, orientierte, flache 2D-Fläche mit freien Fermionen,
und sie zeichnet **keine Zeitrichtung** aus: Rekonstruiert man nach Osterwalder–Schrader entlang
jedes ihrer beiden Basiszyklen ein Quantensystem, entsteht dasselbe orientierte thermische
System. Jede Rekonstruktion soll ein echter KMS-Zustand sein (Fermionen in der Zeit
antiperiodisch) mit invertierbarem Hamiltonoperator (keine Nullmoden). Zeit kommt dann aus dem
Zustand, nicht aus der Geometrie.

**Befund (exakt bis auf die Gittersuche).**
- Rekonstruktion A hat (β/L, Drall/L) = (Im τ, Re τ), die orientierte Vertauschung B hat
  (Im τ/|τ|², −Re τ/|τ|²). Gleichheit erzwingt |τ| = 1 und Re τ = 0, also **τ = i**. Zusammen mit
  KMS und Invertierbarkeit bleibt im ganzen Fundamentalbereich nur **τ = i mit Spinstruktur
  (A, A)**. Lässt man Parität zu, bleibt die ganze rhombische Familie |τ| = 1; die Orientierung
  von P1 ist also tragend. Am hexagonalen Punkt wechselt der Drall das Vorzeichen.
- Die Zeitneutralitätsgruppe ist die Vierteldrehung S, Ordnung 4 = μ₄. Kein Impuls in (ℤ+½)²
  ist unter ihr fix, alle Bahnen haben Länge 4. Die bisher gewählte **Viertelzählung ist damit
  exakt die Zählung zeitneutraler Bahnen**, und c₃ = (Diracrang 2)/(4π · |μ₄|) = **1/(8π)** ohne
  Wahl einer Fläche (s·Bahnzahl → 0,039788735773). Mit Parität wäre es 1/(4π), ohne Neutralität
  1/(2π).
- **Zustand und Zeit.** KMS-Zustand der Diracfermionen bei selbstdualer Temperatur β = L,
  Einteilchenenergien 2π|n + ½|/L, Lücke π/L. Zeit ist jeder der beiden Zyklen; die zwei Wahlen
  sind durch μ₄ äquivalent. Die Zustandssumme am Fixpunkt ist exakt Z(i) = |θ₃(i)/η(i)|² = 2.

**Was neu ist.** Bekannt waren: τ = i als einziger symmetrieerzwungener kritischer Punkt (KERNSYM,
dort ist die Symmetrie τ ↦ −1/τ vorausgesetzt und τ = i nur Sattel), die Wahl (−,−) durch
Viertelinvarianz und fehlende Nullmoden, und 1/(8π) aus einer gewählten Viertelzählung (alles
23.11. oder früher). N ersetzt die zusätzlich genannte quadratische Lage und die gewählte
Viertelzählung durch **eine** physikalische Forderung und liefert dazu einen Zustand mit Zeit.
Energieminimierung leistet das nicht (KERNSYM: Sattel).

**Stufe.** Kandidat, Stufe 1 für P1: Lage der Marken (μ₄, τ = i), Spinstruktur, Viertelzählung
und Normierung c₃, KMS-Zustand. Nicht Stufe 2, weil die Naht selbst (flach, 2D, Fermionen) und
die zwei Blätter mit den vier Marken Eingabe bleiben. Die Blattinvolution der Marken ist eine
π-Drehung um eine Marke; ihr Zentrum liegt um eine halbe Periode neben dem Zentrum der
Vierteldrehung, sie ist also nicht einfach μ₄².

**Basislinie.** Die Rechnung ist elementare Gittergeometrie; die Stärke liegt im Prinzip. Die
Zahlen 4 und 1/(8π) folgen ohne Anpassung; die Gegenproben (Parität zugelassen, keine
Neutralität) geben andere Werte.

**Stärkster Einwand.** „Keine vorab ausgezeichnete Zeit“ ist eine Prämisse (Nähe zur
Thermal-Time-Hypothese von Connes und Rovelli), keine Herleitung. Und ob der KMS-Zustand bei
β = L die operative Saat des Originals ist (Fockzustand der Calderón-Polarisation, Erzeuger
|B_Σ|), ist nicht geprüft.

### N gegen die operative Saat des Originals (`c1_vergleich_zustand.py`)

**Wie das Original den Zustand baut** (KERNSIX20261120 Schritte 1, 3, 4; P1BIND20261120). Die Naht
Σ ist der Raum, die Zeit läuft normal dazu. Kragen D = γ_n(∂_n + B), Polarisation = Calderón-
Projektor P_C, auf den zwei Normalhälften flach. Der Randzustand eines Kragens der Breite ℓ ist
rein; auf einer Fläche ist er thermisch für B bei β = 2ℓ. Breite null gibt den flachen
Kernzustand C = ½, aus dem die aufgezeichneten Zahlen des E₈ Worts stammen. **ℓ ist im Original
frei.** Nullmoden von B bleiben für jede Breite flach; dort lebt das E₈ Wort.

Der Vergleich verlangt, N auf die volle euklidische Geometrie zu heben: thermischer Kreis 2ℓ mal
Naht T²(L, L), Fermionen überall antiperiodisch. N in seiner ersten Form (Zeit innerhalb der Naht)
war eine andere Schnittrichtung derselben Geometrie.

**Befund.**
1. Die Rekonstruktion entlang des thermischen Kreises und die entlang eines Nahtzyklus sind genau
   bei **2ℓ = L** dasselbe System (Spektren β·λ identisch; Abweichung bei ℓ = 0,25…1,0 sonst
   8,7 bis 36). Zeitneutralität **legt die offene Breite des Originals fest**.
2. Damit wird die im Original freie Größe T/λ₁ = 1/(2ℓλ₁) zu **1/(√2 π) = 0,2251** (λ₁ = √2 π/L).
   Die Verschränkung der beiden Kragenflächen ist endlich, **S = 0,5154**, ohne Abschneiden;
   der flache Kernzustand hätte log 2 je Mode.
3. Bei (A, A) hat das geometrische B **keine Nullmoden** (nur (P, P) hat zwei).
4. Lage der Marken: Die mit λ₁ normierte Spektralwirkung hat am Quadrat in der rechteckigen
   Familie ein Minimum, unter Verdrillung einen Sattel (−0,33 bei Re τ = ±0,1), wie KERNSIX
   Schritt 4. N schließt die Verdrillung durch die orientierte Zeitneutralität aus und braucht
   keine Minimierung.

**Urteil: gemischt, mit scharfer Gabel.**
- *Übereinstimmung:* Der zeitneutrale Zustand ist genau der Kragenzustand des Originals (rein,
  auf jeder Fläche thermisch für B), und N liefert die im Original fehlende Breite ℓ = L/2.
- *Abweichung:* Die aufgezeichneten E₈ Zahlen stammen aus dem flachen Zustand der Breite null.
  Den schließt N für jeden Mode aus, der an das geometrische B koppelt, und bei (A, A) gibt es
  keine geometrischen Nullmoden. Die Abbruchregel greift damit für den flachen Sektor, nicht für
  den Zustand als Ganzes.
- *Gabel:* Beide passen nur zusammen, wenn das E₈ Wort auf Nullmoden eines **inneren, nicht an
  die Nahtgeometrie gekoppelten Teils** von B sitzt (B = B_geom ⊕ B_int statt Tensorprodukt).
  Das Original lässt das offen: P1BIND liest B als den endlichen D_F, KERNSIX modelliert B als
  Dirac auf der Naht. Diese Frage entscheidet, ob N trägt.

**Kleinster nächster Test (Kern).** Im Original nachlesen, ob B_Σ die geometrische Naht und den
inneren Teil als direkte Summe oder als Tensorprodukt trägt, und ob die Nullmoden von D_F auch bei
eingeschalteter Nahtgeometrie Nullmoden bleiben. Tensorprodukt mit invertierbarem geometrischem
Teil: N ist mit dem E₈ Wort unverträglich. Direkte Summe: N trägt und fixiert die Breite.
Parallel verfolgt eine andere Sitzung „P1 aus positiver, verlustfreier Nahtverklebung“; beide
Wege zielen auf dieselbe Lücke und sollten vor einer Atlasübernahme abgeglichen werden.

## G: Ein Objekt für P1 und P2 — der Gaußsche Torus ℂ/ℤ[i] (`g1_gauss_torus.py`)

**Prinzip.** Startklasse sind flache 2-Tori, deren Gitter ein Ganzheitsring O_K ist (Tori mit
komplexer Multiplikation, Klassenzahl 1, neun Fälle). Auf einem Torsionsstück des Torus lebt ein
endlicher Oszillator, die Weil-Darstellung auf den Funktionen dieses Stücks. Eine einzige
physikalische Forderung: Auf dem Torsionsstück, auf dem die Symmetrien des Torus als volle
Einheitengruppe eines Primkörpers wirken, hat der Oszillator einen Spin-½-Dublett als
fermionischen Teil.

**Kette.** Spin-½-Dublett ⇒ p = 5 (eindeutig, W) ⇒ Einheitengruppe der Ordnung p − 1 = 4 ⇒ ℤ[i],
der quadratische Torus. Aus diesem einen Torus folgt:
- **P2.** Plätze = Funktionen auf der Torsionsgeraden E[2−i] ≅ 𝔽₅. 3 + 2 = bosonischer Dreier ⊕
  fermionisches Dublett. Y ist die einzige spurfreie Größe im Kommutanten (−1/3, +1/2). Die 16
  trägt die Hyperladungen, jeder Sektor ist unter 2I irreduzibel (W3). SL(2,5) = 2I führt zu E₈
  (Ikosiane, McKay, h = 30).
- **μ₄ und Uhren.** Die komplexe Multiplikation mit i ist auf der 5-Torsion genau das Weyl-Element
  [[0,−1],[1,0]]; ihr Weil-Bild ist die Fourier-Vierteldrehung, also TFPTs R₅ bis auf das
  Vorzeichen. Auf E[2−i] wirkt sie als t ↦ 2t = G, die Torsionstranslation ist C₅, und
  G C₅ G⁻¹ = C₅² gilt exakt.
- **P1.** Randnormierung c = Kovolumen/(2π · |Einheiten|) = 1/(8π) = c₃ (Weylterm der
  Dirac-Wärmespur pro Einheitenbahn, numerisch auf 12 Stellen). Die einzige invariante
  Spinstruktur ohne Nullmoden ist (−,−). Der Modul τ = i ist der Fixpunkt, den KERNSYM als
  symmetrieerzwungen kennt. Die vier Marken sind die Bilder der 2-Torsion unter der
  Doppelüberlagerung y² = z⁴ − 1, mit Windung +1 des quadratischen Differentials (Atlas 23.11.).
  Die Orientierung kommt aus der komplexen Struktur, die Reflexion z ↦ z̄ erhält ℤ[i].
- **Familien.** H¹ von P¹ ohne die Marken: Rang 3, Charaktere {−1, i, −i}, über 𝔽₅
  frobeniusgleich (K), Residuenbild der Trägermoden.

**Auswahltabelle (alle neun CM-Tori der Klassenzahl 1).**
- (A) Volle Einheitenwirkung auf einem Primkörper: p = 5 für ℤ[i], 7 für ℤ[ω], 3 für d = 2 und 11,
  sonst kein p.
- (B) Dublett: nur bei p = 5.
- (C) Invariante invertierbare Spinstruktur: alle außer ℤ[ω].

Nur d = 1 erfüllt alle drei Regeln. Die c-Werte der neun Tori liegen zwischen 0,023 und 0,51;
nur ℤ[i] gibt 1/(8π).

**Neu gegenüber TFPT.** Die Einzelteile waren bekannt: quadratischer Torus mit μ₄ (τ = i,
j = 1728), c₃ aus der Viertelzählung (23.11., bedingt), 5 = (2+i)(2−i) als Zahlenbeobachtung
(ARITH.TRICHOTOMY.01), x ↦ 2x und Y_par, PGL(2,5), 2I und E₈. Neu ist die Brücke: P2s fünf
Plätze sind ein Torsionsstück desselben Torus, dessen Wärmespur P1s c₃ gibt. Und eine einzige
Forderung (Spin-½-Dublett) wählt diesen Torus unter allen CM-Tori aus. Damit beantwortet G beide
Leitfragen der Konsolidierung: **Die Vier ist die Einheitengruppe ℤ[i]^× des Torus**, erzwungen
durch das Dublett. **Die Minimalität ist das erste Torsionsstück**, auf dem die Einheiten voll
wirken.

**Stufe.** Stufe 3 unter benannten Lesarten: P1 (Normierung, Spinstruktur, Modul, Marken,
Orientierung) und P2 (fünf Plätze, 3 + 2, Y) kommen aus demselben Objekt, ohne angepasste
Konstante; die vorhandenen Ergebnisse trägt der TFPT Compiler aus P1 und P2. Nicht Stufe 4:
physische Zeit und Zustand, 3 + 1, Kontinuum und Gravitation liefert G nicht.

**Benannte Lesarten (erste freie Wahlen).**
1. Startklasse: Der Rand ist ein flacher Torus mit Ring-Gitter.
2. Plätze = Schrödinger-Modell der Weil-Darstellung auf dem ersten Torsionsstück mit voller
   Einheitenwirkung.
3. c₃ = Weylterm der Dirac-Wärmespur, gezählt pro Einheitenbahn (in TFPT bisher bedingt).
4. Die Paritätsbasis (Y_par) ist die physische Farb/Schwach-Basis. Das ist eine Identifikations-
   änderung gegenüber der markierten Ladung, siehe W.

**Basislinie.** Die Auswahl ist strukturell, kein Zahlentreffer. Effektiv entscheidet
„Einheitengruppe der Ordnung 4“, und das Dublett erzwingt sie. Der Wert 1/(8π) ist keine
unabhängige Vorhersage, weil die Viertelzählung am quadratischen Torus formuliert wurde.
Unabhängig ist, dass die Regel von der P2-Seite denselben Torus wählt. Gegenproben: ℤ[ω] (sechs
Einheiten) scheitert an zwei Regeln, die Tori mit nur ±1 an der Dublettregel.

**Stärkster Einwand.** Die Lesarten 2 und 3 sind gewählt. SU(3) × SU(2) entsteht nicht, nur 2I
und U(1) × U(1). Die Dynamik fehlt ganz. G erklärt die Eingaben, nicht die Physik dahinter.

**Kleinster nächster Test (Stufe 4, kein Nebenschauplatz).** Zeit und Zustand aus demselben
Torus: Die Vierteldrehung vertauscht die beiden Zyklen (τ ↦ −1/τ hat τ = i als Fixpunkt).
Prüfen, ob der Weil-Oszillator auf E[5] zusammen mit der Torus-Wärmespur eine reflexionspositive
Halbgruppe mit ausgezeichneter Zeitrichtung definiert, deren KMS-Zustand die Normierung c₃
trägt. Abbruch, wenn Zeit oder Zustand eine zusätzliche freie Wahl brauchen.

## W: Weil-Darstellung von SL(2, 𝔽_p) (`w1_weil_darstellung.py`)

**Befund (Gruppen exakt erzeugt, Matrizen numerisch auf 10⁻⁹).**
- Für p = 3, 5, 7, 11, 13 erzeugen Fouriertransformation S und Chirp T genau SL(2, p) der
  Ordnung p(p²−1); die Darstellung ist linear.
- Das zentrale −I wirkt als ± Parität und zerlegt ℂ^p in einen bosonischen (−I ↦ +1) und einen
  fermionischen Teil (−I ↦ −1).
- **Nur bei p = 5** ist der fermionische Teil ein treues SU(2) Dublett (Dimension 2, Kern
  trivial, det = 1, Frobenius–Schur −1) **und** der bosonische Teil eine treue reelle
  3D-Drehung von SL(2,5)/±1 = A₅ (Dimension 3, Kern {±I}, det = 1, Frobenius–Schur +1).
  Bei p = 3 hat das Dublett det ≠ 1 und der bosonische Teil Dimension 1; ab p = 7 ist der
  fermionische Teil größer als 2. Beide Weil-Darstellungen bei p = 5 (additiver Charakter mit
  Rest oder Nichtrest) geben dasselbe Bild.
- SL(2,5) ≅ 2I (binäre Ikosaedergruppe), also dieselbe Gruppe, deren 120 Einheitsquaternionen
  E₈ erzeugen und deren McKay-Graph Ê₈ mit Coxeterzahl 30 ist (in TFPT bekannt).
- **Hyperladung erzwungen.** Der Kommutant der ganzen Darstellung auf ℂ⁵ hat Dimension 2. Die
  einzige spurfreie Größe darin ist Y ∝ I − 5P (P = S² = Parität), mit Skala 1/12 also −1/3 auf
  dem bosonischen Dreier und +1/2 auf dem fermionischen Zweier. Das ist TFPTs Y_par.
- **Vierteldrehung ohne Eingabe.** S ist die Fouriertransformation des diskreten Phasenraums,
  S⁴ = I für jedes p, S² = Parität. Ihr Spektrum ist {1, −1, −1, i, −i}, also bis auf das globale
  Vorzeichen das Spektrum {1, 1, −1, i, −i} von TFPTs kohärentem Vierteloperator R₅
  (Farbe 1, 1, −1; schwach ±i).
- **Kein Spin-1 Vektor des Dubletts.** Der bosonische Dreier ist nicht Sym² des Dubletts, sondern
  Sym² des galoiskonjugierten Dubletts (√5 → −√5), also das andere Ikosaedertriplett 3′.

**Basislinie.** Die Teilung ((p+1)/2, (p−1)/2) ist bei jedem p die balancierteste; das Verhältnis
−2/3 der Hyperladungen folgt aus 3 + 2 und Spurfreiheit. Informativ ist allein, dass eine
einfache physikalische Regel p = 5 auswählt und dass dann der fermionische Teil das schwache
Dublett ist. Zwei natürliche Regeln („Spinorteil ist ein Dublett“, „Bosonteil ist eine reelle
3D-Drehung“) wählen beide p = 5; einfache Gegenregeln („Bosonteil eindimensional“,
„Spinorteil dreidimensional“) wählen p = 3 bzw. keine Primzahl. Ein strukturelles Ergebnis,
kein Zahlentreffer.

**Stärkster Einwand.** Es entstehen endliche Gruppen, keine Eichgruppen: 2I ⊂ SO(3) × SU(2) und
der Kommutant U(1) × U(1), nicht SU(3) × SU(2). Der Dreier ist ein reelles Ikosaedertriplett,
kein komplexes Farbtriplett. Und die Auswahlregel ist eine Wahl: Warum soll der fermionische
Teil eines endlichen Oszillators das schwache Dublett sein?

**Abgleich mit Atlas und Prüfkatalog.** Die Permutation x ↦ 2x auf den fünf Plätzen, Y_par =
(I − 5G²)/12 und PGL(2,5) als Uhrengruppe stehen seit dem 11. November bzw. Oktober in TFPT
(`TFPT-CLOCK-LADUNG-TRANSFER-20261111.txt`, Uhrenbündel 20261028); 2I, E₈ und McKay seit v219.
Neu sind die Weil-Darstellung als gemeinsamer Ursprung, die Auswahlregel für p = 5, der
Kommutant als Begründung von Y und die Zerlegung der 16 (W3). **Benannte Identifikations-
änderung:** W übernimmt Y_par, also die Teilung in gerade und ungerade Funktionen auf 𝔽₅, nicht
die ursprüngliche Ladung, die auf markierten Platzkoordinaten diagonal ist. Laut Befund vom
11. November erhält der dort vorgeschlagene Operator H* zwar Y_par, aber nicht die markierte
Ladung (‖[H*, Y_original]‖² = 19/12 bei der Standardaufteilung, 17/6 bei den mit G²
verträglichen Aufteilungen); der Wechsel ändert also die gemeinsame Wirkung. W behauptet,
dass die physische Farb/Schwach-Basis die Paritätsbasis des Fünfersystems ist; das ist zu
prüfen, nicht zu übergehen.

**Trägt / liefert nicht / bricht.** Trägt: fünf Plätze, 3 + 2 mit kanonischer Markierung, Y,
μ₄ (als Fourier-Vierteldrehung), die 16 mit Hyperladungen, 2I und damit E₈ und h = 30.
Liefert nicht: drei Familien, D₅ ⊕ A₃ mit ℤ₄ Glue, c₃, φ₀, 41, α, Eichgruppen, Raumzeit,
Gravitation. Bricht: nichts, ersetzt aber die markierte Koordinatenladung durch Y_par.

## W2: Clifford-Gruppe eines Fünfersystems (`w2_qudit5_clifford.py`)

- Die Dilatationen x ↦ 2x und x ↦ 3x liegen bis auf das Vorzeichen −1 im Weil-Bild; die
  Parität mit Vorzeichen +1.
- Verschiebung X (TFPT C₅) und Dilatation G erfüllen G X G⁻¹ = X², G⁴ = I und erzeugen
  AGL(1,5) = F₂₀. **TFPTs Uhrenrelation ist die affine Gruppe des Körpers 𝔽₅**, also der
  affine Teil der Clifford-Gruppe. Y kommutiert mit SL(2,5) und G, nicht mit X.
- **Verworfen:** Die sechs gegenseitig unverzerrten Basen eines Fünfersystems werden von der
  Clifford-Gruppe nur als PSL(2,5) ≅ A₅ (Ordnung 60, Stabilisator 10) permutiert, auch mit
  antiunitärer Konjugation. TFPTs sechs Fasern tragen PGL(2,5) mit Stabilisator F₂₀. Die
  Lesart „sechs Fasern = sechs Basen“ fehlt ein Element mit Nichtrest-Multiplikator.

## W3: die 16 unter 2I (`w3_sechzehn_unter_2I.py`)

Auf Λ^gerade ℂ⁵ = 16 ist **jeder Hyperladungssektor eine irreduzible Darstellung von 2I**:
ν^c 1, u^c 3′, Q 6, e^c 1, L 2, d^c 3′. Q ist die sechsdimensionale irreduzible Darstellung
(der Verzweigungsknoten von Ê₈). Das liegt am galoiskonjugierten Dreier: Wäre der Dreier der
Spin-1 Vektor des Dubletts, zerfiele Q als Spin 1 ⊗ Spin ½ = 2 ⊕ 4. Damit spricht W3 gegen die
Lesart „drei Farbplätze = drei Raumrichtungen“ (der stärkste Einwand gegen A1 vom 20.11.) und
für einen inneren Dreier.

## H: minimaler holografischer Code (`h1_holocode.py`)

**Befund.** [[5,1,3]] ist ein perfekter Sechsbein-Tensor (alle 20 Dreierreduktionen = I/8);
vier Qubits erreichen bei der Reinheitssumme nur 1 statt 0,75 (kein perfekter Vierbein-Tensor).
Quanten-Singleton gibt n = 5 für eine logische Einheit mit Abstand 3; das kleinste hyperbolische
q für Fünfecke ist 4. Die Dreiecksgruppe (2,4,5) hat als kleinstes glattes endliches Bild S₅
(Bring-Fläche, Geschlecht 4, 30 Ecken, 60 Kanten, 24 Fünfecke). In allen 120 erzeugenden
Paaren ist a eine Transposition (markiert 3 + 2 der Plätze), b² hat auf ℂ⁵ die Eigenräume 3 + 2.

**Basislinie.** Andere Signaturen liefern ebenso TFPT-nahe Gruppen (kleinste gefundene Bilder
in transitiven Darstellungen bis Grad 8): (2,3,7) → PSL(2,7) ≅ GL(3,2) (Hamming und Fano),
(2,5,5) → A₅, (2,4,6) und (2,3,8) → Ordnung 48. Für (2,4,5) ist S₅ beweisbar minimal: Gruppen der
Ordnung 40 oder 80 haben keine glatte (2,4,5)-Erzeugung (Sylow-Argument). Mit Abstand 2 statt 3 entsteht
(2,3,7) statt (2,4,5). S₅ ist also kein seltener Treffer.

**Stärkster Einwand.** Die Uhrenrelation G C₅ G⁻¹ = C₅² wird in keinem der 120 Paare erfüllt
(b normalisiert ⟨ab⟩ nie). Der Innenraum ist zweidimensional, liefert also 2 + 1, nicht 3 + 1.
Die Gravitationsseite (Ryu–Takayanagi) ist nicht gerechnet.

## P, T, J, R

P und T stehen bereits in TFPT (Poincaré-Sphäre und McKay in v219; Einstein 8π = 1/c₃ mit
Bekenstein 1/4 = 1/|μ₄|). J und R sind nicht gerechnet; ihre kleinsten Tests stehen in IDEAS.md.

## Was jetzt noch Eingabe ist (aus Sicht von W)

1. Die Auswahlregel für p = 5 (fermionischer Teil = treues SU(2) Dublett).
2. Die Lesart „fünf Plätze = Ortsbasis“ samt Übernahme von Y_par statt der markierten Ladung.
3. Familienzahl 3, D₅ ⊕ A₃, c₃, φ₀ und alles Dynamische wie bisher.

## F: Familientest (Nachtrag am selben Tag, `f1_familien.py`)

**A. Mehrere Fünfersysteme (der angekündigte Test).** Weil-Darstellung von Sp(2n, 5) auf
ℂ^(5ⁿ), zerlegt nach dem Howe-Paar aus diagonalem 2I und O(n, 5); O(n,5) vertauscht exakt mit
der Weil-Wirkung.
- n = 1: Dublett einmal, Partner das Vorzeichen von O(1) (das ist W).
- n = 2 (bosonisch 13, fermionisch 12): **kein** Dublett, weder 2 noch 2′.
- n = 3 (bosonisch 63, fermionisch 62): Dublett 2 **fünfmal**, sein Multiplizitätsraum ist eine
  irreduzible Fünferdarstellung von O(3,5) ≅ PGL(2,5) × ℤ₂; 2′ kommt nicht vor.

**Urteil A: verworfen.** Die Abbruchregel greift: Weder n = 2 noch n = 3 liefert drei Kopien
oder ein Triplett, und D₅ ⊕ A₃ erscheint nicht. Nebenbefund: Bei n = 3 trägt das Dublett
wieder eine Fünf. Dass der bosonische Teil bei n = 3 die Dimension 63 hat (TFPT: 68 = 5 + 63
Gedächtnisrichtungen), ist ein ungeprüfter Zufall; 63 = 2⁶ − 1 ist ebenso naheliegend.

**B. Wörterbuch über 𝔽₅.** Exakt geprüft:
- Die Nullstellen von z⁴ − 1 in 𝔽₅ sind **ganz 𝔽₅^×**: die vier Marken sind die Einheiten.
- Polynome vom Grad ≤ 4 sind über 𝔽₅ genau alle Funktionen auf 𝔽₅ (3125 von 3125): die fünf
  Schnitte von O(D) entsprechen den fünf Punkten der affinen Geraden.
- Die Dilatation x ↦ 2x (TFPTs G) hat auf ℂ^𝔽₅ den invarianten Teil span(δ₀, 1_{𝔽₅^×}) und die
  drei nichttrivialen Dirichlet-Charaktere mod 5 mit Eigenwerten −i, −1, +i. Das sind genau die
  μ₄-Charaktere {−1, i, −i}, die der Atlas dem Familienraum H¹ zuschreibt (Befund vom 23.11.).
- Die sechs Fasern sind P¹(𝔽₅) (in TFPT seit Oktober bekannt), ihr Punktstabilisator
  AGL(1,5) = F₂₀ ist die Uhrengruppe ⟨C₅, G⟩.
- Gauss-Summen: |g(χ)|² = 5 für alle drei, die Legendre-Summe ist √5 (goldener Schnitt, Ikosiane).
- Auswahl: „Die Einheitengruppe ist genau μ₄“ gilt nur für p = 5. Das ist eine zweite,
  unabhängige Regel neben „Spinorteil = Dublett“; beide wählen p = 5.

Damit ist das Zahlengerüst 3 Familien, 4 Marken, 5 Plätze, 6 Fasern dasselbe wie p − 2, p − 1,
p, p + 1 auf der projektiven Geraden über 𝔽₅, und die Familiencharaktere sind die nicht-
trivialen Einheitencharaktere.

**Basislinie zu B.** In TFPT sind 3, 4, 5 schon n − 1, n, n + 1 für n = 4 Marken
(Riemann–Roch). Das Wörterbuch ersetzt die Frage „warum vier Marken“ durch „warum 𝔽₅“; neu ist
die Lesart der Marken als Einheiten und der Familiencharaktere als Dirichlet-Charaktere, nicht
die Zahlen selbst.

**C. Familien als getwistete Sektoren.** Induziert man auf den sechs Punkten von P¹(𝔽₅) mit
jedem Charakter von 𝔽₅^×, entstehen: trivial 1 ⊕ 5 (Steinberg), Legendre 3 ⊕ 3′ (bosonisch),
die beiden quartischen Charaktere jeweils dieselbe irreduzible fermionische 6.
Die drei Familiensektoren sind also **nicht gleich**: Unter der Galoiswirkung i ↦ −i zerfallen
die Familiencharaktere in 1 + 2 (Legendre reell, quartisches Paar konjugiert).

**Urteil F.** Kopienweg verworfen. Das Wörterbuch ist **offen**, Stufe 1 für die Zahl 3 und die
Familiencharaktere, aber ohne Mechanismus für drei gleiche Kopien der 16. Stärkster Einwand:
In TFPT sind Familien eine eigene Kohomologie H¹, hier liegen die Familiencharaktere im selben
ℂ⁵ wie die Plätze; und die drei Sektoren sind verschieden statt gleich. Die Zerlegung 1 + 2 der
Familien unter Galois ist eine prüfbare Vorhersage gegen die TFPT-Flavorstruktur, kein Befund.

## K: Familien als étale Kohomologie (Nachtrag, `k1_kohomologie.py`)

Test 1 aus Abschnitt F: H¹ von X = P¹ ohne die vier Marken über 𝔽̄_p, mit μ₄-Wirkung, gegen den
TFPT Familienraum (Rang 3, Charaktere {−1, i, −i}). Abbruch bei Abweichung.

**Befund (Punktzählung in explizit konstruierten Körpern 𝔽_{p⁴}, Lefschetz).**
- |X(𝔽_{p^k})| ergibt tr(Frobᵏ | H¹_c). Für p = 5 und 13 (p ≡ 1 mod 4) hat der Frobenius auf
  H¹_c die Eigenwerte **{1, 1, 1}**, auf H¹ also dreimal p (Gewicht 2). Für p = 3 und 7
  (p ≡ 3 mod 4) sind es {1, 1, −1}.
- Getwistete Zählungen bei p = 5: #Fix(Frob ∘ g_a) = 2, 6, 6, 6 für a = 1, 2, 4, 3. Die
  Fourierzerlegung über μ₄ gibt: trivialer Charakter **0**, jeder der drei nichttrivialen
  Charaktere **genau einmal** mit Frobenius 1. Dasselbe bei p = 13.
- **Residuenbrücke.** Die Residuen der TFPT Trägermoden e_m = z^m/(z⁴ − 1) an den vier Marken
  bilden die fünf Moden modulo Konstanten **surjektiv auf den Familienraum** ab (Rang 3). Der
  Kern ist span(e₃, e₄ − e₀), die Markenverschiebung ζ ↦ iζ hat auf dem Bild die Eigenwerte
  −1, i, −i. Über 𝔽₅ ist das Gegenstück die Einschränkung Fun(𝔽₅) → Fun(𝔽₅^×)/Konstanten mit
  Kern span(δ₀, 1).

**Urteil K: Test bestanden, aber nicht unabhängig.** Rang 3 und die Charaktere {−1, i, −i} folgen
schon aus der komplexen Geometrie von P¹ ohne vier Punkte; über 𝔽̄₅ ändert sich daran nichts.
**Neu** sind zwei Dinge:
1. **Arithmetisch gleiche Familien.** Genau wenn i in 𝔽_p liegt (p ≡ 1 mod 4, kleinstes p = 5),
   wirkt der Frobenius auf alle drei Familienklassen gleich. Die Familien unterscheiden sich
   dann nur durch ihre μ₄-Ladung (i, −1, −i), sonst durch nichts. Bei p ≡ 3 mod 4 trennt der
   Frobenius 2 + 1.
2. **Familien = Trägermoden modulo Invarianten.** Die Residuen machen aus den fünf Trägermoden
   genau den Familienraum; zwei Moden tragen keinen Familienanteil. Das beantwortet den Einwand
   aus F teilweise: Die Familien liegen nicht „im selben ℂ⁵“, sie sind dessen Residuenbild in
   H¹.

Weiter ohne Mechanismus: drei gleiche Kopien der 16. Der Befund legt eine Lesart nahe (Materie =
Familienraum ⊗ 16, die drei Kopien arithmetisch gleich und nur durch eine ℤ₄ Familienladung 1, 2, 3
unterschieden), aber diese Lesart ist gewählt, nicht hergeleitet.

## Kleinste entscheidende Tests für eine Folgesitzung

1. **ℤ₄ Familienladung gegen Flavor.** Wenn die Familien die μ₄-Ladungen i, −1, −i tragen, ist
   die Massenmatrix in Potenzen nach Ladungssumme geordnet. Prüfen, ob TFPTs Transportlängen L
   und Potenzen K (oder die Leptonkoeffizienten 16/7, 4/3, 7/6) mit Ladungen 1, 2, 3 mod 4
   verträglich sind. Abbruch, wenn keine Zuordnung der Ladungen die Ordnung der K trägt.
2. **Galois 1 + 2 gegen Flavor.** Prüfen, ob eine TFPT Flavorgröße (Leptonkoeffizienten
   16/7, 4/3, 7/6, Transportlängen) eine Familie gegen zwei auszeichnet, und ob es die zum
   Legendre-Charakter passende ist.
3. **Sechs Fasern.** Welches natürliche Element mit Nichtrest-Multiplikator (etwa J aus der
   Standardform kombiniert mit einer Dilatation) erweitert A₅ auf PGL(2,5) mit Stabilisator F₂₀?
4. **Y_par gegen die markierte Ladung.** Den Befund ‖[H*, Y]‖² = 19/12 gegen 17/6 mit der Weil-
   Lesart neu bewerten: Welche Dynamik ist mit der Paritätsbasis verträglich?
5. **Fockspur-Anschluss.** Wählt die positive Fockspur über Λ ℂ^(𝔽₅) mit der Weil-induzierten
   quadratischen Quelle (Exploration F1 vom 23.11.) dieselbe Orientierung Σ ∝ −6Y?
