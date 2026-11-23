# Z21: Was einer vollständigen Theorie aus TFPT plus H6 noch fehlt

Stand 23. November 2026, Explorationsmodus. Eine Karte, kein Ergebnis: Nichts hier ändert einen
Ledgerstatus, nichts wird befördert, keine Atlasänderung. Maßgeblich bleiben
`verification/status_ledger.csv` und die Originalquellen. H6-Aussagen stammen aus `VERDICT.md`
(Stand H6h12, z20 eingeschlossen) und `ausgabe_z20.txt`; sie sind Explorationsbefunde und stehen
nicht im Ledger.

Statuswörter wie im Ledger: **[E]** exakt im Geltungsbereich, **[N]** numerische Übereinstimmung,
**[C]** bedingte Identifikation (inkl. [P]), **[O]** offen (inkl. [A]), **[X]** Kill-Test angehängt.
H6-Wörter wie im VERDICT: **P** hergeleitet unter den Prämissen, **L** unter benannter Lesart.

Zählung und Klassenliste: `z21_ledger_zensus.py`, Ausgabe `ausgabe_z21.txt` (14/14 Prüfungen).

## 0. Kurzfassung

1. Der Dachvertrag **TFPT.TOE.COMPLETE.01 [O]** zerlegt „vollständig“ in acht Tore T1–T8. **Alle
   acht sind [O].** Kein Kind eines Tors ist im Ledger [E] geschlossen; die stärksten sind [C]
   (SEAM.EQUIV.MMST.01, DIMENSION.SELECTOR.4D.01).
2. H6 liefert Herleitungskandidaten fast ausschließlich für **T1** (Ursprung von P1) und **T2**
   (Naht = (E8)₁, auf Ebene der Zustandssumme), dazu perturbative Vorarbeit für T3/T5 (Lorentz ohne
   Feinabstimmung, induzierte Yang–Mills-Wirkung) und Zuflussrechnungen für T7. Zu T4 (chirale
   Materie ohne Spiegel), T6 (Kopplungen, Higgs, Neutrinotextur) und T8 (Anfangszustand) trägt H6
   nur Vorbedingungen bei (z8: ℤ₁₆-Bedingung; y1: sin²θ_W = 3/8 an der Gitterskala; z16:
   KMS-Zustand der Naht).
3. Die Bilanz „15 P, 5 L von 20“ des H6-Prüfstands misst überwiegend T1/T2-Struktur (15 der 20
   Kriterien: K1–K11, K14, K17, K19, K20). Sie ist kein Maß für TOE-Vollständigkeit.
4. Am meisten blockiert die Schicht **2D-Naht → lokale wechselwirkende 3+1D-Theorie** (T3/T5);
   direkt dahinter **chirale Materie ohne Spiegelgeneration** (T4). Beide sind Rechen- und
   Konstruktionsaufgaben, aber schwer.
5. Nicht durch Rechnung allein schließbar: die Einheit v_geo (No-Unit-Satz, nur Messung), die
   Rahmenprämissen (TFPT.IRREDUCIBLE.01 [O], H6-Prämissen A\*, C″, D, E, G), die H6-Entscheidung
   **D3 „Paarung nur auf ξ₊“** (exakte Obstruktion, keine Symmetrie trennt die Kanäle), die
   **Anfangszustandswahl T8** (Ledger: braucht eine nichtentropische Auswahl), die
   **Neutrino-Absolutskala** ohne Mechanismus, δ_PMNS (zugewiesener Texturwert) und die
   **Kosmologiegabel** (Daten).
6. Neue Querbefunde dieser Karte (Abschnitt 6), unter anderem: H6s Einzelnaht (z20, Deck intern)
   steht in Spannung zur Linsen-Geometrie S³/ℤ₄ von CCC.SEAM.CROSSOVER.01 (D1), und H6s
   sin²θ_W = 3/8 an der Gitterskala trifft auf den bekannten [X]-Befund, dass reine SM-Läufe nicht
   vereinigen (QFT4D.RGTEST.01, PS.RGTEST.01).

## 1. Ledgerzensus

789 Zeilen, 784 aktiv. Die offizielle Anzeigeklasse ist exakt aus
`verification/v649_discipline_audit.py` (`row_class`) nachgebaut; dort zählt [N] als E und jedes
Statusfeld mit „contract“ als O.

| Klasse (aktiv) | offiziell (v649) | Feinklasse „schwächstes Glied“ |
|---|---|---|
| E exakt | 192 | 166 |
| N numerisch | (in E enthalten) | 35 |
| C bedingt | 216 | 247 |
| O offen | 325 | 323 |
| X Kill-Test führend | 39 | (Merkmal: 77 Zeilen mit [X] oder Kill) |
| Axiom | 2 (AX.P1.01, AX.P2.01) | 2 |
| Regel (Firewall/Entscheidung) | in OTHER/E | 3 |
| sonstig | 10 | 8 |

Die Feinklasse ordnet zusammengesetzte Status nach ihrem schwächsten Teil (O < C < N < E). Sie ist
eine Lesehilfe für diese Karte, keine Statusänderung.

Von den 325 offenen Zeilen sind nur **38 rein offen** (Status beginnt mit „Open“ oder trägt nur
[O]/[A]); 287 sind gemischt: exakte Teilschritte mit einem benannten offenen Rest. Das Atlasinventar
zählt mit eigenem Vokabular (theory_overview: open 26, mixed 546, recorded_established 160,
conditional 34); es ist Navigation.

Zwei mechanische Fehllesungen, von Hand korrigiert (Ledgertext bleibt maßgeblich):
CCC.SEAM.KINEMATICS.01 ist Kinematik [E] mit CCC-Lesart [C]; das „[O]“ im Status gehört zu
WOIT.OS.TWISTOR.01. QFT4D.LATTICE.FUNDAMENTAL.01 ist eine typisierte Entscheidung; das „[E]“ steht
dort verneint und wird vom v649-Klassifikator als E gezählt.

### Die acht Tore (TFPT.TOE.COMPLETE.01) und die F1-Tore

| Tor | Inhalt | Kinder (Ledger) | Stand | F1-Tor (`research/originlab/gates.json`) |
|---|---|---|---|---|
| T1 | Strukturpostulat → P1/P2/Compiler | AX.P1.01, AX.P2.01 (Axiom), TFPT.IRREDUCIBLE.01 [E]/[O], DIMENSION.SELECTOR.4D.01 [C] | O | G0–G2 |
| T2 | Naht = (E8)₁ bewiesen | SEAM.EQUIV.01 [O], SEAM.EQUIV.MMST.01 [C], SEAM.EQUIV.TWISTOR.01 [O], SEAM.SIMPLECURRENT.GENERATOR.01 [O], SEAM.DETLINE.UNIFICATION.01 [O], SEAM.STATE.RPMIXING.01 [E]+[O] | O | G1, G3 |
| T3 | quasilokale unitäre 3+1D-Hamiltonfamilie | QFT4D.LATTICE.FUNDAMENTAL.01 (Entscheidung), TFPT4D.LATTICE.ACTION.01 [O], DYN.UNITARY.DILATION.01 [O], DYN.MARKOV.EMBED.01 [E]+[O] | O | G4 |
| T4 | chirales SM mit lokalem Maß, Spiegellücke | CHIRAL4D.NOMIRROR.01 [O] | O | G5 |
| T5 | IR-Lorentz, Confinement, Clustering, Streuung | QFT4D.OS.RECON.01 [O], SEAM.BULK4D.RECON.01 [O] | O | G4 |
| T6 | drei Kopplungen + Neutrinotextur intern | ALPHA.QUILLEN.EXACT.01 [O], GAUGE.DETLINE.FIXPOINT.01 [O], FLAV.NUSCALE.01–.06 (Kandidaten [C]/[N], Mechanismus [O], .01–.03 mit [X]-Anteil) | O | – |
| T7 | quantisierter masseloser Spin 2, universelle Kopplung | GRAV.SPIN2.EMERGENCE.01 [O], GRAV.NONCIRCULAR.01 (Firewall), GRAV.WEINBERG.WITTEN.01 [O] | O | G6 |
| T8 | eindeutiger Anfangszustand + ein SK-Funktional | FTRANSFER.SK.RHO0.01 [O], FTRANSFER.GENERATING.01 [O], OBS.TRANSDUCTION.01 [O], PRED.JOINTLIKELIHOOD.01 [O] | O | G7 |

Ledgerwortlaut zum Gesamtstand (Runde 4, v1026–v1030): „All T1–T8 remain [O]; no shared complete
3+1D parent.“ RESIDUAL.CERTIFICATION.01 stellt klar, dass das historische Residuenaudit keine
TOE-Vollständigkeit bescheinigt. Für F1 gilt: `f1_native_contract.json` ist `NOT_CONFIGURED`; ein
H6-Kandidatenvertrag nach `research/originlab/F1-CANDIDATE-CONTRACT.md` existiert noch nicht.

## 2. Was H6 zu den Toren beiträgt

H6-Prämissen: A\* (hintergrundfreie Naht), C″ (Familien = chiraler Index der vier Marken),
D (Gravitation als Zustandsgleichung), E (3+2 lokal), G (Naht = verdrillter Defekt eines
Paarungsmassenfelds).

| H6-Befund | Skript | H6-Stufe | trifft Tor | Grenze |
|---|---|---|---|---|
| GSO-Summe aus großen Diffeos, 16 MW, E8 c = 8, (P,P) flach | v1 (H6h), z1 | P | T1, T2 (Zustandssumme) | Netzebene (128 Spinorströme als Netzaussage) offen |
| ⟨T⟩ = 0 ⇒ E₆(τ) = 0 ⇒ τ = i; Harmonizitätsverbot für Minimierung | z3 | P | T1 | Lean: nur E₆(i) = 0, Valenzformel fehlt |
| τ = i stabil; S(0) = log 12; KMS bei β = L | z4 | P | T1 | Naht-Zustand, nicht 4D |
| c₃ = 1/(8π) je Spinorsektor ohne Fixpunktkorrektur | z1 | P | T1 | Viertelzählung bleibt Lesart der Wärmespur |
| Familiencharaktere {−1, i, −i}, 3+2, H⁰(K(D)) = 3 | z1, z10, z16 | P | T1 (P2 aus P1) | physische Kopplung an K^{1/2}⊗O(D) offen |
| Fork c = 8 gegen 24, 3D-Wand unsichtbar bei c ∈ 8ℤ | z5 | P (gegen Operatorlesart) | T2 | c = 24 im Mittel nicht ausgeschlossen |
| Lorentz ohne Feinabstimmung (Dim. 6), RP über T² | z6 | P, perturbativ | T3/T5 | nichtperturbativ, Confinement, Higgs: nicht gezeigt |
| Induzierte YM-Wirkung, sin²θ_W = 3/8 an der Gitterskala | y1 | P im Modell | T6 (teilweise) | Schwellen nicht behandelt; siehe 6.d |
| Naht als verdrillter Defekt (G), 1 chirale Mode je Paar | z7 | P | T1 | freies Modell |
| 16 aus ℤ₁₆, X ≡ 1 mod 4 für die SM-Generation mit ν^c | z8 | P (partnerfreie Fassung) | T4 (Voraussetzung SMG) | Spiegelentkopplung nicht gerechnet |
| Zufluss SO(2)-Ladung; Wilson bricht Phase ⇒ Wand | z9 | P/N | T2, D4 | Grundzustand erzwingt Defekt nicht |
| Einheitswindung = relative Chernzahl; Nähte paarweise unter G | z11, z16 | L / gescheitert als Zwang | T1 | – |
| Elternoperator: ein Multiplett scheitert, zwei als Lesart | z12 | L | T1, T4 | Aufteilung eingebaut |
| Gravitative Anomalie c = 8, h_R = 1, Wand trägt I₄ | z13 | N/E | T7 (nur Anomaliebilanz) | Fluss in Krümmung nicht gerechnet |
| Satzzahl: gleichförmig n = 3 widerspricht α (EM1) | z14 | N/C | T1, T6 | verborgen n = 3 offen, c₋/k entscheidet |
| a/ℓ_P ≈ 10,8 im H6-Schema, Spanne 3,3–14,3; robust ≥ 1,65 | z2, z15, z16 | N | T7 | reglerabhängig |
| Gibbs-Variationsprinzip für den KMS-Zustand | z16 | P | T8 (nur Naht) | siehe 6.e |
| Sektor 1: keine Grundzustandsregel; Lesarten | z17 | L | T1 (D1, D2) | – |
| Gemeinsamer Operator auf dem Nahttorus; Obstruktion 0, ½, ½ | z18 | L mit Obstruktion | T1, T4 (D3) | – |
| GW-Paarung: Wand oder Positivitätsverletzung | z19 | N | T3, T4 (D4) | Lüschers projizierte Weylfelder nicht geprüft |
| Deck-interner Zweig: physische Deckholonomie auf S³ erzwingt ≥ 1 Schlauch; minimal genau einer | z20 | L, gestärkt | T1 (D1) | Eingaben: Deck intern, Holonomie physisch, Reziprokes nicht erzwungen |

Im Deck-internen Zweig (z20, H6h12) wird G selbst hergeleitet: Φ trägt J-Ladung 2, wechselt um den
erzwungenen Flussschlauch das Vorzeichen und muss dort verschwinden. Die Prämisse G ist dort ersetzt
durch drei Eingaben: Deck intern, Deckholonomie physisch, Reziprokes nicht von der Konfiguration
verlangt. Die spontane Wahl eines Pols ist laut H6h12 eine Wahl zwischen den beiden Äonen-Enden,
also eine Zeitrichtung.

**Die vier offenen H6-Entscheidungen** (H6h11 plus H6h12):

| | Entscheidung | Stand | Was sie entscheiden könnte | Art |
|---|---|---|---|---|
| D1 | genau eine Naht | z17: keine Grundzustandsregel; z20: im Deck-internen Zweig von v957 minimal genau eine, Reziprokes spontan gebrochen; im räumlichen Zweig (Linse) unmöglich | Prinzip „Deck wirkt intern, Holonomie physisch“, oder Datum N_ν (LEP, k = 1 bei +2σ) | Prinzip oder Daten |
| D2 | genau ein E8-Block | nur die Forderung „E8 auf der Naht“ wählt; sonst D₂₄⁺ = so(48)₁ (z17); gleichförmig n = 3 durch α unter EM1 ausgeschlossen (z14) | c₋/k der Naht, also ihre eigene T-Phase, in einem wechselwirkenden Modell | Rechnung (nach S2) plus Lesart |
| D3 | Paarung nur Φ ξ₊ξ₊ | exakte Obstruktion (z18): Nahtkanal liegt in den Familienkanälen, Profil eindeutig s ∝ q; mit ψχ aus G Familienmassen 0, ½, ½ × Nahtmasse; keine Symmetrie trennt | eine neue Symmetrie oder Auswahlregel | Prinzip |
| D4 | Wand gegen Positivität | z19: Wilson mit Wand und RP, oder flache GW-Formen mit RP-Verletzung bis −0,23 | die registrierte Entscheidung QFT4D.LATTICE.FUNDAMENTAL.01 (T = e^{−aH} positiv per Konstruktion) wählt als Lesart die Wand | Entscheidung plus Folgerechnung |

## 3. Lückenkarte nach Schichten

Jede Schicht: was existiert, was fehlt, was H6 beiträgt, und **die eine entscheidende nächste
Rechnung oder der nächste Beweis**.

### S0 Ursprung der Eingaben P1, P2 und der Einheit (T1; F1 G0–G2)

- **Existiert:** AX.P1.01 und AX.P2.01 sind deklarierte Axiome mit Reduktion auf Anker
  a = (1,1,2) plus π. P2.PARTITION.01 [E] mit [C]-Gewichtstypisierung: g_car = 5, N_fam = 3 aus
  der Partition von 4. TFPT.IRREDUCIBLE.01 [E]/[O]: irreduzibel bleiben „Rahmen + π“.
  DIMENSION.SELECTOR.4D.01 [C]: d = 4 eindeutig unter A1–A4, deren Herkunft aus dem Compiler offen
  ist. QGEO.REALIZE.01 [C]/[O]: die Naht als P¹ ohne μ₄ ist Realisierungsprämisse.
  Einheit: ANCHOR.VGEO.01 [I]+[A], VGEO.SHARPEN.01 [E]/[O], METROLOGY.CLOSURE.01 [E]/[N]/[O]:
  v_geo ist durch den No-Unit-Satz verbotene Ausgabe.
- **Fehlt:** ein Prinzip, das die Naht als physisches Objekt im Grundzustand erzwingt (K19 bei
  H6 nur L), die primitive Ladung bzw. Sektor 1 (D1), ein Block (D2), der gemeinsame Operator für
  Naht und leichte Familien (D3, K20 L), die Kopplung physischer Fermionen an K^{1/2}⊗O(D), und
  die Eindeutigkeit von τ = i in Lean (Valenzformel).
- **H6:** leitet unter A\*, C″, E, G die Marken, τ = i, c₃, μ₄, E8 (Zustandssumme), 3+2 und die
  Familiencharaktere her (z1, z3, z4, z10, z16). z20 schließt D1 bis auf drei benannte Eingaben.
  H6 liefert außerdem eine Herkunft für die Axiome A3/A4 von DIMENSION.SELECTOR.4D.01 (Λ²₊ =
  Λ²(Familien), drei chirale Mittelformen nur in d = 4), als Lesart unter C″.
- **Nächste entscheidende Rechnung: D3.** Endliche Symmetriesuche am z18-Operator: alle
  unitären und antiunitären Operatoren, die mit dem Dirac an K_E(R), dem geeichten Deck
  σ̃⊗J⁻¹, dem Verflechter ι, der μ₄-Uhr und der Eichgruppe vertauschen; Test, ob einer den
  Nahtkanal q = f₂² − f₀² von den Familienkanälen trennt, so dass Φ ξ₊ξ₊ erlaubt und ψχ verboten
  ist. Abbruch nach vollständiger Aufzählung der endlichen Gruppe. Ergebnis ist entweder eine
  Herleitung von D3 oder ein No-Go, nach dem der Ein-Multiplett-Elternoperator für leichte
  Familien scheitert. Das ist die schärfste Stelle, an der H6 empirisch falsch werden kann
  (Familienmassen 0, ½, ½ × Nahtmasse).

### S1 Naht = (E8)₁ als physischer Randzustand (T2; F1 G1, G3)

- **Existiert:** SEAM.EQUIV.01 [E]+[O] konzentriert alles auf einen benannten Satz.
  SEAM.EQUIV.MMST.01 [C]: geschlossen modulo zitierter MMST/OS-Sätze. SEAM.EQUIV.UNIQ.01 [E]/[C]:
  Holomorphie erzwingt E8 unter den drei c = 8-Kandidaten. SEAM.SIMPLECURRENT.GENERATOR.01 [O]:
  offener Rest C7 (modulare Invarianz) und G2–G4 auf Netzebene, Kern N1: die 128 Spinorströme sind
  keine Bilineare. ALPHA.QUILLEN.EXACT.01 [E]/[O] als Fläche davon (α als Quillen-Stationarität).
- **Fehlt:** die μ₄/GSO-Erweiterung SO(16)₁ → (E8)₁ als Netzaussage; der Twist- bzw.
  Spinfeldsektor (N2); die Quillen-Variationsexaktheit.
- **H6:** leitet die GSO-Summe aus großen Diffeomorphismen her (A\*), also C7 auf Ebene der
  Zustandssumme, und liefert in z7/z13 ein explizites Gittermodell, dessen Defekt c = 8 trägt.
  Wichtig: Nach DIMENSION.UPLIFT.FIREWALL.01 (F1) schließt auch ein vollständiges T2 nichts in 3+1D.
- **Nächste entscheidende Rechnung:** im z7-Defektmodell die Spinfelder (Twistoperatoren) der
  Naht konstruieren und ihre Zweipunktfunktionen messen: Gewicht h = 1 für die 128 (Twistgewicht
  aus z10) und Schließen der Ströme 120 + 128 zu (E8)₁ (OPE-Strukturkonstanten). Das ist ein
  H6-eigener Weg zu N1/N2 auf Gitterebene.

### S2 2D-Naht → lokale wechselwirkende 3+1D-Theorie (T3, T5; F1 G4)

- **Existiert:** DIMENSION.UPLIFT.FIREWALL.01 (Regel): kein 1+1D-Satz schließt eine 3+1D-Aussage;
  Reihenfolge Naht als Randmodul → echte 4D-Gittertheorie → OS → chirale Materie.
  SEAM.BULK4D.RECON.01 [O] (B1–B7, „nichts im Korpus schließt einen Teil“; B7 per v1020 als
  unerfüllbar umformuliert zu BULK.PRIME_COMPLETION.01 [O]). WOIT.OS.TWISTOR.01 [O] mit
  Spielzeug-Kill (RP bricht für g > 0, SEAM.INT.FKTOY.01). QFT4D.OS.RECON.01 [O] (O1–O6;
  perturbative Formulierungsregel). QFT4D.LATTICE.FUNDAMENTAL.01 (Entscheidung): das endliche
  lokale unitäre Gitter ist das physische Kriterium, der OS-Limes mathematische Verstärkung.
  TFPT4D.LATTICE.ACTION.01 [O]: das eine explizite Objekt, an dem fast alles testbar wird.
  DYN.UNITARY.DILATION.01 [O], DYN.MARKOV.EMBED.01 [E]+[O], QFT4D.EG.ALLORDER.01 [C].
  CCC.SEAM.KINEMATICS.01 [E] (Kinematik, ausdrücklich keine Dynamik), CCC.SEAM.CROSSOVER.01 [O].
- **Fehlt:** eine explizite Hamiltonfamilie in 3+1D mit dynamischen Eichfeldern, Fermionen und
  Skalaren, die T3 erfüllt (Gaußgesetz, lokales hermitesches H, positiver Transfer,
  Lieb–Robinson); dann T5: nichtperturbative IR-Lorentz-Invarianz, Confinement, Clustering,
  nichttriviale Streuung.
- **H6:** H6 kehrt die Richtung von SEAM.BULK4D.RECON.01 um. Nicht die Naht erzeugt das Volumen,
  sondern die Naht ist ein Defekt von 4D-Feldern mit Zufluss (z7, z9, z13). Damit folgt H6 der
  Routenordnung der Firewall (Naht als Randmodul, Volumen als eigene Gittertheorie). Beiträge: das
  isotrope Lipschitzgitter ℤ⁴ gibt Lorentz ohne Feinabstimmung bis Dimension 6 (z6, perturbativ);
  RP über Doppelschritte; Wilson r ≤ 1 linkpositiv (z19); induzierte YM-Wirkung (y1). Nicht
  geleistet: dynamische Eichfelder, Wechselwirkung jenseits einer Schleife, Higgs-Sektor,
  Confinement.
- **Nächste entscheidende Konstruktion:** TFPT4D.LATTICE.ACTION.01 mit H6-Inhalt explizit
  aufschreiben, als Hamiltonroute: Lipschitzgitter, Kogut–Susskind-Links für SU(3)×SU(2)×U(1),
  Wilson-Fermionen mit r ≤ 1 (Wandzweig von D4), die 16 schweren Paare mit komplexem Φ
  (ℤ₄-Ladung 2), Defektsektor. Dann die T3-Liste exakt prüfen. Erster T5-Zeuge: gemeinsame
  Lichtgeschwindigkeit von Fermionen und dynamischen Eichanregungen in kleinem Volumen (v1012:
  „common-c needs dynamical gauge excitations“).

### S3 Chirale Eichtheorie und Fermioninhalt (T4; F1 G5)

- **Existiert:** CAR.SM.ANOM.01 [E]+[C]/[O]: Λ^gerade ℂ⁵ = eine anomaliefreie Generation, SM-RG;
  S-Matrix und Yukawa offen. CHIRAL4D.NOMIRROR.01 [O] (C1–C6; SMG-Mechanismus im Spielzeug [E];
  4D-Spin(10)-Satz [O]). CHIRAL4D.MIRROR.DFILTER.BLOCK.01 [E] Blocklemma, GW-Lesart [C].
  CHIRAL4D.MIRROR.FACTORIZED.BOUND.01 [E] mit Obstruktion. PS.NCG.FLUCT.01 [E]/[C]: innere
  Fluktuationen ergeben das SM-Higgsdublett. CONTRACT.QFT4D.DIRAC.01 [O].
- **Fehlt:** ein lokales eichinvariantes chirales Fermionmaß für den TFPT-Inhalt, globale
  Anomaliefreiheit samt Witten-SU(2), Spiegelentkopplung mit gleichmäßiger Lücke im Kontinuum,
  richtiger Index für N_fam = 3.
- **H6:** z8 zeigt die notwendige Bedingung für symmetrische Massenerzeugung der Spiegel:
  X = 5(B−L) − 4Y ≡ 1 mod 4 für alle 16 Weylfelder mit ν^c, also ℤ₁₆-anomaliefrei. z19 zeigt,
  dass die Wahl des chiralen Reglers mit D4 verflochten ist (Overlap bei θ = 0 positiv; flache
  Formen verletzen RP bei gedrehter Phase). Die eigentliche T4-Konstruktion fehlt in H6.
- **Nächste entscheidende Rechnung:** RP für Lüschers projizierte Weyl-Overlapfelder im
  abelschen Fall bei endlichem Gitterabstand (frei, numerisch machbar; z19 hat es offen gelassen).
  Sie entscheidet, ob ein chiraler Regler positiv sein kann, und damit D4 und den Weg zu C1.
  Die fundamentale Folgefrage ist danach die wechselwirkende Spiegellücke ohne Symmetriebruch
  (SMG) für eine 16 in 3+1D.

### S4 Kopplungen, Higgs/EW-Brechung, Yukawa, Flavor, starkes CP (T6)

- **Existiert, Kopplungen:** α-Wurzel global eindeutig (TFPT-STAND, [E] für Wurzel, [C] für
  EM-Lesart); ALPHA.QUILLEN.EXACT.01 [O]; GAUGE.DETLINE.FIXPOINT.01 [O]: α_s(M_Z) ist externe
  Eingabe. FIRSTPRINCIPLES.BOUNDARY.01 [E]/[X]/[C]: α_s extern, weil das SM nicht vereinigt
  (QFT4D.RGTEST.01 [X]: Spreizung Δα⁻¹ ≈ 8,8 bei 2·10¹⁶ GeV); PS.RGTEST.01 [X]/[O] und
  PS.THRESHOLD.02 [E]/[C]: der Pati–Salam-Zweig vereinigt mit M_PS auf der Skalaronskala;
  PS.PROTON.02 [X]/[C]: Zweig B nur mit +(15,1,1), Hyper-K entscheidet.
- **Existiert, Higgs/EW:** PS.NCG.FLUCT.01 [E]/[C] (Inhalt); HIGGS.FREESEAM.01 [I]+[P]:
  λ = β_λ = 0 an der Naht gibt ein m_H-Band 129–134 GeV, gemessen 125,25 GeV liegt darunter
  (Metastabilität); SM.RG.01 [C]. Ein dynamischer Mechanismus der EW-Brechung (Potential,
  Vorzeichen, Hierarchie als Rechnung in einer 4D-Theorie) ist nicht hergeleitet;
  TFPT4D.LATTICE.ACTION.01 setzt EWSB „downstream“. FLAV.UPOINT.01 [I]+[A] führt die absolute
  Flavoramplitude auf v_geo zurück.
- **Existiert, QCD:** QCD.LAMBDA.01 [I]+[P]+[A]; FR.MPME.01 [O]; FR.QCD.BUDGET.01 [C]/[N]/[X]
  (C_p aus Gitter-QCD extern).
- **Existiert, Flavor/Neutrinos:** Winkel FLAV.PMNS.03 [C]/[O], FLAV.PMNS.CLOSE.01 [E]/[C]/[O];
  θ₁₃: FLAV.TH13.01 [E]/[C]/[O], FLAV.TH13.PRESSURE.01 [N]/[X] (+2,0σ); δ_PMNS = 240° und
  sin²θ₂₃ = ½ sind in `predictions_frozen.json` Schicht „assigned“, keine Herleitung; CKM-Phase
  FLAV.CP.01 [N]+[P] (Phasenmechanismus offen). Absolutskala FLAV.NUSCALE.01 [E]/[N]/[C]/[O]/[X]
  (ein irreduzibler UV-Input m₃ = (y_ν v)²/M_R), Kandidaten FLAV.NUSCALE.05/.06 [C]/[N] mit
  Mechanismus [O], FLAV.NU.TEXTURE.MECHANISM.01 [O] (beide strukturellen Nullbefunde bindend),
  FLAV.NU.PARENT.INTERTWINER.01 [E/N] endlich. Koide FR.KOIDE.03 [C] (Vermutung 53/54),
  FR.POLE.SOLVE.01 [E]/[C].
- **Existiert, starkes CP:** FLAV.STRONGCP.02 [F]: θ_eff = 0 als Modellsatz im selbstdualen
  16-Majorana-Nahtmodell (Pfaffian reell, RP schließt θ ≠ 0 aus). Geltungsbereich Nahtmodell; nach
  Firewall F1 nicht auf das 4D-θ̄ = θ − arg det M_q übertragbar ohne T4 und Yukawaphasen.
- **Fehlt:** α_s ohne Messung; EW-Brechung als Dynamik; Yukawas als Ausgabe einer 4D-Theorie;
  CP-Phasen hergeleitet; M_R und absolute Neutrinomassen mit Mechanismus; θ̄ in 4D.
- **H6:** y1 gibt sin²θ_W = 3/8 an der Gitterskala (gleiche Dynkinsumme je Generation). C″ mit
  S₄-Restsymmetrien gibt tribimaximale Mischung in führender Ordnung; θ₁₃ und δ sind aus TFPT
  geerbt (K18 L). z14: α schließt gleichförmig n = 3 unter EM1 aus. H6 enthält kein EW-Higgsfeld:
  Φ ist ein SO(16)-Singulett mit ℤ₄-Ladung 2.
- **Nächste entscheidende Rechnung:** RG-Verträglichkeit von H6. Zwei Schleifen von M_Z bis zur
  H6-Gitterskala (a/ℓ_P im Fenster 1,65 bis 14,3) mit der Randbedingung gleicher Kopplungen aus
  y1, einmal mit reinem SM-Inhalt (Erwartung nach QFT4D.RGTEST.01: Kill), einmal mit dem
  PS-Zweig aus PS.THRESHOLD.02. Ergebnis: entweder benennt H6 den Zwischeninhalt (dann
  Protonzerfall bei Hyper-K als Vorhersage), oder die Lesart „induziert an der Gitterskala“ fällt.
  Machbar mit vorhandenen RGEs.

### S5 Gravitation jenseits Jacobson, Quantengravitation (T7; F1 G6)

- **Existiert:** GRAV.ENTROPY.EQUILIBRIUM.01 und GRAV.NONLINEAR.01 [E]/[C]/[O]: lokale
  Einsteingleichung mit 8π = 1/c₃ und Λ-Vorfaktor; nach GRAV.NONCIRCULAR.01 (Firewall) sind
  BW/CHM dabei Annahmen, bedingt auf QFT4D.OS.RECON.01. GATE.METRIC.01 [O]: nur das globale Maß;
  QG.AMB.01 [O]: allgemeines Konformfaktorproblem, nicht TFPT-spezifisch, von allen Ausgaben durch
  Δ_eff ≈ 1,648 > 0 entkoppelt. GRAV.SPIN2.EMERGENCE.01 [O], GRAV.WEINBERG.WITTEN.01 [O] (bindend
  vor jedem Spin-2-Schritt), GRAV.PARENT.GAP.NO_GO.01 [E]: ein gleichmäßig gelückter Elter schließt
  einen masselosen Spin-2-Pol aus. GRAVITY.COMPLETE.01 [E]/[C]/[O].
- **Fehlt:** ein masseloser transversaler Spin-2-Pol mit zwei Helizitäten, universeller Kopplung
  und Diffeo-Wardidentitäten aus derselben Theorie; die ausdrückliche Angabe, welche
  Weinberg–Witten-Prämisse umgangen wird; die Flächennormierung η = 1/4 aus dem Feldinhalt.
- **H6:** D ist Prämisse (Jacobson). z13: gravitative Anomalie c = 8, h_R = 1, Wand trägt I₄;
  Fluss in Krümmung nicht gerechnet. z2/z15/z19: a/ℓ_P nicht festgelegt, nur ≥ 1,65 robust. Unter
  D ist η = 1/4 die Definition der Planckeinheit (z2), H6 prüft 8π = 1/c₃ also nicht unabhängig.
  H6s Gitter hat keinen fundamentalen Lorentz-kovarianten T_μν; das ist die Ausgangslage, in der
  Weinberg–Witten umgangen werden kann, aber H6 erzeugt unter D gar kein Graviton.
- **Nächste entscheidende Rechnung:** Sakharov gegen Jacobson im selben Regler. Für den H6-Inhalt
  den induzierten Einstein–Hilbert-Koeffizienten (Wärmekern a₁ einer Schleife, inklusive
  Kontaktterme der Eichfelder) und die Verschränkungsflächendichte κ mit demselben Gitterregler
  rechnen. Stimmen 1/(4G_ind) und κ/a² überein, wird a/ℓ_P eine Vorhersage; stimmen sie nicht,
  ist D mit H6-Inhalt widersprüchlich. Danach erst die TT-Spektraldichte (T7).

### S6 Anfangszustand und Kosmologie (T8)

- **Existiert:** FTRANSFER.SK.RHO0.01 [O] (globales ρ₀ und SK-Eindeutigkeit offen);
  FTRANSFER.GENERATING.01 [O]; v1012 typisiert entropische Auswahl als leer („NO_SET_WORKS“,
  verlangt eine nichtentropische Auswahl). Λ: LAMBDA.BRANCH.01 [I]+[N]+[P]
  (ρ_Λ/M̄⁴ = (3/4π²)e^{−2α⁻¹}), GRAV.H0.01 [N]+[C] (H₀ 66,5–67,1, Hubble-Spannung nicht gelöst).
  Inflation: COSMO.NSTAR.01 [P], GATE.NSTAR.01 [C] (N\* ist Eingabe), AS.RECONCILE.01 [X]:
  langsamer Higgs-Kanal gibt A_s bei −11,4σ. Kosmologiegabel CCC.SEAM.CROSSOVER.01 D3 [O]
  (Starobinsky r = 12/N² gegen Crossover gegen Hybrid). Dunkle Materie: FR.DM.01 [O]
  (Reliktskala), FR.DM.02 [C] (post-inflationär, f_a = M_scal/128, m_a ≈ 23,8 μeV),
  FR.RELIC.SOLVE.01 [C]/[X]. Baryogenese: FR.ETAB.01 [C], FR.ETAB.04 [C]/[X],
  FR.BOLTZMANN.SOLVE.01 [C]/[X] (Faktor ≈ 1,1 bei M₁ = M_scal φ₀²/10); Ω_b: COSMO.OMB.01 [N].
- **Fehlt:** ρ₀ mit Auswahlprinzip; der Mechanismus für Λ in einer 4D-Theorie (die
  Zahlübereinstimmung ist [N]/[C]); Wiedererwärmungskanal und damit N\*; Entscheidung der
  Kosmologiegabel; Verträglichkeit von M₁ ≈ 8,65·10⁹ GeV mit dem Seesaw-Spektrum (M_R ≈ 9·10¹³
  GeV in FLAV.NUSCALE.05).
- **H6:** H6 liefert für die Naht den KMS-Zustand als Gibbs-Maximum (z16) und Energie null bei
  τ = i (z3/z4). Beides ist 2D und entropisch; siehe 6.e. Zur 4D-Kosmologie trägt H6 nichts bei.
- **Nächste entscheidende Rechnung:** den Wiedererwärmungskanal des Skalarons aus dem
  TFPT-Inhalt herleiten (welche Zerfallskanäle, N_eff), weil er N\* und damit A_s, n_s und r
  festlegt. Bleibt nur der langsame Higgs-Kanal, steht A_s bei −11,4σ, also ein Kill-Kandidat.
  Entschieden wird die Gabel aber durch Daten (LiteBIRD/CMB-S4, σ(r) ≈ 0,001).

### S7 Datenprüfung und Statistik (T8/Validated; F1 G7)

- **Existiert:** KILL.BOARD.01 [N] (acht Tests: JUNO θ₁₂ und θ₁₃, r, w, δ = 240°, Ordnung und
  m_ββ, n2EDM θ_eff = 0, Doppelbrechung β), PRED.LIVE.01 [N]/[C]/[X], PS.PROTON.02 [X]/[C],
  PRED.JOINTLIKELIHOOD.01 [O] (v1 [N] ausgeführt), OBS.TRANSDUCTION.01 [O].
- **Fehlt:** die gemeinsame Likelihood in Relationsbündeln (die Treffer hängen an denselben Atomen
  φ₀, c₃, g_car, μ₄ und sind korreliert); eigene, eingefrorene Vorhersagen von H6.
- **H6:** H6 erbt TFPTs Vorhersagen. Eigene niederenergetische Vorhersagen hat H6 bisher nicht;
  seine eigenen Zahlen (a/ℓ_P, Lorentzverletzung ab 2–7·10¹⁸ GeV, sin²θ_W = 3/8 an der
  Gitterskala) sind nicht direkt messbar oder reglerabhängig.
- **Nächster Schritt:** PRED.JOINTLIKELIHOOD.01 in Bündeln ausführen und für H6 mindestens eine
  unabhängige, eingefrorene Vorhersage registrieren, bevor neue Daten kommen (JUNO θ₁₃ ist der
  stärkste Druckpunkt: +2,0σ gegen NuFIT 6.0, Daya Bay 2,17σ laut z16).

## 4. Rangfolge

Fundamentalität = wie viele Tore oder Schichten blockiert werden. Machbarkeit = ob der nächste
Schritt mit vorhandenen Mitteln in einer Sitzung oder wenigen erreichbar ist.

| Rang | Lücke | blockiert | Machbarkeit | Art |
|---|---|---|---|---|
| 1 | S2: lokale wechselwirkende 3+1D-Theorie (T3/T5) | T4, T5, T6 (α_s, Higgs, Yukawa), T7 (Firewall verlangt 4D-Netz), T8, 4D-θ̄, Λ-Mechanismus | niedrig für das Ganze; mittel für die explizite H6-Hamiltonfamilie | Konstruktion |
| 2 | S3: chirale Materie ohne Spiegel (T4) | T5-S-Matrix, EW-Brechung, Yukawas, 4D-θ̄, Sphaleronen/Baryogenese | niedrig (offenes Problem der Gittertheorie); RP-Test Lüscher hoch | Konstruktion |
| 3 | S0: Ursprung von P1, D1–D3 (T1) | den Anspruch „fundamental“ und die Familienkopplung; physikalisch kann P1 gesetzt bleiben | hoch für D3-Symmetriesuche und Lean-Valenzformel; D1 durch z20 weitgehend benannt | Rechnung plus Prinzip |
| 4 | S5: Spin 2 und Quantengravitation (T7) | nur die Gravitation, aber vollständig | niedrig; Sakharov/Jacobson-Konsistenz mittel | Prinzip plus Rechnung |
| 5 | S1: Naht = (E8)₁ auf Netzebene (T2) | α- und c₃-Lesarten; keine 4D-Aussage (Firewall) | mittel | Beweis |
| 6 | S4: Kopplungen, Higgs, Flavor (T6) | quantitative Vorhersagen | hoch für den H6-RG-Test; niedrig für Mechanismen | Rechnung, teils Prinzip |
| 7 | S6: Anfangszustand, Kosmologie (T8) | kosmologische Ausgaben | mittel (Wiedererwärmung); sonst datengetrieben | Prinzip plus Daten |
| 8 | S7: Statistik und Daten | Validierung | hoch | Rechnung plus Daten |

Nach Machbarkeit geordnet, mit dem größten Informationsgewinn je Aufwand: (a) D3-Symmetriesuche
(S0), (b) H6-RG-Verträglichkeit (S4), (c) RP für Lüschers Weylfelder (S3, D4), (d)
Sakharov-gegen-Jacobson im selben Regler (S5), (e) explizite H6-Hamiltonfamilie (S2). (a), (b)
und (c) können H6 jeweils widerlegen.

## 5. Was Rechnung allein nicht schließt

| Lücke | warum nicht | was es bräuchte |
|---|---|---|
| Absolute Skala v_geo | No-Unit-Satz [E] (METROLOGY.CLOSURE.01, VGEO.SHARPEN.01): dimensionslose Daten geben keine Masse | eine Messung (jede Theorie braucht eine Einheit) |
| Rahmen bzw. H6-Prämissen A\*, C″, D, E, G | TFPT.IRREDUCIBLE.01 [O]: „diskreter RP-Rand mit endlichem Träger“ ist Grundhaltung, nicht ableitbar | ein Prinzip; prüfbar nur über Folgen |
| D1 genau eine Naht | z17: keine Grundzustandsregel; z20: erzwungen nur unter „Deck intern, Holonomie physisch, Reziprokes nicht verlangt“; k = 1 bisher empirisch über N_ν (LEP) | ein Prinzip für die drei Eingaben oder das Datum N_ν |
| D3 Paarung nur auf ξ₊ | exakte Obstruktion in z18, keine bekannte Symmetrie trennt die Kanäle | eine neue Symmetrie oder Auswahlregel (die D3-Suche entscheidet, ob es eine gibt) |
| D2 genau ein E8-Block | wählt nur „E8 auf der Naht“; c₋/k erst in einem wechselwirkenden Modell rechenbar | Prinzip jetzt, Rechnung nach S2; α ist bereits als Datum benutzt |
| T7 Quantengraviton | unter D gibt es keines; GRAV.WEINBERG.WITTEN.01 verlangt eine ausdrücklich genannte umgangene Prämisse | Entscheidung zwischen D (thermodynamisch, T7 neu typisieren) und mikroskopischem Spin 2 |
| T8 Anfangszustand | v1012: entropische Auswahl ist leer, eine nichtentropische Auswahl fehlt | ein Auswahlprinzip |
| Neutrino-Absolutskala, M_R | FLAV.NUSCALE.01 [E]: ein irreduzibler UV-Input ohne Mechanismus; Mechanismus [O] | Prinzip (Mechanismus) oder Daten (Σm_ν, m_ββ, KATRIN) |
| δ_PMNS = 240°, sin²θ₂₃ = ½ | Schicht „assigned“ in `predictions_frozen.json` | Herleitung (Prinzip) oder DUNE/Hyper-K |
| Kosmologiegabel, N\* | GATE.NSTAR.01: N\* ist Eingabe; CCC D3 offen | Wiedererwärmung herleiten und Daten (r) |
| α_s(M_Z) | FIRSTPRINCIPLES.BOUNDARY.01 [X]: SM vereinigt nicht | Zwischeninhalt als Prinzip (PS-Zweig) plus Daten (Hyper-K) |
| globales Metrikmaß | QG.AMB.01: allgemeines QG-Problem, von den Ausgaben entkoppelt | allgemeine QG, nicht TFPT-spezifisch |
| C_p = m_p/Λ_QCD | nichtperturbativ, Gitter-QCD | Rechnung, aber außerhalb TFPT (kein Prinzipmangel) |

## 6. Querbefunde dieser Karte

Beobachtungen beim Vergleich von Ledger und H6, keine Beweise, nichts befördert.

a. **z20 gegen CCC.SEAM.CROSSOVER.01.** D1 des CCC-Vertrags verlangt die Nahtalgebra auf der
   Linse X = S³/ℤ₄, also das Deck als Raumquotient. z20 zeigt: In genau diesem räumlichen Zweig
   ist eine einzelne Naht unmöglich (k = 1 nicht deckinvariant, kleinste Ladung 2, zwei
   Nahtstränge). Unter C″ (drei Familien je Naht) gäbe das c = 16 und N_fam = 6. H6s Einzelnaht
   verlangt den Deck-internen Zweig und zusätzlich, dass die Konfiguration das Reziproke
   z → −1/z nicht trägt (sonst wieder zwei Nähte, H6h12). CCC.SEAM.KINEMATICS.01 benutzt genau dieses
   Reziproke (Penroses Eichung, Tausch der Äonenpole) und lässt es den letzten konformen Modul
   festlegen. Wenn diese Festlegung verlangt, dass die Konfiguration das Reziproke trägt, gibt es
   zwei Nähte. Entweder wird der CCC-Vertrag auf den internen Zweig mit spontan gebrochenem
   Reziproken umformuliert, oder H6 und CCC schließen sich unter C″ aus. Zu prüfen ist, ob v957 das
   Reziproke als Symmetrie des Zustands oder nur der Kinematik braucht.
b. **D4 gegen FLAV.STRONGCP.02.** Der Modellsatz θ_eff = 0 benutzt Reflexionspositivität. Im
   flachen GW-Zweig von D4 ist RP bei gedrehter Phase verletzt (z19, bis −0,23). Das RP-Argument
   für θ = 0 trägt in H6 also nur im Wandzweig.
c. **D4 und QFT4D.LATTICE.FUNDAMENTAL.01.** Die registrierte Entscheidung (positiver Transfer per
   Konstruktion) wählt als Lesart den Wandzweig. Dann ist die Wand physisch: Ihre Spannung liegt auf
   Gitterskala (z9), und der freie Grundzustand erzwingt keinen Defekt. Zu rechnen bleibt, ob die
   Wand mit einem Grundzustand mit Naht verträglich ist.
d. **y1 gegen die Vereinigungstests.** H6 gibt gleiche Kopplungen (sin²θ_W = 3/8) an der
   Gitterskala 1/a, mit a/ℓ_P etwa zwischen 1,65 und 14 (nur die Untergrenze ist robust), also
   nahe der Planckskala. Mit reinem SM-Inhalt treffen sich die gemessenen Kopplungen nicht
   (QFT4D.RGTEST.01 [X], PS.RGTEST.01 [X]). H6 braucht also Zwischeninhalt (PS-Zweig) oder
   Schwellen, die es nicht herleitet; einen EW-Higgs hat H6 nicht. Das ist der Hintergrund des
   RG-Tests in S4.
e. **z16 gegen v1012 (T8).** H6s Gibbs-Prinzip wählt den KMS-Zustand der Naht als
   Entropiemaximum bei Energie null. v1012 zeigt auf endlicher Stufe, dass jede antwortfähige
   entropische Auswahl genau ρ_KMS ergibt, auch ohne μ₄, und typisiert entropische Auswahl daher als
   leer für T8. H6s Zustandsprinzip ist mit T8 verträglich, liefert aber nicht die geforderte
   nichtentropische Auswahl.
f. **H6 und DIMENSION.SELECTOR.4D.01.** Der Selektor ist [C], weil die Herkunft von A1–A4 offen
   ist. H6 liefert für A3 (selbstduale 2-Formen) und A4 (chirale Weylfelder) eine Herkunft unter
   C″: Λ²₊ = Λ²(Familien), drei chirale Mittelformen nur in d = 4. A1/A2 (propagierende Eichfelder,
   dimensionslose Kopplung) bleiben ohne Herkunft.
g. **SEAM.BULK4D.RECON.01 umgekehrt.** H6 baut die Naht als Defekt eines 4D-Volumens und nicht das
   Volumen aus der Naht. Für H6 ist B1–B7 daher die falsche Richtung; der passende Vertrag ist
   TFPT4D.LATTICE.ACTION.01 mit der Naht als Defektsektor. Die Eindeutigkeitsfrage verschiebt sich
   zu BULK.PRIME_COMPLETION.01: Hat das H6-Volumen einen von der Naht entkoppelten Faktor?
h. **H6s Bilanz und die TOE-Tore.** Von den 20 H6-Kriterien betreffen 15 vor allem T1/T2
   (K1–K11, K14, K17, K19, K20); K12 (3+1) ist Kinematik, K13 (Gravitation) steht auf D, K15/K16 sind perturbativ, K18 ist von
   TFPT geerbt. 15 P heißen deshalb nicht, dass drei Viertel einer TOE stehen.

## 7. Dateien

- `experiments/exploration/2026-11-23-phasenraum-f5-holocode/z21_toe_luecken.md` (diese Karte)
- `experiments/exploration/2026-11-23-phasenraum-f5-holocode/z21_ledger_zensus.py` (Zensus, nur
  lesend; Aufruf mit `/Users/stefanhamann/Projekte/tfpt-theoryv4/.venv/bin/python -B`)
- `experiments/exploration/2026-11-23-phasenraum-f5-holocode/ausgabe_z21.txt` (Ausgabe, 14/14
  Prüfungen)

Gelesen: `AGENTS.md`, `docs/research/TFPT-STAND.md`, `VERDICT.md` (H6h12 bis H6), `ausgabe_z20.txt`,
`z*_ergebnis.json` (Urteile), `verification/status_ledger.csv`,
`verification/v649_discipline_audit.py` (Klassifikator), `verification/v1012_toe_gate_witnesses.py`
(T8-Diagnose), `verification/predictions_frozen.json`, `research/originlab/F1-CANDIDATE-CONTRACT.md`,
`research/originlab/gates.json`, Atlas `theory_overview` (nur Navigation). Diese Karte hat außer den
drei Dateien oben nichts geschrieben; Ledger, Atlas, VERDICT, TFPT-STAND, PLAN und progress hat sie
nicht verändert.
