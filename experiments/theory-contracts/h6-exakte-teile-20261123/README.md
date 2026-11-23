# H6EXAKT20261123: Exakte endliche Teile von Kandidat H6 in Lean 4

Datum: 23. November 2026. Entwurf eines Theorievertrags mit begrenzten Urteilen. Nichts wird
befördert. `verification/status_ledger.csv`, Paper und Atlas bleiben unverändert, ebenso
`VERDICT.md`, `docs/research/TFPT-STAND.md`, `v1_vollhypothesen.py`, die Skripte `z*.py` und
`semantic_explorer/research/program.json`.

Ausgangspunkt ist die Exploration `experiments/exploration/2026-11-23-phasenraum-f5-holocode`
(Stand H6h8, Abschnitte H6h3 bis H6h8) und der Abschnitt „Kandidat H6“ in `TFPT-STAND.md`. Dieser
Vertrag sucht kein neues Modell. Er formalisiert die endlichen, exakt entscheidbaren Aussagen, die
dort numerisch oder symbolisch geprüft wurden, und trennt sie von dem, was Prämisse oder Numerik
bleibt. Er ist keine Variante der letzten Quellverträge (etwa `f1-quellenphase-20261123`,
`kernel-followups-20261120`).

## Reproduktion

```bash
cd experiments/lean4-carrier-rigidity
lake env lean ../theory-contracts/h6-exakte-teile-20261123/H6ExakteTeile.lean
```

Ergebnis am 23.11.2026 auf diesem Rechner:

- Exit-Code 0, etwa 34 bis 40 s Wandzeit, 75 bis 86 s CPU. Die Mathlib-Oleans lagen vor, es
  wurde nichts neu übersetzt.
- Ausgabe: nur die 47 Zeilen von `#print axioms`. 44 der geprüften Sätze hängen von
  `[propext, Classical.choice, Quot.sound]` ab, `E8_counts` nur von `propext`,
  `lambda_even_C5` und `matter_components` von keinem Axiom.
- Kein `sorryAx`, kein eigenes `axiom`, kein `native_decide`. Die Datei enthält 87 Sätze.

Die Build-Konfiguration ist unverändert. Wie frühere Vertragsdateien (zum Beispiel
`seam-source-character-search-20261015/lean/SourceCharacter.lean` oder
`chat-triangulation-20261108/*.lean`) wird die Datei eigenständig mit `lake env lean` im
Mathlib-Projekt `experiments/lean4-carrier-rigidity` geprüft. Sie gehört nicht zu
`lake build TfptCarrier`.

`lake env lean` übernimmt die `leanOptions` des Pakets nicht. Die Datei setzt deshalb
`autoImplicit false` selbst. Ein lokales `set_option maxHeartbeats 2000000 in` gilt nur für
`gOrb_comm_diracOrb`, weil dort 64 symbolische Einträge zu prüfen sind.

**Umgebung.**

- Lean 4.29.1 (arm64-apple-darwin24.6.0, Commit f72c35b3f637c8c6571d353742168ab66cc22c00),
  Lake 5.0.0-src+f72c35b.
- Mathlib v4.29.1, Revision 5e932f97dd25535344f80f9dd8da3aab83df0fe6.
- SHA-256:
  - `H6ExakteTeile.lean`: 0ee44b5923f0126ff929c1ed0caa2ce5ff2cf6cc2769ab00c782b010e79a775d
  - `lean-toolchain`: 7dc000621e0046d1aada809e2b7177e64454645cf4c741e9daaf79c99ec2e7a2
  - `lakefile.lean`: 618ca688bf172a3c6de210592061888fef028223244c2a2ad973c007faeb3530
  - `lake-manifest.json`: 79252b4623ddc4d1555ff884e4332b0b3e426bd7b481eb87843bdb4ebea35ebb

## Kurzfassung

| Nr | Aussage | Quelle | Urteil |
|---|---|---|---|
| A1 | g̃⁴ = −1, tr g̃ᵏ = 0 für k = 1, 2, 3, [g̃, D] = 0, auf jeder Bahn | z1, Punkt 2 | bewiesen in Lean |
| A2 | charakteristisches Polynom (X⁴ + 1)² = ∏ₘ (X − ζ^{2m+1})² | z1 | bewiesen in Lean |
| A3 | jeder der vier Sektoren e^{iπ(2m+1)/4} hat Dimension genau 2 von 8; kein invarianter Sektor; die vier Sektoren sind verschieden | z1 | bewiesen in Lean |
| A4 | jeder Sektor trifft jeden Eigenraum von D in genau einer Dimension (ein Viertel), auf jeder Bahn und für jedes Drehzentrum; Bahnen in (ℤ + ½)² haben Länge 4 | z1 | bewiesen in Lean |
| A5 | c₃ = 1/(8π) je Spinorsektor (Wärmespur) | z1 | numerisch in Skript z1; der exakte Viertelanteil ist A4 |
| B1 | P(z) ↦ i P(iz) auf {1, z, z²} ist diag(i, −1, −i), Spektrum {−1, i, −i} | z10, Punkt 3 | bewiesen in Lean |
| B2 | Q(z) ↦ Q(iz) auf {1, …, z⁴} ist diag(1, i, −1, −i, 1); Halbdrehung diag(1, −1, 1, −1, 1), Teilung 3 + 2 | z10 | bewiesen in Lean |
| B3 | Λ^gerade ℂ⁵: 1 + 10 + 5 = 16; 3 · 16 = 48 | z10 | bewiesen in Lean (Arithmetik) |
| B4 | die Modellräume sind H⁰(K(D)) und H⁰(O(D)); ker D₊ = 3, ker D₋ = 0 | z10 | offen (Standardgeometrie, nicht formalisiert) |
| C1 | c = ν/2 ∈ 8ℤ ⇔ 16 ∣ ν; die Wandphase exp(2πi c σ/8) ist für alle σ trivial ⇔ 16 ∣ ν | z5, z8 | bewiesen in Lean (Arithmetik); die Form der Phase ist Eingabe |
| C2 | h_R = N/16 = 1 ⇔ N = 16; 120 + 128 = 248 und die Zerlegungen | z10, Punkt 4; z8 | bewiesen in Lean |
| C3 | e^{iπN/8} = 1 ⇔ 16 ∣ N; e^{−2πi c/24} = 1 mit c = N/2 ⇔ 48 ∣ N; bei E₈ ist e^{−2πi/3} ≠ 1 und gleich e^{4πi/3}; 24 ∣ 8n ⇔ 3 ∣ n | z1, Punkt 4; z8 | bewiesen in Lean (Arithmetik der Phasen) |
| D1 | X = 5(B−L) − 4Y ≡ 1 mod 4 für alle 16 Felder; 6 + 3 + 3 + 2 + 1 + 1 = 16, ohne ν^c 15 | z8 | bewiesen in Lean |
| D2 | ΣY = ΣY³ = Σ(B−L) = Σ(B−L)³ = 0; ohne ν^c ist Σ(B−L) = −1 | Kontrolle | bewiesen in Lean |
| E1 | jede Modulform vom Gewicht 6 zu SL(2, ℤ) verschwindet bei τ = i, insbesondere E₆(i) = 0 | z3 | bewiesen in Lean |
| E2 | i ist die einzige Nullstelle von E₆ im Fundamentalbereich | z3 | offen; bewiesen ist nur die Arithmetik der Valenzformel |
| E3 | d log χ/dτ = −(2πi/3) E₆/E₄, also ⟨T⟩ = 0 ⇔ E₆(τ) = 0 | z3 | numerisch in Skript z3 (auf 10⁻²¹) |

## Prämissen

Die Lean-Sätze selbst sind Algebra und Arithmetik. Sie brauchen keine der H6-Prämissen. Nur ihre
Lesart als Aussage über die Naht hängt an Prämissen.

| Prämisse | Inhalt (H6h8) | Was sie hier trägt |
|---|---|---|
| A\* | Naht ohne Hintergrundgeometrie, ⟨T⟩ = 0 | τ = i über E3 und damit μ₄ als Stabilisator, also dass g̃ eine Symmetrie der Naht ist (A1–A5); die GSO-Summe hinter C3 |
| C″ | Familien sind der chirale Index der vier Marken | Lesart von B1 als Familiencharaktere |
| D | Gravitation als Zustandsgleichung | berührt keinen Satz dieses Vertrags |
| E | die Teilung 3 + 2 ist lokal | Lesart von B2 als Teilung der Hyperladung |
| G | Naht als verdrillter Defekt eines Paarungsmassenfelds | Lesart von C1 (ν windende Fermionen, c = ν/2) und von D1 („welche 16“) |

Mathematische Voraussetzungen der Lean-Sätze:

- A1 und A2 gelten über jedem kommutativen Ring mit einem Element ζ, für das ζ⁴ = −1 gilt.
- A3 und A4 gelten über jedem Körper der Charakteristik 0 mit ζ⁴ = −1, zum Beispiel ℚ(ζ₈) oder ℂ.
- In ℂ gilt alles mit ζ = e^{iπ/4} ohne weitere Annahme.
- E2 enthält die Valenzformel nur als Hypothese im Satz, nicht als Axiom.

## A. Gelifte Vierteldrehung (z1, Punkt 2)

**Modell.**

- Index (Impulslinie j ∈ Fin 4, Spinorkomponente ∈ Fin 2).
- g̃ = C_p ⊗ S. C_p ist die gewichtete zyklische Verschiebung der vier Linien (Zeile = Ziel,
  Spalte = Quelle, wie `block` in z1). S = diag(ζ⁻¹, ζ) = diag(−ζ³, ζ).
- D ist blockdiagonal mit k·σ, wobei i = ζ².
- Die allgemeine Bahn von (u, v) besteht aus (u, v), (−v, u), (−u, −v), (v, −u). Die Phasen p
  sind beliebig mit p₀p₁p₂p₃ = 1. In ℂ gilt das für exp(−2πi k_j·a) bei jeder Verschiebung a
  (`orbit_phase_product`), also für jedes Drehzentrum.
- Die Bahn aus z1 hat k₀ = (1/2, 3/2), verdoppelt (1, 3), und a = (1/2, 0). Die Phasen sind
  (−i, −i, i, i) (`phases_zeta8`), S = diag(e^{−iπ/4}, e^{iπ/4}) (`spinRot_zeta8`), und es gilt
  `dirac2` = D_z1/π.

**Lean-Sätze** (Namensraum `TFPT.H6Exact`):

- `gOrb_pow_four` und `gQ_pow_four`: g̃⁴ = −1.
- `trace_gOrb_pow` und `trace_gQ_pow`: tr g̃ᵏ = 0 für k = 1, 2, 3.
- `spinRot_intertwines`: S (k·σ) = ((rot k)·σ) S.
- `gOrb_comm_diracOrb` und `gQ_comm_dirac2`: [g̃, D] = 0.
- `diracOrb_sq`: D² = (u² + v²)·1. `trace_gOrb_pow_mul_diracOrb`: tr(g̃ᵏ D) = 0 für k = 0, …, 3.
- `charpoly_gOrb` und `charpoly_gQ`: (X⁴ + 1)².
- `X4_add_one_eq_prod` und `charpoly_gQ_sectors`: ∏ₘ (X − ζ^{2m+1})².
- `four_mul_finrank_eigenspace` ist der allgemeine Satz: Aus G⁴ = −1 und tr G = tr G² = tr G³ = 0
  folgt 4 · dim Eig(G, l) = dim V für jedes l mit l⁴ = −1. Die Folgerungen
  `finrank_sector_gOrb`, `finrank_sector_gQ` und `complex_sector_dim` geben je Sektor
  Dimension 2.
- `eight_mul_finrank_joint` ist der allgemeine Satz für gemeinsame Eigenräume. Die Folgerungen
  `finrank_joint_gOrb`, `complex_every_orbit_quarter` und `complex_z1_orbit_quarter` sagen: Der
  gemeinsame Eigenraum von Sektor m und D = ±s hat Dimension 1, wobei s² = u² + v² ≠ 0 (für z1:
  s = ±√10).
- `eigenspace_one_trivial` und `complex_no_invariant_sector`: Der Eigenwert 1 kommt nicht vor.
- `zeta8_sectors_injective`: Die vier Sektoreigenwerte sind verschieden.
- `rot2_orbit_length_four`: rot⁴ = id. Bei ungeraden verdoppelten Koordinaten haben rot und rot²
  keinen Fixpunkt.

**Abgleich mit z1.** Lean belegt Phasen und Spinordrehung eintragsweise (`phases_zeta8`,
`spinRot_zeta8`). Die Blockanordnung ist per Definition dieselbe wie in z1. Zusätzlich wurde in
dieser Sitzung `gQ` mit ζ = e^{iπ/4} gegen `G` aus `z1_kritik_pruefung.py` verglichen; die Prüfung
ist nicht abgelegt. Ergebnis: max |Δ| = 6·10⁻¹⁶ für G und 1,8·10⁻¹⁵ für π · `dirac2` gegen D.

**Urteil.** A1 bis A4 sind in Lean bewiesen. Die Aussage aus z1, jeder Sektor trage exakt ein
Viertel jedes Eigenraums, ohne Fixpunktkorrektur, ist damit je Bahn exakt. Das gilt für alle
Bahnen in (ℤ + ½)² und unabhängig vom Drehzentrum.

Nicht in Lean:

- die Wärmespur und ihr Grenzwert c₃ = 1/(8π) je Sektor (z1 numerisch: 0,039788735773);
- die Summation über unendlich viele Bahnen;
- dass g̃ die Symmetrie der physischen Naht ist (A\*).

## B. Familiencharaktere aus dem Markenindex (z10, Punkt 3)

**Lean-Sätze.**

- Differentialform ω = P(z) dz/(z⁴ − 1):
  - `quarter_on_differential_coeff`: Die Koeffizienten werden zu i^{k+1} P_k.
  - `quarter_on_differential_degree`: Der Raum deg P ≤ 2 ist invariant.
  - `marks_polynomial_invariant`: z⁴ − 1 ist invariant.
  - `famMatrix_eq`, `famMatrix_represents` und `famMatrix_charpoly`: die Matrix diag(i, −1, −i)
    mit charakteristischem Polynom (X − i)(X + 1)(X + i).
- Funktion f = Q(z)/(z⁴ − 1):
  - `quarter_on_function_coeff` und `carrierMatrix_represents`: Die Koeffizienten werden zu
    i^k Q_k.
  - `carrierMatrix_sq`: Das Quadrat ist diag(1, −1, 1, −1, 1).
  - `carrierMatrix_pow_four`: Die vierte Potenz ist 1.
  - `carrier_half_turn_charpoly`: (X − 1)³ (X + 1)².
  - `carrier_split_three_two`: Die Teilung ist 3 + 2.
- `lambda_even_C5`, `lambda_even_C5_sum` und `matter_components`: 1 + 10 + 5 = 16 und
  3 · 16 = 48.

**Geltungsbereich.** Modellräume sind die Polynome vom Grad ≤ 2 und ≤ 4. Lean beweist die Wirkung
von z ↦ iz auf den Koeffizienten. Der Faktor i stammt aus dz, und der Markennenner z⁴ − 1 ist
invariant.

Nicht in Lean:

- dass diese Räume H⁰(K(D)) und H⁰(O(D)) sind (Regularität bei ∞, Riemann–Roch auf ℙ¹);
- der Diracindex ker D₊ = 3, ker D₋ = 0;
- die Identifikation von Λ^gerade ℂ⁵ mit dem Spinor von Spin(10) (bewiesen ist nur die
  Dimension);
- die Kopplung der physischen Fermionen an K^{1/2} ⊗ O(D) (offen in H6h8).

Die Rücknahme aus Z10 (Punkt 2) bleibt bestehen: Die Charaktere legen die Satzzahl nicht fest.
Das Gegenbeispiel ist nicht formalisiert.

## C. Arithmetik von Zentralladung, Twistgewicht und GSO-Phasen

**Lean-Sätze.**

- `central_charge_in_8Z_iff` und `cut_phase_trivial_iff` (mit dem Hilfssatz
  `exp_two_pi_I_mul_eq_one_iff`): c = ν/2 ∈ 8ℤ ⇔ 16 ∣ ν. Dasselbe gilt für die Wandphase
  exp(2πi c σ/8) bei allen σ ∈ ℤ.
- `twist_weight_one_iff`: N/16 = 1 ⇔ N = 16.
- `T_covariance_iff`: e^{iπN/8} = 1 ⇔ 16 ∣ N.
- `T_invariance_iff` und `central_charge_in_24Z_iff`: e^{−2πi (N/2)/24} = 1 ⇔ 48 ∣ N ⇔ c ∈ 24ℤ.
- `E8_T_phase_nontrivial` und `E8_T_phase_eq`: Bei N = 16 ist die Phase e^{−2πi/3} ≠ 1, und sie
  ist gleich e^{4πi/3}, also 240°. Die Gleichheit mit δ_PMNS ist nur eine Zahl; ein Zusammenhang
  ist nicht hergeleitet.
- `sets_strict_T_iff` und `sets_minimal_strict`: 24 ∣ 8n ⇔ 3 ∣ n; das kleinste solche n ist 3.
- `E8_counts`:
  - C(16, 2) + 2⁷ = 248;
  - 45 + 15 + 60 = 120;
  - 2 · (16 · 4) = 128;
  - 78 + 8 + 2 · 27 · 3 = 248;
  - 10 + 6 = 16.
- `kms_energy_c8`: c/24 = 1/3.

**Geltungsbereich.** Bewiesen ist die reine Arithmetik der Phasen, deren Form vorgegeben ist.

Nicht in Lean:

- dass die Wandphase exp(2πi c σ/8) lautet (Anomaliepolynom (c/24) p₁, ∫p₁ = 3σ, eine
  4-Mannigfaltigkeit mit σ = 1);
- die T-Wirkung auf θ₂ (Faktor e^{iπ/4}) und auf η (Faktor e^{iπ/12});
- die Identität χ₁₆ = E₄/η⁸;
- die ℤ₁₆-Klassifikation (García-Etxebarría und Montero 2019);
- der Callan–Harvey-Zufluss, auf dem die Wahl n = 1 beruht.

## D. Ladungen einer Generation (z8)

**Lean-Sätze.**

- `chargeX_eq_Xint`, `Xint_mod_four` und `chargeX_one_mod_four`: X nimmt die Werte 1, 1, −3, −3,
  1, 5 an, und jeder ist ≡ 1 mod 4.
- `generation_count`, `generation_count_without_nuc` und `z16_condition`: Es sind 16 Felder;
  16 ∣ 16, aber 16 ∤ 15.
- `anomaly_sums` und `BL_sum_without_nuc`: die Kontrollsummen.

**Geltungsbereich.** Bewiesen ist die Ladungstabelle der linkshändigen Weyl-Multipletts, mit
Konjugierten für rechtshändige Felder. Die Werte sind:

| Feld | B−L | Y |
|---|---|---|
| Q | 1/3 | 1/6 |
| u^c | −1/3 | −2/3 |
| d^c | −1/3 | 1/3 |
| L | −1 | −1/2 |
| e^c | 1 | 1 |
| ν^c | 1 | 0 |

Nicht in Lean:

- die Spin-ℤ₄-Struktur und die ℤ₁₆-Anomalie als Bordismusinvariante;
- dass ν^c ν^c der einzige eichinvariante Majorana-Bilinear einer Generation ist (z8);
- Wechselwirkungen.

## E. E₆ bei τ = i (z3)

**Lean-Sätze.**

- `S_smul_I`: S · i = i.
- `weight6_vanishes_at_I`: Für jedes `f : ModularForm Γ(1) 6` gilt f(i) = 0, wegen
  f(−1/τ) = τ⁶ f(τ) und i⁶ = −1.
- `E6_vanishes_at_I`: dasselbe für Mathlibs `ModularForm.E` vom Gewicht 6. Mathlib normiert diese
  Reihe mit dem Faktor 1/2 über teilerfremde Paare; die Nullstelle hängt davon nicht ab.
- `valence_weight6`: Aus 12v_∞ + 6v_i + 4v_ρ + 12Σ = 6 folgt v_i = 1, und alle anderen
  Ordnungen sind null.
- `valence_weight4`: Bei Gewicht 4 ist nur v_ρ = 1. Unter derselben Hypothese verschwinden also
  E₄ und χ_E₈ = E₄/η⁸ nur bei ρ, am Sechseck aus z3.

**Urteil.**

- E1 ist ein vollständiger Beweis.
- E2 ist offen. Die Valenzformel steht in Mathlib v4.29.1 nicht zur Verfügung; sie geht nur als
  Hypothese `hval` in den Satz ein. Ebenfalls nicht formalisiert ist E₆ ≢ 0 (konstanter
  q-Koeffizient 1).
- Das Argument für E2: Eine von null verschiedene Form vom Gewicht 6 erfüllt
  6/12 = v_∞ + v_i/2 + v_ρ/3 + Σ_p v_p. Die einzige Lösung in natürlichen Zahlen ist v_i = 1, alle
  anderen null (`valence_weight6`). Also hat E₆ im Fundamentalbereich genau die einfache
  Nullstelle bei i.
- E3 ist nur in z3 numerisch geprüft.

## Was dieser Vertrag nicht zeigt

- Keine physikalische Identifikation:
  - g̃ als Symmetrie der Naht (A\*);
  - Familien als Index (C″);
  - 3 + 2 als Hyperladung (E);
  - die 16 als windende Paare (G);
  - Gravitation (D).
- Keinen der offenen Punkte aus H6h8:
  - das Massenbündel oder die Defektsumme;
  - den gravitativen Zufluss;
  - die Kopplung an K^{1/2} ⊗ O(D);
  - die Satzzahl;
  - den gemeinsamen Elternoperator.
- Keinen Ledgerstatus.

## Weg zu einer späteren Beförderung (Entwurf)

1. **Unabhängige Wiederholung.** Auf einem sauberen Checkout mit derselben Toolchain neu prüfen.
   Danach eine Quittung als JSON anlegen, im Stil von
   `seam-source-character-search-20261015/lean-receipt.json` (Hash, Kommando, Exit-Code, Axiome).
   Bewusst noch nicht angelegt, ebenso wenig ein `contract_index.json`.
2. **Offene Formalisierungen schließen oder als Prämisse benennen.**
   - Valenzformel und E₆ ≢ 0, oder dim M₆(SL(2, ℤ)) = 1 mit einer Nullstellenanalyse: schließt
     E2.
   - d log(E₄/η⁸)/dτ = −(2πi/3) E₆/E₄: schließt E3. Mathlib hat `ModularForms/Derivative.lean`
     und `DedekindEta.lean`; der Weg ist zu prüfen.
   - T-Transformation von θ₂ und η für die GSO-Phasen: schließt C3 vollständig (Mathlib
     `JacobiTheta`).
   - Die Modellräume als H⁰ auf ℙ¹, über eine explizite Entwicklung bei ∞ statt Riemann–Roch:
     schließt B4.
   - Optional: das Gegenbeispiel aus Z10 (Punkt 2) als negativen Satz formalisieren.
3. **Was dann gilt.** Die endlichen Aussagen wären exakte Konsequenzen in ihrem Modell. Die
   physikalischen Lesarten bleiben unter A\*, C″, E und G, für Gravitation unter D. Lean kann sie
   nicht befördern.
4. **Akzeptanzkriterien.**
   - Kompilieren mit Exit-Code 0.
   - Keine Axiome außer `propext`, `Classical.choice` und `Quot.sound`.
   - Kein `sorry`, und der Hash stimmt.
   - Jede beförderte Aussage nennt ihren Lean-Satz und ihren Geltungsbereich.
   - Jede physikalische Lesart nennt ihre Prämissen.
5. **Ledger.** Welche Ledgerzeile mit welchem Status betroffen wäre, entscheidet erst ein Review.
   Dieser Entwurf schlägt keine Ledger-ID vor.

## Quellen

- `experiments/exploration/2026-11-23-phasenraum-f5-holocode/VERDICT.md`, Abschnitte H6h3 bis H6h8.
- `docs/research/TFPT-STAND.md`, Abschnitt „Kandidat H6“.
- Skripte und Ausgaben im selben Explorationsordner:

  | Skript | Ausgaben | Inhalt |
  |---|---|---|
  | `z1_kritik_pruefung.py` | `ausgabe_z1.txt`, `z1_ergebnis.json` | Punkt 2: gelifte Vierteldrehung; Punkt 4: T-Phase |
  | `z3_hintergrundfreiheit.py` | `ausgabe_z3.txt`, `z3_ergebnis.json` | E₆ und τ = i |
  | `z5_gabel.py` | `ausgabe_z5.txt`, `z5_ergebnis.json` | Wand bei c ∈ 8ℤ; S und T auf den Sektoren |
  | `z8_sechzehn_und_saetze.py` | `ausgabe_z8.txt`, `z8_ergebnis.json` | X mod 4, Satztabelle, E₈-Zerlegung |
  | `z10_kritik_zwei.py` | `ausgabe_z10.txt`, `z10_ergebnis.json` | Markenindex, 3 + 2, 16, h_R |

- Lean: `H6ExakteTeile.lean` in diesem Ordner.
- Verwendete Mathlib-Bausteine: `ModularForm.E`, `SlashInvariantForm.slash_action_eqn_SL''`,
  `LinearMap.IsProj.trace`, `Matrix.kronecker_diagonal`, `Matrix.det_blockDiagonal`,
  `Polynomial.comp_C_mul_X_coeff`.

## Ledger

`verification/status_ledger.csv` ist unverändert. Dieser Vertrag ändert keinen Status, befördert
nichts und legt keinen Ledgereintrag an. In den Atlas wurde nicht geschrieben.
