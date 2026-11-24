# TFPT aus der Vogelperspektive (Stand H6h35, 24. November 2026)

## Neueste Zusammenführung (H6h35, `z72`–`z74`, `PRUEFUNG-BESTAND-H6h35.md`)

Die Lösung, so weit sie trägt, als eine Kette aus der originalen Naht:
1. **Naht und Zustand.** Die quasifreie Calderón/Szegő-Kovarianz ist KMS bei β = 4π. Daraus folgen die TFPT-Uhr L₀
   (elliptisch, periodisch), die Leiter und die 2D-Kanäle; am KMS-Punkt gilt ⟨T⟩ = 0. [E]
2. **Holomorphe Twistorgerade.** π_*(η ⊗ π*O(2)) = O(1) ⊕ O(1) mit dem kanonischen η = (A,A) − (P,P). [E]
3. **3+1-Kinematik.**
   - Der Dublettzustand ρ₁ (Windungszahl m = 1 nach der Uhrenquotient-Regel) und seine Tomita-Standardform liefern J = †,
     Herm(2) = R^{1,3}, die Zeitorientierung, V ⊕ V̄* und die Twistorform (2,2). [E]
   - Die eine verbleibende Prämisse ist I_pos, gleichwertig zu T′: Die Lagen der Naht sind die dualen GNS-Vektoren. [C]
4. **Felder.** Hardy- und Stromauslesung erzwingen zusammen den Kern det(I − ZW†)⁻¹ (ν = 1). Daraus folgen die
   SU(2,2)-Leiterdarstellung und die masselosen Gleichungen. [E unter bi-U(2) = T′]
5. **Zeit.** Die kontinuierliche, positive Zeittranslation ist die parabolische Wiesbrock-Differenz der Modulgeneratoren
   geschachtelter Nahtbögen, exakt im Vakuumsektor der Uhr. Im thermischen Original gibt es keine halbseitigen Inklusionen.
   [E, eindimensional]
6. **Teilcheninhalt in 4D.** Die Zählung Σm²qᵐ (konformer 4D-Skalar) ist nativ. Die Gruppe (su(2), SU(2,2)) ist eine
   Lesart, denn die nativen Heber sind Toeplitz statt Oszillatoren. [E Zählung, L Gruppe]
7. **Gravitation.** Linear gibt es beide Helizitäten, nichtlinear gilt die googly-Dualität; aus dem eigenen Spektrum folgt
   keine Skala. [E/O]
8. **Gelöste Bausteine** (Teilchen, Flavour, Massen, α, E₈) bleiben unverändert.

Offen, genau benannt:
- **(a) I_pos aus den Antworten.**
- **(b) Oszillatornormierung.**
  - Es fehlt ein Prinzip, das Toeplitz → CCR erzwingt, analog zu „Leiter statt Band“ (z58).
  - Dieselbe Frage ist „bosonische Paare für ν = 1“, denn fermionische Paare geben ein Determinantenwort ≠ 0.
- **(c) Vakuum gegen thermisch.** Die Zeit lebt im Vakuumsektor, das Original ist thermisch. Wie beide zusammenhängen, ist
  offen; nicht über „modular ⇒ BW“ (GRAV.NONCIRCULAR.01).
- **(d) Nichtlineare volle Gravitation und Flächennormierung.**
- **(e) Mehrpunktantworten und Spin 2.** Der Agentenlauf brach am Nutzungslimit ab.
- **(f) Restflavor.**

**Firewalls:** 3+1-Aussagen bleiben nach DIMENSION.UPLIFT.FIREWALL.01 typisiert. P₀ gilt als Vakuumgröße, nicht als
Folge des Modulflusses.

## Stand H6h34 (Einzelheiten)


Gesamtbild nach den Rechnungen `z58`–`z70` in `experiments/exploration/2026-11-23-phasenraum-f5-holocode/`, dem
Bestandsbericht `PRUEFUNG-BESTAND-H6h34.md` dort und dem Atlas. Maßgeblich bleiben Ledger und Originalquellen.
Statuskürzel:
- **E:** exakt unter genannten Voraussetzungen;
- **C:** bedingt;
- **L:** Lesart;
- **O:** offen.

## Kurzfassung

- **Gelöste Bausteine:** TFPT hat viel bereits gelöst (TFPT-STAND §3): 16 je Generation, E₈ und Clocks, φ₀ und die
  Massenformel, α, Mischungswinkel und CKM, die Readouts des Ledgers. Diese Bausteine werden übernommen und nicht neu
  aufgemacht.
- **Was hier offen ist:** der Ursprung von P1 und die physikalische Brücke, also Zustand und Zeit, 3+1 und Gravitation.
- **Das Bild dazu:** eine holomorphe Naht mit ihren Antworten. Es ist jetzt dichter:
  - Holomorph ist die Naht kanonisch eine Twistorgerade, über die Spinstruktur ihrer Doppelüberlagerung.
  - Lorentzisch wird sie nur über das konjugierte Partnerdublett, das ist die eine Prämisse T′.
  - Der Dublettzustand folgt aus dem Original bis auf die Windungszahl m.
  - Beide Gravitonhelizitäten gibt es linear.
- **Offen:** T′ selbst, m, die physikalische Zeit, die nichtlineare volle Gravitation und die Flächennormierung. Eine
  vollständige Theorie ist das nicht.

## 1. Vorbestand und Firewalls

**Nicht neu.** Die Twistor- und Himmelssphären-Lesart steht im Bestand:
- CELEST.SEAM.01, WOIT.OS.TWISTOR.01 (sieben Kill-Tests) und SEAM.EQUIV.TWISTOR.01, alle [O]; PT′ = Tot(O(1) ⊕ O(1))
  gesetzt in v514;
- Naht als Himmelssphäre (v180);
- τ = i ⇔ J = 0 (v493, CELEST.WP2.01, Pauli-Quadrupol);
- Herm(2)-Kegel (compiler-cone-object-audit, VERDICT_CAUSAL);
- E₈-Karte = Standardform von M₂(ℂ) (E8.STANDARDFORM.CHART.01);
- Hardy-Kern und Wellengleichung (z40, z41).

Neu sind die Rekonstruktion aus Antworten (z58–z60), die kanonische Herkunft von O(1) ⊕ O(1) (z64, z70), der
Dublettzustand (z63), die Reduktion von T auf T′ und die googly-Dualität der Blätter (z65).

**Bindend:**
- DIMENSION.UPLIFT.FIREWALL.01: Kinematik schließt keine 3+1-Aussage.
- GRAV.NONCIRCULAR.01: nie „modular ⇒ Lorentz ⇒ BW“.
- GRAV.WEINBERG.WITTEN.01: Das Graviton ist ein Eichfeld ohne kovarianten Energietensor.

## 2. Wörterbuch (korrigiert)

| TFPT | Bedeutung | Status |
|---|---|---|
| Naht CP¹, Variable z = S | Raum der reinen Zustände des Spinors V = H⁰(O(1))^*; Himmelssphäre | E |
| H⁰(O(1)) = 2 | Spinor V, **nicht** das schwache Dublett (das ist {1 − z⁴, z²} ⊂ H⁰(O(4)), Z34) | E |
| Träger H⁰(O(4)) = 5 | Spin-2-Raum, mit D₄-Twist χ_q (quaternion-frames) | E/L |
| Familien H⁰(K(D)) = 3 | Λ²₊ ⊗ χ_q: selbstduale Zweiformen mit Charakter-Twist | E/L |
| vier Marken | Majorana-Sterne des Spin-2-Zustands (\|2⟩ − \|−2⟩)/√2; J = 0 ⇔ τ = i, dazu das Ausrichtungsbit | E |
| Spinstrukturen (A,A), (P,P) | ihre Differenz η = (1+i)/2 ist die einzige vierteldrehungsinvariante 2-Torsion | E |
| π_*(η ⊗ π*O(2)) | O(1) ⊕ O(1), also die holomorphe Twistorgerade | E |
| Partnerdublett V̄* (Tomita) | einzige Quelle der Signatur (1,3) | E (Satz), C (T′) |
| Originalzustand | ρ_m = diag(1, e^{−πm})/(1 + e^{−πm}) auf End(V) | E je m, m offen |
| Modulfluss | thermische Zeit der Naht, periodisch; keine Raumzeit-Zeittranslation | E |
| zwei Blätter | PT und PT*: linear beide Helizitäten, nichtlinear googly-dual | E linear, O nichtlinear |

## 3. Die Kette

1. **Naht.** P1: Sphäre mit Einheitswindung, vier Marken bei τ = i, Doppelüberlagerung E mit den Spinstrukturen
   (A,A) und (P,P). [E, Bestand]
2. **Zustand und thermische Zeit.** Die quasifreie Calderón/Szegő-Kovarianz ist KMS bei β = 4π. Daraus folgen der
   Modulgenerator 4πL₀, gaußsche Ströme, die Leiter und der 2D-Spur- und Tensorkanal; am KMS-Punkt gilt ⟨T⟩ = 0.
   [E, z58, z59, z68]
3. **Holomorphe Raumzeit.** π_*(η ⊗ π*O(2)) = O(1) ⊕ O(1) mit dem kanonischen η. Die Naht bewegt sich in einer
   vollständigen 4-dimensionalen komplexen Familie. [E, z64, z70]
4. **Lorentzform.** Die Signatur (1,3) kommt nur über das konjugierte Partnerdublett: T′, Naht = P(V) ⊂ P(V ⊕ V̄*). Der
   Dublettzustand ρ_m (z63) gibt den Zeitpfeil im Kegel. [C: T′; m offen]
5. **Felder.** Holomorphie gibt die masselosen Gleichungen; positive Frequenz ist die Szegő-Bedingung. [E unter T′, z61]
6. **Gravitation.**
   - Linear tragen die beiden Blätter (PT, PT*) beide Helizitäten.
   - Nichtlinear ist jedes Blatt das Graviton des anderen (googly); beide zusammen verlangen Daten um die
     Ambitwistorquadrik. [E linear, O nichtlinear, z65]
7. **Gelöste Bausteine.** Teilchen, Flavour, Massen, α und E₈ werden unverändert übernommen. [Bestand]

## 4. Was wirklich offen ist

- **(1) T′:** Warum geht das Partnerdublett konjugiert ein (V̄*)? Kandidat ist die Tomita-Konjugation des Originalzustands;
  dann wäre T′ eine Folge von (2).
- **(2) Die Windungszahl m:**
  - Neutral ist nur m = 4, getrennte Marken verlangen ungerades m, und kein neutrales Dublett trennt die Marken.
  - Die externe Analyse arbeitet mit ρ₁.
- **(3) Physikalische Zeit:**
  - Der globale Modulfluss ist periodisch, also keine Zeittranslation (z70).
  - Die lokale Modulzeit desselben Zustands auf Nahtbögen ist kontinuierlich (z71).
  - Es fehlen die lokalen Translationsbeziehungen (halbseitige Inklusionen bzw. Kozyklen).
- **(4) Nichtlineare volle Gravitation:** Gebraucht werden Daten in der formalen Umgebung der Ambitwistorquadrik (für
  Yang–Mills dritte Ordnung).
- **(5) Flächennormierung:** G und dimensionslose Skalen aus Antworten. Das eigene Spektrum ist dilatationsinvariant.
- **(6) Restflavor:** absolute Quarkamplitude, Neutrinozweig, Phasen (TFPT-STAND §5).
- **Ursprungskandidat für P1 (z71, externe Analyse):** zwei maximal inkompatible binäre Beobachtungen. Sie geben die vier
  Marken mit τ = i (genau bei p = ½), die Vierteldrehung als Clifford-Drehung und das kanonische η als Paarung nach
  Beobachtung. Der Zustand bleibt offen: Das Gleichgewicht ist der Spurzustand.
- **4D-Fortsetzung (z71):** Die S³-Fortsetzung stimmt auf der Linie, ist aber nicht erzwungen. Längs des Nahtwinkels gilt die
  2D-Identität.

## 5. Prüfstand und Suchprogramm

- **Prüfstand:** 13 hergeleitet, 7 unter Lesart.
  - K12 ist Lesart unter T′ mit ρ_m.
  - K6 ist als thermische Zeit der Naht typisiert.
  - K13 zählt linear beide Helizitäten.
- **Suchprogramm:** OriginLab `reconstruction_search` (Neunter Lauf) zählt Normalenbündel und Dublettzustände mit
  berechneten Klauseln auf. Ledger unverändert.

## 6. Nächster entscheidender Schritt

(1) und (2) zusammen:
- **T′ aus Tomita:** Gibt die Tomita-Konjugation des Originalzustands genau das konjugierte Partnerdublett V̄*?
- **m aus T′:** Welches m macht ρ_m mit der Spinbrücke (vier getrennte Marken, η antipodal) verträglich?

Wenn beides gelingt, sind 3+1-Kegel, Orte und Feldgleichungen aus dem Original rekonstruiert. Dann bleiben die
physikalische Zeit, die nichtlineare Gravitation und die Flächennormierung.
