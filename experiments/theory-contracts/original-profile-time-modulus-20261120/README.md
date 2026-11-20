# KERNORIG20261120: Profil, Zeit und Modul aus dem Original selbst

Datum: 20. November 2026. Theorievertrag mit begrenzten Urteilen. Nichts wird befördert, kein
Ledger und kein Paper wird geändert. Anschluss an KERNSIX20261120
(`../kernel-followups-20261120`), dessen offene Punkte „Profil f" und „Breite ℓ" hier am Original
geprüft werden.

Reproduktion (19 Prüfungen, alle bestanden):

```bash
cd experiments/theory-contracts/original-profile-time-modulus-20261120
for s in a_profile b_os c_modulus; do ../../../.venv/bin/python -B check_$s.py; done
```

## Was im Original steht

Pfade relativ zu `documents/archive--tfpt-45--source_extracts--`.

- `05_metrology_source/sections/001-…tex`, Zeilen 98–127: Die geometrische Schließung wird aus
  der Familie Z_rel(s) = Tr_rel(|D|^(−2s) − |D_ref|^(−2s)) + (iπ/2) Δη gelesen, mit dem Residuum
  bei s = 1 und dem endlichen Teil bei s = 0, ausdrücklich „rather than from a freely chosen
  spectral profile".
- `02_carrier_source/sections/009-source.tex`, Zeilen 334–350: Die projizierte modulare freie
  Energie enthält −log Z_rel + Re log det_ζ(D² + μ²) + ½ η.
- `02_carrier_source/sections/010-…tex`, Zeilen 1–40: Die stetigen Schließungsgleichungen sind
  die kritischen Punkte dieser einen freien Energie.
- `04_qft_source/sections/002-…tex`, Zeilen 256–312 und `008-…tex`, Zeilen 140–253: relatives
  erzeugendes Funktional mit Fermionquellen, Form (F,G) = ⟨ΘF,G⟩, Quotient nach dem Nullraum,
  Zeithalbgruppe exp(−tH) mit H ≥ 0, rekonstruierte CAR Theorie mit Vakuum Ω. Die
  Reflexionspositivität des fermionischen Kerns ist dort Voraussetzung des Lemmas.
- `01_boundary_kernel_source/sections/001-…tex`, Zeilen 345–351: Ein Bernsteinprofil ist nur
  eine analytische Realisierung, kein Teil des Datums.

## A: Profil (`check_a_profile.py`, 4/4)

- Die halbe Ableitung der Zeta Familie bei s = 0 ist die primitive Wirkung mit f = −log|x|;
  zusammen mit dem Phasenterm des Originals ist sie der negative Logarithmus der relativen
  Fermiondeterminante.
- Der Wert bei s = 0 trägt nur den Phasenterm.
- Der Term Re log det(D² + μ²) der freien Energie ist bei μ = 0 das Profil +2 log|x|: gleiches
  Profil, umgekehrtes Vorzeichen, doppeltes Gewicht.

Urteil: Das Original legt das Profil fest. Es ist die Zeta Familie, deren endlicher Teil der
Logarithmus der Determinante ist. Die Frage aus KERNSIX „legt das Original f fest?" ist damit mit
ja beantwortet, was die Form angeht. Offen bleibt das Vorzeichen: Ein Fermionfeld gibt
−log|det D|, die freie Energie schreibt +log det(D² + μ²) neben −log Z_rel. Ob das ein
Regulatorterm zu einem fermionischen Z_rel ist, sagt der gelesene Text nicht.

## B: Zustand und Zeit (`check_b_os.py`, 8/8)

Die Rekonstruktion des Originals, angewandt auf seinen eigenen Kragenoperator
D = γ_n(∂_τ + B):

- Der reflektierte Kern auf positiven Zeiten ist exp(−(t + t')|B|) auf einer Spektralhälfte und
  positiv.
- Sein Nullraum ist die andere Spektralhälfte. Der Quotient ist die Polarisation.
- Die Zeithalbgruppe ist exp(−t|B|). Der Erzeuger |B| hat die Lücke λ₁, die Saatskala.
- Der Gleichzeitgrenzwert der Zweipunktfunktion ist der Projektor; im Fockraum ist das der
  Grundzustand von dΓ(B).
- Der Flächenzustand eines Kragens der Breite ℓ aus KERNSIX läuft für ℓ → ∞ gegen dieses
  Vakuum. Das Original sitzt bei Temperatur null.
- Planck Koeffizient des Originals M² = χ²/(2π²) = 32 c₃² χ² gegen das Flächengesetz aus
  KERNSIX: Beide stimmen nur für ℓ χ = (3/4)√ζ(3) = 0,8223 überein. Das ist eine Bedingung,
  kein Ergebnis.

Urteil: Unter den Voraussetzungen, die das Original selbst nennt (Reflexionspositivität des
fermionischen Kerns, Kern = Greenfunktion des Kragenoperators), ist der Randzustand der
Fockzustand der Polarisation, und die Zeit ist die rekonstruierte Zeit mit Erzeuger |B_Σ| und
Lücke λ₁. Die Zeit hat damit eine Länge in Saateinheiten, ohne freie Breite. Die thermische
Lesart mit endlicher Breite aus KERNSIX Schritt 3 wird für die Zeit nicht gebraucht; sie bleibt
eine Familie, deren Grenzfall das Original ist. Gezeigt auf dem Zylinder und auf der Stufe einer
Teilchenzahl plus kleinem Fockraum; dass die Spiegelung Θ des Originals die modulare Konjugation
aus P1BIND ist, ist nicht gezeigt.

## C: Modul der vier Marken (`check_c_modulus.py`, 7/7)

Der endliche Teil bei s = 0 ist für den Dirac Operator des flachen Torus mit antiperiodischer
Spinstruktur die Zeta Determinante |θ₃(τ)/η(τ)|². Sie braucht weder ein Testprofil noch die
Normierung durch λ₁.

- Das Modenprodukt mit der regularisierten Nullpunktsumme gibt θ₃/η; invariant unter τ → −1/τ.
- τ = i ist kritischer Punkt in beiden Richtungen: Minimum entlang der Spiegelfamilie, Maximum
  quer dazu, mit entgegengesetzt gleichen Krümmungen (±1,5964).
- Im Gebiet |Re τ| ≤ 0,9, 0,45 ≤ Im τ ≤ 3 hat die Ableitung genau eine Nullstelle, die bei τ = i.
- Die beiden anderen geraden Spinstrukturen sind am Quadrat nicht kritisch.
- Das Original verlangt einen kritischen Punkt, kein Minimum.

Urteil: Mit dem profilfreien Objekt des Originals ist das Quadrat (Marken bei μ₄) der einzige
kritische Punkt im abgesuchten Gebiet. Das Modellprofil und die Normierung durch λ₁ aus KERNSIX
Schritt 4 entfallen, ebenso die Voraussetzung „Minimum"; die Spiegelbedingung wird für die
Stationarität nicht mehr gebraucht. Es bleiben zwei Voraussetzungen: Der Modul geht über den
Dirac Operator des elliptischen Überlagerers der vier Marken ein, und die Spinstruktur ist
antiperiodisch auf beiden Zyklen. Der Residuenterm bei s = 1 und die Referenz D_ref sind nicht
nachgebildet; das abgesuchte Gebiet deckt den Fundamentalbereich nicht bis in die Spitzen ab.

## Stand

| Punkt aus KERNSIX | vorher | jetzt |
|---|---|---|
| Profil f | vom Original offen gelassen | Zeta Familie des Originals; Vorzeichen offen |
| Randzustand | bedingt auf f = −log\|x\| | bedingt auf die Reflexionspositivität, die das Original voraussetzt |
| Zeit mit Länge | freie Breite ℓ | rekonstruierte Zeit, Erzeuger \|B_Σ\|, Lücke λ₁ |
| Lage der Marken | Minimum eines Modells auf der Spiegelfamilie | einziger kritischer Punkt der Zeta Determinante im abgesuchten Gebiet |
| Kubik | nicht hergeleitet | unverändert |
| Gravitation | nicht hergeleitet | unverändert; Planck Koeffizient des Originals gibt nur eine Bedingung an ℓ |

F1 bleibt in Bearbeitung.
