# Positivitätsvermerk für die Spin-2-Zeilen mit ganzem Formfaktor (24.11.2026)

**Umgesetzt am 24.11.2026 mit Freigabe des Nutzers.** Sieben Ledgerzeilen tragen jetzt einen
Positivitätsvermerk; ihre Zahlen und Polaussagen bleiben unverändert.

## Anlass

In der H6-Exploration (`experiments/exploration/2026-11-23-phasenraum-f5-holocode/z40_gemeinsamer_ursprung.py`,
Teil 10) fiel auf: Die Zeilen GRAV.SPIN2.UNITARITY.01, GRAV.KMS.HESSIAN.01, QGAMB.IDG.01,
GRAV.AMPLITUDE.01, GRAV.LOOP.FINITE.01, GRAV.NONLOCAL.01 und GRAV.SPIN2.EMERGENCE.01 stützen den
Spin-2-Sektor auf den ganzen Formfaktor P(p²) = e^{−p²/M²}/p². v386 liest das als „Spektraldichte =
ein einziger Pol mit positivem Residuum, also Im A ≥ 0“.

## Was v1037 zeigt (`verification/v1037_spin2_positivity.py`, 9/9)

- Ein Källén–Lehmann-positiver Propagator C = Z/Q² + ∫dρ(s)/(Q² + s), ρ ≥ 0, erfüllt Q²C ≥ Z, und der
  Q⁴-Koeffizient von C⁻¹ ist ≤ 0 (50 Zufallsmaße, symbolische Entwicklung).
- Ist die Spektraldichte nur der Pol bei 0, so ist C exakt Z/Q². Der ganze Formfaktor gibt
  C − Z/Q² = (e^{−u} − 1)/Q² < 0 und Q²P = e^{−Q²/M²} < Z; im inversen Propagator steht +Q⁴/(M²Z).
- Kontrollen: Die Stelle-Trunkierung 1/Q² − 1/(Q² + M²) (Geist) verletzt den Test ebenso; der positive
  Zweipol Z/Q² + r/(Q² + M²) besteht ihn, mit C⁻¹ = Q²/Z − rQ⁴/(M²Z²) + O(Q⁶).
- Bei zeitartigem p² = −s wächst der ganze Propagator wie e^{s/M²}/s; er hat keinen Schnitt und damit
  keine Spektraldichte.

Ausgabe: [`ausgabe_v1037.txt`](ausgabe_v1037.txt).

## Was der Vermerk sagt

- [E] Die Polaussagen (ein Pol bei p² = 0, Residuum 1, kein Stelle-Geist) und das Power Counting
  bleiben. v304, v370, v380, v386, v389 und v400 laufen unverändert durch.
- [X] Als physischer positiver Spin-2-Propagator gelesen ist der ganze Formfaktor ausgeschlossen;
  „ein Pol ⇒ Im A ≥ 0 ⇒ Unitarität“ gilt nur für den Pol.
- [O] Ein positiver Spin-2-Sektor braucht ein positives Spektralmaß.

## Geänderte Stellen

- `verification/status_ledger.csv`: Anmerkung an `claim`, Zusatz an `status`, `v1037_spin2_positivity.py` an `script`
  der sieben Zeilen; `canonical_status` unverändert.
- `verification/run_all.py`, `verification/script_registry.csv`, Website-Spiegel, generierte Indizes.
- `documents/changelog/sections/001-source.tex` (Eintrag vom 24.11.) und
  `documents/tfpt_4_frontier/sections/004-source.tex` (Vermerk an der Spin-2-Passage).

Kein T1–T8- oder TOE-Marker bewegt sich.
