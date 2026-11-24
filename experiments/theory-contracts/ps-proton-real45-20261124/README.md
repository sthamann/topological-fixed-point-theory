# PS.THRESHOLD.02 / PS.PROTON.02: Nachzählung mit der einen reellen 45 (2026-11-24)

Anlass: Der H6-RG-Test z25 (`experiments/exploration/2026-11-23-phasenraum-f5-holocode/z25_rg_test.py`,
Commit 4ac3cc23f) hat nebenbei gefunden, dass `v266_ps_threshold_proton.py` die SU(4)-Adjungierte
(15,1,1) doppelt zählt. Dieser Bericht prüft das unabhängig und schlägt neue Ledgertexte vor.
**Das Ledger ist nicht geändert; die Entscheidung liegt beim Nutzer.**

Prüfskript: [`v1036_ps_proton_real45.py`](v1036_ps_proton_real45.py) in diesem Ordner, Ausgabe
[`ausgabe_v1036.txt`](ausgabe_v1036.txt); eigene Gruppentheorie (Brüche) und eigene Löser, keine
Übernahme aus z25. Das Skript liegt bis zur Ledgerentscheidung bewusst nicht unter `verification/`
(siehe Abschnitt 7); es lädt `tfpt_constants` von dort.

```bash
.venv/bin/python -B experiments/theory-contracts/ps-proton-real45-20261124/v1036_ps_proton_real45.py   # 9 passed, 0 failed, ~11 s
```

## Urteil

**Der Zählfehler ist bestätigt, und er entscheidet das Urteil.** v266 zählt die (15,1,1) als
komplexen Skalar (Δb₄ = 4/3), also als zwei reelle Adjungierte. Die E₈-Hülle liefert genau eine
reelle 45 und damit eine reelle (15,1,1) (Δb₄ = 2/3). v266 sagt in seiner eigenen Prüfung 5 selbst
„the E8 hull supplies only ONE 45“.

Mit der richtigen Zählung gilt:

- **eine Schleife:** τ_p = 2,2·10³⁴ a, zentral unter Super-K (2,4·10³⁴ a). Das ×3-Band
  [7,2·10³³; 6,5·10³⁴] überdeckt die Schranke. Nach v266s eigener Regel (τ/3 > Super-K) ist der
  Inhalt also nicht protonsicher, sondern nur im Band verträglich.
- **zwei Schleifen:** τ_p = 1,1·10³³ a. Auch das obere Bandende (3,3·10³³ a) liegt etwa
  siebenfach unter Super-K, der Inhalt ist ausgeschlossen.
- **Unabhängig von der Zählung** fällt v266s Überlebender schon in zwei Schleifen: Mit der
  komplexen Zählung kommt τ_p = 7,2·10³³ a heraus, das obere Bandende 2,2·10³⁴ a liegt knapp unter
  Super-K. v266 rechnet nur in einer Schleife. Das Repo-Experiment `experiments/gauge-unification`
  gab im Hybrid (SM in zwei Schleifen, PS in einer) schon 7,7·10³³ a.

Die Aussage „Zweig B mit +(15,1,1) überlebt Super-K, Hyper-K entscheidet“ ist damit nicht haltbar.
Zweig B überlebt den Protonzerfall nur noch mit anderen E₈-Inhalten, und zwar in einer Ecke mit
gespaltenen Multipletts (siehe Abschnitt 3).

## 1. Zählung (Konventionen von v249/v266)

GUT-Normierung α₁ = (5/3) α_Y. Einschleifenkoeffizient b = −11/3 C₂(G) + 2/3 Σ T(Weyl) + 1/3 Σ T(komplexer
Skalar) + 1/6 Σ T(reeller Skalar). T(15) = C₂(SU(4)) = 4.

| (15,1,1) als | Δb₄ | b(PS) |
|---|---:|---|
| reeller Skalar (eine reelle 45) | 2/3 | (−29/3, −3, −7/3) |
| komplexer Skalar = zwei reelle (v266 `B_45`) | 4/3 | (−9, −3, −7/3) |

Geprüft: Die allgemeine Formel gibt die SM-Werte b_i = (41/10, −19/6, −7) und die SM-Matrix b_ij
exakt wieder, außerdem v266s `B_MIN` und `B_45`. Die PS-Zerlegung der 248
[(45,1)+(1,15)+(10,6)+(16,4)+(16̄,4̄)] enthält genau eine (15,1,1), und zwar in der reellen (45,1);
45 → (15,1,1)+(1,3,1)+(1,1,3)+(6,2,2) besteht nur aus reellen Darstellungen.

v266s Minimalinhalt ist E₈-verträglich: Das komplexe Bidublett entspricht zwei der sechs reellen 10,
das komplexe (4̄,1,2) einer der vier komplexen 16. Doppelt gezählt ist nur die (15,1,1).

## 2. Neue Zahlen (gleiche Eingaben, Matching und τ_p-Formel wie v266)

Eingaben α_i⁻¹(M_Z) = (59,01; 29,59; 8,47), Matching auf Baumniveau bei M_PS,
τ_p = 10³⁶ (M_GUT/10¹⁶)⁴ (α_U⁻¹/40)² a, Band ×3, M_s = c₃^{7/2} M̄ = 3,06·10¹³ GeV. Zwei Schleifen:
SM-Eichterme in zwei Schleifen plus Top-Yukawa (y_t(M_t) = 0,9369, führende Ordnung), PS-Eichterme in
zwei Schleifen. Der PS-Top-Yukawa (Bidublett) hat in allen drei PS-Faktoren denselben Koeffizienten
(4, 4, 4). Er verschiebt M_PS und M_GUT daher nicht, und α_U⁻¹ ändert sich um weniger als 0,004.

| Inhalt | Schleifen | M_PS [GeV] (κ) | M_GUT [GeV] | α_U⁻¹ | τ_p [a] | ×3-Band | Urteil |
|---|---|---|---|---:|---:|---|---|
| minimal (v266) | 1 | 4,19·10¹³ (1,37) | 2,43·10¹⁵ | 45,06 | 4,4·10³³ | 1,5·10³³–1,3·10³⁴ | ausgeschlossen |
| minimal (v266) | 2 | 3,69·10¹³ (1,21) | 1,19·10¹⁵ | 44,44 | 2,5·10³² | 8,3·10³¹–7,5·10³² | ausgeschlossen |
| +(15,1,1) komplex (v266) | 1 | 4,02·10¹³ (1,31) | 5,89·10¹⁵ | 45,48 | 1,55·10³⁵ | 5,2·10³⁴–4,7·10³⁵ | sicher (nur so in v266) |
| +(15,1,1) komplex (v266) | 2 | 3,54·10¹³ (1,16) | 2,75·10¹⁵ | 44,83 | 7,2·10³³ | 2,4·10³³–2,2·10³⁴ | ausgeschlossen |
| **+(15,1,1) reell (E₈)** | **1** | **4,11·10¹³ (1,34)** | **3,61·10¹⁵** | **45,25** | **2,2·10³⁴** | 7,2·10³³–6,5·10³⁴ | **nur im Band** |
| **+(15,1,1) reell (E₈)** | **2** | **3,62·10¹³ (1,18)** | **1,73·10¹⁵** | **44,61** | **1,1·10³³** | 3,7·10³²–3,3·10³³ | **ausgeschlossen** |

M_PS bleibt in allen Fällen auf der Skalaronskala (κ = 1,2–1,4) und unter M_GUT; die κ-Aussage von
PS.THRESHOLD.02 und PS.KAPPA.01 ändert sich nicht.

α_s-Fehler (±0,065 in α₃⁻¹), reelle 45: τ_p = 1,65–2,86·10³⁴ a in einer Schleife,
0,86–1,44·10³³ a in zwei Schleifen. Die Kalibrierungen: Der Hybrid (SM in zwei Schleifen, PS in
einer, Eingaben des Experiments) gibt `experiments/gauge-unification/results/pati_salam.json` auf
0,02 % wieder. RK4 mit Sehnen-Newton und scipy DOP853 mit fsolve stimmen auf 1,6·10⁻¹⁰ überein.

## 3. Andere E₈-Inhalte (Scan, zwei Schleifen)

Inhaltsraum wie z25: reelle Bidubletts 1–6 und reelle (6,1,1) 0–6 (aus sechs reellen 10),
komplexe (4̄,1,2) 1–4 und komplexe (4,2,1) 0–4 (aus vier komplexen 16), dazu je 0–1 der reellen
Stücke (15,1,1), (1,3,1), (1,1,3), (6,2,2) der 45. Das sind 13 440 Inhalte; alle konvergieren.

- Super-K-sicher insgesamt: 193 zentral, 556 im ×3-Band.
- **Mit κ ∈ [1; 1,7] (PS.KAPPA.01): 11 zentral, 101 im Band**, τ_p = 8,7·10³³–3,5·10³⁴ a,
  κ = 1,00–1,43. Das obere Bandende liegt bei höchstens 1,06·10³⁵ a, also unter der Hyper-K-Reichweite
  von 1,4·10³⁵ a. Bei einer Reichweite von 1·10³⁵ a (Wert in `experiments/gauge-unification`) liegen
  4 Bandenden darüber.
- **Alle 101 brauchen die reelle (15,1,1) und zusätzlich mindestens vier der sechs reellen (6,1,1)
  bei M_PS.** Das sind Farbtriplettpaare, abgespalten von ihren Bidublett-Partnern. Ohne die
  (15,1,1) liegt das Maximum im κ-Band bei 3,9·10³³ a, das obere Bandende bei 1,2·10³⁴ a.
- Die elf zentral sicheren Inhalte haben alle bd_r = 1 und t6_r = 6, z. B.
  {1 Bidublett, 6 (6,1,1), 1 (4̄,1,2), 1 (15,1,1)} mit κ = 1,13, M_GUT = 4,1·10¹⁵ GeV,
  τ_p = 3,55·10³⁴ a.
- Minimal plus eine reelle (15,1,1) ist nicht darunter.
- Gegenprobe zur Schleifenordnung: In einer Schleife hätte dasselbe κ-Band 1331 zentral und 2902 im
  Band sichere Inhalte. Die Schleifenordnung entscheidet.

## 4. Abgleich mit den Referenzzahlen

Alle in der Aufgabe genannten z25-Zahlen stimmen auf alle gedruckten Stellen. Das liegt an
derselben Physik (gleiche Eingaben, Baum-Matching, y_t in führender Ordnung), nicht an gemeinsamem
Code:

- eine Schleife, reelle 45: M_GUT = 3,61·10¹⁵ GeV, τ_p = 2,17·10³⁴ a
- zwei Schleifen: 2,50·10³² a (minimal) und 7,21·10³³ a (komplex)
- Scan: 193/556 sicher, im κ-Band 11/101, τ_p = 8,74·10³³–3,55·10³⁴ a

Vorbefund: `experiments/theory-contracts/universalraum-fable-kernfragen-20260910/fable-runde3/BERICHT.md`
§11 (11. 9. 2026) fand schon, dass die (15,1,1)-Rettung in zwei Schleifen nicht hält (7,4·10³³ a). Dort
wurde die (15,1,1) aber ebenfalls komplex gezählt („komplexe Skalare“, v249 `B_45`), und das Ergebnis
als „disfavoured, nicht ausgeschlossen“ gewertet. Ins Ledger kam der Befund nicht.

## 5. Vorschlag für die Ledgerzeilen

Die Texte sind im Stil des Ledgers (englisch) gehalten; `location`, `dependencies`, `external_data`
und `active` bleiben unverändert.

### PS.THRESHOLD.02

- **status:** `Numerical (two-step RG) [E]/[C]/[X]`
- **script:** `v266_ps_threshold_proton.py; v1036_ps_proton_real45.py`
- **claim:**
  > EXPLICIT PATI-SALAM THRESHOLD MODEL (2026-06-18; recount 2026-11-24, v1036): the two-step
  > carrier-native PS chain made explicit for the E8-allowed contents, at one and two loops. Minimal
  > 16-Higgs (complex bidoublet + (4bar,1,2)): M_PS~4.2e13 (kappa 1.37), M_GUT~2.4e15 (1L);
  > M_PS~3.7e13, M_GUT~1.2e15 (2L). +(15,1,1) from the single REAL 45 of the E8 hull (Delta b4 =
  > T(15)/6 = 2/3): M_PS~4.1e13 (kappa 1.34), M_GUT~3.6e15, alpha_U^-1~45.3 (1L); M_PS~3.6e13
  > (kappa 1.18), M_GUT~1.7e15 (2L). The v266 value M_GUT~5.9e15 counted the (15,1,1) complex
  > (Delta b4 = 4/3, two real adjoints), which the 248 does not supply; it is superseded. M_PS stays on
  > the scalaron scale in all cases. Decides branch B together with PS.PROTON.02.
- **canonical_status:**
  > [E] explicit two-step solutions (M_PS, M_GUT, alpha_U) at 1L and 2L for minimal and +real
  > (15,1,1), M_PS on the scalaron scale (kappa 1.2-1.4); [X] v266's +(15,1,1) value M_GUT~5.9e15
  > rests on a complex (doubled) count and is superseded by v1036; [C] realises the threshold policy.

### PS.PROTON.02

- **status:** `Kill-test (proton decay) [X]/[N]/[O]`
- **script:** `v266_ps_threshold_proton.py; v1036_ps_proton_real45.py`
- **claim:**
  > PROTON DECAY EXCLUDES THE v266 CONTENTS OF BRANCH B (2026-06-18; corrected 2026-11-24, v1036):
  > d=6 tau_p(p->e+pi0) = 1e36 (M_GUT/1e16)^4 (alpha_U^-1/40)^2 yr with the x3 hadronic band, against
  > Super-K 2.4e34 yr. Minimal 16-Higgs: 4.4e33 (1L), 2.5e32 (2L) -> KILLED. +(15,1,1) from the
  > single real 45: 2.2e34 (1L, band 7e33-6.5e34 straddles Super-K -> only band-compatible) and
  > 1.1e33 (2L, band top 3.3e33) -> KILLED. v266's survivor (tau_p~1.6e35) counted two real adjoints
  > and one loop; even that count gives 7.2e33 at two loops (band top 2.2e34 < Super-K). A two-loop
  > scan of 13440 E8-hull contents leaves at kappa in [1,1.7] 101 band-compatible contents (11 safe
  > centrally) with tau_p = 8.7e33-3.5e34 yr; every one needs the real (15,1,1) plus >=4 of the six
  > real (6,1,1) light at M_PS (split multiplets); band top <= 1.06e35 < Hyper-K reach ~1.4e35.
  > Branch B is not content-selected by +(15,1,1); it survives only in this split-multiplet corner,
  > which Hyper-K decides completely. A (boundary-only) stays canonical.
- **canonical_status:**
  > [X] minimal and +real (15,1,1) excluded by Super-K at two loops (the latter only band-compatible
  > at one loop); [X] v266's survival verdict rested on counting the (15,1,1) complex and on one loop
  > -- superseded by v1036; [N] 2L E8 scan: 101 kappa-band contents with tau_p 8.7e33-3.5e34 yr, all
  > needing the real 15 + >=4 light (6,1,1), all within Hyper-K reach; [O] tree matching, step
  > thresholds, split multiplets, x3 hadronic band, gauging fork (QFT4D.FORK.01).

## 6. Abhängige Zeilen und Stellen

Ledger (keine Änderung vorgenommen):

| Zeile | Bezug | Vorschlag |
|---|---|---|
| **PS.UVBRANCH.01** | Zweig B als bedingter UV-Kandidat (v265) | Text um einen Satz ergänzen: protonsicher nur in der Ecke mit gespaltenen Multipletts (reelle (15,1,1) und ≥4 leichte (6,1,1)), τ_p ~ 1·10³⁴ a, Hyper-K entscheidet; Status [C] bleibt. |
| **QFT4D.FORK.01** | eingefrorene Regel: B fällt bei „fail proton safety“ | Regel unverändert; ihre Anwendung vermerken: Für die v266-Inhalte greift die Abstufung. B lebt nur noch über die gespaltenen Inhalte, die keine verbotene Darstellung brauchen und im κ-Band liegen. |
| **PRED.LIVE.01** | nennt „tau_p=1.55e35 yr … +(15,1,1) (safe)“, Skript v283 setzt 1.55e35 fest ein | [X]-Eintrag zum Protonzerfall ersetzen (v266-Inhalte ausgeschlossen; κ-Band-Überlebende 0,9–3,5·10³⁴ a, Hyper-K entscheidet); `v283_live_killtest_scorecard.py` Zeile mit `tau_sel = 1.55e35` anpassen. |
| PS.RGTEST.01 | v249 benutzt dasselbe `B_45` (komplex) | Aussage (5) „+(15,1,1) raises M_GUT proton-safer“ bleibt qualitativ richtig (2,4→3,6·10¹⁵); Hinweis auf die reelle Zählung und auf zwei Schleifen ergänzen. Die v249-Prüfung besteht auch mit reeller Zählung. |
| PS.KAPPA.01 | v253 scannt κ über `B_MIN`, `B_45` (komplex) | κ ändert sich kaum (1,31 → 1,34 in einer Schleife), das Intervall bleibt; nur Hinweis. |
| PS.E8BRANCH.01 | „one 45 … proton-MARGINAL“ | Wort „marginal“ ist nach v1036 zu schwach: Mit der einen 45 ist der Inhalt in zwei Schleifen ausgeschlossen. Hinweis ergänzen. |
| PS.PROTON.KILL.01 | Regel | unverändert; wird jetzt auf die v266-Inhalte angewandt. |

Weitere Stellen ohne Ledger, die die alten Zahlen tragen:
- `verification/ps_threshold_proton.json` (v266-Ausgabe, „proton-safe only with +(15,1,1)“)
- die v266-Beschreibung in `run_all.py` und `script_registry.csv`
- `documents/tfpt_research_contracts/sections/004-source.tex`, `documents/changelog/sections/018-source.tex`
- `experiments/proton-decay/`, `experiments/gauge-unification/README.md` und die generierten Register und Webseiten
- `docs/research/TFPT-STAND.md` (Z. 223–225, „Prüfung offen“): kann nach der Entscheidung auf diesen
  Bericht verweisen.

## 7. Einbindung von v1036 (nach der Ledgerentscheidung)

Unter `verification/` meldete `audit_sync.py` v1036 als „nicht in run_all.py registriert“, und
`website/lib/discipline.ts` wurde veraltet. Die Registrierung setzt voraus, dass eine Ledgerzeile
v1036 nennt; `make_script_index.py` prüft das. Deshalb liegt das Skript bis zur Freigabe in diesem
Ordner (Nummer v1036 reserviert), und das Audit von main bleibt bei den zwei älteren
`F.wolfram`-Meldungen. Nach der Freigabe:

0. Skript nach `verification/v1036_ps_proton_real45.py` zurückverschieben und die
   `sys.path`-Zeile entfernen.
1. Ledgerspalte `script` von PS.THRESHOLD.02 und PS.PROTON.02 auf
   `v266_ps_threshold_proton.py; v1036_ps_proton_real45.py` setzen.
2. Eintrag in `verification/run_all.py` (`MODULES`) und Zeile in `verification/script_registry.csv`.
3. `python3 verification/make_script_index.py` und `make_discipline_stats.py` neu erzeugen, dann
   `python3 verification/audit_sync.py`. Das Manifest (`make_manifest.py`) kommt wie üblich als letzter
   Schritt; es ist schon vor dieser Änderung veraltet (vier Dateien).

Schon vorher rot und ohne Bezug zu dieser Änderung: zwei `F.wolfram`-Meldungen im Audit.

## 8. Offene Punkte

- Matching auf Baumniveau, ohne MS-bar-Schwellenkonstanten bei M_PS und M_GUT; Schwellen als
  Stufenfunktion.
- Der Scan nimmt gespaltene Multipletts an (Überlebenshypothese verletzt); die Überlebenden brauchen
  vier bis sechs leichte Farbtriplettpaare (6,1,1) bei M_PS.
- Hadronisches Band ×3 und der SU(5)-kalibrierte Vorfaktor sind aus v266 übernommen; eine
  kanalaufgelöste Rechnung fehlt. Für die reelle 45 in zwei Schleifen bräuchte es einen Faktor von
  etwa 20 statt 3, damit der Inhalt überlebt.
- Die Schranke gilt nur, wenn SO(10) oberhalb M_GUT geeicht ist (Gauging-Gabel, QFT4D.FORK.01).
