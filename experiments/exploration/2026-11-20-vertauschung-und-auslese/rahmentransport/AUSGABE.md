# CODE und tatsächliche Ausgabe

Prüfung mit Python 3.14.3 und SymPy 1.14.0. TFPT_REPO bezeichnet den Checkout der ursprünglichen Quellen auf Stand 5b3f26d7befd7550127fe5d60fa94c899ed81ed3. Alle Rechnungen der vier Skripte sind exakt. Die Skripte importieren die vorigen Flaggenproben und den ursprünglichen Cliffordvergleich.

## interferenz.py

```sh
python3 interferenz.py --repo "$TFPT_REPO"
```

```text
Six primitive axes + binary phase erasure + idle 1/6 reproduce native channel: 256/256
Antisymmetric port at relative phases 0, pi/2, pi in the tracial state: [3/5, 1/2, 2/5]
Mean fixed-frame trine contrast at the second mark and half turn: 7/18 2/5
Pair weights: 1/10 per collinear pair, 7/120 per other pair; branch marginal 1/6
Antisymmetric port norm squared 3/5 and symmetric port 2/5 for 15/15 pairs
Antisymmetric axis response: {7/40: 2, 13/60: 3}
Symmetric axis response: {19/80: 2, 7/40: 3}
Uniform-pair antisymmetric control: {1/5: 5}
```

## rahmen.py

```sh
python3 rahmen.py --repo "$TFPT_REPO"
```

```text
Original 16x16 transfer and transported edge covariance: 256/256
Based label holonomy: S6, 720 elements; fixed traceless dimension 0
Even holonomy A6 already has fixed traceless dimension 0
```

## referenz.py

```sh
python3 referenz.py --repo "$TFPT_REPO"
```

```text
Blind density after a common blank: cos(alpha)^2 on blank, sin(alpha)^2/5 on the other 5
Internal two-particle interferometer reads relative phase z^2 and cancels a common phase
Common unit spectator absent; best possible worst-branch overlap probability: 1/2
Reference Gram spectrum 6 : {3: 1, 3/5: 5}
Common reference with two four-mode banks at fixed total occupation four: phase z^2 in one arm, 1 in the other
Universal reference trial saturating the bound: normalized bivector of the native B
```

## windung.py

```sh
python3 windung.py --repo "$TFPT_REPO"
```

```text
Digon native Pfaffian ratio: z**2 ; rank of positive generator: 2
First opposite-line contact at acos(1/6)/pi of one turn
Triangle length between pi and 3*pi/2; min positive odd winding = 3
Phase-normalized stationary current density in this model: 1/(8*pi)
Positive lift condition for either orientation: phase_rate >= abs(speed)
Native signed quarter reached as a primitive endpoint; generator spectrum {3/2: 1, 1/2: 1, 0: 2}
Midpath quarter is excluded: its quaternionic normal is already back at the reference line.
```

## Unveränderter ursprünglicher Septemberprüfer

```sh
python3 "$TFPT_REPO/experiments/theory-contracts/marked-seam-cp-lift-20260921/checker.py" --repo "$TFPT_REPO" --out /tmp/original-marked-check.json
```

```text
{"verdict": "PARTIAL", "count": 65, "weights": ["7/12", "1/12", "0"]}
```
