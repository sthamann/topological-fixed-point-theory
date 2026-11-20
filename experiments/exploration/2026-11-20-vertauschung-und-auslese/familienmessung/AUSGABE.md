# CODE und Ausgabe

Alle sieben Skripte und der ursprüngliche Prüfer wurden mit Python 3.14.3 und SymPy 1.14.0 erfolgreich ausgeführt. TFPT_REPO bezeichnet den Originalcheckout, geprüfter Stand 5b3f26d7befd7550127fe5d60fa94c899ed81ed3. Die Python Dateien verwenden die mitgelieferte vorherige Fortsetzung über relative Importe. Alle hier berichteten Werte sind exakt.

## ereigniszeit.py

```sh
python3 ereigniszeit.py
```

```text
Ordered allowed transitions: one orbit of 168
Exact event correlations, zero to six intervening jumps: [1/28, -1/168, 1/224, -1/1344, 1/1792, -1/10752, 1/14336]
Continuous unit-rate positive-time correlation: (-sqrt(2)*sinh(sqrt(2)*t/4)/84 + cosh(sqrt(2)*t/4)/28)*exp(-t)
Slope relative to zero-lag density: -7/6
Spectral gap of the flag generator, unsigned lab time: 1 - sqrt(2)/4
```

## flaggen.py

```sh
python3 flaggen.py
```

```text
repeat: independently resampled 1/7, shared endpoint 1/4 , frozen 1
next label after 1 1 [1/4, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8]
next label after 1 2 [1/12, 1/4, 1/6, 1/8, 1/8, 1/8, 1/8]
labeled Fano geometries matching all distinct triples: 1 / 30
prefix/suffix span starting at length zero [1, 7, 21, 21]
nonzero Gram determinant for independent columns: 10625324586456701730816
```

## holonomie.py

```sh
python3 holonomie.py
```

```text
Tangent axis equals difference / projected sum: True False
Line triangles' axis sum spectrum: {1: 3, 0: 4}
All closing flag triangles, exact mean tangent axis projector: {7/40: 2, 13/60: 3, 0: 2}
Fifteen equally drawn directions, exact mean: {1/5: 5, 0: 2}
Closed triangle spin adjoint = one Clifford conjugation; SO minus-rank 2
```

## pfaffkreis.py

```sh
python3 pfaffkreis.py
```

```text
Shortest projective triangle holonomy: -1 on all 35/35 triples, all 64 gauges
Fano line triples / other triples: 7 28
Four-event closed triangles, line and nonline probability: 3/128 7/128
Probability conditional on a triangle that it is incident: 3/10
```

## spintransport.py

```sh
python3 spintransport.py --repo "$TFPT_REPO"
```

```text
Native reconstructed form, modulo scalar: [[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]]
Intertwiner dimension: 1; reversed products: 256/256; antiunitary square -I
Six native skew forms: full Clifford identities; Pfaffian quadratic form is sum of squares
Two-spinor-chart overlap from native matrices: 10/9
Unrecorded shared-flag tick 1 exact spectrum {5/9: 15}
Unrecorded shared-flag tick 2 exact spectrum {101/324: 15}
Unrecorded shared-flag tick 3 exact spectrum {24427/139968: 7, 195479/1119744: 8}
Different full three-tick maps on 30 labeled incidence tables: 30
Rank-seven overlap traces: [7/3, 3, 13/3, 7] ; gap: 7/124416
Closest maps, squared HS superoperator distance: 49/2902376448
Unrecorded shared-flag tick 4 exact spectrum {3945217/40310784: 7, 3946351/40310784: 8}
Seven rank-five channel projectors: span dimension 7; distinct trace overlaps 10/9
Traceless responses: one tick 5/9 ; two ticks {'fresh': 25/81, 'shared': 101/324, 'frozen': 1/3}
```

## trine.py

```sh
python3 trine.py --repo "$TFPT_REPO"
```

```text
Original F index and second gamma of A, zero based: 3 4
Native response deficit: 0 ; wrong-vertex minimum: 35/108
Uniform missing label deficit: 1/3
Original E3 probability levels with tangent contrast (1, 1, 1, -1, -1, -1) : 2 [2/9, 11/81, 11/81, 11/81, 11/81, 11/81, 11/81]
Original E3 probability levels with tangent contrast (-9, -4, -1, 2, 5, 7) : 7 [2/9, 4957/28512, 449/3564, 3277/28512, 835/7128, 3781/28512, 4285/28512]
```

## zyklen.py

```sh
python3 zyklen.py
```

```text
whole incidence 168 cycle character multiplicity 0
one mark 24 cycle character multiplicity 0
one mark, sign 24 cycle character multiplicity 4
unordered apartment 6 cycle character multiplicity 7
apartment sign 6 cycle character multiplicity 15
one flag 8 cycle character multiplicity 3
cycle rank 64 ; character values with frequencies Counter({1: 104, 0: 42, -8: 21, 64: 1})
```

## Unveränderter Septemberprüfer

```sh
python3 "$TFPT_REPO/experiments/theory-contracts/marked-seam-cp-lift-20260921/checker.py" --repo "$TFPT_REPO" --out /tmp/marked-response.json
```

```text
{"verdict": "PARTIAL", "count": 65, "weights": ["7/12", "1/12", "0"]}
```
