# Structure scenarios, first run of 19 November 2026

Atlas revision `1823e3ca9e529dc86c8e8c0a`, seed 1, 1000 random fingerprints.
Reading aid and search tool. No proof, no promotion, no ledger change.

## Input

- Fingerprint: 93 integers and fractions that established findings of at
  least three of the 58 Atlas problems state (weight 482.5).
- Pool: 731 objects in 23 families, 47 of them closures of an object with
  the objects it contains by construction.

## Result

| Object | Score | Share | Random median | p search-wide |
|---|---|---|---|---|
| 8 modes split (3 + 2) + 3 with its charge, with its parts | 170.5 | 35 % | 69.7 | 0.001 |
| 8 modes split (4 + 1) + 3 | 154.2 | 32 % | | 0.001 |
| 8 modes split (5 + 1) + 2 | 141.5 | 29 % | | 0.001 |
| 8 modes split (2 + 1) + 5 | 133.8 | 28 % | | 0.001 |

The closure of the first object contains 8, 5 and 3 fermion modes, the
split 5 + 3 with D5 + A3 in E8, the split 3 + 2 with its traceless charge,
and the disjointness graph of pairs of five points (Petersen). It carries
45 of the 93 numbers, among them 16, 32, 128, 240, 248, 256 (Fock spaces
and the even unimodular closure), 45, 15, 40, 12, 5/4, 3/4, 5/8, 3/8 (D5,
A3 and their glue), 1/3, 2/3, 1/9, 4/9 (normalised Petersen spectrum),
5/6, 5/3, 3/2, 10/3, 41/10, 48, 96 (charge traces, three copies).

- Layers: after removing these 45 numbers no further object beats the
  search-wide null (second layer p = 0.43). One layer, not several.
- Withheld test: selected on two thirds of the fingerprint (it is chosen
  in 90 % of the splits), the object predicts the withheld third with
  p = 0.006.
- Not carried by this object: 1/6, 1/8, 27, 30, 36, 1/16, 4/5, 1/5, 25,
  4/3, 9/8, 1/729, 64/729, 75, 81, 192, 33, 56 and 30 further numbers.
  No generated object carries 480, 7560, 65537 or 92.

## What this does and does not show

1. The recurring numbers of the record are not a loose collection. A
   third of their weight follows from one object, and the search finds
   that object without being told.
2. The fingerprint is the record of the theory, not a measurement. The
   theory was built on E8, D5 and A3, so a high score of this object is
   expected. The run confirms the internal organisation of the record;
   it is no evidence about nature and no closure of F1.
3. Selection after inspection: the families fermion modes, subset graphs,
   repeated steps, matter content and nested splits were added after the
   first run had shown which numbers stayed unexplained. The p-values do
   not account for that choice. Two features written after reading an
   Atlas statement (sixth powers, operator split 1 + 24 + 75) were
   removed again.
4. The margin between (3 + 2) + 3 and (4 + 1) + 3 is 16 points; its
   significance is not computed.
5. Reals (c3 = 1/(8 pi), alpha) are outside the fingerprint.

## Next

The unexplained numbers are the work list. 1/729, 64/729, 75, 81 and 27
are powers and operator multiplicities of the Petersen channel and would
follow from the same object with one more generic feature rule. 4/5, 1/5,
25, 30 and 6/5 belong to the cut of order five. 4/3, 9/8, 33, 56, 192
have no generated home yet.

## Second run: the unexplained list worked through

Same Atlas revision, seed 1, 1000 random fingerprints.

Three changes, all made after reading the 48 unexplained numbers of the
first run. The p-values do not account for that choice.

1. Bookkeeping removed from the fingerprint: error sizes (3.2e-13),
   decimals without a leading digit (.35), counts of checks ("26 exakte
   Kontrollen"), digits of agreement, moduli of exact arithmetic (65537)
   and line references. 93 numbers become 86.
2. Glue classes computed from the 240 roots for one explicit embedding of
   each of the 15 full-rank subsystems of E8 (D5 + A3: 60, 64, 64;
   E6 + A2: 81, 81; A4 + A4: 50 four times; D4 + D4: 64 three times).
   The first version of this code took a highest root from a glue class;
   the sums of all class sizes now equal the roots outside.
3. Generic features: operators on the modes and those that keep a split,
   charges of the even Fock states with their multiplicities, the share
   of one place among n, and the two operator decompositions of pairs of
   n points.

| | first run | second run |
|---|---|---|
| fingerprint | 93 numbers | 86 numbers |
| carried by (3 + 2) + 3 | 45 numbers, 35 % | 53 numbers, 48 % |
| score against random median | 170.5 / 69.7 | 200.2 / 84.0 |
| p search-wide | 0.001 | 0.001 |
| withheld p | 0.006 | 0.002 |
| selected in the withheld splits | 90 % | 97 % |
| runner-up | (4 + 1) + 3, 154.2 | (5 + 1) + 2, 172.2 |
| ahead of the runner-up when the fingerprint is resampled | not computed | 94 % |
| second layer | p 0.43 | p 0.32 |
| unexplained | 48 | 33 |

The 33 numbers that remain, by cause:

- Other cuts of the same E8, not facts of the split 5 + 3: 27 and 81
  (A2 four times, E6 + A2), 192 (D4 + D4), 56 (A7 + A1), 36
  (A5 + A2 + A1), 6/5 (A4 + A4), 7560 and 480 (numbers of cuts and of
  sign choices), 33 and 54 (root census of an element of order six), 576.
- Six steps of the Petersen channel: 1/729 and 64/729 are the sixth
  powers of 1/3 and 2/3. The object carries the rates; nothing in it
  fixes the number six.
- Powers of one half and 1/12: 1/16, 1/32, 1/64.
- Settings of single computations: 1/10, 1/20, 2/9, 5/16 and the small
  integers 11, 14, 17, 18, 21, 26, 30, 49, which mix counts and content.
- Derivable but not generated: 4/3 as difference of the charges 2/3 and
  -2/3, 5/2.
- Without a home in the pool: 9/8, 92.

Reading: the split (3 + 2) + 3 of eight modes organises about half of the
recurring numbers of the record and predicts withheld ones. What it does
not give is a reason for six steps, the other cuts of E8, and 9/8. The
limits stated above stay: the fingerprint is the record of a theory built
on these objects, and no real quantity is in it.

## Third run: move rules as an object family

`gen_moves` adds `moves:{d}x{q}` for d = 2..6 places and q = 2..5 values: one of the d*q moves "set place i to value
v" is drawn uniformly per step, and one hand is d*q steps. Each group of n fermion modes (n = 2..6) declares its cube
`moves:{n}x2` as a part. On Atlas revision `3d4edf28f96138a983aa8968` (88 numbers, 1000 fakes, seed 1) the best object
stays the closure of `nested:(3+2)+3`: score 213.90 against a null median of 88.46 (200.15 against 85.17 without the
family), share 51.6 % (46.4 %), withheld p 0.0005 (0.0008), 31 numbers unexplained (35). Newly carried: 1/729, 64/729,
1/10, 1/16. The family was written after reading the remainder list, so the p values do not account for it. Exact
checks of the rule against the original modules: `experiments/theory-contracts/petersen-six-20261119/check_cube.py`.
