# Harmonic P1 seam and transfer search

Use `t = tan(delta/2)` in `(0,1]` and mark order `(1,u,-1,-u)` with
`u = (1+i t)/(1-i t)`. In the ordered cross-ratio convention
`(1,-1,u,-u)`, the ratio is `lambda = -t²`. The six anharmonic images
are distinct generically. The only negative harmonic value in `[-1,0)`
is `-1`, and the equianharmonic values are nonreal. Thus a projective
four-cycle on the marks exists exactly at `t=1`. This is a consequence
of a *geometric* four-cycle, a premise which a formal operator grading
does not supply. The finite search independently checks all 24
permutations with homogeneous integer determinants at every grid point.
The [independent geometry audit](P1-ORIGINAL-REAUDIT-CROSSRATIO-TIME-20261115.md)
also checks the normalized Vandermonde volume and square time candidate.

The verification [v512](../../verification/v512_seam_tau_flag.py)
already computes the side-sensitive odd Green contrast
`Delta = (2/pi) log(1/t)` on the free antipodal seam. Its exponential
`exp(pi Delta/2) = 1/t` is consumed exactly by `odd_green_ratio` when
an independently bound observation is available. Bare `V4` transitivity
holds throughout the family. The cycle rank three of the punctured
sphere and dimension five of the divisor section space also hold
throughout. The graph cycle space is `wedge² U = U tensor sign` for the
reduced four-point permutation space `U`. The four-cycle has traces
`+1` on graph cycles and `-1` on puncture cycles. The identification
requires an orientation line with the appropriate action.

## Square transfer and deformed source

In cyclic mark order, the Hermitian `D4` commutant has diagonal `d`,
adjacent entry `e`, and opposite entry `o`. Its eigenvalues are
`d+2e+o`, `d-2e+o`, and `d-o` with multiplicities `1,1,2`.
With eigenvalue one on the stationary singlet, the last two eigenvalues
remain independent. The condition `e=o` adds the tetrahedral degeneracy.
`V4` instead allows distinct adjacent entries `e1,e2`. In the two odd
characters, the difference of the form eigenvalues is `2(e1-e2)`.
An exact square witness uses the transfer spectrum `(1,1/2,3/4,3/4)`.

Claude Fable 5.2 independently suggested comparing the odd and even
Green contrasts away from the square, and testing reciprocal rational
models with withheld points. For the unit point form with off-diagonal
kernel `-log|z-w|`, the normalized character eigenvalue differences
are `2 log(1/t)` and `log(t/(1+t²))` (alternating minus mean odd),
up to a common kernel scale. The ratio of these contrasts cancels a
common scale; this requires the same charge normalization and field
map. In particular the odd contrasts at `1/3` and `1/2` have ratio
`log(3)/log(2)`. It is irrational by unique factorization and
transcendental by Gelfond–Schneider. Hence models with algebraic raw
splittings at both points and a common nonzero scale cannot match both.
This exclusion does not apply to algebraic *exponentiated* splittings.

The reciprocal model test enumerates `t^k prod_a (1+a t)/(a+t)` with
the configured integer exponent and bounded ordered multiset of factors.
Distinct parameter multisets have distinct poles, so the enumeration
does not count equivalent rational denotations twice. Each model obeys
`R(1/t)=1/R(t)` and `R(1)=1`. The fit and withheld observations are
disjoint inputs. For an unrestricted rational class of degree at most
`d` in numerator and denominator, agreement at `2d+1` distinct
nonpoles forces equality by cross multiplication. The fixture uses
three points and its uniqueness conclusion is for the enumerated class.

## Reproduce the bounded run

```sh
python3 -B research/originlab/main.py --out /tmp/p1-harmonic-search \
  --marked-search-input research/originlab/fixtures/p1-harmonic-search.json \
  --campaign marked_search --campaign r113 --campaign f1
```

The fixture evaluates 79,852 reduced rational moduli, 263,169 square
transfer spectra, and 191,620 reciprocal bridge models. The only modulus
with a projective four-cycle is `t=1`. Of the square transfers, 513 have
the additional tetrahedral degeneracy. The bridge observations are
synthetic controls. Two models fit `R(1/2)=16/25`; only the squared
factor with `a=2` gives both withheld values `25/49` at `1/3` and
`49/64` at `2/3`. The competing `a=34/7` factor gives `55/109` and
`89/116`. Stream hashes appear in the output.

The physical next test needs the original unnormalized odd and even
reflected response forms on the same transported affine carrier at
`t=1/2,1/3,2/3`. Calibrate at `1/2` and withhold the latter two values.
This requires the original fibre weight, time and absolute empty reference.
Even an exact match of ratios leaves an overall exponent and other common
source data to determine. The native R113 diagnostic remains blocked and
no R113 values are assigned by this search.
