# Native E8 quarter clock, colour and its cyclic obstruction

The native clock on the affine E8 algebra is the lattice shift
`h=(0,0,0,0,0,1/2,1/2,1/2)`. Its phase on a root of charge `a` is
`exp(2 pi i (h,a))`. The exact root spectrum including the eight Cartan
currents is `(60,64,60,64)`. The native colour convention is the A2 on
the last three coordinates. This note uses the chiral level-one lattice
VOA and its stated toral action. The physical P1 ensemble is still
needed to identify this action with a physical clock.

## Exhaustive replacement test and the twisted ground spaces

For `L=E8`, the toral actions of order dividing four are indexed by
`(1/4)L/L`. The native simple-root Gram matrix is even unimodular, so
the eight coefficients modulo four in this basis enumerate all 65,536
distinct shifts. For an order-four shift `x`, the twisted charge lattice
is `L+x`, and the twisted lowest conformal weight is
`rho = min ||L+x||²/2`. Writing `v=4x` gives the cyclic type
`16 rho = ||v||²/2 mod 4`, since replacing `x` by `x+l` changes
the last expression by a multiple of four.

| Exact requirement | Candidates | Cyclic type |
| --- | ---: | ---: |
| Native phase multiplicities | 15120 | All type 2 |
| Native multiplicities and both native colour roots fixed | 594 | All type 2 |
| Both native colour roots fixed, exact order four, type zero | 864 | The only spectrum is `(64,64,56,64)` |

For the native representative, the four ground spaces are

| Twist | Lowest conformal weight | Dimension | Colour |
| --- | ---: | ---: | --- |
| Identity | 0 | 1 | Singlet |
| Quarter turn | 3/8 | 4 | Singlet plus fundamental triplet |
| Half turn | 1/2 | 16 | Ten singlets, triplet, conjugate triplet |
| Inverse quarter turn | 3/8 | 4 | Singlet plus conjugate triplet |

At the first twist the integer part of the lattice gives the four sign
patterns on the last three entries of `h` with an even number of minus
signs. The half-integer part has larger norm. The three nonconstant
patterns differ by native colour roots, and their colour zero modes give
the triplet. The triplet provides an exact three-dimensional colour
corner of the twisted ground space.

For the half turn, `2h` is an odd integral vector. The fixed lattice
of that involution in the bosonic E8 extension is D8. The lowest
twisted charges are the 16 vectors `±e_i`, of weight 1/2. The vector
extension of D8 is the odd integral lattice `Z^8`, whose vertex
operators give the usual 16 Majorana fields. A vector charge and a
spinor charge have half-integral pairing. Their mutual monodromy is
`−1`. This permits a common description using different sectors and
rules out a single local extension carrying both nonzero charge
families with the same lattice operator products.

## The two cyclic operations on the fixed lattice

The native quarter fixed lattice is `A3 ⊕ D5`. Its discriminant group
has the quadratic form
`q(a,b)=(3a²+5b²)/8 mod 1` on `Z4 × Z4`.
The E8 extension comes from an isotropic order-four glue subgroup.
There are exactly two isotropic subgroups of size four, the diagonal
and the anti-diagonal. They meet in `{(0,0),(2,2)}`. The order-four
glue subgroup has no transverse isotropic complement. The shared
order-two glue is compatible with the D8 intermediate lattice. This
explains why the existence of the even E8 gluing and the gauging of the
dual quarter clock have different conditions.

For a cyclic action of order `n`, take the representative cocycle
`omega(a,b,c)=exp(2 pi i t a floor((b+c)/n)/n)`. The product of
`omega(1,j,1)` over a complete cycle is invariant under rephasing the
binary junctions, because all factors in the coboundary telescope.
For this quarter clock the product is `−1`. Its square has type zero
as an order-two symmetry. The geometric quarter representative also
has type 2. The colour centre shift
`c=(0,0,0,0,0,2/3,−1/3,−1/3)` has type zero. Their product has
order 12 and type 6, so adding the colour centre leaves the obstruction.

## One explicit compatible bulk response

The previously identified E8 Chern–Simons bulk has Cartan matrix `K`.
With actual native simple roots `alpha_i`, choose integral background
charge vectors `(alpha_i,4h)` and `(alpha_i,3c)`. In the reviewed
basis they are `(-1,0,0,0,0,0,2,0)` and
`(-1,0,0,0,0,0,2,-3)`. The response matrix
`Q^T K^(-1) Q` is `diag(12,6)` up to orientation. Restricting the
two bosonic background fields to the clock and colour centre gives
classes `12/2=2 mod 4` and `6/2=0 mod 3` with no mixed entry.
With the bulk orientation matched to the edge, this one standard
symmetry-coupled bulk accounts for both the gravitational response and
the cyclic symmetry obstruction. Its appearance in P1 is a physical
source condition. The original seam must supply this inflow or an
equivalent inverse class. Its regulated junction amplitudes, reflected
Gram form and common absolute R113 reference have yet to be derived.

The [exact audit](f1_native_clock_orbifold.mjs) checks the census,
shortest shifted lattice vectors, colour steps, discriminant groups
and cyclic cocycle equation using integer arithmetic. Run it with
`node research/originlab/f1_native_clock_orbifold.mjs`.

[Native colour and clock source](F1-QUTRIT-CAR-E8-COLOR-SOURCE-GATE-20261115.md).
[Cyclic orbifold convention and conformal weight](https://arxiv.org/abs/1507.08142).
[Symmetry coupling in the bosonic K-matrix description](https://arxiv.org/abs/1205.3156).
