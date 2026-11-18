# MATH.SEAM.EQUIVARIANCE.SEARCH.01

24 October 2026. Exact results in declared mathematical categories.

## 1. Fixed source and definitions

Retain the original integral E8 basis, its ordered sign cocycle, and J from
`marked-seam-monodromy-20261023` and its upstream sources. Put
T0=diag(1,-1,1,-1,1,-1,1,-1). Write j and t0 for the generator normalized
field lifts. The original cocycle verifies j^4=t0^2=(t0*j)^2=1 on all lattice
generators. Equality on these generators determines these automorphisms of
the central lattice extension. These are specific lifts of this subgroup;
no splitting of the original global G31 extension follows.

Let V=Hom(E8,F2), with its fixed basis, and D_c the output sign character.
Let A=(J^-1)^T and K=(T0^-1)^T on V. In the fixed coordinates their columns,
encoded as binary masks with coordinate zero in the least significant bit,
are

```
A = [1,253,114,9,185,32,96,129]
K = [1,253,250,9,16,32,65,129]
```

They obey A^2=K^2=I and AK=KA. The field automorphisms generated here have
faithful normal forms (c,k,e)=D_c j^k t0^e, with c in V, k modulo4, e modulo2.
Their product and inverse are

```
(c,k,e)(d,l,f) = (c+A^k K^e d, k+(-1)^e l mod4, e+f mod2),
(c,k,e)^-1    = (A^k K^e c, (-1)^(e+1) k mod4, e).
```

The character action is reduced modulo2 in the exponent of A. Projection
to the eight distinct dihedral lattice isometries and faithfulness of the
characters prove uniqueness. The group has order2048. Its subgroup Gamma
over <J> has order1024 and is the previously studied group.

Let X=CP1 minus(1,i,-1,-i), with basepoint0 and positive radial meridians
gamma_i in cyclic order, product1. Every lifted meridian is M_i=D_ci j.
For a reflection tau_s(z)=i^s conjugate(z), choose the fibre marking
T=D_t j^s t0 and require T^2=1 and

```
T M_i T^-1 = M_((s-i) mod4)^-1.
```

For rotation z->iz choose R=D_r j and require R^4=1 and
R M_i R^-1=M_(i+1 mod4). With both markings at s=0 also require
T R T^-1=R^-1. These are extra declared structures. They do not derive a
physical space, state, preferred dimension or physical dynamics from P1.

Morphisms in the first count fix the projected operators and are simultaneous
conjugations by D_g. The second count forgets the additional lifts and counts
the image in the established meridian moduli modulo all Gamma conjugations.
These counts must not be interchanged.

## 2. Complete equations and finite classification

Put H=A^s K. Normalizing the group words gives precisely these equations:

```
all cases:   c0+A*c1+c2+A*c3 = 0;
reflection: (I+H)*t = 0;
             (I+A)*t+H*ci+A*c_((s-i) mod4) = 0, for every i;
rotation:   (I+A)*r+A*ci+c_(i+1 mod4) = 0, for every i;
dihedral:   (I+A)*t+(K+A)*r = 0, in addition to s=0 reflection and rotation.
```

R^4 imposes no extra character equation because I+A+A^2+A^3=0.
Necessity and sufficiency follow by the unique normal forms: each original
relator has identity lattice projection and is exactly the character displayed.
There is no replacement by a dimension or charge count.

| Marking | Binary variables | Constraint rank | Solution dimension | Character gauge rank | Marked classes modulo D_g |
|---|---:|---:|---:|---:|---:|
| Reflection s=0 | 40 | 20 | 20 | 4 | 65536 |
| Reflection s=1 | 40 | 20 | 20 | 4 | 65536 |
| Rotation | 40 | 24 | 16 | 4 | 4096 |
| Joint rotation and reflection s=0 | 48 | 30 | 18 | 4 | 16384 |

For a projected operator p the gauge translation is (I+p)*g. Simultaneous
translations have rank4 in every row of the table and lie in the solution
space. Their stabilizer has order16. Thus every orbit has16 elements and
the number of classes is 2^(solution dimension-4). This is a complete finite
classification with normal forms obtained by reducing a solution by this
explicit translation subspace.

On the image obtained by forgetting r,t, Gamma acts by uniform translations
(I+A)g on the four ci and the simultaneous involution A. The image has
dimensions16,16,12,10 in the table order. After translation reduction the
dimensions are12,12,8,6, and the A fixed dimensions are8,10,8,6. Burnside gives

```
reflection s=0: (2^12+2^8)/2  = 2176 underlying Gamma classes;
reflection s=1: (2^12+2^10)/2 = 2560 underlying Gamma classes;
rotation:      (2^8+2^8)/2    = 256;
joint:         (2^6+2^6)/2    = 64.
```

The full projected image and every Gamma orbit were also enumerated separately.
This verifies the fixed space calculation without using its Burnside formula.
These are images in the existing526336 classes; they are not counts of all
possible outer isometries, reflections or physical realizations.

## 3. Pure mapping class quotient of the original unmarked stratum

This section adds a different equivalence to the original four meridian model.
It does not quotient fixed reflection markings by incompatible base changes.
The standard pure mapping class generators are the two adjacent Dehn twists,
represented by squared Hurwitz moves b0^2,b1^2. The two twist generation of
PMod(S0,4) is recalled in Section2 of Ferihe Atalan, MPIM24-26:
https://archive.mpim-bonn.mpg.de/5138/1/mpim-preprint_2024-26.pdf

Use b_i(x,y)=(x*y*x^-1,x). With N=I+A, one multiplication gives

```
b_i(ci,cj) = (N*ci+A*cj, ci), j=i+1;
b_i^2: ci -> ci+N*(ci+cj), cj -> cj+N*(ci+cj).
```

N^2=0 proves these squared moves are involutions and commute on the character
tuples. The meridian product is preserved. On the product1 space b2^2=b0^2
modulo uniform character gauge: N*(c0+c1)=N*(c2+c3), so their difference is
the same N*(c0+c1) in all four positions. They preserve each peripheral
Gamma conjugacy class. Hence the full pure mapping class image on the original
Gamma moduli is C2 x C2, generated by P0=b0^2 and P1=b1^2.

The product1 character space has dimension24. Uniform character gauge has
dimension4, leaving a vector space Q of dimension20. A,P0,P1 induce commuting
involutions on Q. For masks0..7, with bits representing A,P0,P1, the respective
fixed point counts are

```
[1048576,4096,65536,4096,65536,4096,4096,65536].
```

Burnside gives 157696 orbits. Independently every one of the1048576 points
of Q was enumerated. Taking the least encoded integer among its eight images
gives precisely157696 different normal forms. Their orbit census, measured in
the previously counted Gamma classes, is

```
size1:   4096 orbits;
size2:  46080 orbits;
size4: 107520 orbits.
```

The counts sum to157696; the size weighted sum is526336. This is the complete
pure mapping class quotient in the original ordered quarter clock stratum.
It includes no puncture permutations or arbitrary W(E8) transformations.

## 4. Independent checks and research system comparison

`independent_check.py` recomputes all168 character columns by literal products
of the original cocycle field lifts on every lattice basis generator. Inverse
operators use their third powers, independently of normal form inversion.
Ranks20,20,24,30 are also checked with SymPy DomainMatrix over GF(2), separately
from the prototype elimination implementation.

The active learner receives only evaluations of individual original group
words, not the target coefficient matrix. Its stated hypothesis grammar is
XOR sums of I,A,K,AK. Closure under composition follows from A^2=K^2=I and
AK=KA, so every requested coefficient belongs to this grammar. Semantically
there are eight distinct maps in this source. The learner selects observations
to distinguish them, saves17 verified word lemmas and composes all four tasks.

Discovery needs43 observations with worst case elimination selection,129 in
basis order, and65..78 with random choice over20 seeds, mean72.55. Every mode
receives the same separate344 basis evaluations for the distinct lemma
certificates. The counts43 and129 are discovery counts, not total oracle cost
or a measured speedup. Restoring saved maps in a new process reproduces the
actual downstream classification without new discovery observations.

This learner is an exact symbolic hypothesis system. No neural weights were
trained in these experiments. The tested local language model variants did
not close any complete structural task. None of the finite results in this
document prove that a proposed neural architecture will solve TFPT as a whole.
