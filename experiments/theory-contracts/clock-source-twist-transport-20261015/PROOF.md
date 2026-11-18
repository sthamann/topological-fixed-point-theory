# The actual J action on the existing source twist family

15 October 2026. Research ID `UR.CLOCK.SOURCE_TWIST_TRANSPORT.01`.
Overall verdict **PARTIAL**. This is an exact conditional operator result in
an already stipulated E8 lattice Fock representation, with an explicit
fixed-sector exclusion. No physical T1–T8 gate or ledger claim is promoted.

## 1. The requested identification and the retained inputs

The question is whether the actual Gaussian compiler complex structure J can
act on the same source state, charge, time and uncompensated quarter transition.
We test the first operator identification, rather than inferring physical
selection from matching dimensions or clock orders.

The supplied representation is

    H = M(1) tensor C_eps[Ls],   Ls = D8 union (D8+s),   s=(1,...,1)/2.

It is the even E8 completion already used in
`seam-joint-ward-state-20261014`. Its existence and physical inclusion are
premises here, not derived from the raw seam or P1/P2. The source charge is
`Q_q=H(q)_0`, q=(1,...,1). On the finite-energy core the retained source time is

    K_delta(p,N)=|p|²/2 + delta.p + N,
    delta_k=(k+alpha)/8 - 1/2,  k=0,...,7.

The original physical scaling H_source=8 K_delta, where used in the folded
source, is unchanged by the covariance statement. The actual quarter cases
are alpha=1/4 and 3/4. As in the original fold audit, the sign relating this
alpha to the source-circle holonomy depends on the retained orientation; no
orientation or mode reindexing is changed in this calculation.

The actual J is read from `verification/v634_st31_structure.py:J_vec`:

    J(x0,x1,...,x6,x7)=(-x1,x0,...,-x7,x6).

Thus J²=-I and J⁴=I. It is an isometry of the original Ls lattice. The checker
tests this on its integral basis, as well as checking the 240 original roots.
The basis test, with integral inverse, proves lattice preservation; the root
census is a separate finite regression.

This is the Gaussian root-lattice J, not the distinct cyclotomic Frobenius
operator G in `G C5 G^-1=C5²`. No identification of those representations is
made by this contract.

## 2. A fixed-source affine lift fails for an exact reason

Let A be a lattice isometry and suppose an affine action p -> Ap+r preserves
the *same* displayed K_delta for every p, with no oscillator or Hamiltonian
change. Comparing linear coefficients in p gives necessarily

    A^T(r+delta)=delta,  hence r=A delta-delta.

The constant then vanishes because |A delta|=|delta|. To map Ls onto itself,
r must belong to Ls. For an isometry with A²=-I this entails

    (I+A)r = -2 delta in Ls.

But adjacent coordinates of 2 delta differ by 1/4, whereas coordinate
differences of every Ls vector are integers. Consequently 2 delta is not in
Ls, for any alpha. There is **no affine charge-lattice lift of any such A**
which preserves this same K_delta. This includes every lattice-conjugate of
the actual J when both objects are compared in the same source frame.

For alpha=1/4 the checker obtains explicitly

    delta=(-15,-11,-7,-3,1,5,9,13)/32,
    J delta-delta=(13,-2,5,-2,-3,-2,-11,-2)/16, not in Ls.

This excludes the specified fixed-sector affine identification. It does not
exclude a clock which transports backgrounds, an operator with another
specified charge action, or a larger source construction. A clock action
does not in general have to commute with one fixed Hamiltonian.

Keeping the original fixed q does not rescue that identification:
Jq=(-1,1,-1,1,-1,1,-1,1) is not q. No affine constant can remove the
nonconstant linear difference q.Ap-q.p. Thus the actual J is not a symmetry
of the original fixed charge operator either.

## 3. The positive result: one covariant operator family

Retain the entire covariant tuple, rather than identifying its members:

    delta_a=J^a delta,   q_a=J^a q,   a in Z4.

All four delta_a are distinct in either quarter case; delta_2=-delta. For
the original source offset j define

    nu_j=((k+j)/8)_k,
    tau_(alpha,j)=nu_j+alpha q/8,
    t_(alpha,j)=exp(2 pi i H(tau_(alpha,j))_0).

This is the full uncompensated torus transition already retained in the
source contract: `exp(2 pi i alpha Q/8) g0`. The alpha term is not discarded
or replaced by its square. Put tau_(a,alpha,j)=J^a tau_(alpha,j).

For the phase-corrected unitary U_J constructed below,

    K_(delta_(a+1)) U_J = U_J K_(delta_a),
    Q_(q_(a+1)) U_J = U_J Q_(q_a),
    t_(a+1,alpha,j) U_J = U_J t_(a,alpha,j),
    U_J^4=I.

These equations hold on the full finite-energy core and extend to the
corresponding diagonal selfadjoint operators and their unitary groups. The
first equality is the exact quadratic identity

    |Jp|²/2+(J delta_a).Jp+N = |p|²/2+delta_a.p+N.

The others follow from equality of the *full* cocharacter pairings, not just
their mod-1 or adjoint shadows. On the Heisenberg oscillators U_J acts by J,
which preserves the total level N. On charged root currents it acts as

    U_J Y(e^p,z) U_J^-1 = eta_J(p) Y(e^(Jp),z).

The phase eta_J(p) is fixed by the actual source cocycle and the explicitly
stated basis convention below. Current central terms and adjoints are
preserved. In particular the quarter-transition character of a charged
current is transported verbatim. The checker covers both quarter cases,
four frames, all eight j offsets, and the lattice basis. Basis equality
proves the character relation for every lattice charge; the 240-root
comparison is not substituted for an all-charge argument.

The coefficient q must be transported too. If q is held fixed while J is
applied, this is not the same tuple and the Q equality fails. Transporting
q is a real limitation, not an omitted gauge compensation.

## 4. The actual ordered CAR/E8 cocycle, not a new sign model

The checker reads the original functions `eps_car`, `split_even_lattice` and
`eps_e8` from `seam-joint-ward-state-20261014/operator_transport_check.py`.
No replacement Hamming/Pauli cocycle is inserted. The original construction
uses the ordered CAR cocycle on D8 and adjoins X=E_s by

    X E_d X^-1=(-1)^(s.d) E_d,   X²=-E_q.

Every p has the unique form d+c s, d in D8 and c in {0,1}. The resulting
epsilon restricts exactly to the source ordered CAR cocycle on D8. Its
associativity and all-charge commutator were proved in the cited source;
here its code and source proof are pinned by SHA-256.

Let b1,...,b8 be the exact original integral E8 basis from that same file.
The lifted central group has generators z, t1,...,t8, where z²=1 is central,
the ti have unrestricted integer powers, and

    [ti,tj]=z^(bi.bj).

There are no extra finite-order relations on the ti, because Ls is free
abelian. The assignment

    z -> z,   E_bi -> E_(Jbi)

preserves all these commutators: J preserves the inner product. Since J is
an integral lattice automorphism, this defines an automorphism of the
*same* central extension. Negative powers use their actual inverses,
`E_b^-1=epsilon(b,-b) E_-b`; ignoring that sign would change the lift.

For n in Z8, multiply the ordered generator powers in the actual source
cocycle. Write

    product_i E_bi^ni = c_B(n) E_(Bn),
    product_i E_(Jbi)^ni = c_JB(n) E_(JBn).

Then

    eta_J(Bn)=c_JB(n)/c_B(n)

is defined for every lattice charge. The induced operator
`U_J|p>=eta_J(p)|Jp>` is unitary on the charge Hilbert space, with the ordinary
orthogonal second quantization J on the Heisenberg Fock factor. It fixes the
vacuum and preserves adjoints. Its all-charge cocycle identity is

    eta_J(u+v) epsilon(u,v)
       = eta_J(u) eta_J(v) epsilon(Ju,Jv).

The identity follows from the generator-extension construction; the 1,024
generator-pair checks are regressions of its implementation.

Since J⁴=I, the fourth power of this extension automorphism is a lattice
character. The exact calculation gives +1 on each of the eight original
basis generators. It is therefore +1 on **every** lattice vector. The
oscillator action has fourth power I as well. The vacuum normalization
removes a global scalar: U_J⁴=I on the complete stipulated Fock space.

The displayed choice E_bi -> E_(Jbi) is a phase convention, not a physical
selection rule. Other cocycle-compatible lifts differ by a torus character
chi. Their fourth-power correction is the norm character
`product_(a=0)^3 chi(J^a p)`, which is one because
`I+J+J²+J³=0`. Thus the order-four result is not achieved by secretly tuning
a new physical coupling or replacing the original signs. It does not
canonically select the lift's torus character.

## 5. State, positivity and what the four frames mean

The earlier unbounded-lattice source theorem proves K_delta>=0 for
0<alpha<1, with precisely two ground vectors, at p=0 and
r0=(1,1,1,1,-1,-1,-1,-1)/2, and oscillator vacuum. It is not a finite-root
scan. Unitary covariance transports that exact statement to every frame:
the two ground labels are 0 and J^a r0. All spectra, multiplicities, positivity
and domains are unitarily equivalent. No degeneracy is removed.

The same vacuum |0> is fixed by U_J. Within the specified normal ground
doublet, the existing full source Cartan Ward values select its density
|0><0| by positivity. This result, already proved in the source contract,
transports with the marked Cartan tuple. It is not a new selection from
G31 symmetry, a proof of one-dimensional Hamiltonian ground space, or a
state theorem on arbitrary representations. The finite irreducible C4
Clifford representation and this infinite E8 Fock representation must not
be substituted for one another.

The quartet delta_a is at present a covariance orbit of supplied backgrounds
and charge markings. It has **not** been derived as four physically existing
superselection sectors of the primitive source. One can package the family
as a direct sum with a cyclic U_J action, but that is bookkeeping; it does
not select a physical sector weight, introduce a physical register, or
derive a sector-switch rate. If J is merely a change of source frame, the
result is precisely a phase-corrected covariance statement under that frame
change. It remains useful because it states exactly which time, charge and
quarter-transition data have to move together.

## 6. What is reused and the remaining decisive interface

The common even-sector `W4=U_s R_c` factorization, all its Klein signs and
current central terms are already established in the preceding contract.
The present result neither reopens it nor derives its E8 completion.
Likewise a known complete positive word functional already has a conditional
GNS reconstruction; repeating that theorem cannot select its input source.

The new tested distinction is exact: actual J fails as a fixed K/Q symmetry
in this affine source class, but has a genuine order-four lift transporting
the complete stipulated tuple. A further physical identification must state
whether these four tuples are allowed source frame descriptions or actual
physical sectors and must supply its marked operator/state map. This
covariance result supplies no instrument or local charged-field insertion,
no ordinary-local inclusion of odd vector CAR fields in E8, no physical
event law, no raw Xi, and no selection of the E8 representation from P1/P2.

## 7. Validation and original sources

`python3 experiments/theory-contracts/clock-source-twist-transport-20261015/checker.py`
passed **123/123 exact groups**, with zero failures. `results.json` records
the original source hashes, cocycle basis character, actual quarter data
and every check. Numerical spectral fits, GPU checks, new Lean builds and
the full verification suite were not executed or claimed.

The source pins are in `source_manifest.json`:

- `verification/v634_st31_structure.py`, J_vec at lines 109–113 and its E8
  construction: the actual Gaussian complex structure.
- `seam-joint-ward-state-20261014/checker.py`: original Ls roots and full
  lattice membership predicate.
- `seam-joint-ward-state-20261014/operator_transport_check.py`, lines 31–43,
  125–160: source delta/q and the original basis/ordered cocycle.
- `seam-joint-ward-state-20261014/PROOF.txt`, lines 55–80, 152–171,
  1003–1080: full-lattice positivity and ground count, conditional Ward-state
  selection, and actual common even-sector cocycle/factorization.
- `v16-operator-pullback-20261014/contract_index.json`: retained source
  orientation, holonomy and folded physical-time conventions.
