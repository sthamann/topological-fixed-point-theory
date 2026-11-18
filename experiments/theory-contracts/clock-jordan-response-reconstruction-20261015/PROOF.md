# Joint clock probes close neutral quadratic response reconstruction

15 October 2026. Research result, not verification-ledger promotion.
Overall verdict **PARTIAL**. Algebraic verdict: complete reconstruction in the
specified real symmetric quadratic current sector. No complete TFPT solution.

## Question and reused source

The original `neutral-response-tomography-20261013/PROOF.txt` solves the
three-block scalar problem, but explicitly leaves arbitrary off-diagonal
quadratic matrices open. Its self-correlators cannot distinguish isospectral
matrices with different carrier alignment. Separately,
`native-mark-memory-20261014/clock_probe.py` proves a scalar joint commutant for
four actual clock transports of the same carrier projector. We combine these
two results through the native free-current algebra, without changing a
Hamiltonian or fitting a new source.

Use `tfpt_explorer/neutral_source_response.py`, its original vacuum contraction
`h_i(z) h_j(w) ~ delta_ij/(z-w)^2`, and

    P = diag(1,1,1,1,1,0,0,0),
    Q0=P, Q1=C P C^T, Q2=J P J^T, Q3=C^2 P (C^2)^T.

C and J are the actual `clock_matrices.*.vector` in
`compiler-integral-triality-20260918/certificate.json`. The checker verifies
the source fingerprint, imports the native P, and checks that C and J preserve
the same 240 E8 roots in this basis. Eight indices are not eight physical
quantum states; these clocks are not assumed to be time symmetries.

## 1. Actual current products supply Jordan multiplication

For real symmetric A, define the existing local field

    T_A(z) = (1/2) sum_ij A_ij :h_i(z) h_j(z):.

The standard vacuum Wick rule gives

    T_A(z) T_B(w) = (Tr(AB)/2)/(z-w)^4
                   + 2 T_(A o B)(w)/(z-w)^2 + O((z-w)^(-1)),
    A o B = (AB+BA)/2.

Proof: two complete pairings give Tr(AB)/2. The four single contractions
cancel the prefactor 1/4, leaving h(z)^T AB h(w). At the double pole the fields
are evaluated at w and their symmetric normal product kills the antisymmetric
part of AB. Thus h^T(AB+BA)h/2 = 2 T_(A o B). This holds for arbitrary A,B,
not only for the checked projectors. The simple pole is intentionally not
identified with just a derivative: for noncommuting matrices it can contain
additional antisymmetric weight-three terms.

This uses standard free-boson Wick algebra, not a new general theorem. A
primary reference for the normalizations is Paul Ginsparg, *Applied Conformal
Field Theory*, section 2.3, equations (2.16)-(2.17):
https://arxiv.org/pdf/hep-th/9108028.

Consequently repeated extraction of the weight-two double-pole coefficient
generates exactly the corresponding iterated Jordan products as local fields
in this stipulated source algebra. This is a field/OPE identity, not the claim
that multiplying scalar measurement outputs implements an OPE.

## 2. Four source fields generate all 36 symmetric directions

Start only with Q0,Q1,Q2,Q3. Iteratively adjoin B o Qj when linearly independent.
The explicit rational certificate contains 36 matrices, with degree histogram

    degree 1: 4; degree 2: 6; degree 3: 16; degree 4: 10.

Every matrix is symmetric. Their coordinate matrix is nonsingular over Q:
the checker gives an exact rational two-sided inverse. Full rank modulo 101
is a second algebraic certificate, not a floating-point singular-value test.
An independent reconstruction without the identity seed also attained rank
36 modulo 103; the traceless projections attained rank 35 modulo 107.

The vector space of real symmetric 8 by 8 matrices has dimension 36. These
36 independent Jordan words therefore span the entire space. Identity is
already in their span and is not required as a fifth generator. No theorem
that every irreducible associative algebra has full Jordan closure is assumed;
the explicit certificate resolves this question for the actual four Q here.

The displayed degree-four construction is sufficient. We make no universal
minimality claim about probe count, correlator order, or experimental cost.

## 3. Explicit inverse from mixed responses

Let B_l be the certificate basis. The native mixed two-point coefficient is

    m_l = coefficient of <T_(B_l)(z) T_A(w)> = Tr(B_l A)/2.

Writing a_(ij)=A_ij for i<=j, we obtain m=W a, where

    W_(l,ii)=B_(l,ii)/2;  W_(l,ij)=B_(l,ij) for i<j.

The stored `response_inverse` gives a=W^(-1)m. The checker verifies both
inverse products and reconstructs all 36 canonical symmetric input matrices,
so linearity proves recovery for EVERY real symmetric A, not just a sampled
matrix. The original R is also recovered exactly. If Tr A is prescribed,
the unknown subspace has dimension 35, and 35 independent functionals on that
subspace can be selected. The certificate uses 36 for the general case.

This closes the earlier algebraic off-diagonal reconstruction problem with
an explicit source-defined probe family. Four generator fields are not four
scalar measurements: the inverse uses the generated family of mixed responses.

## 4. Implication for the common-source search

A candidate claiming the SAME marked neutral quadratic source can now be
compared against a complete response vector instead of just norm, spectrum,
or carrier trace. Equal vectors force equal A within this precise class.
This removes an avoidable ambiguity when linking Clock, carrier and neutral
response calculations. It does not create a D5 component in the selected R:
Tr(P5 R)=0 remains exactly unchanged. It does not derive the physical 5/4
susceptibility or repair an excluded scalar P0-only identification.

The extra premise is access to these mixed vacuum/OPE coefficients. Original
P1/P2 have not selected a detector, its instrument, state-transfer map, or
this complete menu. Interacting charged fields, odd CAR/spin sectors,
multitime instruments, unnormalized geometric response and arbitrary
operators are outside the reconstruction theorem. No complete source,
initial state, physical dynamics, or T1-T8 gate is claimed.

## Reproduction and independent review

Run `python3 experiments/theory-contracts/clock-jordan-response-reconstruction-20261015/checker.py`.
`results.json` records original source hashes and exact checks;
`certificate.json` stores the explicit expressions, matrices and inverse.
Independent review checked the OPE coefficient, exact inverse, separate
prime-field ranks, native root-coordinate agreement and the access premise.
The analytic Wick identity and finite certificate have distinct roles.
