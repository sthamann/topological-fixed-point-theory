# F1: a genuinely different P1/quarter origin from infinite exchangeable twisted source algebras

**16 November 2026. New literature-grounded exact finite bridge, F1 OPEN.**
This route starts neither with four original TFPT punctures nor an
assumed E8 or physical spacetime. It asks whether the algebraic
consistency of an INFINITE quantum-event source can select a natural
four-sector commutation law.

## External mathematical classification (do NOT misstate the premises)

Francesco Fidaleo and Elia Vincenzi, *Infinite twisted C*-tensor
product and symmetric states*, arXiv:2609.09494 (8 Sep 2026),
https://arxiv.org/abs/2609.09494.

Theorem 5.1 addresses nondegenerate Hermitian bicharacters
u:Ghat x Ghat -> U(1) for which the self-isotropic subset

    Delta_+ = {g in Ghat : u(g,g)=1}

has trivial u on Delta_+ x Delta_+. Under these STRONG AXIOMS
there are only three classes up to equivalence: ordinary
trivial statistics, Fermi Z2, and a unique Klein four-group
G=Z2 x Z2 twist. The theorem provides the foundations
for an infinite exchangeable (de Finetti) twisted C*-source.
It is a CLASSIFICATION CONDITIONAL on the chosen axioms,
not a claim nature requires exactly this source.

## Exact Klein tensor source and a crucial original form mismatch

For G=F2^2 and x=(x1,x2), y=(y1,y2), the Klein twist is

    u_K(x,y)=(-1)^(x1*y1+x2*y2).

Its 4x4 pairing matrix is a real Hadamard with H H^T=4I,
so the four characters are nondegenerate, not an arbitrary
choice of weights. The self-isotropic subgroup is exactly
Delta_+={00,11}, and u_K=1 on that subgroup squared.
This passes the theorem's unusually strong positivity/
permutation-compatibility restriction.

Compare the nondegenerate alternating pairing

    u_symp(x,y)=(-1)^(x1*y2+x2*y1).

It has self-isotropic subset ALL of G, but is nontrivial
on G x G. It FAILS the additional isotropy-triviality
criterion. Both have H H^T=4I and both are symmetric
bicharacters, yet they define different physical source
structures. In particular u_K(10,10)=-1 but every
alternating bilinear/symplectic form has self-pairing +1
in this exponentiated convention. No pairing-preserving
linear identification of Klein u_K with the ORIGINAL
TFPT alternating F2^4 symplectic bilinear form is possible.
A genuinely new quadratic/polarization/spin interface is
needed. Calling them both "V4" is NOT a derivation.

Cyclic Z4 with nondegenerate u(m,n)=i^(+mn) or i^(-mn)
fails the Hermitian flip criterion already at (1,1).
Thus an abstract quarter-phase cyclic source cannot simply
replace the selected Klein tensor source under these axioms.

## Exact source-to-spin bridge, four possible outcomes

Solve for a quadratic refinement q:G -> mu4,

    q(0)=1,
    q(x+y)=q(x)q(y)u_K(x,y).

Because u_K(e1,e1)=u_K(e2,e2)=-1 and 2e_j=0,
we must have q(e1)^2=q(e2)^2=-1. Therefore

    q(e1) in {+i,-i},  q(e2) in {+i,-i},
    q(e1+e2)=q(e1)q(e2).

All 4 solutions were exhaustively checked against all
16 polarization equations. The four possibilities are:

    (q(0),q(e1),q(e2),q(e1+e2)) | Gauss sum | c_- mod 8
    (1,+i,-i,1)                       1          0
    (1,-i,+i,1)                       1          0
    (1,+i,+i,-1)                     +i          2
    (1,-i,-i,-1)                     -i          6

In particular BOTH physically known double-semion twist
orientations arise from the unique Klein u_K by choosing
one of the two opposite refinements. The exchange
e1 <-> e2 swaps these two choices. Hence the source
cannot pick an orientation without original chiral/
spin information. The other two refinements have
chiral central charge 2 and 6 modulo 8, respectively.

These are valid finite pointed quadratic forms; a
physical association of each to a standalone
anyon topological order additionally requires its
consistent modular/edge realization. The double-semion
example supplies that realization for the c=0 pair.

## Why this is closer to fundamental source selection but still not F1

Under the declared constraints, the four-group/commutation
form no longer needs to be inserted arbitrarily: the
INFINITE algebra classification singles it out among the
nontrivial new options. Quadratic spin then supplies
possible quarter phases from source algebraic relations,
not from guessing six 4-cycle coefficient phases.
This addresses an important missing source TYPE in
the previous finite GHZ/combinatorial P1 compiler.

However the following remain separate:
- The actual TFPT q0 is a refinement of an alternating
  four-dimensional F2 bilinear pairing; a Klein
  nonalternating two-dimensional twist is NOT that source.
- Selecting the quadratic refinement still leaves four
  choices, including two mirror double-semion sectors.
- Anyon topological spins fix the gravitational central
  charge only modulo 8. Adding an E8 invertible
  2+1D layer changes c_- by 8 without changing the
  anyon/spin data. Physical absolute Ref/Eta/framing
  and source Z0 must select the lift.
- An infinite exchangeable C*-algebra does not by itself
  imply a positive-energy half-sided modular inclusion,
  actual physical OS-time, a local chiral 3+1D fermion
  theory, or spin-2 gravity.
- The original R113 whole_source_bundle is still null.
  No actual joint full-field source functional, 
  P_adm, history/field Gram, absolute empty reference
  or first/second weighted cutoff moments were evaluated.
- Original D4 2+4 orbit weights and six coherent permutation
  histories do NOT follow automatically from u_K.

## One precise F1 synthesis experiment

1. Determine from the ORIGINAL boundary/fermion R113 source
   its actual exchange/graded bicharacter between elementary
   separated history words. Check Hermiticity,
   nondegeneracy, Delta_+ and de Finetti isotropy
   on TWO genuinely blocked levels.
2. Test whether the original q0 symplectic code and
   physical P_adm supply an independently defined
   polarized spin lift mapping its physical exchange
   relations to Klein u_K. An algebraic V4 name
   equivalence is INSUFFICIENT.
3. Compute in the same original all-field source the
   signs q(e1), q(e2), their mirror action, the
   original absolute framing phase/eta and Z0,
   as well as ZU,ZY,ZUY and the complete positive
   reflected history Gram. These must select the
   same source state and oriented time.
4. Verify whether the resulting infinite algebra's
   cyclic/separating state and actual half-sided
   modular inclusion yield positive translations.
   Only after this, test 3+1D chirality,
   E8_1 Ward currents, full Higgs Dirichlet form,
   stress/Spin2 and original two-level source pushing.

F1 remains OPEN. This is a substantial alternative
selection program for a four-sector pregeometric
quantum source, not a final theory.

Exact code and regressions:
research/originlab/f1_klein_twisted_source.py
research/originlab/test_f1_klein_twisted_source.py.
