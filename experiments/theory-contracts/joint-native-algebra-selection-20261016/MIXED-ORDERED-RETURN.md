# Actual mixed ordered return on the original code and erasure carrier

2026-10-16. UR.NATIVE.JOINT_OPERATOR.01, mathematical/scoped PARTIAL.
No physical source, instrument, time law, coupling or ledger promotion.
This continues the ACTUAL mixed source response, rather than combining
the separate inner-E8 and neutral-quadratic exclusions into a false no-go.

## 1. Use the isometry already supplied by the original source

`cartan_source._glue_deck_character_certificate` supplies both the native
Construction-A simple-root matrix S and standard doubled-root matrix T.
Its actual map is

    O=T S^(-1)/sqrt(2),  O^T O=I8,  OJ=JO.

All240 original roots are checked: O(r_native/sqrt2)=r_source/2, and the
result is exactly the Chevalley root set used by `_actual_current_source`.
No fitted Cartan identification or chosen Dynkin chamber is introduced.
Its physical selection as a P1 source map remains an existing extra premise.
Return beta and v_beta to the native frame with O^T. The original positive
code coordinates remain u_j=(e_2j-i e_(2j+1))/sqrt2.

## 2. Full first operator action, then its actual adjoint return

Reuse W_beta=(v_beta.h(-1))e^beta with
v_beta=-AR beta+(beta^T AR beta)beta/2, beta.v_beta=0.
Let b=U^T beta, c=U^T v_beta; B_n=beta(-n), V_n=v_beta(-n).
For a homogeneous input f of degree d=4 (code) or d=3 (trimer), set

    F_k=(-1)^k D_b^k f/k!,
    G_k=(-1)^k D_b^k(D_c f)/k!,
    E_0=1, E_1=B1,
    E_2=(B1^2+B2)/2,
    E_3=B1^3/6+B1 B2/2+B3/3.

The complete energy-preserving source action is

    o(W_beta)f
      =sum_(m=1..d-1) sum_(k=m+1..d) V_m F_k E_(k-m-1)
       +sum_(k=0..d-1) G_k E_k.

Derivation: normal-order :v.h(z)Y(e^beta,z):. A creation term V_m z^(m-1)
requires k=m+l+1 to contribute to z^-2. The annihilation term is
z^-2 D_c f, requiring k=l. Transversality removes an extra root contraction.
For d=4 all output vectors lie in charge beta, oscillator grade3, total
grade4. The complete PBW carrier has dimension8+64+120=192. For d=3 it
has oscillator grade2, dimension8+36=44, total grade3.

Define the ordered return matrix by the actual positive Fock adjoint:

    (K_beta)_ij=<o(W_beta)f_i, o(W_beta)f_j>.

There is NO code projection inserted between the two source operators.
This is a source mode-insertion amplitude between specified original code
states. Restricting its endpoints to the code does not define a physical
postselection, repeated projection, measurement or transition channel.
Each original root-field cocycle phase cancels between an amplitude and
its adjoint. Different root charges are orthogonal.

## 3. Actual five-code return: nontrivial, but wrong marking and clock

The exact sparse calculation uses the ORIGINAL five code polynomials with
normalization squares(96,8,8,8,1). It keeps original v_beta coefficients.
All194 nonzero first-response directions give positive definite rank5
return matrices; none is scalar or commutes with the original P2. The
remaining46 first responses and return matrices are zero.

The sum over the original unit-norm root-current INPUT basis is

    K_sum=diag(3581,4191,3419,3419,3790).

The W_beta response states have NOT been separately normalized. This is
an original-R-weighted derivative-response Gram, not a selected physical
equal-rate event ensemble. The eight Cartan input derivatives are zero.
Its exact characteristic polynomial is

    (x-4191)(x-3790)(x-3581)(x-3419)^2.

It commutes with original tau, but not original sigma or P2. Exact witnesses:

    [K_sum,sigma]_(1,2)=772,
    [K_sum,P2]_(0,1)=(488 sqrt2-61 sqrt3)/5 !=0.

Indices are zero-based. The original P2 has entry
(P2)_(0,1)=(sqrt3-8sqrt2)/50. This compares ACTUAL operators, not only
degeneracy counts. The response has spectral multiplicities1+1+1+2;
collecting selected eigenvectors would add a new rule, and its unique
two-dimensional eigenspace is not the required original marked plane.

All FOUR original `omitted_family` choices give EXACTLY THE SAME K_sum in
the held native chart. Choosing another supplied family representative
therefore does not repair this marking/clock mismatch. This says nothing
about a different original process with different derived contraction data.

## 4. A directly checked original channel and its erasure continuation

For source doubled root r=(0,0,0,0,0,0,-2,-2), the actual map gives

    beta_native=sqrt2 e5,
    v_native=4sqrt2(e6+e7),
    b2=-i, c3=4-4i (all other components zero).

Its five-return matrix is

    [ 188/3, -40sqrt3 i/3, 0, 0, 0 ]
    [ 40sqrt3 i/3, 80,    0, 0, 0 ]
    [ 0,             0, 48, 0, 0 ]
    [ 0,             0,  0,48, 0 ]
    [ 0,             0,  0, 0,32 ].

The actual native erasure twenty is Sym^3(H+), NOT the differently charged
E8 F20 current block. Its original isometry has Fock columns

    g_(i,a)=partial_a f_i, i=0..4, a=0..3.

All20 columns are checked orthonormal. The factor is fixed: annihilating
one oscillator in degree4 equals sqrt4 times tensor contraction, exactly
the existing W20 decoder's factor2. Thus this is the original intertwiner,
not a freely selected twenty-space of the same dimension.

Apply the SAME source operator o(W_beta) and its adjoint to these original
trimer columns. The computed return K20 has

    rank K5=5,  rank K20=13,
    operator-Schmidt rank K20 (5 x 4 factors)=4,
    ||K20-K5 tensor I4||_HS^2=124096/3.

An analytic witness needs no large matrix: g_(0,0)=x0^3/sqrt6 is killed
because D_b and D_c both vanish on it. Therefore K20 g_(0,0)=0, whereas

    <g_(0,0), (K5 tensor I4) g_(0,0)>=188/3.

This tests the direct decoder continuation of THIS return response. Even
the stronger allowance of an arbitrary register unitary cannot repair its
complete branchwise erasure action with given invertible K5. If

    T(rho tensor I4)T*=K5 rho K5* tensor I4  for every rho,

put S=(K5^(-1) tensor I4)T. Taking rho=I makes S unitary; the full matrix
algebra then makes S=I5 tensor V4. Hence T=K5 tensor V4 has rank20, in
contradiction to actual K20 rank13. Operator-Schmidt rank4 also rules out
any factorization with a register matrix, independent of an overall scale.
No unnormalized amplitude above is a detector probability.

The statement excludes this ordered response branch as a common full-erasure
lift. It does NOT classify arbitrary sums of histories, nonselective
instruments, other source terms, or all combined root/quadratic evolution.
Nor is K_beta itself the unitary quarter R5; constructing a desired
function of it would require an additional source-selected rule.

## 5. Positive root reweighting cannot repair this two-insertion class

The individual mismatch is not just an unfortunate channel choice. For
ANY original root beta and ANY REAL v perpendicular to beta, define
K5(beta,v), K20(beta,v) by this same normalized adjoint-return calculation.
Resolve v into its component parallel to Jbeta and its component orthogonal
to both beta and Jbeta. Exact trace forms give

    tr K5=(45/16)||v_parallel||^2+(203/48)||v_orthogonal||^2,
    tr K20=(23/3)||v_parallel||^2+(32/3)||v_orthogonal||^2,

    4 tr K5-tr K20
      =(43/12)||v_parallel||^2+(25/4)||v_orthogonal||^2 >0

whenever v is nonzero. This is for real transverse v; it is not a claimed
classification of arbitrary complex weight-two fields.

Proof: the actual native G31 acts transitively on all240 roots, preserves J
and the code, and acts unitarily on the entire Sym3 trimer. Therefore it
suffices to use beta=sqrt2 e0. The action is linear in the seven real
coordinates of v, so both traces are real quadratic forms. Computing ALL
entries for v=sqrt2 e_j, j=1..7, gives trace matrices

    C5=diag(45/8,203/24,203/24,203/24,203/24,203/24,203/24),
    C20=diag(46/3,64/3,64/3,64/3,64/3,64/3,64/3).

The exact difference is positive definite; each tested basis vector has
squared norm2. The imaginary antisymmetric Gram part is retained, not
claimed zero, and drops out for REAL v. Root transitivity is checked by
an exact240-element orbit under the original five native generators; code
invariance uses the existing native representation theorem. As a separate
cross-check, these formulas reproduce the traces of ALL240 complete code
return matrices. For the original aggregate they give

    tr K5_sum=18400, tr K20_sum=46848, gap=26752.

Now let W(a)=sum_beta a_beta W_beta, with arbitrary complex input profile
coefficients and the ORIGINAL real transverse v_beta. Orthogonal charge
sectors force Kd(a)=sum_beta |a_beta|^2 Kd,beta, d=5,20. Every nontrivial
profile has the same strictly positive trace gap. It therefore cannot
obey K20(a)=K5(a) tensor I4. Since every nonzero original K5,beta is positive
definite, so is K5(a); positivity plus the preceding commutant argument
also excludes complete unrenormalized branchwise erasure with an arbitrary
register unitary. No positive root weights can repair this response class.

The SAME original source normalization is essential. Inserting an extra
relative normalization on the two carriers would change the tested map.
The theorem covers these positive ADJOINT two-insertion profiles, not
arbitrary nonadjoint outgoing/return pairs, new complex transverse fields,
longer charge-balanced histories, or all physical channels.

## 6. Time and the now decisive next source object

Between the quartic outgoing and return insertions, the existing L0 gives
only e^(-4it)K_beta, or e^(-4t)K_beta for Euclidean propagation. On the
trimer the corresponding factor is e^(-3it) or e^(-3t). These are common
sector amplitudes, not relative logical clock rates. In density-matrix
conjugation scalar phases cancel. The erased degree-one participant and
its physical role cannot be inferred from this grading alone.

The source now demonstrably has nontrivial RETURN effects on the actual
code. The missing claim is no longer merely existence of an operator or
of a code return. The original process must specify a JOINT history or
inserted kernel on code and erasure carrier whose contractions and time
rule preserve their actual intertwiner and select the original mark.
Do not fit root weights, normalize response directions individually,
replace sigma, or insert a code-reprojection rule to make these matrices
match. Even arbitrary positive weights of the same original response
directions are now excluded by the trace theorem. Reopen this class only
with a concrete changed original contraction/history premise. In particular
a nonadjoint return field or a longer charge-balanced history is outside
this theorem; derive it from the original process before evaluating it.
Different combined source histories remain unclassified, not disproved.

## Reproduction and independent check

- check_mixed_return.py: actual original O, all240 code returns, exact ranks,
  sum, original marking/clock comparisons and positive Fock Gram.
- check_mixed_return_extensions.py: all four original family omissions,
  complete20-column erasure isometry and actual20x20 return for the witness.
- check_mixed_return_trace_gap.py: all seven-dimensional quadratic trace
  forms, actual240-root orbit, every original code trace and the strict
  positive-profile erasure obstruction with fixed source normalization.
- Independent analytic review checked the complete mode formula, original
  isometry/normalization, adjoint/phase and time scopes. An independent
  Sympy Laurent expansion (not sparse QQI arithmetic) reproduces the witness
  unnormalized Gram with diagonal(6016,640,384,384,32), entries01=-640i,
  10=640i; it also independently checks the displayed commutator witnesses.
  The240-root aggregate is not claimed to have two independent full runs.
- Independent Laurent/Fock checks also reproduce the two trace-form
  representatives: beta=sqrt2 e0, v=sqrt2 e1 gives traces45/8 and46/3;
  v=sqrt2 e2 gives203/24 and64/3. The actual root stabilizer's permutations
  of the other three complex coordinates and paired i-rotations force
  the transverse isotropy and remove cross terms. Together with the
  original root transitivity this independently supports the universal
  two-coefficient trace-gap formula. No independent full-profile search
  or arbitrary physical normalization theorem is claimed.
- Source: original cartan_source.py, current_block_geometry.py,
  neutral_source_response.py, process.py and prior decoder/code proofs.
- Standard vertex and Fock normalization: Dong and Nagatomo, Section2.2,
  https://arxiv.org/pdf/math/9808088#page=5 . The source-specific return,
  actual map, marking and erasure comparisons are derived here.

The main physical goal stays OPEN/PARTIAL. No unconditional P1-to-source,
3+2, Higgs/fermion coupling, physical clock, alpha or Flavor derivation is
promoted by this source response calculation.
