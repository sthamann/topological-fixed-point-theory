# Original source/erasure Ward identity and full return code closure

2026-10-16. Contract `UR.NATIVE.JOINT_OPERATOR.01`, experiments only.
Verdict: EXACT SCOPED IDENTITIES AND EXCLUSIONS; physical selection PARTIAL.
This addendum retains every result of `MIXED-ORDERED-RETURN.md`.

## 1. The first necessary identification

The previous positive matrix is

    K5_ij = <S_beta f_i, S_beta f_j>,   S_beta = o(W_beta).

It is the endpoint compression Pi L_beta Pi, with L_beta=S_beta^dagger S_beta.
Neither positivity nor rank5 implies that L_beta preserves the code.
The necessary full-operation test is (I-Pi)L_beta Pi=0. We test this
before treating K5 as a logical generator or searching longer histories.
No intermediate projection, physical measurement or new coupling is used.

Premises retained together: given E8_1 lattice vacuum/current realization;
original J-positive Cartan/Fock embedding and native quartic code;
original neutral moment branch R_neut and existing root isometry
O=T S^-1/sqrt2. The choice of this physical realization from P1 remains open.

## 2. Fixed source action and erasure

For beta^2=2 and the original real vector v with beta.v=0,

    W_beta = (v.h(-1)) e^beta,  E_beta = e^beta,
    b_a = u_a.beta, c_a = u_a.v, a_a = bar(u_a)(1).

Here u_a=(e_(2a)-i e_(2a+1))/sqrt2, and the actual original erasure
acts as a_a f=partial_a f on the positive Fock polynomials. The derivative
columns partial_a f_i are the previously checked original20-column isometry.
We do not replace this neutral Sym3 carrier by a charged E8 current block.

With Y(A,z)=sum_n A_n z^(-n-1), o(W_beta)=W_beta,1 since wt(W_beta)=2.
The original Heisenberg action on the STATE W_beta gives

    a_a(0) W_beta = bar(b_a) W_beta,
    a_a(1) W_beta = bar(c_a) E_beta,
    a_a(j) W_beta = 0 for j>=2.

The VOA commutator identity at mode1 therefore gives, without a free coefficient,

    [a_a,Y(W_beta,z)] = bar(c_a)Y(E_beta,z)
                       +bar(b_a)zY(W_beta,z),
    [a_a,o(W_beta)] = bar(c_a)E_beta,1 + bar(b_a)W_beta,2.       (1)

These are algebraic Ward terms of the given lattice source. They are not
an additional physical contact/measurement prescription or a P1 time law.
The general commutator identity is Eq.(12) in Tuite and Zuevsky,
arXiv:1106.6149v2, p4 of the PDF (printed p4):
https://arxiv.org/pdf/1106.6149#page=4 . The identification of its states
with the actual native code and erasure is checked here, not asserted by
that external reference.

An explicit independent polynomial comparison checks (1) on every one of
the actual5 code inputs and4 erasures for the original witness below.
If F_k=(-1)^k D_b^k f/k!, G_k=(-1)^k D_b^k D_c f/k!,
B_m=beta(-m), V_m=v(-m), E2=(B1^2+B2)/2, then on quartic inputs:

    E_beta,1 f = F2 + B1 F3 + E2 F4,
    W_beta,2 f = V1(F3+B1 F4) + V2 F4 + G1+B1 G2+E2 G3.

Each term has total conformal grade3, charge beta and oscillator grade2.
The checker differentiates charged output mode1 variables without discarding
other modes. All20 Ward identities hold exactly, including7 nonzero terms.
The universal identity (1) follows analytically, not from finite sampling.

## 3. Analytic full-return leakage witness

Use the already present doubled source root (0,0,0,0,0,0,-2,-2).
The original O and R_neut transport it to

    beta_native=sqrt2 e5, v_native=4sqrt2(e6+e7),
    b2=-i, c3=4-4i, all other b,c zero.

Write P0=sum_j x_j^4 and f0=P0/sqrt96. Define

    h=x0^4/sqrt24,
    q=(2/sqrt3)(h-f0/2)=(x0^4-P0/4)/sqrt18.

The original Fock product gives ||h||=||q||=1, <h,f0>=1/2,
and q is orthogonal to ALL5 code vectors. Because b,c have no0 component,
the FULL original vertex action gives S_beta h=0. Thus

    <h,L_beta f0>=<S_beta h,S_beta f0>=0.

The independently reproduced earlier K00=188/3, however, implies

    <h,Pi L_beta f0>=K00/2=94/3.

Consequently the exact full-operation matrix element is

    <q,L_beta f0>=-K00/sqrt3=-188/(3sqrt3) !=0,               (2)
    ||(I-Pi)L_beta f0||^2 >= 35344/27.                       (3)

This proves failure of code closure without constructing a large full matrix.
All norms are fixed original source/Fock norms; the bound is not a probability.
The sparse check also evaluates <S_beta q_raw,S_beta P0>=-1504 with
||q_raw||^2=18 and ||P0||^2=96, reproducing (2) separately.

For this witness (1) is zero for a0,a1; it is iW_beta,2 for a2 and
(4+4i)E_beta,1 for a3. In particular there is no missing a0 Ward coefficient
that can repair the old decoder comparison. The adjoint also commutes with
a0: beta and v have no0 component. Therefore a0 L_beta f0=0, whereas the
compressed logical decoder has its nonzero188/3 diagonal entry. Compression
is a valid mathematical response operation; physical permission to impose it
has not been derived.

## 4. The unchanged all-root sum also fails full code closure

This is NOT an arbitrary fitted linear combination. Use the same original
root-current INPUT basis and unrenormalized R_neut responses as before:

    L_sum = sum_beta S_beta^dagger S_beta.

Direct exact evaluation on the same q_raw and P0 gives

    sum_beta <S_beta q_raw,S_beta P0> = -50208,
    <q,L_sum f0> = -2092/sqrt3 !=0,                         (4)
    ||(I-Pi)L_sum f0||^2 >= 4376464/3.                      (5)

All4 ORIGINAL family omissions yield exactly the same (4). The computation
rebuilds each original A16 and A_R=4I-A16 before transport by O; it does not
rotate the code to fit a branch. Thus summing all original root responses or
changing among these4 supplied branches does not restore code closure.

As full operators, L_sum and S^dagger S for S=sum_beta a_beta S_beta
are different: the latter also has cross terms with nonzero final lattice
charge gamma-beta. For the NEUTRAL q and f0 used in (4), charge orthogonality
eliminates all beta!=gamma cross terms. The same witness therefore also
excludes the full coherent S^dagger S with |a_beta|=1, whatever the phases.
No arbitrary positive-weight code-closure theorem is claimed here. The
separate previously proved strict trace-gap theorem excludes every such
nontrivial positive profile as the fixed common erasure lift, with its stated
real-transverse/common-normalization premises.

## 5. Consequence for the requested joint derivation

The actual source supplies nontrivial coherent amplitudes and the exact
source/erasure interchange terms (1); these parts are no longer missing.
The paired source return tested here and its original all-root sum are NOT
operations confined to the actual5 code states. Their5x5 matrices therefore
cannot be identified with physical code dynamics without deriving the
additional operation that removes or transports the off-code part.

The next source-derived candidate must specify the actual insertion/history
and explain its code closure jointly with erasure. It must pass

    (I-Pi) U_source Pi = 0

on the FULL source operation before its compressed5x5 spectrum is used to
claim a clock, marking or3+2 selection. If the original process instead
supplies a recovery/protection step, that step and its permission must be
derived; do not silently add Pi, postselection, a detector/reset rule or
independent20/5 normalization. Do not launch arbitrary longer words as a
substitute for identifying this original process.

This is a scoped exclusion of the displayed adjoint pair and original sum,
not a universal E8/P1 no-go. Other nonadjoint/higher/composite source histories,
joint collective processes and different originally derived realizations are
not classified here. Original physical source selection, shared time/clock,
polar selection and physical3+2 remain OPEN. The native common4/5/20
representation, decoder product theorem and prior scoped results remain.

## 6. Reproduction and review

    /opt/homebrew/bin/python3 -OO check_source_erasure_ward.py

Uses existing check_mixed_return.py and original repository constructors;
all validation uses explicit exceptions and survives optimized Python.
Certificate: source_erasure_ward_results.json. Original source pins are in
source_manifest.json. Independent analytic review confirmed (1)-(3), signs,
normalizations and the a0 argument. An independent targeted original-data
check for the original branch gives sum<S_beta x0^4,S_beta P0>=35736;
the stored original sum||S_beta P0||^2=343776 then gives -50208.
This separately confirms the aggregate witness; the other3 family omissions
were not independently rerun. No ledger or physical gate is promoted.
Exact check success is not an unconditional physical derivation.
